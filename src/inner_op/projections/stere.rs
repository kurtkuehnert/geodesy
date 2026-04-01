//! Stereographic projection, following PROJ 9.8.0:
//! - https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/stere.cpp
//! - https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/stere.rst
//!
//! The stere, sterec, and ups entry points all share the same core setup and
//! pointwise math, with only the wrapper-specific parameters layered on top.
use crate::authoring::*;
use crate::projection::{
    ConformalLatitude, ProjectionAspect, ProjectionFrame, projection_gamut,
    spherical_inverse_equatorial, spherical_inverse_oblique,
};

const DENOMINATOR_TOLERANCE: f64 = 1e-10;
const POLAR_TOLERANCE: f64 = 1e-8;
const POLAR_SOURCE_TOLERANCE: f64 = 1e-15;

#[rustfmt::skip]
pub const GAMUT: &[OpParameter] = projection_gamut!(
    OpParameter::Real { key: "lat_ts", default: Some(90_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
);

#[derive(Clone, Copy, Debug)]
pub(crate) struct Stere {
    frame: ProjectionFrame,
    conformal: ConformalLatitude,
    aspect: ProjectionAspect,
    oblique: Option<(f64, f64)>,
    akm1: f64,
}

impl Stere {
    pub(crate) fn build_core(
        params: &ParsedParameters,
        mut frame: ProjectionFrame,
        lat_0: f64,
        mut lat_ts: f64,
        variant_c: bool,
    ) -> Result<Self, Error> {
        let def = &params.name;
        if lat_ts.abs() > FRAC_PI_2 + DENOMINATOR_TOLERANCE {
            return Err(Error::BadParam("lat_ts".to_string(), def.clone()));
        }

        let aspect = ProjectionAspect::classify(lat_0, POLAR_TOLERANCE);
        if aspect.is_polar() {
            lat_ts = lat_ts.copysign(lat_0);
        }

        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let conformal = ConformalLatitude::new(ellps);
        let e = ellps.eccentricity();
        let k_0 = frame.k_0;

        if variant_c && !aspect.is_polar() {
            return Err(Error::Unsupported(
                "sterec is only supported for polar stereographic aspects".into(),
            ));
        }

        let oblique = if aspect.is_oblique() {
            if e != 0.0 {
                let xang = conformal.geographic_to_conformal(lat_0);
                Some((xang.sin(), xang.cos()))
            } else {
                Some((lat_0.sin(), lat_0.cos()))
            }
        } else {
            None
        };

        let akm1 = if e != 0.0 {
            if aspect.is_polar() {
                let lat_ts_abs = lat_ts.abs();
                if (lat_ts_abs - FRAC_PI_2).abs() < POLAR_TOLERANCE {
                    let numerator = 2.0 * k_0;
                    let denominator = ((1.0 + e).powf(1.0 + e) * (1.0 - e).powf(1.0 - e)).sqrt();
                    a * numerator / denominator
                } else {
                    let sin_ts = lat_ts_abs.sin();
                    let factor = lat_ts_abs.cos() / conformal.ts_from_latitude(lat_ts_abs);
                    a * k_0 * factor / (1.0 - (e * sin_ts).powi(2)).sqrt()
                }
            } else {
                let te = e * lat_0.sin();
                2.0 * a * k_0 * lat_0.cos() / (1.0 - te * te).sqrt()
            }
        } else if !aspect.is_polar() {
            2.0 * a * k_0
        } else if (lat_ts.abs() - FRAC_PI_2).abs() >= DENOMINATOR_TOLERANCE {
            a * lat_ts.abs().cos() / (FRAC_PI_4 - 0.5 * lat_ts.abs()).tan()
        } else {
            2.0 * a * k_0
        };

        if variant_c {
            let lat_ts_abs = lat_ts.abs();
            let sin_ts = lat_ts_abs.sin();
            let cos_ts = lat_ts_abs.cos();
            let rho_f = a * cos_ts / (1.0 - (e * sin_ts).powi(2)).sqrt();
            frame.y_0 += if aspect.is_south_polar() {
                -rho_f
            } else {
                rho_f
            };
        }

        frame.lat_0 = lat_0;

        Ok(Self {
            frame,
            conformal,
            aspect,
            oblique,
            akm1,
        })
    }

    pub(crate) fn build_with_variant_c(
        params: &ParsedParameters,
        variant_c: bool,
    ) -> Result<Self, Error> {
        let lat_0 = params.lat(0);
        let lat_ts = params.real("lat_ts").unwrap_or(FRAC_PI_2);
        let frame = ProjectionFrame::from_params(params);
        Self::build_core(params, frame, lat_0, lat_ts, variant_c)
    }

    fn fwd_spherical(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let sin_lam = lam.sin();
        let mut cos_lam = lam.cos();
        let sin_phi = lat.sin();
        let cos_phi = lat.cos();

        if self.aspect.is_equatorial() {
            let denominator = 1.0 + cos_phi * cos_lam;
            if denominator <= DENOMINATOR_TOLERANCE {
                return None;
            }

            let scale = self.akm1 / denominator;
            return Some((scale * cos_phi * sin_lam, scale * sin_phi));
        }

        if let Some((sin_x1, cos_x1)) = self.oblique {
            let denominator = 1.0 + sin_x1 * sin_phi + cos_x1 * cos_phi * cos_lam;
            if denominator <= DENOMINATOR_TOLERANCE {
                return None;
            }

            let scale = self.akm1 / denominator;
            return Some((
                scale * cos_phi * sin_lam,
                scale * (cos_x1 * sin_phi - sin_x1 * cos_phi * cos_lam),
            ));
        }

        let mut phi = lat;
        if self.aspect.is_north_polar() {
            cos_lam = -cos_lam;
            phi = -phi;
        }
        if (phi - FRAC_PI_2).abs() < POLAR_TOLERANCE {
            return None;
        }

        let scale = self.akm1 * (FRAC_PI_4 + 0.5 * phi).tan();
        Some((scale * sin_lam, scale * cos_lam))
    }

    fn fwd_ellipsoidal(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let sin_lam = lam.sin();
        let mut cos_lam = lam.cos();
        let mut phi = lat;
        let mut sin_x = 0.0;
        let mut cos_x = 0.0;
        if !self.aspect.is_polar() {
            let xang = self.conformal.geographic_to_conformal(phi);
            sin_x = xang.sin();
            cos_x = xang.cos();
        }

        let (mut x_local, y_local) = if let Some((sin_x1, cos_x1)) = self.oblique {
            let denominator = cos_x1 * (1.0 + sin_x1 * sin_x + cos_x1 * cos_x * cos_lam);
            if denominator == 0.0 {
                return None;
            }

            let scale = self.akm1 / denominator;
            (
                scale * cos_x,
                scale * (cos_x1 * sin_x - sin_x1 * cos_x * cos_lam),
            )
        } else if self.aspect.is_equatorial() {
            let denominator = 1.0 + cos_x * cos_lam;
            if denominator == 0.0 {
                return None;
            }

            let scale = self.akm1 / denominator;
            (scale * cos_x, scale * sin_x)
        } else {
            if self.aspect.is_south_polar() {
                phi = -phi;
                cos_lam = -cos_lam;
            }

            let scale = if (phi - FRAC_PI_2).abs() < POLAR_SOURCE_TOLERANCE {
                0.0
            } else {
                self.akm1 * self.conformal.ts_from_latitude(phi)
            };
            (scale, -scale * cos_lam)
        };

        x_local *= sin_lam;
        Some((x_local, y_local))
    }

    fn inv_spherical(&self, x_local: f64, y_local: f64) -> Option<(f64, f64)> {
        let rho = x_local.hypot(y_local);
        let c = 2.0 * (rho / self.akm1).atan();

        if self.aspect.is_equatorial() {
            return Some(spherical_inverse_equatorial(
                x_local,
                y_local,
                rho,
                c,
                DENOMINATOR_TOLERANCE,
                0.0,
            ));
        }

        if let Some((sin_x1, cos_x1)) = self.oblique {
            return Some(spherical_inverse_oblique(
                x_local,
                y_local,
                rho,
                c,
                DENOMINATOR_TOLERANCE,
                self.frame.lat_0,
                sin_x1,
                cos_x1,
            ));
        }

        let mut y_local = y_local;
        let cos_c = c.cos();
        if self.aspect.is_north_polar() {
            y_local = -y_local;
        }

        let lat = if rho.abs() <= DENOMINATOR_TOLERANCE {
            self.frame.lat_0
        } else if self.aspect.is_south_polar() {
            (-cos_c).asin()
        } else {
            cos_c.asin()
        };

        let lam = if x_local == 0.0 && y_local == 0.0 {
            0.0
        } else {
            x_local.atan2(y_local)
        };

        Some((lam, lat))
    }

    fn inv_ellipsoidal(&self, x_local: f64, y_local: f64) -> Option<(f64, f64)> {
        let rho = x_local.hypot(y_local);

        let (chi, sin_lam, cos_lam) = if let Some((sin_x1, cos_x1)) = self.oblique {
            let tp = 2.0 * (rho * cos_x1).atan2(self.akm1);
            let cos_phi = tp.cos();
            let sin_phi = tp.sin();
            let chi = if rho == 0.0 {
                (cos_phi * sin_x1).asin()
            } else {
                (cos_phi * sin_x1 + y_local * sin_phi * cos_x1 / rho).asin()
            };
            (
                chi,
                x_local * sin_phi,
                rho * cos_x1 * cos_phi - y_local * sin_x1 * sin_phi,
            )
        } else if self.aspect.is_equatorial() {
            let tp = 2.0 * rho.atan2(self.akm1);
            let cos_phi = tp.cos();
            let sin_phi = tp.sin();
            let chi = if rho == 0.0 {
                0.0
            } else {
                (y_local * sin_phi / rho).asin()
            };
            (chi, x_local * sin_phi, rho * cos_phi)
        } else {
            let mut y_local = y_local;
            if self.aspect.is_north_polar() {
                y_local = -y_local;
            }
            (
                FRAC_PI_2 - 2.0 * (-rho / self.akm1).atan(),
                x_local,
                y_local,
            )
        };

        let lat = if self.aspect.is_polar() {
            let lat = self.conformal.conformal_to_geographic(PI - chi);
            if self.aspect.is_south_polar() {
                -lat
            } else {
                lat
            }
        } else {
            self.conformal.conformal_to_geographic(chi)
        };

        let lam = if sin_lam == 0.0 && cos_lam == 0.0 {
            0.0
        } else {
            sin_lam.atan2(cos_lam)
        };

        Some((lam, lat))
    }
}

impl PointOp for Stere {
    const NAME: &'static str = "stere";
    const TITLE: &'static str = "Stereographic";
    const GAMUT: &'static [OpParameter] = GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Self::build_with_variant_c(params, false)
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = self.frame.remove_central_meridian_raw(lon);

        let (x_local, y_local) = if self.conformal.spherical() {
            self.fwd_spherical(lam, lat)?
        } else {
            self.fwd_ellipsoidal(lam, lat)?
        };

        let (x, y) = self.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = self.frame.remove_false_origin(coord[0], coord[1]);

        let (lam, lat) = if self.conformal.spherical() {
            self.inv_spherical(x_local, y_local)?
        } else {
            self.inv_ellipsoidal(x_local, y_local)?
        };

        let lon = self.frame.apply_central_meridian(lam);
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::{assert_forward_and_roundtrip, assert_inverse};

    #[test]
    fn south_polar_inverse_matches_proj() -> Result<(), Error> {
        assert_inverse(
            "stere lat_0=-90 lat_ts=-71 lon_0=0 ellps=WGS84",
            Coor4D::raw(819_391.619_181_387, 1_419_227.915_756_798, 0., 0.),
            Coor4D::geo(-75., 30., 0., 0.),
            1e-6,
        )
    }

    #[test]
    fn custom_ellipsoid_inverse_matches_proj() -> Result<(), Error> {
        assert_inverse(
            "stere lat_0=90 lat_ts=70 lon_0=-45 ellps=6378273,298.279411123064",
            Coor4D::raw(1_085_943.187_924_962, 0.0, 0., 0.),
            Coor4D::geo(80., 45., 0., 0.),
            1e-11,
        )
    }

    #[test]
    fn stere_equatorial_matches_proj_sample() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "stere ellps=GRS80",
            Coor4D::geo(1., 2., 0., 0.),
            Coor4D::raw(222_644.854_550_117, 110_610.883_474_174, 0., 0.),
            1e-6,
            1e-10,
        )
    }
}
