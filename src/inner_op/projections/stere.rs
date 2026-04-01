//! Stereographic.
//!
//! Attribution:
//! - PROJ 9.8.0 `stere.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/stere.cpp>
//! - PROJ 9.8.0 `stere` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/stere.rst>
//! - Snyder, 1987: *Map Projections: A Working Manual*, pp. 154-163:
//!   <https://pubs.usgs.gov/publication/pp1395>
//! - EPSG method 9829, Polar Stereographic (variant B):
//!   <https://epsg.io/9829-method>
//! - EPSG method 9830, Polar Stereographic (variant C):
//!   <https://epsg.io/9830-method>

use crate::authoring::*;
use crate::projection::{spherical_inverse_equatorial, spherical_inverse_oblique};

const DENOMINATOR_TOLERANCE: f64 = 1e-10;
const POLAR_TOLERANCE: f64 = 1e-8;
const POLAR_SOURCE_TOLERANCE: f64 = 1e-15;

#[rustfmt::skip]
pub const GAMUT: &[OpParameter] = projection_gamut!(
    OpParameter::Real { key: "lat_ts", default: Some(90_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
);

#[derive(Clone, Copy, Debug)]
enum StereMode {
    SouthPole,
    NorthPole,
    Oblique { sin_x1: f64, cos_x1: f64 },
    Equatorial,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Stere {
    frame: ProjectionFrame,
    conformal: ConformalLatitude,
    mode: StereMode,
    akm1: f64,
}

impl Stere {
    fn mode_from_lat_0(conformal: ConformalLatitude, lat_0: f64) -> StereMode {
        let aspect = ProjectionAspect::classify(lat_0, POLAR_TOLERANCE);
        if aspect.is_north_polar() {
            StereMode::NorthPole
        } else if aspect.is_south_polar() {
            StereMode::SouthPole
        } else if aspect.is_equatorial() {
            StereMode::Equatorial
        } else {
            debug_assert!(aspect.is_oblique());
            let xang = conformal.geographic_to_conformal(lat_0);
            StereMode::Oblique {
                sin_x1: xang.sin(),
                cos_x1: xang.cos(),
            }
        }
    }

    fn akm1(
        a: f64,
        k_0: f64,
        lat_0: f64,
        lat_ts: f64,
        conformal: ConformalLatitude,
        mode: StereMode,
        e: f64,
    ) -> f64 {
        if e != 0.0 {
            match mode {
                StereMode::NorthPole | StereMode::SouthPole => {
                    let lat_ts_abs = lat_ts.abs();
                    if (lat_ts_abs - FRAC_PI_2).abs() < POLAR_TOLERANCE {
                        let denominator =
                            ((1.0 + e).powf(1.0 + e) * (1.0 - e).powf(1.0 - e)).sqrt();
                        2.0 * a * k_0 / denominator
                    } else {
                        let sin_ts = lat_ts_abs.sin();
                        let factor = lat_ts_abs.cos() / conformal.ts_from_latitude(lat_ts_abs);
                        a * k_0 * factor / (1.0 - (e * sin_ts).powi(2)).sqrt()
                    }
                }
                StereMode::Equatorial | StereMode::Oblique { .. } => {
                    let te = e * lat_0.sin();
                    2.0 * a * k_0 * lat_0.cos() / (1.0 - te * te).sqrt()
                }
            }
        } else {
            match mode {
                StereMode::Equatorial | StereMode::Oblique { .. } => 2.0 * a * k_0,
                StereMode::NorthPole | StereMode::SouthPole => {
                    if (lat_ts.abs() - FRAC_PI_2).abs() >= DENOMINATOR_TOLERANCE {
                        a * lat_ts.abs().cos() / (FRAC_PI_4 - 0.5 * lat_ts.abs()).tan()
                    } else {
                        2.0 * a * k_0
                    }
                }
            }
        }
    }

    fn variant_c_false_northing(ellps: Ellipsoid, lat_ts: f64, south_pole: bool) -> f64 {
        let a = ellps.semimajor_axis();
        let e = ellps.eccentricity();
        let lat_ts_abs = lat_ts.abs();
        let sin_ts = lat_ts_abs.sin();
        let cos_ts = lat_ts_abs.cos();
        let rho_f = a * cos_ts / (1.0 - (e * sin_ts).powi(2)).sqrt();
        if south_pole { -rho_f } else { rho_f }
    }

    pub(crate) fn build_core(
        params: &ParsedParameters,
        mut frame: ProjectionFrame,
        lat_0: f64,
        lat_ts: f64,
    ) -> Result<Self, Error> {
        let def = &params.name;
        if lat_ts.abs() > FRAC_PI_2 + DENOMINATOR_TOLERANCE {
            return Err(Error::BadParam("lat_ts".to_string(), def.clone()));
        }

        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let conformal = ellps.conformal();
        let e = ellps.eccentricity();
        let mode = Self::mode_from_lat_0(conformal, lat_0);
        let lat_ts = match mode {
            StereMode::NorthPole | StereMode::SouthPole => lat_ts.copysign(lat_0),
            StereMode::Equatorial | StereMode::Oblique { .. } => lat_ts,
        };

        frame.lat_0 = lat_0;

        Ok(Self {
            frame,
            conformal,
            mode,
            akm1: Self::akm1(a, frame.k_0, lat_0, lat_ts, conformal, mode, e),
        })
    }

    pub(crate) fn build_stere(params: &ParsedParameters) -> Result<Self, Error> {
        let lat_0 = params.lat(0);
        let lat_ts = params.real("lat_ts").unwrap_or(FRAC_PI_2);
        let frame = ProjectionFrame::from_params(params);
        Self::build_core(params, frame, lat_0, lat_ts)
    }

    pub(crate) fn build_variant_c(params: &ParsedParameters) -> Result<Self, Error> {
        let lat_0 = params.lat(0);
        let lat_ts = params.real("lat_ts").unwrap_or(FRAC_PI_2);
        let mut frame = ProjectionFrame::from_params(params);
        let aspect = ProjectionAspect::classify(lat_0, POLAR_TOLERANCE);
        if !aspect.is_polar() {
            return Err(Error::Unsupported(
                "sterec is only supported for polar stereographic aspects".into(),
            ));
        }
        let south_pole = aspect.is_south_polar();

        frame.y_0 += Self::variant_c_false_northing(params.ellps(0), lat_ts, south_pole);
        Self::build_core(params, frame, lat_0, lat_ts)
    }

    pub(crate) fn build_ups(params: &ParsedParameters) -> Result<Self, Error> {
        let south = params.boolean("south");
        let mut frame = ProjectionFrame::from_params(params);
        frame.lat_0 = if south { -FRAC_PI_2 } else { FRAC_PI_2 };
        Self::build_core(params, frame, frame.lat_0, frame.lat_0)
    }

    fn fwd_spherical(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let sin_lam = lam.sin();
        let cos_lam = lam.cos();
        let sin_phi = lat.sin();
        let cos_phi = lat.cos();

        match self.mode {
            StereMode::Equatorial => {
                let denominator = 1.0 + cos_phi * cos_lam;
                if denominator <= DENOMINATOR_TOLERANCE {
                    return None;
                }

                let scale = self.akm1 / denominator;
                Some((scale * cos_phi * sin_lam, scale * sin_phi))
            }
            StereMode::Oblique { sin_x1, cos_x1 } => {
                let denominator = 1.0 + sin_x1 * sin_phi + cos_x1 * cos_phi * cos_lam;
                if denominator <= DENOMINATOR_TOLERANCE {
                    return None;
                }

                let scale = self.akm1 / denominator;
                Some((
                    scale * cos_phi * sin_lam,
                    scale * (cos_x1 * sin_phi - sin_x1 * cos_phi * cos_lam),
                ))
            }
            StereMode::NorthPole => {
                let phi = -lat;
                let cos_lam = -cos_lam;
                if (phi - FRAC_PI_2).abs() < POLAR_TOLERANCE {
                    return None;
                }

                let scale = self.akm1 * (FRAC_PI_4 + 0.5 * phi).tan();
                Some((scale * sin_lam, scale * cos_lam))
            }
            StereMode::SouthPole => {
                if (lat - FRAC_PI_2).abs() < POLAR_TOLERANCE {
                    return None;
                }

                let scale = self.akm1 * (FRAC_PI_4 + 0.5 * lat).tan();
                Some((scale * sin_lam, scale * cos_lam))
            }
        }
    }

    fn fwd_ellipsoidal(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let sin_lam = lam.sin();
        let cos_lam = lam.cos();

        match self.mode {
            StereMode::Oblique { sin_x1, cos_x1 } => {
                let xang = self.conformal.geographic_to_conformal(lat);
                let sin_x = xang.sin();
                let cos_x = xang.cos();
                let denominator = cos_x1 * (1.0 + sin_x1 * sin_x + cos_x1 * cos_x * cos_lam);
                if denominator == 0.0 {
                    return None;
                }

                let scale = self.akm1 / denominator;
                let x = scale * cos_x * sin_lam;
                let y = scale * (cos_x1 * sin_x - sin_x1 * cos_x * cos_lam);
                Some((x, y))
            }
            StereMode::Equatorial => {
                let xang = self.conformal.geographic_to_conformal(lat);
                let sin_x = xang.sin();
                let cos_x = xang.cos();
                let denominator = 1.0 + cos_x * cos_lam;
                if denominator == 0.0 {
                    return None;
                }

                let scale = self.akm1 / denominator;
                Some((scale * cos_x * sin_lam, scale * sin_x))
            }
            StereMode::NorthPole => {
                let phi = lat;
                if (phi - FRAC_PI_2).abs() < POLAR_SOURCE_TOLERANCE {
                    return Some((0.0, 0.0));
                }

                let scale = self.akm1 * self.conformal.ts_from_latitude(phi);
                Some((scale * sin_lam, scale * -cos_lam))
            }
            StereMode::SouthPole => {
                let phi = -lat;
                let cos_lam = -cos_lam;
                if (phi - FRAC_PI_2).abs() < POLAR_SOURCE_TOLERANCE {
                    return Some((0.0, 0.0));
                }

                let scale = self.akm1 * self.conformal.ts_from_latitude(phi);
                Some((scale * sin_lam, scale * -cos_lam))
            }
        }
    }

    fn inv_spherical(&self, x_local: f64, y_local: f64) -> Option<(f64, f64)> {
        let rho = x_local.hypot(y_local);
        let c = 2.0 * (rho / self.akm1).atan();

        match self.mode {
            StereMode::Equatorial => Some(spherical_inverse_equatorial(
                x_local,
                y_local,
                rho,
                c,
                DENOMINATOR_TOLERANCE,
                0.0,
            )),
            StereMode::Oblique { sin_x1, cos_x1 } => Some(spherical_inverse_oblique(
                x_local,
                y_local,
                rho,
                c,
                DENOMINATOR_TOLERANCE,
                self.frame.lat_0,
                sin_x1,
                cos_x1,
            )),
            StereMode::NorthPole => {
                let y_local = -y_local;
                let cos_c = c.cos();
                let lat = if rho.abs() <= DENOMINATOR_TOLERANCE {
                    self.frame.lat_0
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
            StereMode::SouthPole => {
                let cos_c = c.cos();
                let lat = if rho.abs() <= DENOMINATOR_TOLERANCE {
                    self.frame.lat_0
                } else {
                    (-cos_c).asin()
                };
                let lam = if x_local == 0.0 && y_local == 0.0 {
                    0.0
                } else {
                    x_local.atan2(y_local)
                };
                Some((lam, lat))
            }
        }
    }

    fn inv_ellipsoidal(&self, x_local: f64, y_local: f64) -> Option<(f64, f64)> {
        let rho = x_local.hypot(y_local);
        match self.mode {
            StereMode::Oblique { sin_x1, cos_x1 } => {
                let tp = 2.0 * (rho * cos_x1).atan2(self.akm1);
                let cos_phi = tp.cos();
                let sin_phi = tp.sin();
                let chi = if rho == 0.0 {
                    (cos_phi * sin_x1).asin()
                } else {
                    (cos_phi * sin_x1 + y_local * sin_phi * cos_x1 / rho).asin()
                };
                let sin_lam = x_local * sin_phi;
                let cos_lam = rho * cos_x1 * cos_phi - y_local * sin_x1 * sin_phi;
                let lam = if sin_lam == 0.0 && cos_lam == 0.0 {
                    0.0
                } else {
                    sin_lam.atan2(cos_lam)
                };
                let lat = self.conformal.conformal_to_geographic(chi);
                Some((lam, lat))
            }
            StereMode::Equatorial => {
                let tp = 2.0 * rho.atan2(self.akm1);
                let cos_phi = tp.cos();
                let sin_phi = tp.sin();
                let chi = if rho == 0.0 {
                    0.0
                } else {
                    (y_local * sin_phi / rho).asin()
                };
                let sin_lam = x_local * sin_phi;
                let cos_lam = rho * cos_phi;
                let lam = if sin_lam == 0.0 && cos_lam == 0.0 {
                    0.0
                } else {
                    sin_lam.atan2(cos_lam)
                };
                let lat = self.conformal.conformal_to_geographic(chi);
                Some((lam, lat))
            }
            StereMode::NorthPole => {
                let lat = self.conformal.latitude_from_ts(rho / self.akm1);
                let y_local = -y_local;
                let lam = if x_local == 0.0 && y_local == 0.0 {
                    0.0
                } else {
                    x_local.atan2(y_local)
                };
                Some((lam, lat))
            }
            StereMode::SouthPole => {
                let lat = -self.conformal.latitude_from_ts(rho / self.akm1);
                let lam = if x_local == 0.0 && y_local == 0.0 {
                    0.0
                } else {
                    x_local.atan2(y_local)
                };
                Some((lam, lat))
            }
        }
    }
}

impl PointOp for Stere {
    const NAME: &'static str = "stere";
    const TITLE: &'static str = "Stereographic";
    const GAMUT: &'static [OpParameter] = GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Self::build_stere(params)
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
    use crate::projection::{assert_op_err, assert_proj_match};

    #[test]
    fn south_polar_inverse_matches_proj() -> Result<(), Error> {
        assert_proj_match(
            "stere lat_0=-90 lat_ts=-71 lon_0=0 ellps=WGS84",
            Coor4D::geo(-75., 30., 0., 0.),
            Coor4D::raw(819_391.619_203_618_1, 1_419_227.915_756_797_2, 0., 0.),
        )
    }

    #[test]
    fn custom_ellipsoid_inverse_matches_proj() -> Result<(), Error> {
        assert_proj_match(
            "stere lat_0=90 lat_ts=70 lon_0=-45 ellps=6378273,298.279411123064",
            Coor4D::geo(80., 45., 0., 0.),
            Coor4D::raw(1_085_943.187_870_963_5, 0.0, 0., 0.),
        )
    }

    #[test]
    fn stere_equatorial_matches_proj_sample() -> Result<(), Error> {
        assert_proj_match(
            "stere ellps=GRS80",
            Coor4D::geo(1., 2., 0., 0.),
            Coor4D::raw(222_644.854_550_117, 110_610.883_474_174, 0., 0.),
        )
    }

    #[test]
    fn stere_oblique_matches_proj_sample() -> Result<(), Error> {
        assert_proj_match(
            "stere lat_0=45 lon_0=10 k_0=0.9996 ellps=GRS80",
            Coor4D::geo(46., 12., 0., 0.),
            Coor4D::raw(154_877.361_182_499_2, 113_025.062_719_405_45, 0., 0.),
        )
    }

    #[test]
    fn sterec_rejects_non_polar_aspects() {
        assert_op_err("sterec lat_0=45 lat_ts=70");
    }
}
