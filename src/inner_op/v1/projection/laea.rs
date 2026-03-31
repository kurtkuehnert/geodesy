//! Lambert Azimuthal Equal Area (EPSG method 9820, IOGP 2019).
//!
//! Attribution:
//! - PROJ 9.8.0 `laea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/laea.cpp>
//! - PROJ 9.8.0 `laea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/laea.rst>
use crate::authoring::*;
use crate::projection::{AuthalicLatitude, ProjectionAspect, ProjectionFrame};

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const EPS10: f64 = 1e-10;
const POLAR_DOMAIN_TOLERANCE: f64 = 1e-15;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug)]
struct LaeaState {
    frame: ProjectionFrame,
    authalic: AuthalicLatitude,
    aspect: ProjectionAspect,
    qp: f64,
    rq: f64,
    dd: f64,
    xmf: f64,
    ymf: f64,
    sinb1: f64,
    cosb1: f64,
}

#[derive(Clone, Copy, Debug)]
enum PolarInput {
    Spherical { lat: f64 },
    Ellipsoidal,
}

impl LaeaState {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let frame = ProjectionFrame::from_params(params);
        let lat_0 = frame.lat_0;
        if lat_0.is_nan() || lat_0.abs() > FRAC_PI_2 + EPS10 {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = ProjectionAspect::classify(lat_0, EPS10);
        let ellps = params.ellps(0);
        let authalic = AuthalicLatitude::new(ellps);
        let e = ellps.eccentricity();
        let es = ellps.eccentricity_squared();
        let qp = authalic.qp();
        let rq = (0.5 * qp).sqrt();

        let (dd, xmf, ymf, sinb1, cosb1) = match aspect {
            ProjectionAspect::NorthPolar | ProjectionAspect::SouthPolar => {
                (1.0, 1.0, 1.0, 0.0, 1.0)
            }
            ProjectionAspect::Equatorial => (rq.recip(), 1.0, 0.5 * qp, 0.0, 1.0),
            ProjectionAspect::Oblique => {
                let (sinphi_0, cosphi_0) = lat_0.sin_cos();
                let b1 = (ancillary::qs(sinphi_0, e) / qp).asin();
                let (sinb1, cosb1) = b1.sin_cos();
                let dd = cosphi_0 / ((1.0 - es * sinphi_0 * sinphi_0).sqrt() * rq * cosb1);
                (dd, rq * dd, rq / dd, sinb1, cosb1)
            }
        };

        Ok(Self {
            frame,
            authalic,
            aspect,
            qp,
            rq,
            dd,
            xmf,
            ymf,
            sinb1,
            cosb1,
        })
    }

    fn fwd_spherical(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let (sinlam, coslam) = lam.sin_cos();
        let (sinphi, cosphi) = lat.sin_cos();

        self.fwd_common(
            sinlam,
            coslam,
            sinphi,
            cosphi,
            PolarInput::Spherical { lat },
        )
    }

    fn fwd_ellipsoidal(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let (sinlam, coslam) = lam.sin_cos();
        let xi = self.authalic.beta_from_phi(lat);
        let (sinb, cosb) = xi.sin_cos();

        self.fwd_common(sinlam, coslam, sinb, cosb, PolarInput::Ellipsoidal)
    }

    fn fwd_common(
        &self,
        sinlam: f64,
        coslam: f64,
        sinv: f64,
        cosv: f64,
        polar_input: PolarInput,
    ) -> Option<(f64, f64)> {
        match self.aspect {
            ProjectionAspect::Equatorial => {
                let b = 1.0 + cosv * coslam;
                if b <= EPS10 {
                    return None;
                }
                let b = (2.0 / b).sqrt();
                Some((self.xmf * b * cosv * sinlam, self.ymf * b * sinv))
            }
            ProjectionAspect::Oblique => {
                let b = 1.0 + self.sinb1 * sinv + self.cosb1 * cosv * coslam;
                if b <= EPS10 {
                    return None;
                }
                let b = (2.0 / b).sqrt();
                Some((
                    self.xmf * b * cosv * sinlam,
                    self.ymf * b * (self.cosb1 * sinv - self.sinb1 * cosv * coslam),
                ))
            }
            ProjectionAspect::NorthPolar | ProjectionAspect::SouthPolar => {
                let b = match polar_input {
                    PolarInput::Spherical { lat } => {
                        if (lat + self.frame.lat_0).abs() < EPS10 {
                            return None;
                        }

                        2.0 * if self.aspect.is_south_polar() {
                            (FRAC_PI_4 - 0.5 * lat).cos()
                        } else {
                            (FRAC_PI_4 - 0.5 * lat).sin()
                        }
                    }
                    PolarInput::Ellipsoidal => {
                        if (self.frame.lat_0 + self.authalic.phi_from_beta(sinv.asin())).abs()
                            < EPS10
                        {
                            return None;
                        }

                        if self.aspect.is_north_polar() {
                            self.qp - sinv * self.qp
                        } else {
                            self.qp + sinv * self.qp
                        }
                        .sqrt()
                    }
                };

                let b = if b < POLAR_DOMAIN_TOLERANCE { 0.0 } else { b };

                let y = if self.aspect.is_north_polar() {
                    -b * coslam
                } else {
                    b * coslam
                };

                Some((b * sinlam, y))
            }
        }
    }

    fn inv_spherical(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let rho = x.hypot(y);
        let half_rho = 0.5 * rho;
        if half_rho > 1.0 {
            return None;
        }

        let c = 2.0 * half_rho.asin();
        let (sin_c, cos_c) = c.sin_cos();

        let (x, y, lat) = match self.aspect {
            ProjectionAspect::Equatorial => {
                let lat = if rho <= EPS10 {
                    0.0
                } else {
                    (y * sin_c / rho).asin()
                };
                let x = x * sin_c;
                let y = cos_c * rho;
                (x, y, lat)
            }
            ProjectionAspect::Oblique => {
                let lat = if rho <= EPS10 {
                    self.frame.lat_0
                } else {
                    (cos_c * self.sinb1 + y * sin_c * self.cosb1 / rho).asin()
                };
                let x = x * sin_c * self.cosb1;
                let y = (cos_c - lat.sin() * self.sinb1) * rho;
                (x, y, lat)
            }
            ProjectionAspect::NorthPolar | ProjectionAspect::SouthPolar => {
                let (y, lat) = if self.aspect.is_north_polar() {
                    (-y, FRAC_PI_2 - c)
                } else {
                    (y, c - FRAC_PI_2)
                };

                (x, y, lat)
            }
        };

        let lam = if y == 0.0 && x == 0.0 {
            0.0
        } else {
            x.atan2(y)
        };

        Some((self.frame.lon_0 + lam, lat))
    }

    fn inv_ellipsoidal(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let (x, y, rho, sin_c, cos_c) = if self.aspect.is_equatorial() || self.aspect.is_oblique() {
            let x = x / self.dd;
            let y = y * self.dd;
            let rho = x.hypot(y);
            if rho < EPS10 {
                return Some((self.frame.lon_0, self.frame.lat_0));
            }

            let asin_argument = 0.5 * rho / self.rq;
            if asin_argument > 1.0 {
                return None;
            }

            let c = 2.0 * asin_argument.asin();
            let (sin_c, cos_c) = c.sin_cos();
            (x, y, rho, sin_c, cos_c)
        } else {
            (x, y, 0.0, 0.0, 0.0)
        };

        let (x, y, ab) = match self.aspect {
            ProjectionAspect::Equatorial => {
                let x = x * sin_c;
                let ab = y * sin_c / rho;
                let y = rho * cos_c;
                (x, y, ab)
            }
            ProjectionAspect::Oblique => {
                let x = x * sin_c;
                let ab = cos_c * self.sinb1 + y * sin_c * self.cosb1 / rho;
                let y = rho * self.cosb1 * cos_c - y * self.sinb1 * sin_c;
                (x, y, ab)
            }
            ProjectionAspect::NorthPolar | ProjectionAspect::SouthPolar => {
                let y = if self.aspect.is_north_polar() { -y } else { y };
                let q = x * x + y * y;
                if q == 0.0 {
                    return Some((self.frame.lon_0, self.frame.lat_0));
                }

                let mut ab = 1.0 - q / self.qp;
                if self.aspect.is_south_polar() {
                    ab = -ab;
                }
                (x, y, ab)
            }
        };

        if ab.abs() > 1.0 + EPS10 {
            return None;
        }

        let lam = x.atan2(y);
        let lat = self.authalic.phi_from_beta(ab.clamp(-1.0, 1.0).asin());
        Some((self.frame.lon_0 + lam, lat))
    }
}

struct Laea;

impl PointOp for Laea {
    type State = LaeaState;
    const GAMUT: &'static [OpParameter] = &GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        LaeaState::new(params)
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = state.frame.lon_delta(lon);

        let (x, y) = if state.authalic.spherical() {
            state.fwd_spherical(lam, lat)?
        } else {
            state.fwd_ellipsoidal(lam, lat)?
        };

        let (x, y) = state
            .frame
            .apply_false_origin(state.frame.a * x, state.frame.a * y);

        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (x, y) = state.frame.remove_false_origin(coord[0], coord[1]);
        let x = x / state.frame.a;
        let y = y / state.frame.a;

        let (lon, lat) = if state.authalic.spherical() {
            state.inv_spherical(x, y)?
        } else {
            state.inv_ellipsoidal(x, y)?
        };

        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

pub fn new(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    Op::point::<Laea>(parameters, ctx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::{
        assert_forward_and_roundtrip, assert_inverse_rejects, assert_roundtrip,
    };

    #[test]
    fn laea_matches_proj_reference() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000",
            Coor4D::geo(50.0, 5.0, 0.0, 0.0),
            Coor4D::raw(3_962_799.45, 2_999_718.85, 0.0, 0.0),
            0.01,
            1e-11,
        )
    }

    #[test]
    fn laea_roundtrips_origin() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "laea lon_0=10 lat_0=52 x_0=4321000 y_0=3210000",
            Coor4D::geo(52.0, 10.0, 0.0, 0.0),
            Coor4D::raw(4_321_000.0, 3_210_000.0, 0.0, 0.0),
            1e-8,
            1e-11,
        )
    }

    #[test]
    fn laea_polar_roundtrip() -> Result<(), Error> {
        assert_roundtrip(
            "laea ellps=GRS80 lat_0=90 lon_0=0 x_0=0 y_0=0",
            Coor4D::geo(20.0, 80.0, 0.0, 0.0),
            1e-11,
        )?;
        assert_roundtrip(
            "laea ellps=GRS80 lat_0=-90 lon_0=0 x_0=0 y_0=0",
            Coor4D::geo(-45.0, -70.0, 0.0, 0.0),
            1e-11,
        )?;
        Ok(())
    }

    #[test]
    fn laea_inverse_rejects_distant_points() -> Result<(), Error> {
        assert_inverse_rejects(
            "laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000",
            Coor4D::raw(1e30, 1e30, 0.0, 0.0),
        )
    }
}
