//! Lambert Azimuthal Equal Area (EPSG method 9820, IOGP 2019).
//!
//! Attribution:
//! - PROJ 9.8.0 `laea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/laea.cpp>
//! - PROJ 9.8.0 `laea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/laea.rst>
use crate::authoring::*;
use crate::projection::{AuthalicLatitude, AzimuthalAspect, ProjectionFrame, projection_gamut};

use std::f64::consts::FRAC_PI_2;

const ANGULAR_TOLERANCE: f64 = 1e-10;
const POLAR_DOMAIN_TOLERANCE: f64 = 1e-15;

#[rustfmt::skip]
pub const GAMUT: &[OpParameter] = projection_gamut!();

#[derive(Clone, Copy, Debug)]
struct Laea {
    frame: ProjectionFrame,
    authalic: AuthalicLatitude,
    aspect: AzimuthalAspect,
    d: f64,
    x_factor: f64,
    y_factor: f64,
    sin_beta1: f64,
    cos_beta1: f64,
}

impl Laea {}

impl PointOp for Laea {
    type State = Self;
    const GAMUT: &'static [OpParameter] = GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        let frame = ProjectionFrame::from_params(params);

        if frame.lat_0.is_nan() || frame.lat_0.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = AzimuthalAspect::classify(frame.lat_0, ANGULAR_TOLERANCE);
        let ellps = params.ellps(0);
        let authalic = AuthalicLatitude::new(ellps);

        let (d, x_factor, y_factor, sin_beta1, cos_beta1) = match aspect {
            AzimuthalAspect::Polar { .. } => (1.0, 1.0, 1.0, 0.0, 1.0),
            AzimuthalAspect::Oblique => {
                let beta1 = authalic.beta_from_phi(frame.lat_0);
                let (sinphi_0, cosphi_0) = frame.lat_0.sin_cos();
                let (sin_beta1, cos_beta1) = beta1.sin_cos();
                let es = ellps.eccentricity_squared();
                let rq = (0.5 * authalic.q_pole()).sqrt();
                let d = cosphi_0 / ((1.0 - es * sinphi_0 * sinphi_0).sqrt() * rq * cos_beta1);
                (d, rq * d, rq / d, sin_beta1, cos_beta1)
            }
        };

        Ok(Self {
            frame,
            authalic,
            aspect,
            d,
            x_factor,
            y_factor,
            sin_beta1,
            cos_beta1,
        })
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = state.frame.remove_central_meridian(lon);

        let (sin_lam, cos_lam) = lam.sin_cos();
        let beta = state.authalic.beta_from_phi(lat);
        let (sin_beta, cos_beta) = beta.sin_cos();

        let (x_unit, y_unit) = match state.aspect {
            AzimuthalAspect::Oblique => {
                let scale = 1.0 + state.sin_beta1 * sin_beta + state.cos_beta1 * cos_beta * cos_lam;

                if scale <= ANGULAR_TOLERANCE {
                    return None;
                }

                let scale = (2.0 / scale).sqrt();
                (
                    state.x_factor * scale * cos_beta * sin_lam,
                    state.y_factor
                        * scale
                        * (state.cos_beta1 * sin_beta - state.sin_beta1 * cos_beta * cos_lam),
                )
            }
            AzimuthalAspect::Polar { pole_sign } => {
                if state.frame.apply_lat_origin(lat).abs() < ANGULAR_TOLERANCE {
                    return None;
                }

                let scale = (state.authalic.q_pole()
                    - pole_sign * sin_beta * state.authalic.q_pole())
                .sqrt();

                if scale < POLAR_DOMAIN_TOLERANCE {
                    (0.0, 0.0)
                } else {
                    (scale * sin_lam, -pole_sign * scale * cos_lam)
                }
            }
        };

        let x_local = state.frame.a * x_unit;
        let y_local = state.frame.a * y_unit;
        let (x, y) = state.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = state.frame.remove_false_origin(coord[0], coord[1]);
        let x_unit = x_local / state.frame.a / state.d;
        let y_unit = y_local / state.frame.a * state.d;

        let q_pole = state.authalic.q_pole();
        let q = x_unit * x_unit + y_unit * y_unit;
        let rho = x_unit.hypot(y_unit);

        if rho < ANGULAR_TOLERANCE {
            return Some(Coor4D::raw(
                state.frame.lon_0,
                state.frame.lat_0,
                coord[2],
                coord[3],
            ));
        }

        let (lam, sin_beta) = match state.aspect {
            AzimuthalAspect::Oblique => {
                if q > 2.0 * q_pole {
                    return None;
                }

                let cos_c = 1.0 - q / q_pole;
                let sin_c = rho * (2.0 * q_pole - q).sqrt() / q_pole;
                let sin_beta = cos_c * state.sin_beta1 + sin_c * state.cos_beta1 * y_unit / rho;
                let sin_lam = x_unit * sin_c;
                let cos_lam = rho * state.cos_beta1 * cos_c - y_unit * state.sin_beta1 * sin_c;
                let lam = sin_lam.atan2(cos_lam);

                (lam, sin_beta)
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let sin_beta = pole_sign * (1.0 - q / q_pole);
                let lam = x_unit.atan2(-pole_sign * y_unit);

                (lam, sin_beta)
            }
        };

        let lon = state.frame.apply_central_meridian(lam);
        let lat = state.authalic.phi_from_sin_beta(sin_beta)?;
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
    fn laea_spherical_polar_matches_proj() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "laea ellps=6378136.6,0 lat_0=90 lon_0=0 x_0=0 y_0=0",
            Coor4D::geo(20.0, 80.0, 0.0, 0.0),
            Coor4D::raw(7_205_540.644_230_844, -1_270_531.226_169_353, 0.0, 0.0),
            1e-8,
            1e-11,
        )
    }

    #[test]
    fn laea_inverse_rejects_distant_points() -> Result<(), Error> {
        assert_inverse_rejects(
            "laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000",
            Coor4D::raw(1e30, 1e30, 0.0, 0.0),
        )
    }
}
