//! Albers Equal Area.
//!
//! Attribution:
//! - PROJ 9.8.0 `aea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/aea.cpp>
//! - PROJ 9.8.0 `aea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/aea.rst>

use crate::authoring::*;
use crate::math::sqrt_checked;
use crate::projection::{AuthalicLatitude, ProjectionFrame, projection_gamut};
use std::f64::consts::FRAC_PI_2;

const STANDARD_PARALLEL_TOLERANCE: f64 = 1e-10;

#[rustfmt::skip]
pub const GAMUT: &[OpParameter] = projection_gamut!(
    OpParameter::Real { key: "lat_1", default: None },
    OpParameter::Real { key: "lat_2", default: None },
);

#[derive(Clone, Copy, Debug)]
pub(crate) struct Aea {
    frame: ProjectionFrame,
    authalic: AuthalicLatitude,
    n: f64,
    c: f64,
    rho0: f64,
}

impl Aea {
    pub(crate) fn with_standard_parallels(
        params: &ParsedParameters,
        phi0: f64,
        phi1: f64,
        phi2: f64,
    ) -> Result<Self, Error> {
        if phi1.abs() > FRAC_PI_2 || phi2.abs() > FRAC_PI_2 {
            return Err(Error::BadParam(
                "lat_1/lat_2".to_string(),
                params.name.clone(),
            ));
        }
        if (phi1 + phi2).abs() < STANDARD_PARALLEL_TOLERANCE {
            return Err(Error::General(
                "Aea: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
            ));
        }

        let ellps = params.ellps(0);
        let frame = ProjectionFrame::from_params(params);
        let authalic = AuthalicLatitude::new(ellps);

        let m1 = ellps.meridional_scale(phi1);
        let q0 = authalic.q_from_phi(phi0);
        let q1 = authalic.q_from_phi(phi1);

        let n = if (phi1 - phi2).abs() >= STANDARD_PARALLEL_TOLERANCE {
            let m2 = ellps.meridional_scale(phi2);
            let q2 = authalic.q_from_phi(phi2);

            (m1 * m1 - m2 * m2) / (q2 - q1)
        } else {
            phi1.sin()
        };
        let c = m1 * m1 + n * q1;
        let rho0 = (c - n * q0).sqrt() / n;

        Ok(Self {
            frame,
            authalic,
            n,
            c,
            rho0,
        })
    }
}

impl PointOp for Aea {
    type State = Self;
    const GAMUT: &'static [OpParameter] = GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        Self::with_standard_parallels(params, params.lat(0), params.lat(1), params.lat(2))
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = state.frame.remove_central_meridian(lon);
        let q = state.authalic.q_from_phi(lat);
        let rho = sqrt_checked(state.c - state.n * q)? / state.n;
        let theta = lam * state.n;
        let (sin_theta, cos_theta) = theta.sin_cos();
        let x_local = state.frame.a * rho * sin_theta;
        let y_local = state.frame.a * (state.rho0 - rho * cos_theta);
        let (x, y) = state.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = state.frame.remove_false_origin(coord[0], coord[1]);
        let sign = state.n.signum();
        let rho_sin = sign * x_local / state.frame.a;
        let rho_cos = sign * (state.rho0 - y_local / state.frame.a);
        let rho = rho_sin.hypot(rho_cos);
        if rho == 0.0 {
            return Some(Coor4D::raw(
                state.frame.lon_0,
                sign * FRAC_PI_2,
                coord[2],
                coord[3],
            ));
        }
        let theta = rho_sin.atan2(rho_cos);
        let lam = theta / state.n;
        let q = (state.c - (rho * state.n).powi(2)) / state.n;
        let lon = state.frame.apply_central_meridian(lam);
        let lat = state.authalic.phi_from_q_saturating(q)?;
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

pub fn new(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    Op::point::<Aea>(parameters, ctx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::{assert_forward_and_roundtrip, assert_inverse};

    #[test]
    fn aea_origin_roundtrip() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "aea lat_0=23 lon_0=-96 lat_1=29.5 lat_2=45.5 ellps=GRS80",
            Coor4D::geo(23., -96., 0., 0.),
            Coor4D::raw(0.0, 0.0, 0., 0.),
            1e-8,
            1e-10,
        )
    }

    #[test]
    fn aea_inverse_matches_proj_metric_coordinates() -> Result<(), Error> {
        assert_inverse(
            "aea lat_0=0 lon_0=-120 lat_1=34 lat_2=40.5 x_0=0 y_0=-4000000 ellps=GRS80",
            Coor4D::raw(0.0, -112_982.409_1, 0.0, 0.0),
            Coor4D::geo(37.0, -120.0, 0.0, 0.0),
            1e-8,
        )
    }

    #[test]
    fn aea_spherical_inverse_matches_proj() -> Result<(), Error> {
        assert_inverse(
            "aea lat_0=40 lon_0=0 lat_1=20 lat_2=60 ellps=6378136.6,0",
            Coor4D::raw(0.0, 0.0, 0.0, 0.0),
            Coor4D::geo(40.0, 0.0, 0.0, 0.0),
            1e-10,
        )?;

        assert_inverse(
            "aea lat_0=40 lon_0=0 lat_1=20 lat_2=60 ellps=6378136.6,0",
            Coor4D::raw(10_000.0, 20_000.0, 0.0, 0.0),
            Coor4D::geo(40.169_004_441_322_194, 0.124_940_293_483_244, 0.0, 0.0),
            1e-10,
        )
    }

    #[test]
    fn aea_spherical_roundtrip() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "aea lat_0=40 lon_0=0 lat_1=20 lat_2=60 ellps=6378136.6,0",
            Coor4D::geo(35.0, 10.0, 0.0, 0.0),
            Coor4D::raw(863_038.380_142_686_1, -543_990.840_159_122_7, 0.0, 0.0),
            1e-8,
            1e-10,
        )
    }

    #[test]
    fn wraps_longitude_difference_across_dateline() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aea lat_0=50 lon_0=-154 lat_1=55 lat_2=65 ellps=clrk66")?;

        let geo = [Coor4D::geo(60., 179., 0., 0.)];
        let projected = [Coor4D::raw(
            -1_459_959.150_054_334_7,
            1_413_239.948_137_74,
            0.,
            0.,
        )];
        let ellps = Ellipsoid::named("clrk66")?;

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(ellps.distance(&operands[0], &geo[0]) < 1e-8);
        Ok(())
    }
}
