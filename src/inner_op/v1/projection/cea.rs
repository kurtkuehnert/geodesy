//! Equal Area Cylindrical (CEA).
//!
//! Attribution:
//! - PROJ 9.8.0 `cea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/cea.cpp>
//! - PROJ 9.8.0 `cea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/cea.rst>
use crate::authoring::*;
use crate::projection::{AuthalicLatitude, ProjectionFrame, projection_gamut};

use std::f64::consts::FRAC_PI_2;

pub const GAMUT: &[OpParameter] = projection_gamut!(
    OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    OpParameter::Real { key: "k_0",    default: Some(1_f64) },
);

#[derive(Clone, Copy, Debug)]
struct CeaState {
    frame: ProjectionFrame,
    authalic: AuthalicLatitude,
}

impl CeaState {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let mut frame = ProjectionFrame::from_params(params);
        let authalic = AuthalicLatitude::new(ellps);

        if let Some(lat_ts) = params.given_real("lat_ts") {
            if lat_ts.abs() > FRAC_PI_2 {
                return Err(Error::General(
                    "CEA: Invalid value for lat_ts: |lat_ts| should be <= 90°",
                ));
            }

            let (sin_ts, cos_ts) = lat_ts.sin_cos();

            if ellps.is_spherical() {
                frame.k_0 = cos_ts;
            } else {
                frame.k_0 = cos_ts / (1.0 - ellps.eccentricity_squared() * sin_ts * sin_ts).sqrt();
            }
        }

        Ok(Self { frame, authalic })
    }
}

struct Cea;

impl PointOp for Cea {
    type State = CeaState;
    const GAMUT: &'static [OpParameter] = GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        CeaState::new(params)
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = state.frame.lon_delta(lon);
        let q = state.authalic.q_from_phi(lat);

        let x_local = state.frame.a * state.frame.k_0 * lam;
        let y_local = state.frame.a / state.frame.k_0 * 0.5 * q;

        let (x, y) = state.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = state.frame.remove_false_origin(coord[0], coord[1]);

        let lam = x_local / (state.frame.a * state.frame.k_0);
        let q = 2.0 * y_local * state.frame.k_0 / state.frame.a;

        let lon = state.frame.apply_lon_delta(lam);
        let lat = state.authalic.phi_from_q(q)?;
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

pub fn new(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    Op::point::<Cea>(parameters, ctx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_forward_and_roundtrip;

    #[test]
    fn cea_matches_proj_gie_ellipsoidal_case() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "cea ellps=GRS80",
            Coor4D::geo(1.0, 2.0, 0.0, 0.0),
            Coor4D::raw(222_638.981_586_547, 110_568.812_396_267, 0.0, 0.0),
            1e-6,
            1e-10,
        )
    }

    #[test]
    fn cea_lat_ts_overrides_k0_when_present() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "cea a=6400000 k_0=3 lat_ts=0",
            Coor4D::geo(1.0, 2.0, 0.0, 0.0),
            Coor4D::raw(223_402.144_255_274, 111_695.401_198_614, 0.0, 0.0),
            1e-6,
            1e-10,
        )
    }
}
