//! Lambert Equal Area Conic as a thin wrapper over Albers Equal Area.

use super::aea::Aea;
use crate::authoring::*;
use crate::projection::projection_gamut;
use std::f64::consts::FRAC_PI_2;

#[rustfmt::skip]
pub const GAMUT: &[OpParameter] = projection_gamut!(
    OpParameter::Flag { key: "south" },
    OpParameter::Real { key: "lat_1", default: None },
);

#[derive(Clone, Copy, Debug)]
struct Leac(Aea);

impl PointOp for Leac {
    type State = Self;
    const GAMUT: &'static [OpParameter] = GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        let phi0 = params.lat(0);
        let phi1 = if params.boolean("south") {
            -FRAC_PI_2
        } else {
            FRAC_PI_2
        };
        let phi2 = params.lat(1);
        Ok(Self(Aea::with_standard_parallels(
            params, phi0, phi1, phi2,
        )?))
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        Aea::fwd(&state.0, coord)
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        Aea::inv(&state.0, coord)
    }
}

pub fn new(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    Op::point::<Leac>(parameters, ctx)
}
