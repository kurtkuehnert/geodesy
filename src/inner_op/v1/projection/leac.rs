//! Lambert Equal Area Conic as a thin wrapper over Albers Equal Area.

use super::aea::Aea;
use crate::authoring::*;
use crate::projection::projection_gamut;
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Leac(Aea);

impl PointOp for Leac {
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = projection_gamut!(
        OpParameter::Flag { key: "south" },
        OpParameter::Real { key: "lat_1", default: None },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
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

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.fwd(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.inv(coord)
    }
}
