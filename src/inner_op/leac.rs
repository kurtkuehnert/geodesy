//! Lambert Equal Area Conic as a thin wrapper over Albers Equal Area.

use super::aea::AeaInner;
use crate::authoring::*;
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug)]
pub(crate) struct LeacInner(AeaInner);

pub(crate) type Leac = Framed<LeacInner>;

impl FramedProjection for LeacInner {
    const NAME: &'static str = "leac";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lat_1", default: None },
        OpParameter::Flag { key: "south" },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let phi0 = params.lat(0);
        let phi1 = if params.boolean("south") {
            -FRAC_PI_2
        } else {
            FRAC_PI_2
        };
        let phi2 = params.lat(1);
        Ok(Self(AeaInner::with_standard_parallels(
            params, phi0, phi1, phi2,
        )?))
    }

    fn fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        FramedProjection::fwd(&self.0, lam, phi)
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        FramedProjection::inv(&self.0, x, y)
    }
}
