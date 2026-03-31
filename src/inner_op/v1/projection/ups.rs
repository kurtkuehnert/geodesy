//! Universal Polar Stereographic as a thin wrapper over stereographic.

use super::stere::Stere;
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::FRAC_PI_2;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "south" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0.0) },
    OpParameter::Real { key: "lon_0", default: Some(0.0) },
    OpParameter::Real { key: "x_0", default: Some(2_000_000.0) },
    OpParameter::Real { key: "y_0", default: Some(2_000_000.0) },
    OpParameter::Real { key: "k_0", default: Some(0.994) },
];

#[derive(Clone, Copy, Debug)]
pub(crate) struct Ups(Stere);

impl PointOp for Ups {
    type State = Self;
    const GAMUT: &'static [OpParameter] = &GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        let south = params.boolean("south");
        let mut frame = ProjectionFrame::from_params(params);
        frame.lat_0 = if south { -FRAC_PI_2 } else { FRAC_PI_2 };

        Ok(Self(Stere::build_core(params, frame, frame.lat_0, frame.lat_0, false)?))
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        Stere::fwd(&state.0, coord)
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        Stere::inv(&state.0, coord)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_forward_and_roundtrip;

    #[test]
    fn ups_north_roundtrip() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "ups ellps=WGS84",
            Coor4D::geo(85., 0., 0., 0.),
            Coor4D::raw(2_000_000.0, 1_444_542.608_617_322_5, 0., 0.),
            1e-6,
            1e-8,
        )
    }
}
