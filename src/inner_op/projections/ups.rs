//! Universal Polar Stereographic as a thin wrapper over stereographic.

use super::stere::Stere;
use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Ups(Stere);

impl PointOp for Ups {
    const NAME: &'static str = "ups";
    const TITLE: &'static str = "Universal Polar Stereographic";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Flag { key: "south" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0.0) },
        OpParameter::Real { key: "lon_0", default: Some(0.0) },
        OpParameter::Real { key: "x_0", default: Some(2_000_000.0) },
        OpParameter::Real { key: "y_0", default: Some(2_000_000.0) },
        OpParameter::Real { key: "k_0", default: Some(0.994) },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self(Stere::build_ups(params)?))
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.fwd(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.inv(coord)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_proj_match;

    #[test]
    fn ups_north_roundtrip() -> Result<(), Error> {
        assert_proj_match(
            "ups ellps=WGS84",
            Coor4D::geo(85., 0., 0., 0.),
            Coor4D::raw(2_000_000.0, 1_444_542.608_617_322_5, 0., 0.),
        )
    }

    #[test]
    fn ups_south_roundtrip() -> Result<(), Error> {
        assert_proj_match(
            "ups south ellps=WGS84",
            Coor4D::geo(-85., 0., 0., 0.),
            Coor4D::raw(2_000_000.0, 2_555_457.391_382_677_6, 0., 0.),
        )
    }
}
