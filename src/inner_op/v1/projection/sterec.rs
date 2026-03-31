//! Polar stereographic variant C as a thin wrapper over stereographic.

use super::stere::Stere;
use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Sterec(Stere);

impl PointOp for Sterec {
    const GAMUT: &'static [OpParameter] = super::stere::GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self(Stere::build_with_variant_c(params, true)?))
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
    use crate::projection::assert_forward_and_roundtrip;

    #[test]
    fn sterec_matches_epsg_example() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "sterec lat_0=-90 lat_ts=-67 lon_0=140 x_0=300000 y_0=200000 ellps=intl",
            Coor4D::geo(-66.60522777777778, 140.0714, 0., 0.),
            Coor4D::raw(303_169.522, 244_055.721, 0., 0.),
            0.02,
            1e-8,
        )
    }
}
