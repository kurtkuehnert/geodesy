//! Polar stereographic variant C as a thin wrapper over stereographic.

use super::stere::Stere;
use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Sterec(Stere);

impl PointOp for Sterec {
    const NAME: &'static str = "sterec";
    const TITLE: &'static str = "Polar Stereographic Variant C";
    const GAMUT: &'static [OpParameter] = super::stere::GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self(Stere::build_variant_c(params)?))
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
    fn sterec_matches_epsg_example() -> Result<(), Error> {
        assert_proj_match(
            "sterec lat_0=-90 lat_ts=-67 lon_0=140 x_0=300000 y_0=200000 ellps=intl",
            Coor4D::geo(-66.60522777777778, 140.0714, 0., 0.),
            Coor4D::raw(303_169.521_856_970_5, 244_055.720_500_675_5, 0., 0.),
        )
    }
}
