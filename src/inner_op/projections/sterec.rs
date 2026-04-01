//! Polar stereographic variant C as a thin wrapper over stereographic.

use crate::authoring::*;

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Sterec;

impl PointOp for Sterec {
    const NAME: &'static str = "sterec";
    const TITLE: &'static str = "Polar Stereographic Variant C";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lon_0", default: Some(0_f64) },
        OpParameter::Real { key: "x_0", default: Some(0_f64) },
        OpParameter::Real { key: "y_0", default: Some(0_f64) },
        OpParameter::Real { key: "k_0", default: Some(1_f64) },
        OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    ];

    fn build(_params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Err(Error::Unsupported(
            "sterec is temporarily decoupled from stere".into(),
        ))
    }

    fn fwd(&self, _coord: Coor4D) -> Option<Coor4D> {
        None
    }

    fn inv(&self, _coord: Coor4D) -> Option<Coor4D> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_proj_match;

    #[test]
    #[ignore = "sterec is temporarily decoupled from stere during stere rewrite"]
    fn sterec_matches_epsg_example() -> Result<(), Error> {
        assert_proj_match(
            "sterec lat_0=-90 lat_ts=-67 lon_0=140 x_0=300000 y_0=200000 ellps=intl",
            Coor4D::geo(-66.60522777777778, 140.0714, 0., 0.),
            Coor4D::raw(303_169.521_856_970_5, 244_055.720_500_675_5, 0., 0.),
        )
    }
}
