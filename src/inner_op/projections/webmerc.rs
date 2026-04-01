//! Web Mercator / Pseudo Mercator.
//!
//! Attribution:
//! - PROJ 9.8.0 `merc.cpp` (`webmerc` entry point):
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/merc.cpp>
//! - PROJ 9.8.0 `webmerc` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/webmerc.rst>

use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct WebMercInner {
    a: f64,
}

pub(crate) type WebMerc = Framed<WebMercInner>;

impl FramedProjection for WebMercInner {
    const NAME: &'static str = "webmerc";
    const TITLE: &'static str = "Web Mercator / Pseudo Mercator";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!();

    fn build(_params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self {
            a: Ellipsoid::named("WGS84")?.semimajor_axis(),
        })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let x = self.a * lam;
        let y = self.a * (FRAC_PI_4 + lat / 2.0).tan().ln();
        Some((x, y))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let lam = x / self.a;
        let lat = FRAC_PI_2 - 2.0 * (-y / self.a).exp().atan();
        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_proj_match;

    #[test]
    fn webmerc_matches_proj_reference() -> Result<(), Error> {
        assert_proj_match(
            "webmerc",
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::raw(1_335_833.889_519_282_8, 7_361_866.113_051_188, 0.0, 0.0),
        )
    }

    #[test]
    fn webmerc_supports_central_meridian_and_false_origin() -> Result<(), Error> {
        assert_proj_match(
            "webmerc lon_0=10 x_0=100 y_0=200",
            Coor4D::geo(49.0, 12.0, 0.0, 0.0),
            Coor4D::raw(222_738.981_586_547_1, 6_275_061.394_006_575, 0.0, 0.0),
        )
    }
}
