//! Mercator.
//!
//! Attribution:
//! - PROJ 9.8.0 `merc.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/merc.cpp>
//! - PROJ 9.8.0 `merc` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/merc.rst>

use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct MercInner {
    a: f64,
    k_0: f64,
    ellps: Ellipsoid,
}

pub(crate) type Merc = Framed<MercInner>;

impl FramedProjection for MercInner {
    const NAME: &'static str = "merc";
    const TITLE: &'static str = "Mercator";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps",  default: Some("GRS80") },
        OpParameter::Real { key: "k_0",    default: Some(1_f64) },
        OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let mut k_0 = params.k(0);

        if let Some(lat_ts) = params.given_real("lat_ts") {
            if lat_ts.abs() > FRAC_PI_2 {
                return Err(Error::General(
                    "Merc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
                ));
            }

            k_0 = ellps.meridional_scale(lat_ts);
        }

        Ok(Self { a, k_0, ellps })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let psi = self.ellps.latitude_geographic_to_isometric(lat);
        Some((self.a * self.k_0 * lam, self.a * self.k_0 * psi))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let lam = x / (self.a * self.k_0);
        let psi = y / (self.a * self.k_0);
        let lat = self.ellps.latitude_isometric_to_geographic(psi);
        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_proj_match;

    #[test]
    fn merc_matches_proj_reference() -> Result<(), Error> {
        assert_proj_match(
            "merc",
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::raw(1_335_833.889_519_282_8, 7_326_837.714_873_877, 0.0, 0.0),
        )
    }

    #[test]
    fn merc_lat_ts_overrides_k0_when_present() -> Result<(), Error> {
        assert_proj_match(
            "merc lat_ts=56",
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::raw(748_713.257_925_886_8, 4_106_573.862_841_270_4, 0.0, 0.0),
        )
    }

    #[test]
    fn merc_with_false_origin_and_central_meridian() -> Result<(), Error> {
        assert_proj_match(
            "merc lat_ts=-41 lon_0=100 x_0=1234 y_0=5678 ellps=WGS84",
            Coor4D::geo(-33.0, 141.0, 0.0, 0.0),
            Coor4D::raw(3_450_776.589_676_943, -2_920_802.103_023_605, 0.0, 0.0),
        )
    }

    #[test]
    fn merc_with_parser_normalized_prime_meridian() -> Result<(), Error> {
        assert_proj_match(
            "merc lon_0=13.5",
            Coor4D::geo(45.0, 13.5, 0.0, 0.0),
            Coor4D::raw(0.0, 5_591_295.918_405_316, 0.0, 0.0),
        )
    }

    #[test]
    fn merc_spherical_inverse_matches_proj() -> Result<(), Error> {
        assert_proj_match(
            "merc ellps=6356751.9,0",
            Coor4D::geo(0.180267173822635, 0.090133735615956, 0.0, 0.0),
            Coor4D::raw(10_000.000_000_000_05, 19_999.999_999_999_98, 0.0, 0.0),
        )
    }
}
