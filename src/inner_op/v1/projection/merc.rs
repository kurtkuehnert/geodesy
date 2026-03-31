//! Mercator.
//!
//! Attribution:
//! - PROJ 9.8.0 `merc.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/merc.cpp>
//! - PROJ 9.8.0 `merc` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/merc.rst>

use crate::authoring::*;
use crate::projection::{ProjectionFrame, projection_gamut};
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Merc {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
}

impl PointOp for Merc {
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = projection_gamut!(
        OpParameter::Real { key: "k_0",    default: Some(1_f64) },
        OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let mut frame = ProjectionFrame::from_params(params);

        if let Some(lat_ts) = params.given_real("lat_ts") {
            if lat_ts.abs() > FRAC_PI_2 {
                return Err(Error::General(
                    "Merc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
                ));
            }

            frame.k_0 = ellps.meridional_scale(lat_ts);
        }

        Ok(Self { frame, ellps })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = self.frame.remove_central_meridian(lon);
        let psi = self.ellps.latitude_geographic_to_isometric(lat);

        let x_local = self.frame.a * self.frame.k_0 * lam;
        let y_local = self.frame.a * self.frame.k_0 * psi;

        let (x, y) = self.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = self.frame.remove_false_origin(coord[0], coord[1]);

        let lam = x_local / (self.frame.a * self.frame.k_0);
        let psi = y_local / (self.frame.a * self.frame.k_0);

        let lon = self.frame.apply_central_meridian(lam);
        let lat = self.ellps.latitude_isometric_to_geographic(psi);
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::{assert_forward_and_roundtrip, assert_inverse, assert_roundtrip};

    #[test]
    fn merc_matches_proj_reference() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "merc",
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::raw(1_335_833.889_519_282_8, 7_326_837.714_873_877, 0.0, 0.0),
            20e-9,
            20e-9,
        )
    }

    #[test]
    fn merc_lat_ts_overrides_k0_when_present() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "merc lat_ts=56",
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::raw(748_713.257_925_886_8, 4_106_573.862_841_270_4, 0.0, 0.0),
            20e-9,
            20e-9,
        )
    }

    #[test]
    fn merc_with_false_origin_and_central_meridian() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "merc lat_ts=-41 lon_0=100 x_0=1234 y_0=5678 ellps=WGS84",
            Coor4D::geo(-33.0, 141.0, 0.0, 0.0),
            Coor4D::raw(3_450_776.589_667_497, -2_920_802.103_023_224_5, 0.0, 0.0),
            20e-6,
            20e-9,
        )
    }

    #[test]
    fn merc_with_parser_normalized_prime_meridian() -> Result<(), Error> {
        assert_roundtrip("merc lon_0=13.5", Coor4D::geo(45.0, 13.5, 0.0, 0.0), 20e-9)
    }

    #[test]
    fn merc_spherical_inverse_matches_proj() -> Result<(), Error> {
        assert_inverse(
            "merc ellps=6356751.9,0",
            Coor4D::raw(10_000.0, 20_000.0, 0.0, 0.0),
            Coor4D::geo(0.180267173822635, 0.090133735615956, 0.0, 0.0),
            1e-12,
        )
    }
}
