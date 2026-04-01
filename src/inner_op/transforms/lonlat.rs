//! Geographic lon/lat with optional prime meridian offset.
//!
//! Attribution:
//! - PROJ 9.8.0 `latlong.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/latlong.cpp>
//! - PROJ 9.8.0 `latlon` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/conversions/latlon.rst>

use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct LonLat {
    lon_0: f64,
}

impl PointOp for LonLat {
    const NAME: &'static str = "lonlat";
    const TITLE: &'static str = "Lat/long (Geodetic alias)";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self {
            lon_0: params.lon(0),
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D::raw(
            coord[0] - self.lon_0,
            coord[1],
            coord[2],
            coord[3],
        ))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D::raw(
            coord[0] + self.lon_0,
            coord[1],
            coord[2],
            coord[3],
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_forward_and_roundtrip;

    #[test]
    fn lonlat_applies_prime_meridian_offset() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "lonlat lon_0=0.00289027777777778",
            Coor4D::geo(38.0, 125.00289027777778, 0.0, 0.0),
            Coor4D::geo(38.0, 125.0, 0.0, 0.0),
            1e-12,
            1e-12,
        )
    }
}
