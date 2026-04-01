//! Universal Polar Stereographic as a fixed-parameter wrapper over stereographic.

use super::stere::StereInner;
use crate::authoring::*;

const UPS_FALSE_ORIGIN: f64 = 2_000_000.0;
const UPS_SCALE_FACTOR: f64 = 0.994;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Ups {
    inner: StereInner,
}

impl PointOp for Ups {
    const NAME: &'static str = "ups";
    const TITLE: &'static str = "Universal Polar Stereographic";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Flag { key: "south" },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let lat_0 = params.polar_lat();
        Ok(Self {
            inner: StereInner::build_with(params.ellps(0), lat_0, lat_0, UPS_SCALE_FACTOR)?,
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let (x_local, y_local) = self.inner.fwd(lon, lat)?;
        Some(Coor4D::raw(
            x_local + UPS_FALSE_ORIGIN,
            y_local + UPS_FALSE_ORIGIN,
            coord[2],
            coord[3],
        ))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x, y) = (coord[0] - UPS_FALSE_ORIGIN, coord[1] - UPS_FALSE_ORIGIN);
        let (lon, lat) = self.inner.inv(x, y)?;
        Some(Coor4D::raw(
            angular::normalize_symmetric(lon),
            lat,
            coord[2],
            coord[3],
        ))
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

    #[test]
    fn ups_supports_non_default_ellipsoid() -> Result<(), Error> {
        assert_proj_match(
            "ups ellps=intl",
            Coor4D::geo(85.0, 5.0, 0.0, 0.0),
            Coor4D::raw(2_048_413.890_802_832_6, 1_446_626.695_943_447_3, 0.0, 0.0),
        )
    }
}
