//! Web Mercator / Pseudo Mercator.
//!
//! Attribution:
//! - PROJ 9.8.0 `merc.cpp` (`webmerc` entry point):
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/merc.cpp>
//! - PROJ 9.8.0 `webmerc` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/webmerc.rst>

use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

#[derive(Clone, Copy, Debug)]
pub(crate) struct WebMerc {
    frame: ProjectionFrame,
}

impl PointOp for WebMerc {
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Real { key: "lon_0", default: Some(0.0) },
        OpParameter::Real { key: "x_0", default: Some(0.0) },
        OpParameter::Real { key: "y_0", default: Some(0.0) },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let frame = ProjectionFrame {
            lon_0: params.lon(0),
            lat_0: 0.0,
            x_0: params.x(0),
            y_0: params.y(0),
            k_0: 1.0,
            a: Ellipsoid::named("WGS84")?.semimajor_axis(),
        };
        Ok(Self { frame })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = self.frame.remove_central_meridian(lon);
        let x_local = self.frame.a * lam;
        let y_local = self.frame.a * (FRAC_PI_4 + lat / 2.0).tan().ln();
        let (x, y) = self.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = self.frame.remove_false_origin(coord[0], coord[1]);
        let lam = x_local / self.frame.a;
        let lon = self.frame.apply_central_meridian(lam);
        let lat = FRAC_PI_2 - 2.0 * (-y_local / self.frame.a).exp().atan();
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_forward_and_roundtrip;

    #[test]
    fn webmerc_matches_proj_reference() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "webmerc",
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::raw(1_335_833.889_519_282_8, 7_361_866.113_051_188, 0.0, 0.0),
            1e-8,
            2e-9,
        )
    }

    #[test]
    fn webmerc_supports_central_meridian_and_false_origin() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "webmerc lon_0=10 x_0=100 y_0=200",
            Coor4D::geo(49.0, 12.0, 0.0, 0.0),
            Coor4D::raw(222_738.981_586_547_13, 6_275_061.394_006_576, 0.0, 0.0),
            1e-8,
            2e-9,
        )
    }
}
