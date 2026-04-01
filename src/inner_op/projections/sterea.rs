//! Oblique Stereographic Alternative.
//!
//! Attribution:
//! - PROJ 9.8.0 `sterea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/sterea.cpp>
//! - PROJ 9.8.0 `gauss.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/gauss.cpp>
//! - PROJ 9.8.0 `sterea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/sterea.rst>
use crate::authoring::*;
use crate::projection::{Gauss, ProjectionFrame};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Sterea {
    frame: ProjectionFrame,
    gauss: Gauss,
}

impl Sterea {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let frame = ProjectionFrame::from_params(params);
        let ellps = params.ellps(0);
        let Some(gauss) = Gauss::new(ellps.eccentricity(), frame.lat_0) else {
            return Err(Error::Unsupported("sterea gauss setup failed".into()));
        };

        Ok(Self { frame, gauss })
    }
}

impl PointOp for Sterea {
    const NAME: &'static str = "sterea";
    const TITLE: &'static str = "Oblique Stereographic Alternative";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lon_0", default: Some(0_f64) },
        OpParameter::Real { key: "k_0",   default: Some(1_f64) },
        OpParameter::Real { key: "x_0",   default: Some(0_f64) },
        OpParameter::Real { key: "y_0",   default: Some(0_f64) },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Self::new(params)
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let (lam, phi) = self
            .gauss
            .forward(self.frame.remove_central_meridian_raw(lon), lat);
        let sinc = phi.sin();
        let cosc = phi.cos();
        let cosl = lam.cos();
        let denom = 1.0 + self.gauss.sinc0 * sinc + self.gauss.cosc0 * cosc * cosl;
        if denom == 0.0 {
            return None;
        }

        let k = self.frame.a * self.frame.k_0 * self.gauss.r2 / denom;
        let x = k * cosc * lam.sin();
        let y = k * (self.gauss.cosc0 * sinc - self.gauss.sinc0 * cosc * cosl);
        let (x, y) = self.frame.apply_false_origin(x, y);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x, y) = coord.xy();
        let (x_local, y_local) = self.frame.remove_false_origin(x, y);
        let x = x_local / (self.frame.a * self.frame.k_0);
        let y = y_local / (self.frame.a * self.frame.k_0);

        let rho = x.hypot(y);
        let (lon, lat) = if rho != 0.0 {
            let c = 2.0 * rho.atan2(self.gauss.r2);
            let sinc = c.sin();
            let cosc = c.cos();
            let lat = (cosc * self.gauss.sinc0 + y * sinc * self.gauss.cosc0 / rho).asin();
            let lon = (x * sinc).atan2(rho * self.gauss.cosc0 * cosc - y * self.gauss.sinc0 * sinc);
            (lon, lat)
        } else {
            (0.0, self.gauss.phic0)
        };

        let (lon, lat) = self.gauss.inverse(lon, lat)?;
        Some(Coor4D::raw(self.frame.lon_0 + lon, lat, coord[2], coord[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_forward_and_roundtrip;

    const DEFINITION: &str = "sterea lat_0=52.1666666666667 lon_0=19.1666666666667 k_0=0.999714 x_0=500000 y_0=500000 ellps=krass";

    #[test]
    fn sterea_roundtrip_origin() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            DEFINITION,
            Coor4D::geo(52.1666666666667, 19.1666666666667, 0., 0.),
            Coor4D::raw(500000.0, 500000.0, 0., 0.),
            1e-8,
            1e-8,
        )
    }

    #[test]
    fn sterea_matches_proj_sample() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            DEFINITION,
            Coor4D::geo(49.734897261834, 19.151997685019, 0., 0.),
            Coor4D::raw(498942.3343743173, 229504.5906664739, 0., 0.),
            1e-6,
            1e-11,
        )
    }
}
