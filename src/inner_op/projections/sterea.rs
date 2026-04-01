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

#[derive(Clone, Copy, Debug)]
pub(crate) struct StereaInner {
    a: f64,
    k_0: f64,
    gauss: Gauss,
}

pub(crate) type Sterea = Framed<StereaInner>;

impl FramedProjection for StereaInner {
    const NAME: &'static str = "sterea";
    const TITLE: &'static str = "Oblique Stereographic Alternative";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "k_0", default: Some(1_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let lat_0 = params.real("lat_0")?;
        let Some(gauss) = Gauss::new(ellps.eccentricity(), lat_0) else {
            return Err(Error::Unsupported("sterea gauss setup failed".into()));
        };

        Ok(Self {
            a: ellps.semimajor_axis(),
            k_0: params.k(0),
            gauss,
        })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let (lam, phi) = self.gauss.forward(lam, lat);
        let (sin_phi, cos_phi) = phi.sin_cos();
        let cos_lam = lam.cos();
        let denom = 1.0 + self.gauss.sinc0 * sin_phi + self.gauss.cosc0 * cos_phi * cos_lam;
        if denom == 0.0 {
            return None;
        }

        let k = self.a * self.k_0 * self.gauss.r2 / denom;
        let x = k * cos_phi * lam.sin();
        let y = k * (self.gauss.cosc0 * sin_phi - self.gauss.sinc0 * cos_phi * cos_lam);
        Some((x, y))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let x = x / (self.a * self.k_0);
        let y = y / (self.a * self.k_0);
        let rho = x.hypot(y);
        let (lam, lat) = if rho == 0.0 {
            (0.0, self.gauss.phic0)
        } else {
            let c = 2.0 * rho.atan2(self.gauss.r2);
            let (sin_c, cos_c) = c.sin_cos();
            let lat = (cos_c * self.gauss.sinc0 + y * sin_c * self.gauss.cosc0 / rho).asin();
            let lon =
                (x * sin_c).atan2(rho * self.gauss.cosc0 * cos_c - y * self.gauss.sinc0 * sin_c);
            (lon, lat)
        };

        self.gauss.inverse(lam, lat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DEFINITION: &str = "sterea lat_0=52.1666666666667 lon_0=19.1666666666667 k_0=0.999714 x_0=500000 y_0=500000 ellps=krass";

    #[test]
    fn sterea_roundtrip_origin() -> Result<(), Error> {
        assert_proj_match(
            DEFINITION,
            Coor4D::geo(52.1666666666667, 19.1666666666667, 0., 0.),
            Coor4D::raw(500000.0, 500000.0, 0., 0.),
        )
    }

    #[test]
    fn sterea_matches_proj_sample() -> Result<(), Error> {
        assert_proj_match(
            DEFINITION,
            Coor4D::geo(49.734897261834, 19.151997685019, 0., 0.),
            Coor4D::raw(498942.3343743173, 229504.5906664739, 0., 0.),
        )
    }
}
