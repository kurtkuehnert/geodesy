//! Equal Area Cylindrical (CEA).
//!
//! Attribution:
//! - PROJ 9.8.0 `cea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/cea.cpp>
//! - PROJ 9.8.0 `cea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/cea.rst>

use crate::authoring::*;
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Cea {
    a: f64,
    k_0: f64,
    authalic: AuthalicLatitude,
}

impl FramedProjection for Cea {
    const NAME: &'static str = "cea";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = projection_gamut!(
        OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
        OpParameter::Real { key: "k_0",    default: Some(1_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let authalic = ellps.authalic();

        let k_0 = match params.given_real("lat_ts") {
            Some(lat_ts) => {
                if lat_ts.abs() > FRAC_PI_2 {
                    return Err(Error::General(
                        "CEA: Invalid value for lat_ts: |lat_ts| should be <= 90°",
                    ));
                }

                let (sin_ts, cos_ts) = lat_ts.sin_cos();
                cos_ts / (1.0 - ellps.eccentricity_squared() * sin_ts * sin_ts).sqrt()
            }
            None => params.k(0),
        };

        Ok(Self { a, k_0, authalic })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let q = self.authalic.q_from_geographic(lat);

        Some((self.a * self.k_0 * lam, self.a / self.k_0 * 0.5 * q))
    }

    fn inv(&self, x_local: f64, y_local: f64) -> Option<(f64, f64)> {
        let lam = x_local / (self.a * self.k_0);
        let q = 2.0 * y_local * self.k_0 / self.a;

        let lat = self.authalic.geographic_from_q(q)?;
        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cea_matches_proj_ellipsoidal_case() -> Result<(), Error> {
        assert_proj_match(
            "cea ellps=GRS80 lon_0=7 x_0=1234 y_0=5678",
            Coor4D::geo(1.0, 2.0, 0.0, 0.0),
            Coor4D::raw(-555_363.453_966_367_9, 116_246.812_396_267_53, 0.0, 0.0),
        )
    }

    #[test]
    fn cea_matches_proj_spherical_lat_ts_override_case() -> Result<(), Error> {
        assert_proj_match(
            "cea R=6400000 k_0=3 lat_ts=0 lon_0=7 x_0=1234 y_0=5678",
            Coor4D::geo(1.0, 2.0, 0.0, 0.0),
            Coor4D::raw(-557_271.360_638_185_5, 117_373.401_198_614_48, 0.0, 0.0),
        )
    }

    #[test]
    fn cea_rejects_invalid_lat_ts() {
        let mut ctx = Minimal::default();
        assert!(ctx.op("cea lat_ts=91").is_err());
    }
}
