//! Stereographic.
//!
//! Attribution:
//! - PROJ 9.8.0 `stere.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/stere.cpp>
//! - PROJ 9.8.0 `stere` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/stere.rst>

use crate::authoring::*;

const ANGULAR_TOLERANCE: f64 = 1e-10;
const POLAR_SOURCE_TOLERANCE: f64 = 1e-15;

#[derive(Clone, Copy, Debug)]
pub(crate) struct StereInner {
    aspect: AzimuthalAspect,
    conformal: ConformalLatitude,
    lat_0: f64,
    sin_chi_0: f64,
    cos_chi_0: f64,
    rho_factor: f64,
}

pub(crate) type Stere = Framed<StereInner>;

impl StereInner {
    fn build_with(ellps: Ellipsoid, lat_0: f64, lat_ts: f64, k_0: f64) -> Result<Self, Error> {
        if lat_0.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam(
                "lat_0".to_string(),
                "out of range".to_string(),
            ));
        }

        if lat_ts.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam(
                "lat_ts".to_string(),
                "out of range".to_string(),
            ));
        }

        let a = ellps.semimajor_axis();
        let aspect = AzimuthalAspect::classify(lat_0, ANGULAR_TOLERANCE);
        let conformal = ellps.conformal();
        let chi_0 = conformal.geographic_to_conformal(lat_0);
        let (sin_chi_0, cos_chi_0) = chi_0.sin_cos();

        let rho_factor = match aspect {
            AzimuthalAspect::Oblique => {
                let m_0 = ellps.meridional_scale(lat_0);
                2.0 * a * k_0 * m_0 / cos_chi_0
            }
            AzimuthalAspect::Polar { .. } => {
                let signed_lat_ts = lat_ts.abs();
                let m_ts = ellps.meridional_scale(signed_lat_ts);
                let ts = conformal.ts_from_latitude(signed_lat_ts);
                a * k_0 * m_ts / ts
            }
        };

        Ok(Self {
            lat_0,
            aspect,
            conformal,
            sin_chi_0,
            cos_chi_0,
            rho_factor,
        })
    }
}

impl FramedProjection for StereInner {
    const NAME: &'static str = "stere";
    const TITLE: &'static str = "Stereographic";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Flag { key: "south" },
        OpParameter::Text { key: "ellps",  default: Some("GRS80") },
        OpParameter::Real { key: "k_0",    default: Some(1_f64) },
        OpParameter::Real { key: "lat_0",  default: Some(0_f64) },
        OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let lat_0 = params.lat(0);
        let lat_ts = params.given_real("lat_ts").unwrap_or(lat_0);
        let k_0 = params.k(0);
        Self::build_with(ellps, lat_0, lat_ts, k_0)
    }

    fn fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        match self.aspect {
            AzimuthalAspect::Oblique => {
                let chi = self.conformal.geographic_to_conformal(phi);
                let (sin_chi, cos_chi) = chi.sin_cos();
                let (sin_lam, cos_lam) = lam.sin_cos();
                let denom = 1.0 + self.sin_chi_0 * sin_chi + self.cos_chi_0 * cos_chi * cos_lam;
                if denom <= ANGULAR_TOLERANCE {
                    return None;
                }

                let rho = self.rho_factor / denom;
                let x = rho * cos_chi * sin_lam;
                let y = rho * (self.cos_chi_0 * sin_chi - self.sin_chi_0 * cos_chi * cos_lam);
                Some((x, y))
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let signed_phi = pole_sign * phi;
                if (signed_phi - FRAC_PI_2).abs() < POLAR_SOURCE_TOLERANCE {
                    return Some((0.0, 0.0));
                }

                let rho = self.rho_factor * self.conformal.ts_from_latitude(signed_phi);
                Some(self.aspect.polar_xy(lam, rho))
            }
        }
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let rho = x.hypot(y);
        if rho < ANGULAR_TOLERANCE {
            return Some((0.0, self.lat_0));
        }

        match self.aspect {
            AzimuthalAspect::Oblique => {
                let c = 2.0 * (rho / self.rho_factor).atan();
                let (sin_c, cos_c) = c.sin_cos();
                let sin_chi = cos_c * self.sin_chi_0 + y * sin_c * self.cos_chi_0 / rho;
                let chi = sin_chi.clamp(-1.0, 1.0).asin();
                let lam =
                    (x * sin_c).atan2(rho * self.cos_chi_0 * cos_c - y * self.sin_chi_0 * sin_c);
                let lat = self.conformal.conformal_to_geographic(chi);
                Some((lam, lat))
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let lam = self.aspect.polar_lam(x, y);
                let lat = pole_sign * self.conformal.latitude_from_ts(rho / self.rho_factor);
                Some((lam, lat))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stere_matches_proj_oblique_case() -> Result<(), Error> {
        assert_proj_match(
            "stere ellps=GRS80 lat_0=45 lon_0=10 k_0=0.9999",
            Coor4D::geo(50.0, 12.0, 0.0, 0.0),
            Coor4D::raw(143_683.518_553_269_36, 558_126.863_912_303_7, 0.0, 0.0),
        )
    }

    #[test]
    fn stere_matches_proj_equatorial_case() -> Result<(), Error> {
        assert_proj_match(
            "stere ellps=GRS80 lat_0=0 lon_0=5 k_0=1 x_0=100 y_0=200",
            Coor4D::geo(15.0, 10.0, 0.0, 0.0),
            Coor4D::raw(547_504.282_665_646_9, 1_671_859.208_378_212_8, 0.0, 0.0),
        )
    }

    #[test]
    fn stere_matches_proj_north_polar_case() -> Result<(), Error> {
        assert_proj_match(
            "stere ellps=GRS80 lat_0=90 lat_ts=70 lon_0=0",
            Coor4D::geo(80.0, 20.0, 0.0, 0.0),
            Coor4D::raw(371_406.615_760_600_14, -1_020_431.290_238_308_3, 0.0, 0.0),
        )
    }

    #[test]
    fn stere_matches_proj_south_polar_case() -> Result<(), Error> {
        assert_proj_match(
            "stere ellps=GRS80 lat_0=-90 lat_ts=-70 lon_0=0",
            Coor4D::geo(-80.0, -20.0, 0.0, 0.0),
            Coor4D::raw(-371_406.615_760_600_14, 1_020_431.290_238_308_3, 0.0, 0.0),
        )
    }

    #[test]
    fn stere_matches_proj_spherical_oblique_case() -> Result<(), Error> {
        assert_proj_match(
            "stere R=6378136.6 lat_0=45 lon_0=10 k_0=0.9999",
            Coor4D::geo(50.0, 12.0, 0.0, 0.0),
            Coor4D::raw(143_358.809_536_352_6, 558_741.895_045_014_3, 0.0, 0.0),
        )
    }

    #[test]
    fn stere_matches_proj_spherical_polar_case() -> Result<(), Error> {
        assert_proj_match(
            "stere R=6378136.6 lat_0=90 lat_ts=70 lon_0=0",
            Coor4D::geo(80.0, 20.0, 0.0, 0.0),
            Coor4D::raw(370_194.700_049_148_06, -1_017_101.579_186_811_2, 0.0, 0.0),
        )
    }

    #[test]
    fn stere_defaults_lat_ts_to_lat_0_in_polar_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let implicit = ctx.op("stere ellps=GRS80 lat_0=90 lon_0=0")?;
        let explicit = ctx.op("stere ellps=GRS80 lat_0=90 lat_ts=90 lon_0=0")?;
        let mut lhs = [Coor4D::geo(80.0, 20.0, 0.0, 0.0)];
        let mut rhs = lhs;

        ctx.apply(implicit, Fwd, &mut lhs)?;
        ctx.apply(explicit, Fwd, &mut rhs)?;

        assert!(lhs[0].hypot2(&rhs[0]) < 1e-9);
        Ok(())
    }

    #[test]
    fn stere_rejects_invalid_lat_0() {
        assert_op_err("stere lat_0=91");
    }

    #[test]
    fn stere_rejects_invalid_lat_ts() {
        assert_op_err("stere lat_0=90 lat_ts=91");
    }
}
