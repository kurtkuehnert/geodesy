//! Lambert Azimuthal Equal Area (EPSG method 9820, IOGP 2019).
//!
//! Attribution:
//! - PROJ 9.8.0 `laea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/laea.cpp>
//! - PROJ 9.8.0 `laea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/laea.rst>
use crate::authoring::*;

const ANGULAR_TOLERANCE: f64 = 1e-10;
const POLAR_DOMAIN_TOLERANCE: f64 = 1e-15;

#[derive(Clone, Copy, Debug)]
pub(crate) struct LaeaInner {
    aspect: AzimuthalAspect,
    authalic: AuthalicLatitude,
    lat_0: f64,
    a: f64,
    q_pole: f64,
    d: f64,
    sin_beta1: f64,
    cos_beta1: f64,
}

pub(crate) type Laea = Framed<LaeaInner>;

impl FramedProjection for LaeaInner {
    const NAME: &'static str = "laea";
    const TITLE: &'static str = "Lambert Azimuthal Equal Area";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let lat_0 = params.lat(0);

        if lat_0.is_nan() || lat_0.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = AzimuthalAspect::classify(lat_0, ANGULAR_TOLERANCE);
        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let authalic = AuthalicLatitude::new(ellps);
        let q_pole = authalic.q_pole();
        let rq = (0.5 * q_pole).sqrt();

        let (d, sin_beta1, cos_beta1) = match aspect {
            AzimuthalAspect::Polar { .. } => (1.0, 0.0, 1.0),
            AzimuthalAspect::Oblique => {
                let beta1 = authalic.geographic_to_authalic(lat_0);
                let (sin_beta1, cos_beta1) = beta1.sin_cos();
                let (sin_phi0, cos_phi0) = lat_0.sin_cos();
                let d = cos_phi0
                    / ((1.0 - ellps.eccentricity_squared() * sin_phi0 * sin_phi0).sqrt()
                        * rq
                        * cos_beta1);
                (d, sin_beta1, cos_beta1)
            }
        };

        Ok(Self {
            a,
            lat_0,
            authalic,
            aspect,
            q_pole,
            d,
            sin_beta1,
            cos_beta1,
        })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let (sin_lam, cos_lam) = lam.sin_cos();
        let q_authalic = self.authalic.q_from_geographic(lat);

        match self.aspect {
            AzimuthalAspect::Polar { pole_sign } => {
                let q_radius = self.q_pole - pole_sign * q_authalic;
                if q_radius >= 2.0 * self.q_pole - ANGULAR_TOLERANCE {
                    return None;
                }
                let rho = self.a * q_radius.sqrt();
                Some(self.aspect.polar_xy(lam, rho))
            }
            AzimuthalAspect::Oblique => {
                let beta = (q_authalic / self.q_pole).asin();
                let (sin_beta, cos_beta) = beta.sin_cos();

                let cos_c = self.sin_beta1 * sin_beta + (self.cos_beta1 * cos_beta * cos_lam);
                if 1.0 + cos_c <= ANGULAR_TOLERANCE {
                    return None;
                }
                let rho = self.a * (self.q_pole / (1.0 + cos_c)).sqrt();
                let x = (rho * self.d) * (cos_beta * sin_lam);
                let y = (rho / self.d)
                    * (self.cos_beta1 * sin_beta - self.sin_beta1 * cos_beta * cos_lam);
                Some((x, y))
            }
        }
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let projected_radius = x.hypot(y);
        if projected_radius < ANGULAR_TOLERANCE {
            return Some((0.0, self.lat_0));
        }

        match self.aspect {
            AzimuthalAspect::Polar { pole_sign } => {
                let q_radius = (x * x + y * y) / (self.a * self.a);

                if q_radius > self.q_pole * (2.0 + POLAR_DOMAIN_TOLERANCE) {
                    return None;
                }

                let q_authalic = pole_sign * (self.q_pole - q_radius);
                let lam = x.atan2(pole_sign * -y);
                let lat = self.authalic.geographic_from_q(q_authalic)?;
                Some((lam, lat))
            }
            AzimuthalAspect::Oblique => {
                let rho = (x / self.d).hypot(self.d * y);
                let radius_ratio = rho / (self.a * (2.0 * self.q_pole).sqrt());
                if radius_ratio > 1.0 + POLAR_DOMAIN_TOLERANCE {
                    return None;
                }

                let c = 2.0 * radius_ratio.clamp(0.0, 1.0).asin();
                let (sin_c, cos_c) = c.sin_cos();
                let sin_beta = cos_c * self.sin_beta1
                    + self.d * y * sin_c * self.cos_beta1 / rho;
                let q_authalic = sin_beta * self.q_pole;

                let lam = (x * sin_c).atan2(
                    self.d * rho * self.cos_beta1 * cos_c
                        - self.d * self.d * y * self.sin_beta1 * sin_c,
                );
                let lat = self.authalic.geographic_from_q(q_authalic)?;
                Some((lam, lat))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn laea_matches_proj_reference() -> Result<(), Error> {
        assert_proj_match(
            "laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000",
            Coor4D::geo(50.0, 5.0, 0.0, 0.0),
            Coor4D::raw(3_962_799.450_955_067_8, 2_999_718.853_159_564, 0.0, 0.0),
        )
    }

    #[test]
    fn laea_roundtrips_origin() -> Result<(), Error> {
        assert_proj_match(
            "laea lon_0=10 lat_0=52 x_0=4321000 y_0=3210000",
            Coor4D::geo(52.0, 10.0, 0.0, 0.0),
            Coor4D::raw(4_321_000.0, 3_210_000.0, 0.0, 0.0),
        )
    }

    #[test]
    fn laea_polar_roundtrip() -> Result<(), Error> {
        assert_proj_match(
            "laea ellps=GRS80 lat_0=90 lon_0=0 x_0=0 y_0=0",
            Coor4D::geo(20.0, 80.0, 0.0, 0.0),
            Coor4D::raw(7_204_871.251_665_172, -1_270_413.194_199_339_2, 0.0, 0.0),
        )?;
        assert_proj_match(
            "laea ellps=GRS80 lat_0=-90 lon_0=0 x_0=0 y_0=0",
            Coor4D::geo(-45.0, -70.0, 0.0, 0.0),
            Coor4D::raw(-4_594_471.834_923_999, 1_672_250.990_086_757, 0.0, 0.0),
        )?;
        Ok(())
    }

    #[test]
    fn laea_spherical_polar_matches_proj() -> Result<(), Error> {
        assert_proj_match(
            "laea ellps=6378136.6,0 lat_0=90 lon_0=0 x_0=0 y_0=0",
            Coor4D::geo(20.0, 80.0, 0.0, 0.0),
            Coor4D::raw(7_205_540.644_230_844, -1_270_531.226_169_353, 0.0, 0.0),
        )
    }

    #[test]
    fn laea_inverse_rejects_distant_points() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000")?;
        let mut operands = [Coor4D::raw(1e30, 1e30, 0.0, 0.0)];

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0][0].is_nan());
        Ok(())
    }

    #[test]
    fn laea_matches_proj_ellipsoidal_equatorial_case() -> Result<(), Error> {
        assert_proj_match(
            "laea ellps=GRS80",
            Coor4D::geo(1.0, 2.0, 0.0, 0.0),
            Coor4D::raw(222_602.471_450_095_18, 110_589.827_224_410_46, 0.0, 0.0),
        )
    }

    #[test]
    fn laea_matches_proj_spherical_oblique_case() -> Result<(), Error> {
        assert_proj_match(
            "laea R=1 lat_0=45",
            Coor4D::geo(45.0, 45.0, 0.0, 0.0),
            Coor4D::raw(0.519_376_687_376_145_1, 0.152_121_909_742_267_33, 0.0, 0.0),
        )
    }

    #[test]
    fn laea_matches_proj_ellipsoidal_polar_case() -> Result<(), Error> {
        assert_proj_match(
            "laea ellps=GRS80 lat_0=90",
            Coor4D::geo(45.0, 0.0, 0.0, 0.0),
            Coor4D::raw(0.0, -4_889_334.802_992_742_5, 0.0, 0.0),
        )
    }

    #[test]
    fn laea_rejects_invalid_lat_0() {
        assert_op_err("laea lat_0=91");
    }
}
