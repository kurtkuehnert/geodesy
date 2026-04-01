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
    d: f64,
    x_factor: f64,
    y_factor: f64,
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

        let (d, x_factor, y_factor, sin_beta1, cos_beta1) = match aspect {
            AzimuthalAspect::Polar { .. } => (1.0, 1.0, 1.0, 0.0, 1.0),
            AzimuthalAspect::Oblique => {
                let beta1 = authalic.geographic_to_authalic(lat_0);
                let (sinphi_0, cosphi_0) = lat_0.sin_cos();
                let (sin_beta1, cos_beta1) = beta1.sin_cos();
                let es = ellps.eccentricity_squared();
                let rq = (0.5 * authalic.q_pole()).sqrt();
                let d = cosphi_0 / ((1.0 - es * sinphi_0 * sinphi_0).sqrt() * rq * cos_beta1);
                (d, rq * d, rq / d, sin_beta1, cos_beta1)
            }
        };

        Ok(Self {
            a,
            lat_0,
            authalic,
            aspect,
            d,
            x_factor,
            y_factor,
            sin_beta1,
            cos_beta1,
        })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let (sin_lam, cos_lam) = lam.sin_cos();
        let beta = self.authalic.geographic_to_authalic(lat);
        let (sin_beta, cos_beta) = beta.sin_cos();

        let (x_unit, y_unit) = match self.aspect {
            AzimuthalAspect::Oblique => {
                let scale = 1.0 + self.sin_beta1 * sin_beta + self.cos_beta1 * cos_beta * cos_lam;

                if scale <= ANGULAR_TOLERANCE {
                    return None;
                }

                let scale = (2.0 / scale).sqrt();
                (
                    self.x_factor * scale * cos_beta * sin_lam,
                    self.y_factor
                        * scale
                        * (self.cos_beta1 * sin_beta - self.sin_beta1 * cos_beta * cos_lam),
                )
            }
            AzimuthalAspect::Polar { pole_sign } => {
                if (lat + self.lat_0).abs() < ANGULAR_TOLERANCE {
                    return None;
                }

                let scale =
                    (self.authalic.q_pole() - pole_sign * sin_beta * self.authalic.q_pole()).sqrt();

                if scale < POLAR_DOMAIN_TOLERANCE {
                    (0.0, 0.0)
                } else {
                    self.aspect.polar_xy(lam, scale)
                }
            }
        };

        Some((self.a * x_unit, self.a * y_unit))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let x_unit = x / self.a / self.d;
        let y_unit = y / self.a * self.d;

        let q_pole = self.authalic.q_pole();
        let q = x_unit * x_unit + y_unit * y_unit;

        let rho = x_unit.hypot(y_unit);
        if rho < ANGULAR_TOLERANCE {
            return Some((0.0, self.lat_0));
        }

        if q > 2.0 * q_pole + ANGULAR_TOLERANCE {
            return None;
        }

        let (lam, sin_beta) = match self.aspect {
            AzimuthalAspect::Oblique => {
                let cos_c = 1.0 - q / q_pole;
                let sin_c = rho * (2.0 * q_pole - q).sqrt() / q_pole;
                let sin_beta = cos_c * self.sin_beta1 + sin_c * self.cos_beta1 * y_unit / rho;
                let sin_lam = x_unit * sin_c;
                let cos_lam = rho * self.cos_beta1 * cos_c - y_unit * self.sin_beta1 * sin_c;
                let lam = sin_lam.atan2(cos_lam);

                (lam, sin_beta)
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let lam = self.aspect.polar_lam(x_unit, y_unit);
                let sin_beta = pole_sign * (1.0 - q / q_pole);

                (lam, sin_beta)
            }
        };

        let q = sin_beta.clamp(-1.0, 1.0) * q_pole;
        let lat = self.authalic.geographic_from_q(q)?;
        Some((lam, lat))
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
}
