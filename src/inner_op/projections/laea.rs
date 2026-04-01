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
pub(crate) struct Laea {
    frame: ProjectionFrame,
    authalic: AuthalicLatitude,
    aspect: AzimuthalAspect,
    d: f64,
    x_factor: f64,
    y_factor: f64,
    sin_beta1: f64,
    cos_beta1: f64,
}

impl Laea {}

impl PointOp for Laea {
    const NAME: &'static str = "laea";
    const TITLE: &'static str = "Lambert Azimuthal Equal Area";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = projection_gamut!();

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let frame = ProjectionFrame::from_params(params);

        if frame.lat_0.is_nan() || frame.lat_0.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = AzimuthalAspect::classify(frame.lat_0, ANGULAR_TOLERANCE);
        let ellps = params.ellps(0);
        let authalic = AuthalicLatitude::new(ellps);

        let (d, x_factor, y_factor, sin_beta1, cos_beta1) = match aspect {
            AzimuthalAspect::Polar { .. } => (1.0, 1.0, 1.0, 0.0, 1.0),
            AzimuthalAspect::Oblique => {
                let beta1 = authalic.geographic_to_authalic(frame.lat_0);
                let (sinphi_0, cosphi_0) = frame.lat_0.sin_cos();
                let (sin_beta1, cos_beta1) = beta1.sin_cos();
                let es = ellps.eccentricity_squared();
                let rq = (0.5 * authalic.q_pole()).sqrt();
                let d = cosphi_0 / ((1.0 - es * sinphi_0 * sinphi_0).sqrt() * rq * cos_beta1);
                (d, rq * d, rq / d, sin_beta1, cos_beta1)
            }
        };

        Ok(Self {
            frame,
            authalic,
            aspect,
            d,
            x_factor,
            y_factor,
            sin_beta1,
            cos_beta1,
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = self.frame.remove_central_meridian(lon);

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
                if self.frame.apply_lat_origin(lat).abs() < ANGULAR_TOLERANCE {
                    return None;
                }

                let scale =
                    (self.authalic.q_pole() - pole_sign * sin_beta * self.authalic.q_pole()).sqrt();

                if scale < POLAR_DOMAIN_TOLERANCE {
                    (0.0, 0.0)
                } else {
                    (scale * sin_lam, -pole_sign * scale * cos_lam)
                }
            }
        };

        let x_local = self.frame.a * x_unit;
        let y_local = self.frame.a * y_unit;
        let (x, y) = self.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = self.frame.remove_false_origin(coord[0], coord[1]);
        let x_unit = x_local / self.frame.a / self.d;
        let y_unit = y_local / self.frame.a * self.d;

        let q_pole = self.authalic.q_pole();
        let q = x_unit * x_unit + y_unit * y_unit;
        let rho = x_unit.hypot(y_unit);

        if rho < ANGULAR_TOLERANCE {
            return Some(Coor4D::raw(
                self.frame.lon_0,
                self.frame.lat_0,
                coord[2],
                coord[3],
            ));
        }

        let (lam, sin_beta) = match self.aspect {
            AzimuthalAspect::Oblique => {
                if q > 2.0 * q_pole {
                    return None;
                }

                let cos_c = 1.0 - q / q_pole;
                let sin_c = rho * (2.0 * q_pole - q).sqrt() / q_pole;
                let sin_beta = cos_c * self.sin_beta1 + sin_c * self.cos_beta1 * y_unit / rho;
                let sin_lam = x_unit * sin_c;
                let cos_lam = rho * self.cos_beta1 * cos_c - y_unit * self.sin_beta1 * sin_c;
                let lam = sin_lam.atan2(cos_lam);

                (lam, sin_beta)
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let sin_beta = pole_sign * (1.0 - q / q_pole);
                let lam = x_unit.atan2(-pole_sign * y_unit);

                (lam, sin_beta)
            }
        };

        if sin_beta.abs() > 1.0 + ANGULAR_TOLERANCE {
            return None;
        }

        let beta = sin_beta.clamp(-1.0, 1.0).asin();
        let lon = self.frame.apply_central_meridian(lam);
        let lat = self.authalic.authalic_to_geographic(beta);
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_proj_match;

    #[test]
    fn laea_matches_proj_reference() -> Result<(), Error> {
        assert_proj_match(
            "laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000",
            Coor4D::geo(50.0, 5.0, 0.0, 0.0),
            Coor4D::raw(3_962_799.450_955_0678, 2_999_718.853_159_564, 0.0, 0.0),
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
            Coor4D::raw(-4_594_471.834_923_999, 1_672_250.990_086_757_1, 0.0, 0.0),
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
