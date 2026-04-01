//! Albers Equal Area.
//!
//! Attribution:
//! - PROJ 9.8.0 `aea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/aea.cpp>
//! - PROJ 9.8.0 `aea` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/aea.rst>

use crate::authoring::*;
use crate::math::sqrt_checked;
use std::f64::consts::FRAC_PI_2;

const STANDARD_PARALLEL_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Copy, Debug)]
pub(crate) struct AeaInner {
    authalic: AuthalicLatitude,
    a: f64,
    n: f64,
    c: f64,
    rho0: f64,
}

pub(crate) type Aea = Framed<AeaInner>;

impl AeaInner {
    pub(crate) fn with_standard_parallels(
        params: &ParsedParameters,
        phi0: f64,
        phi1: f64,
        phi2: f64,
    ) -> Result<Self, Error> {
        if phi1.abs() > FRAC_PI_2 || phi2.abs() > FRAC_PI_2 {
            return Err(Error::BadParam(
                "lat_1/lat_2".to_string(),
                params.name.clone(),
            ));
        }
        if (phi1 + phi2).abs() < STANDARD_PARALLEL_TOLERANCE {
            return Err(Error::General(
                "Aea: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
            ));
        }

        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let authalic = ellps.authalic();

        let m1 = ellps.meridional_scale(phi1);
        let q0 = authalic.q_from_geographic(phi0);
        let q1 = authalic.q_from_geographic(phi1);

        let n = if (phi1 - phi2).abs() >= STANDARD_PARALLEL_TOLERANCE {
            let m2 = ellps.meridional_scale(phi2);
            let q2 = authalic.q_from_geographic(phi2);

            (m1 * m1 - m2 * m2) / (q2 - q1)
        } else {
            phi1.sin()
        };
        let c = m1 * m1 + n * q1;
        let rho0 = (c - n * q0).sqrt() / n;

        Ok(Self {
            authalic,
            a,
            n,
            c,
            rho0,
        })
    }
}

impl FramedProjection for AeaInner {
    const NAME: &'static str = "aea";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lat_1", default: None },
        OpParameter::Real { key: "lat_2", default: None },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Self::with_standard_parallels(params, params.lat(0), params.lat(1), params.lat(2))
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let q = self.authalic.q_from_geographic(lat);
        let rho = sqrt_checked(self.c - self.n * q)? / self.n;
        let theta = lam * self.n;
        let (sin_theta, cos_theta) = theta.sin_cos();
        let x = self.a * rho * sin_theta;
        let y = self.a * (self.rho0 - rho * cos_theta);
        Some((x, y))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let sign = self.n.signum();
        let rho_sin = sign * x / self.a;
        let rho_cos = sign * (self.rho0 - y / self.a);
        let rho = rho_sin.hypot(rho_cos);
        if rho == 0.0 {
            return Some((0.0, sign * FRAC_PI_2));
        }
        let theta = rho_sin.atan2(rho_cos);
        let lam = theta / self.n;
        let q = (self.c - (rho * self.n).powi(2)) / self.n;
        let lat = self.authalic.geographic_from_q(q)?;
        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aea_matches_proj_ellipsoidal_case() -> Result<(), Error> {
        assert_proj_match(
            "aea lat_0=23 lon_0=-96 lat_1=29.5 lat_2=45.5 x_0=1234 y_0=5678 ellps=GRS80",
            Coor4D::geo(35.0, -75.0, 0.0, 0.0),
            Coor4D::raw(1_886_662.390_542_839_4, 1_541_647.285_801_267_7, 0.0, 0.0),
        )
    }

    #[test]
    fn aea_matches_proj_spherical_case() -> Result<(), Error> {
        assert_proj_match(
            "aea lat_0=40 lon_0=0 lat_1=20 lat_2=60 ellps=6378136.6,0",
            Coor4D::geo(35.0, 10.0, 0.0, 0.0),
            Coor4D::raw(863_038.380_142_686_1, -543_990.840_159_122_7, 0.0, 0.0),
        )
    }

    #[test]
    fn wraps_longitude_difference_across_dateline() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aea lat_0=50 lon_0=-154 lat_1=55 lat_2=65 ellps=clrk66")?;

        let geo = [Coor4D::geo(60., 179., 0., 0.)];
        let projected = [Coor4D::raw(
            -1_459_959.150_054_334_7,
            1_413_239.948_137_74,
            0.,
            0.,
        )];
        let ellps = Ellipsoid::named("clrk66")?;

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(ellps.distance(&operands[0], &geo[0]) < 1e-8);
        Ok(())
    }

    #[test]
    fn aea_rejects_opposite_standard_parallels() {
        assert_op_err("aea lat_1=30 lat_2=-30");
    }
}
