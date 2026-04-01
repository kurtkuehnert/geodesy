//! Lambert Conformal Conic.
//!
//! Attribution:
//! - PROJ 9.8.0 `lcc.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/lcc.cpp>
//! - PROJ 9.8.0 `lcc` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/lcc.rst>

use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct LccInner {
    conformal: ConformalLatitude,
    conic: Conic,
    c: f64,
}

pub(crate) type Lcc = Framed<LccInner>;

impl LccInner {
    fn new(params: &ParsedParameters, phi0: f64, phi1: f64, phi2: f64) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let k_0 = params.k(0);
        let conformal = ellps.conformal();

        let m1 = ellps.meridional_scale(phi1);
        let ts0 = conformal.ts_from_latitude(phi0);
        let ts1 = conformal.ts_from_latitude(phi1);

        let n = Conic::cone_constant(params, phi1, phi2, true, || {
            let m2 = ellps.meridional_scale(phi2);
            let ts2 = conformal.ts_from_latitude(phi2);
            (m1 / m2).ln() / (ts1 / ts2).ln()
        })?;

        let c = m1 * ts1.powf(-n) / n;
        let rho0 = if (phi0.abs() - FRAC_PI_2).abs() < Conic::STANDARD_PARALLEL_TOLERANCE {
            0.0
        } else {
            c * ts0.powf(n)
        };

        Ok(Self {
            conformal,
            conic: Conic::new(a * k_0, n, rho0),
            c,
        })
    }
}

impl FramedProjection for LccInner {
    const NAME: &'static str = "lcc";
    const TITLE: &'static str = "Lambert Conformal Conic";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lat_1", default: None },
        OpParameter::Real { key: "lat_2", default: Some(f64::NAN) },
        OpParameter::Real { key: "k_0", default: Some(1_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let phi1 = params.real("lat_1")?;
        let phi2 = params.given_real("lat_2");
        let phi0 = params
            .given_real("lat_0")
            .unwrap_or_else(|| phi2.map_or(phi1, |_| params.real("lat_0").unwrap()));

        Self::new(params, phi0, phi1, phi2.unwrap_or(phi1))
    }

    fn fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        if (phi.abs() - FRAC_PI_2).abs() < Conic::STANDARD_PARALLEL_TOLERANCE {
            if phi * self.conic.n() <= 0.0 {
                return None;
            }
            return Some(self.conic.polar_point());
        }

        let rho = self.c * self.conformal.ts_from_latitude(phi).powf(self.conic.n());
        Some(self.conic.project(lam, rho))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let Some((lam, rho)) = self.conic.inverse(x, y) else {
            return Some(self.conic.pole());
        };

        let ts = (rho / self.c).powf(1.0 / self.conic.n());
        let lat = self.conformal.latitude_from_ts(ts);

        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lcc_matches_proj_one_standard_parallel_case() -> Result<(), Error> {
        assert_proj_match(
            "lcc lat_1=57 lon_0=12",
            Coor4D::geo(55.0, 10.0, 0.0, 0.0),
            Coor4D::raw(-128_046.472_438_652_24, -220_853.700_160_506_4, 0.0, 0.0),
        )
    }

    #[test]
    fn lcc_matches_proj_two_standard_parallel_case() -> Result<(), Error> {
        assert_proj_match(
            "lcc lat_1=33 lat_2=45 lon_0=10",
            Coor4D::geo(40.0, 12.0, 0.0, 0.0),
            Coor4D::raw(169_863.026_093_938_3, 4_735_925.219_292_451, 0.0, 0.0),
        )
    }

    #[test]
    fn lcc_matches_proj_one_standard_parallel_with_lat_offset() -> Result<(), Error> {
        assert_proj_match(
            "lcc lat_1=39 lat_0=35 lon_0=10",
            Coor4D::geo(40.0, 12.0, 0.0, 0.0),
            Coor4D::raw(170_800.011_728_740_65, 557_172.361_112_929_4, 0.0, 0.0),
        )
    }

    #[test]
    fn lcc_matches_proj_two_standard_parallels_with_lat_offset() -> Result<(), Error> {
        assert_proj_match(
            "lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10",
            Coor4D::geo(40.0, 12.0, 0.0, 0.0),
            Coor4D::raw(169_863.026_093_938_36, 554_155.440_793_916_6, 0.0, 0.0),
        )
    }

    #[test]
    fn lcc_matches_proj_with_false_origin() -> Result<(), Error> {
        assert_proj_match(
            "lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10 x_0=12345 y_0=67890",
            Coor4D::geo(40.0, 12.0, 0.0, 0.0),
            Coor4D::raw(182_208.026_093_938_3, 622_045.440_793_916_6, 0.0, 0.0),
        )
    }

    #[test]
    fn lcc_matches_proj_with_false_origin_and_scaling() -> Result<(), Error> {
        assert_proj_match(
            "lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10 x_0=12345 y_0=67890 k_0=0.99",
            Coor4D::geo(40.0, 12.0, 0.0, 0.0),
            Coor4D::raw(180_509.395_832_998_9, 616_503.886_385_977_5, 0.0, 0.0),
        )
    }

    #[test]
    fn lcc_rejects_opposite_standard_parallels() {
        assert_op_err("lcc lat_1=30 lat_2=-30");
    }

    #[test]
    fn lcc_requires_lat_1() {
        assert_op_err("lcc");
    }

    #[test]
    fn lcc_treats_missing_lat_2_as_tangent_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let tangent = ctx.op("lcc lat_1=30 ellps=GRS80")?;
        let explicit = ctx.op("lcc lat_0=30 lat_1=30 lat_2=30 ellps=GRS80")?;
        let mut lhs = [Coor4D::geo(35.0, 10.0, 0.0, 0.0)];
        let mut rhs = lhs;

        ctx.apply(tangent, Fwd, &mut lhs)?;
        ctx.apply(explicit, Fwd, &mut rhs)?;

        assert!(lhs[0].hypot2(&rhs[0]) < 1e-9);
        Ok(())
    }

    #[test]
    fn lcc_defaults_lat_0_to_lat_1_in_tangent_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let implicit = ctx.op("lcc lat_1=30 ellps=GRS80")?;
        let explicit = ctx.op("lcc lat_0=30 lat_1=30 ellps=GRS80")?;
        let mut lhs = [Coor4D::geo(35.0, 10.0, 0.0, 0.0)];
        let mut rhs = lhs;

        ctx.apply(implicit, Fwd, &mut lhs)?;
        ctx.apply(explicit, Fwd, &mut rhs)?;

        assert!(lhs[0].hypot2(&rhs[0]) < 1e-9);
        Ok(())
    }

    #[test]
    fn lcc_rejects_near_polar_standard_parallel() {
        assert_op_err("lcc lat_1=89.9999999945");
    }
}
