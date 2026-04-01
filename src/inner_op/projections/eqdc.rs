//! Equidistant Conic.
//!
//! Attribution:
//! - PROJ 9.8.0 `eqdc.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/eqdc.cpp>
//! - PROJ 9.8.0 `eqdc` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/eqdc.rst>

use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct EqdcInner {
    meridian: MeridianLatitude,
    conic: Conic,
    c: f64,
}

pub(crate) type Eqdc = Framed<EqdcInner>;

impl EqdcInner {
    fn new(params: &ParsedParameters, phi0: f64, phi1: f64, phi2: f64) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let meridian = MeridianLatitude::new(ellps);
        let a = meridian.semimajor_axis();

        let m1 = ellps.meridional_scale(phi1);
        let ml0 = meridian.distance_from_geographic(phi0) / a;
        let ml1 = meridian.distance_from_geographic(phi1) / a;

        let n = Conic::cone_constant(params, phi1, phi2, false, || {
            let m2 = ellps.meridional_scale(phi2);
            let ml2 = meridian.distance_from_geographic(phi2) / a;
            (m1 - m2) / (ml2 - ml1)
        })?;

        let c = ml1 + m1 / n;
        let rho0 = c - ml0;

        Ok(Self {
            meridian,
            conic: Conic::new(a, n, rho0),
            c,
        })
    }
}

impl FramedProjection for EqdcInner {
    const NAME: &'static str = "eqdc";
    const TITLE: &'static str = "Equidistant Conic";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lat_1", default: None },
        OpParameter::Real { key: "lat_2", default: Some(f64::NAN) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let phi0 = params.real("lat_0")?;
        let phi1 = params.real("lat_1")?;
        let phi2 = params.given_real("lat_2").unwrap_or(phi1);
        Self::new(params, phi0, phi1, phi2)
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let rho = self.c - self.meridian.distance_from_geographic(lat) / self.meridian.semimajor_axis();
        Some(self.conic.project(lam, rho))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let Some((lam, rho)) = self.conic.inverse(x, y) else {
            return Some(self.conic.pole());
        };

        let lat = self
            .meridian
            .geographic_from_distance((self.c - rho) * self.meridian.semimajor_axis());
        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eqdc_matches_proj_ellipsoidal_case() -> Result<(), Error> {
        assert_proj_match(
            "eqdc ellps=GRS80 lat_1=0.5 lat_2=2",
            Coor4D::geo(1.0, 2.0, 0.0, 0.0),
            Coor4D::raw(222_588.440_269_286, 110_659.134_907_347, 0.0, 0.0),
        )
    }

    #[test]
    fn eqdc_matches_proj_spherical_case() -> Result<(), Error> {
        assert_proj_match(
            "eqdc ellps=6378136.6,0 lat_0=10 lat_1=20 lat_2=60 lon_0=5",
            Coor4D::geo(35.0, 15.0, 0.0, 0.0),
            Coor4D::raw(860_776.220_460_427_5, 2_830_344.409_203_822_7, 0.0, 0.0),
        )
    }

    #[test]
    fn eqdc_rejects_opposite_standard_parallels() {
        assert_op_err("eqdc lat_1=30 lat_2=-30");
    }

    #[test]
    fn eqdc_rejects_zero_n_tangent_case() {
        assert_op_err("eqdc lat_1=0 lat_2=0");
    }

    #[test]
    fn eqdc_requires_lat_1() {
        assert_op_err("eqdc");
    }

    #[test]
    fn eqdc_treats_missing_lat_2_as_tangent_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let tangent = ctx.op("eqdc lat_1=30 ellps=GRS80")?;
        let explicit = ctx.op("eqdc lat_1=30 lat_2=30 ellps=GRS80")?;
        let mut lhs = [Coor4D::geo(35.0, 10.0, 0.0, 0.0)];
        let mut rhs = lhs;

        ctx.apply(tangent, Fwd, &mut lhs)?;
        ctx.apply(explicit, Fwd, &mut rhs)?;

        assert!(lhs[0].hypot2(&rhs[0]) < 1e-9);
        Ok(())
    }

    #[test]
    fn eqdc_accepts_near_polar_standard_parallel_like_proj() -> Result<(), Error> {
        assert_proj_match(
            "eqdc ellps=GRS80 lat_1=89.9999999945 lat_2=60",
            Coor4D::geo(10.0, 10.0, 0.0, 0.0),
            Coor4D::raw(1_475_877.204_801_14, 1_229_134.056_910_798_4, 0.0, 0.0),
        )
    }
}
