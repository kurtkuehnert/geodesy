//! Azimuthal Equidistant.
//!
//! Attribution:
//! - PROJ 9.8.0 `aeqd.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/aeqd.cpp>
//! - PROJ 9.8.0 `aeqd` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/aeqd.rst>
use crate::authoring::*;
use crate::projection::{AzimuthalAspect, GeodesicPath, MeridianLatitude};

const ANGULAR_TOLERANCE: f64 = 1e-10;
const LONGITUDE_CANONICALIZATION_TOLERANCE: f64 = 1e-12;

#[derive(Clone, Copy, Debug)]
pub(crate) struct AeqdInner {
    a: f64,
    lat_0: f64,
    aspect: AzimuthalAspect,
    meridian: MeridianLatitude,
    geodesic: GeodesicPath,
    mp: f64,
}

pub(crate) type Aeqd = Framed<AeqdInner>;

impl FramedProjection for AeqdInner {
    const NAME: &'static str = "aeqd";
    const TITLE: &'static str = "Azimuthal Equidistant";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let lat_0 = params.lat(0);
        if lat_0.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = AzimuthalAspect::classify(lat_0, ANGULAR_TOLERANCE);
        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let origin = Coor4D::raw(0.0, lat_0, 0.0, 0.0);
        let meridian = ellps.meridian();
        let geodesic = GeodesicPath::new(ellps, lat_0, origin);
        let mp = match aspect {
            AzimuthalAspect::Polar { pole_sign } => {
                meridian.distance_from_geographic(pole_sign * FRAC_PI_2)
            }
            AzimuthalAspect::Oblique => 0.0,
        };

        Ok(Self {
            a,
            lat_0,
            aspect,
            meridian,
            geodesic,
            mp,
        })
    }

    fn fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        match self.aspect {
            AzimuthalAspect::Oblique => {
                if lam.abs() < ANGULAR_TOLERANCE && (phi - self.lat_0).abs() < ANGULAR_TOLERANCE {
                    return Some((0.0, 0.0));
                }
                if phi.abs() < ANGULAR_TOLERANCE && self.lat_0.abs() < ANGULAR_TOLERANCE {
                    return Some((self.a * lam, 0.0));
                }

                let (distance, azimuth) = self.geodesic.distance_and_azimuth(lam, phi)?;
                Some((distance * azimuth.sin(), distance * azimuth.cos()))
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let coslam = -pole_sign * lam.cos();
                let rho = (self.mp - self.meridian.distance_from_geographic(phi)).abs();
                Some((rho * lam.sin(), rho * coslam))
            }
        }
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let distance = x.hypot(y);
        if distance < ANGULAR_TOLERANCE {
            return Some((0.0, self.lat_0));
        }

        let (lam, lat) = match self.aspect {
            AzimuthalAspect::Oblique => {
                let azimuth = x.atan2(y);
                self.geodesic.destination(azimuth, distance)?
            }
            AzimuthalAspect::Polar { pole_sign } => {
                let lam = x.atan2(-pole_sign * y);
                let lat = self
                    .meridian
                    .geographic_from_distance(self.mp - pole_sign * distance);
                (lam, lat)
            }
        };

        let lam = angular::normalize_symmetric(lam);
        let lam = if (lam + PI).abs() < LONGITUDE_CANONICALIZATION_TOLERANCE {
            PI
        } else {
            lam
        };
        Some((lam, lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aeqd_matches_proj_ellipsoidal_case() -> Result<(), Error> {
        assert_proj_match(
            "aeqd lat_0=52 lon_0=-97.5 x_0=8264722.177 y_0=4867518.353 ellps=WGS84",
            Coor4D::geo(61.390407158824, -101.971128034161, 0., 0.),
            Coor4D::raw(8_024_875.974_381_589, 5_920_866.113_963_526, 0., 0.),
        )
    }

    #[test]
    fn aeqd_matches_proj_spherical_case() -> Result<(), Error> {
        assert_proj_match(
            "aeqd R=6400000 lat_0=30 lon_0=20 x_0=1000 y_0=2000",
            Coor4D::geo(40.0, 25.0, 0.0, 0.0),
            Coor4D::raw(430_838.560_011_896, 1_129_341.938_051_475_4, 0.0, 0.0),
        )
    }

    #[test]
    fn aeqd_matches_proj_ellipsoidal_north_polar_case() -> Result<(), Error> {
        assert_proj_match(
            "aeqd ellps=WGS84 lat_0=90 lon_0=15 x_0=3000 y_0=4000",
            Coor4D::geo(80.0, 20.0, 0.0, 0.0),
            Coor4D::raw(100_337.787_119_382_44, -1_108_575.997_809_590_5, 0.0, 0.0),
        )
    }

    #[test]
    fn aeqd_matches_proj_spherical_north_polar_case() -> Result<(), Error> {
        assert_proj_match(
            "aeqd R=6400000 lat_0=90 lon_0=15 x_0=3000 y_0=4000",
            Coor4D::geo(80.0, 20.0, 0.0, 0.0),
            Coor4D::raw(100_353.899_069_939_5, -1_108_760.158_247_157_5, 0.0, 0.0),
        )
    }

    #[test]
    fn aeqd_matches_proj_spherical_far_side_case() -> Result<(), Error> {
        assert_proj_match(
            "aeqd R=6400000 lat_0=0 lon_0=0",
            Coor4D::geo(0.0, 175.0, 0.0, 0.0),
            Coor4D::raw(19_547_687.622_336_593, 0.0, 0.0, 0.0),
        )
    }

    #[test]
    fn aeqd_matches_proj_equatorial_case() -> Result<(), Error> {
        assert_proj_match(
            "aeqd ellps=WGS84 lat_0=0 lon_0=0 x_0=100 y_0=200",
            Coor4D::geo(0.0, 10.0, 0.0, 0.0),
            Coor4D::raw(1_113_294.907_932_735_7, 200.0, 0.0, 0.0),
        )
    }
}
