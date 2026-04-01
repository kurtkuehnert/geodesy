//! Azimuthal Equidistant.
//!
//! Attribution:
//! - PROJ 9.8.0 `aeqd.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/aeqd.cpp>
//! - PROJ 9.8.0 `aeqd` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/aeqd.rst>
use crate::authoring::*;
use crate::projection::{ProjectionAspect, ProjectionFrame, SphericalGeodesic, projection_gamut};
use std::f64::consts::{FRAC_PI_2, PI};

const ANGULAR_TOLERANCE: f64 = 1e-10;
const LONGITUDE_CANONICALIZATION_TOLERANCE: f64 = 1e-12;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Aeqd {
    pub(crate) frame: ProjectionFrame,
    pub(crate) ellps: Ellipsoid,
    pub(crate) origin: Coor4D,
    pub(crate) aspect: ProjectionAspect,
    pub(crate) spherical: bool,
    pub(crate) sphere: SphericalGeodesic,
    pub(crate) mp: Option<f64>,
}

impl Aeqd {
    pub(crate) fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let frame = ProjectionFrame::from_params(params);
        if frame.lat_0.abs() > FRAC_PI_2 + ANGULAR_TOLERANCE {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = ProjectionAspect::classify(frame.lat_0, ANGULAR_TOLERANCE);

        let ellps = params.ellps(0);
        let spherical = ellps.flattening() == 0.0;
        let sphere = SphericalGeodesic::new(frame.a, frame.lon_0, frame.lat_0);
        let mp = if spherical {
            None
        } else {
            Some(match aspect {
                ProjectionAspect::NorthPolar => ellps.meridian_latitude_to_distance(FRAC_PI_2),
                ProjectionAspect::SouthPolar => ellps.meridian_latitude_to_distance(-FRAC_PI_2),
                ProjectionAspect::Equatorial | ProjectionAspect::Oblique => 0.0,
            })
        };

        Ok(Self {
            frame,
            ellps,
            origin: Coor4D::raw(frame.lon_0, frame.lat_0, 0.0, 0.0),
            aspect,
            spherical,
            sphere,
            mp,
        })
    }

    pub(crate) fn spherical_fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        let lon = self.frame.lon_0 + lam;
        let (azimuth, distance) = self.sphere.geodesic_inv(lon, phi)?;
        Some((distance * azimuth.sin(), distance * azimuth.cos()))
    }

    pub(crate) fn ellipsoidal_fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        match self.aspect {
            ProjectionAspect::NorthPolar | ProjectionAspect::SouthPolar => {
                let coslam = if self.aspect.is_north_polar() {
                    -lam.cos()
                } else {
                    lam.cos()
                };
                let rho =
                    (self.mp.unwrap_or(0.0) - self.ellps.meridian_latitude_to_distance(phi)).abs();
                Some((rho * lam.sin(), rho * coslam))
            }
            ProjectionAspect::Equatorial | ProjectionAspect::Oblique => {
                if lam.abs() < ANGULAR_TOLERANCE
                    && (phi - self.frame.lat_0).abs() < ANGULAR_TOLERANCE
                {
                    return Some((0.0, 0.0));
                }
                if self.aspect.is_equatorial()
                    && phi.abs() < ANGULAR_TOLERANCE
                    && self.frame.lat_0.abs() < ANGULAR_TOLERANCE
                {
                    return Some((self.frame.a * lam, 0.0));
                }

                let target = Coor4D::raw(self.frame.lon_0 + lam, phi, 0.0, 0.0);
                let inv = self.ellps.geodesic_inv(&self.origin, &target);
                let azimuth = inv[0];
                let distance = inv[2];
                (azimuth.is_finite() && distance.is_finite())
                    .then_some((distance * azimuth.sin(), distance * azimuth.cos()))
            }
        }
    }

    pub(crate) fn spherical_inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let distance = x.hypot(y);
        if distance < ANGULAR_TOLERANCE {
            return Some((self.frame.lon_0, self.frame.lat_0));
        }

        let sigma = self.sphere.central_angle(distance)?;

        let (lon, lat) = match self.aspect {
            ProjectionAspect::Equatorial | ProjectionAspect::Oblique => {
                let azimuth = x.atan2(y);
                self.sphere.geodesic_fwd(azimuth, distance)?
            }
            ProjectionAspect::NorthPolar => (self.frame.lon_0 + x.atan2(-y), FRAC_PI_2 - sigma),
            ProjectionAspect::SouthPolar => (self.frame.lon_0 + x.atan2(y), sigma - FRAC_PI_2),
        };

        Some((self.finalize_inverse_lon(lon), lat))
    }

    pub(crate) fn ellipsoidal_inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let distance = x.hypot(y);
        if distance < ANGULAR_TOLERANCE {
            return Some((self.frame.lon_0, self.frame.lat_0));
        }

        let (lon, lat) = match self.aspect {
            ProjectionAspect::Equatorial | ProjectionAspect::Oblique => {
                let azimuth = x.atan2(y);
                let dest = self.ellps.geodesic_fwd(&self.origin, azimuth, distance);
                (dest[0].is_finite() && dest[1].is_finite()).then_some((dest[0], dest[1]))?
            }
            ProjectionAspect::NorthPolar => (
                self.frame.lon_0 + x.atan2(-y),
                self.ellps
                    .meridian_distance_to_latitude(self.mp.unwrap_or(0.0) - distance),
            ),
            ProjectionAspect::SouthPolar => (
                self.frame.lon_0 + x.atan2(y),
                self.ellps
                    .meridian_distance_to_latitude(self.mp.unwrap_or(0.0) + distance),
            ),
        };

        Some((self.finalize_inverse_lon(lon), lat))
    }

    fn finalize_inverse_lon(&self, lon: f64) -> f64 {
        let lon = angular::normalize_symmetric(lon);
        if (lon + PI).abs() < LONGITUDE_CANONICALIZATION_TOLERANCE {
            PI
        } else {
            lon
        }
    }
}

impl PointOp for Aeqd {
    const NAME: &'static str = "aeqd";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = projection_gamut!();

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Self::new(params)
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = self.frame.remove_central_meridian(lon);

        let (x, y) = if self.spherical {
            self.spherical_fwd(lam, lat)?
        } else {
            self.ellipsoidal_fwd(lam, lat)?
        };

        let (x, y) = self.frame.apply_false_origin(x, y);

        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x, y) = self.frame.remove_false_origin(coord[0], coord[1]);

        let (lon, lat) = if self.spherical {
            self.spherical_inv(x, y)?
        } else {
            self.ellipsoidal_inv(x, y)?
        };

        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
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
}
