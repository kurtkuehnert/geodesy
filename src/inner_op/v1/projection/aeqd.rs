//! Azimuthal Equidistant.
//!
//! Attribution:
//! - PROJ 9.8.0 `aeqd.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/aeqd.cpp>
//! - PROJ 9.8.0 `aeqd` documentation:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/docs/source/operations/projections/aeqd.rst>
use crate::authoring::*;
use crate::projection::{ProjectionAspect, ProjectionFrame, SphericalGeodesic};
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
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
        OpParameter::Real { key: "lon_0", default: Some(0_f64) },
        OpParameter::Real { key: "x_0",   default: Some(0_f64) },
        OpParameter::Real { key: "y_0",   default: Some(0_f64) },
    ];

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
    use crate::projection::{assert_forward, assert_forward_and_roundtrip};

    const DEFINITION: &str =
        "aeqd lat_0=52 lon_0=-97.5 x_0=8264722.177 y_0=4867518.353 ellps=WGS84";

    #[test]
    fn aeqd_roundtrip_origin() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            DEFINITION,
            Coor4D::geo(52., -97.5, 0., 0.),
            Coor4D::raw(8_264_722.177, 4_867_518.353, 0., 0.),
            1e-8,
            1e-10,
        )
    }

    #[test]
    fn aeqd_forward_reference() -> Result<(), Error> {
        assert_forward(
            DEFINITION,
            Coor4D::geo(61.390407158824, -101.971128034161, 0., 0.),
            Coor4D::raw(8_024_875.974_4, 5_920_866.114_0, 0., 0.),
            1e-3,
        )
    }

    #[test]
    fn aeqd_rejects_invalid_lat_0() {
        let mut ctx = Minimal::default();
        assert!(matches!(
            ctx.op("aeqd R=1 lat_0=91"),
            Err(Error::BadParam(_, _))
        ));
    }
}
