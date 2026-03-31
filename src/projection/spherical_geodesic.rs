use crate::authoring::*;
use std::f64::consts::PI;

const ANGULAR_TOLERANCE: f64 = 1e-10;
const SINGULARITY_TOLERANCE: f64 = 1e-14;

#[derive(Clone, Copy, Debug)]
pub(crate) struct SphericalGeodesic {
    radius: f64,
    lon_0: f64,
    lat_0: f64,
    sinph0: f64,
    cosph0: f64,
}

impl SphericalGeodesic {
    pub fn new(radius: f64, lon_0: f64, lat_0: f64) -> Self {
        let (sinph0, cosph0) = lat_0.sin_cos();
        Self {
            radius,
            lon_0,
            lat_0,
            sinph0,
            cosph0,
        }
    }

    pub fn geodesic_inv(&self, lon: f64, lat: f64) -> Option<(f64, f64)> {
        let lam = angular::normalize_symmetric(lon - self.lon_0);
        let (sinphi, cosphi) = lat.sin_cos();
        let (sinlam, coslam) = lam.sin_cos();
        let cosc = self.sinph0 * sinphi + self.cosph0 * cosphi * coslam;

        if (cosc.abs() - 1.0).abs() < SINGULARITY_TOLERANCE {
            return (cosc >= 0.0).then_some((0.0, 0.0));
        }

        let azimuth = (cosphi * sinlam).atan2(self.cosph0 * sinphi - self.sinph0 * cosphi * coslam);
        Some((azimuth, self.radius * cosc.acos()))
    }

    pub fn geodesic_fwd(&self, azimuth: f64, distance: f64) -> Option<(f64, f64)> {
        let sigma = self.central_angle(distance)?;

        if sigma < ANGULAR_TOLERANCE {
            return Some((self.lon_0, self.lat_0));
        }

        let (sinsig, cossig) = sigma.sin_cos();
        let lat = (self.sinph0 * cossig + self.cosph0 * sinsig * azimuth.cos()).asin();
        let lon = self.lon_0
            + (azimuth.sin() * sinsig * self.cosph0).atan2(cossig - self.sinph0 * lat.sin());
        Some((lon, lat))
    }

    pub fn central_angle(&self, distance: f64) -> Option<f64> {
        let sigma = distance / self.radius;
        if sigma > PI && sigma - ANGULAR_TOLERANCE > PI {
            return None;
        }
        Some(sigma.min(PI))
    }
}
