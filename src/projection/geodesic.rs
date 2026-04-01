use crate::authoring::*;
use std::f64::consts::PI;

const ANGULAR_TOLERANCE: f64 = 1e-10;
const SINGULARITY_TOLERANCE: f64 = 1e-14;

#[derive(Clone, Copy, Debug)]
pub(crate) enum GeodesicPath {
    Spherical {
        radius: f64,
        lat_0: f64,
        sinph0: f64,
        cosph0: f64,
    },
    Ellipsoidal {
        ellps: Ellipsoid,
        origin: Coor4D,
    },
}

impl GeodesicPath {
    pub(crate) fn new(ellps: Ellipsoid, lat_0: f64, origin: Coor4D) -> Self {
        if ellps.flattening() == 0.0 {
            let (sinph0, cosph0) = lat_0.sin_cos();
            Self::Spherical {
                radius: ellps.semimajor_axis(),
                lat_0,
                sinph0,
                cosph0,
            }
        } else {
            Self::Ellipsoidal { ellps, origin }
        }
    }

    pub(crate) fn destination(self, azimuth: f64, distance: f64) -> Option<(f64, f64)> {
        match self {
            Self::Spherical {
                radius,
                lat_0,
                sinph0,
                cosph0,
            } => {
                let sigma = distance / radius;
                if sigma > PI && sigma - ANGULAR_TOLERANCE > PI {
                    return None;
                }
                let sigma = sigma.min(PI);
                if sigma < ANGULAR_TOLERANCE {
                    return Some((0.0, lat_0));
                }

                let (sinsig, cossig) = sigma.sin_cos();
                let lat = (sinph0 * cossig + cosph0 * sinsig * azimuth.cos()).asin();
                let lam = (azimuth.sin() * sinsig * cosph0).atan2(cossig - sinph0 * lat.sin());
                Some((lam, lat))
            }
            Self::Ellipsoidal { ellps, origin } => {
                let dest = ellps.geodesic_fwd(&origin, azimuth, distance);
                (dest[0].is_finite() && dest[1].is_finite()).then_some((dest[0], dest[1]))
            }
        }
    }

    pub(crate) fn distance_and_azimuth(self, lam: f64, phi: f64) -> Option<(f64, f64)> {
        match self {
            Self::Spherical {
                radius,
                sinph0,
                cosph0,
                ..
            } => {
                let (sinphi, cosphi) = phi.sin_cos();
                let (sinlam, coslam) = lam.sin_cos();
                let cosc = sinph0 * sinphi + cosph0 * cosphi * coslam;

                if (cosc.abs() - 1.0).abs() < SINGULARITY_TOLERANCE {
                    return (cosc >= 0.0).then_some((0.0, 0.0));
                }

                let azimuth = (cosphi * sinlam).atan2(cosph0 * sinphi - sinph0 * cosphi * coslam);
                Some((radius * cosc.acos(), azimuth))
            }
            Self::Ellipsoidal { ellps, origin } => {
                let target = Coor4D::raw(lam, phi, 0.0, 0.0);
                let inv = ellps.geodesic_inv(&origin, &target);
                let azimuth = inv[0];
                let distance = inv[2];
                (azimuth.is_finite() && distance.is_finite()).then_some((distance, azimuth))
            }
        }
    }
}
