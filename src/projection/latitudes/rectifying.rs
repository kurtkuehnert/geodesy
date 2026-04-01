use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub enum RectifyingLatitude {
    Spherical {
        a: f64,
    },
    Ellipsoidal {
        a: f64,
        coefficients: FourierCoefficients,
    },
}

impl RectifyingLatitude {
    pub(crate) fn new(ellps: Ellipsoid) -> Self {
        if ellps.flattening() == 0.0 {
            Self::Spherical {
                a: ellps.semimajor_axis(),
            }
        } else {
            Self::Ellipsoidal {
                a: ellps.semimajor_axis(),
                coefficients: ellps.coefficients_for_rectifying_latitude_computations(),
            }
        }
    }

    pub fn semimajor_axis(self) -> f64 {
        match self {
            Self::Spherical { a } | Self::Ellipsoidal { a, .. } => a,
        }
    }

    pub fn distance_from_latitude(self, lat: f64) -> f64 {
        match self {
            Self::Spherical { a } => a * lat,
            Self::Ellipsoidal { a, coefficients } => {
                a * coefficients.etc[0] * (lat + fourier::sin(2.0 * lat, &coefficients.fwd))
            }
        }
    }

    pub fn latitude_from_distance(self, distance: f64) -> f64 {
        match self {
            Self::Spherical { a } => distance / a,
            Self::Ellipsoidal { a, coefficients } => {
                let rectifying = distance / (a * coefficients.etc[0]);
                rectifying + fourier::sin(2.0 * rectifying, &coefficients.inv)
            }
        }
    }
}
