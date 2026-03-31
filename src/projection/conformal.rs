use crate::authoring::Ellipsoid;
use crate::ellipsoid::EllipsoidBase;
use crate::ellipsoid::latitudes::Latitudes;
use crate::math::{FourierCoefficients, ancillary, fourier};

#[derive(Clone, Copy, Debug)]
pub(crate) enum ConformalLatitude {
    Spherical,
    Ellipsoidal {
        ellps: Ellipsoid,
        coefficients: FourierCoefficients,
    },
}

impl ConformalLatitude {
    pub fn new(ellps: Ellipsoid) -> Self {
        if ellps.flattening() == 0.0 {
            return Self::Spherical;
        }

        Self::Ellipsoidal {
            ellps,
            coefficients: ellps.coefficients_for_conformal_latitude_computations(),
        }
    }

    pub fn spherical(self) -> bool {
        matches!(self, Self::Spherical)
    }

    pub fn eccentricity(self) -> f64 {
        match self {
            Self::Spherical => 0.0,
            Self::Ellipsoidal { ellps, .. } => ellps.eccentricity(),
        }
    }

    pub fn reduced(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => phi,
            Self::Ellipsoidal {
                ellps,
                coefficients,
            } => ellps.latitude_geographic_to_conformal(phi, &coefficients),
        }
    }

    pub fn series_reduced(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => phi,
            Self::Ellipsoidal { coefficients, .. } => {
                phi + fourier::sin(2.0 * phi, &coefficients.fwd)
            }
        }
    }

    pub fn geographic(self, chi: f64) -> f64 {
        match self {
            Self::Spherical => chi,
            Self::Ellipsoidal {
                ellps,
                coefficients,
            } => ellps.latitude_conformal_to_geographic(chi, &coefficients),
        }
    }

    pub fn series_geographic(self, chi: f64) -> f64 {
        match self {
            Self::Spherical => chi,
            Self::Ellipsoidal { coefficients, .. } => {
                chi + fourier::sin(2.0 * chi, &coefficients.inv)
            }
        }
    }

    pub fn ts(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => (std::f64::consts::FRAC_PI_4 - 0.5 * phi).tan(),
            Self::Ellipsoidal { ellps, .. } => {
                ancillary::ts((phi.sin(), phi.cos()), ellps.eccentricity())
            }
        }
    }
}
