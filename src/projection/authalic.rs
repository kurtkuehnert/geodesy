use crate::authoring::Ellipsoid;
use crate::ellipsoid::EllipsoidBase;
use crate::ellipsoid::latitudes::Latitudes;
use crate::math::{FourierCoefficients, ancillary};

#[derive(Clone, Copy, Debug)]
pub(crate) enum AuthalicLatitude {
    Spherical,
    Ellipsoidal {
        ellps: Ellipsoid,
        coefficients: FourierCoefficients,
        qp: f64,
    },
}

impl AuthalicLatitude {
    pub fn new(ellps: Ellipsoid) -> Self {
        if ellps.flattening() == 0.0 {
            return Self::Spherical;
        }

        Self::Ellipsoidal {
            qp: ancillary::qs(1.0, ellps.eccentricity()),
            coefficients: ellps.coefficients_for_authalic_latitude_computations(),
            ellps,
        }
    }

    pub fn spherical(self) -> bool {
        matches!(self, Self::Spherical)
    }

    pub fn qp(self) -> f64 {
        match self {
            Self::Spherical => 2.0,
            Self::Ellipsoidal { qp, .. } => qp,
        }
    }

    pub fn qs(self, sinphi: f64) -> f64 {
        match self {
            Self::Spherical => 2.0 * sinphi,
            Self::Ellipsoidal { ellps, .. } => ancillary::qs(sinphi, ellps.eccentricity()),
        }
    }

    pub fn beta_from_phi(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => phi,
            Self::Ellipsoidal { .. } => (self.qs(phi.sin()) / self.qp()).clamp(-1.0, 1.0).asin(),
        }
    }

    pub fn phi_from_beta(self, beta: f64) -> f64 {
        match self {
            Self::Spherical => beta,
            Self::Ellipsoidal {
                ellps,
                coefficients,
                ..
            } => ellps.latitude_authalic_to_geographic(beta, &coefficients),
        }
    }
}
