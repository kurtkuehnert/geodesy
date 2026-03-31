use crate::authoring::Ellipsoid;
use crate::ellipsoid::EllipsoidBase;
use crate::ellipsoid::latitudes::Latitudes;
use crate::math::{FourierCoefficients, ancillary};
use std::f64::consts::FRAC_PI_2;

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
    const DOMAIN_TOLERANCE: f64 = 1e-10;
    const SATURATION_TOLERANCE: f64 = 1e-7;

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

    pub fn beta_from_phi(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => phi,
            Self::Ellipsoidal { .. } => (self.q_from_phi(phi) / self.qp()).clamp(-1.0, 1.0).asin(),
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

    pub fn phi_from_q(self, q: f64) -> Option<f64> {
        let normalized_q = q / self.qp();
        if normalized_q.abs() > 1.0 + Self::DOMAIN_TOLERANCE {
            return None;
        }
        Some(self.phi_from_beta(normalized_q.clamp(-1.0, 1.0).asin()))
    }

    pub fn q_from_phi(self, phi: f64) -> f64 {
        let sinphi = phi.sin();
        match self {
            Self::Spherical => 2.0 * sinphi,
            Self::Ellipsoidal { ellps, .. } => ancillary::qs(sinphi, ellps.eccentricity()),
        }
    }

    pub fn phi_from_q_saturating(self, q: f64) -> Option<f64> {
        if (self.qp() - q.abs()).abs() <= Self::SATURATION_TOLERANCE {
            return Some(if q < 0.0 { -FRAC_PI_2 } else { FRAC_PI_2 });
        }
        self.phi_from_q(q)
    }
}
