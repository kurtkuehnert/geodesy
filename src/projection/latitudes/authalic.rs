use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub enum AuthalicLatitude {
    Spherical,
    Ellipsoidal {
        e: f64,
        q_pole: f64,
        ellps: Ellipsoid,
        coefficients: FourierCoefficients,
    },
}

impl AuthalicLatitude {
    const DOMAIN_TOLERANCE: f64 = 1e-10;

    pub fn new(ellps: Ellipsoid) -> Self {
        if ellps.flattening() == 0.0 {
            return Self::Spherical;
        }

        let e = ellps.eccentricity();
        let q_pole = 1.0 + (1.0 - e * e) * e.atanh() / e;

        Self::Ellipsoidal {
            e,
            q_pole,
            ellps,
            coefficients: ellps.coefficients_for_authalic_latitude_computations(),
        }
    }

    pub fn q_pole(self) -> f64 {
        match self {
            Self::Spherical => 2.0,
            Self::Ellipsoidal { q_pole, .. } => q_pole,
        }
    }

    pub fn geographic_to_authalic(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => phi,
            Self::Ellipsoidal { .. } => self.beta_from_q(self.q_from_geographic(phi)),
        }
    }

    pub fn authalic_to_geographic(self, beta: f64) -> f64 {
        match self {
            Self::Spherical => beta,
            Self::Ellipsoidal {
                e,
                q_pole,
                ellps,
                coefficients,
            } => {
                let mut geographic = beta + fourier::sin(2.0 * beta, &coefficients.inv);
                if ellps.third_flattening().abs() < 0.01 {
                    return geographic;
                }

                let es = e * e;
                let one_es = 1.0 - es;
                let q = beta.sin() * q_pole;
                let q_div_one_minus_es = q / one_es;

                for _ in 0..10 {
                    let sinphi = geographic.sin();
                    let cosphi = geographic.cos();
                    let one_minus_es_sin2phi = 1.0 - es * sinphi * sinphi;
                    let dphi = (one_minus_es_sin2phi * one_minus_es_sin2phi) / (2.0 * cosphi)
                        * (q_div_one_minus_es
                            - sinphi / one_minus_es_sin2phi
                            - (e * sinphi).atanh() / e);
                    if dphi.abs() < 1e-15 {
                        break;
                    }
                    geographic += dphi;
                }

                geographic
            }
        }
    }

    pub fn q_from_geographic(self, phi: f64) -> f64 {
        let sinphi = phi.sin();
        match self {
            Self::Spherical => 2.0 * sinphi,
            Self::Ellipsoidal { e, .. } => ancillary::qs(sinphi, e),
        }
    }

    pub fn geographic_from_q(self, q: f64) -> Option<f64> {
        let q_pole = self.q_pole();
        if q.abs() > q_pole + Self::DOMAIN_TOLERANCE {
            return None;
        }

        let q = q.clamp(-q_pole, q_pole);
        let beta = self.beta_from_q(q);

        Some(self.authalic_to_geographic(beta))
    }

    fn beta_from_q(self, q: f64) -> f64 {
        (q / self.q_pole()).clamp(-1.0, 1.0).asin()
    }
}
