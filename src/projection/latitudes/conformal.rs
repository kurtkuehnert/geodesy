use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub enum ConformalLatitude {
    Spherical,
    Ellipsoidal {
        ellps: Ellipsoid,
        eccentricity: f64,
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
            eccentricity: ellps.eccentricity(),
            coefficients: ellps.coefficients_for_conformal_latitude_computations(),
        }
    }

    pub fn spherical(self) -> bool {
        matches!(self, Self::Spherical)
    }

    pub fn geographic_to_conformal(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => phi,
            Self::Ellipsoidal { eccentricity, .. } => {
                let (sinphi, cosphi) = phi.sin_cos();
                ((sinphi / cosphi).asinh() - eccentricity * (eccentricity * sinphi).atanh())
                    .sinh()
                    .atan()
            }
        }
    }

    pub fn conformal_to_geographic(self, chi: f64) -> f64 {
        match self {
            Self::Spherical => chi,
            Self::Ellipsoidal { eccentricity, .. } => {
                ancillary::sinhpsi_to_tanphi(chi.tan(), eccentricity).atan()
            }
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

    pub fn series_geographic(self, chi: f64) -> f64 {
        match self {
            Self::Spherical => chi,
            Self::Ellipsoidal { coefficients, .. } => {
                chi + fourier::sin(2.0 * chi, &coefficients.inv)
            }
        }
    }

    pub fn ts_from_latitude(self, phi: f64) -> f64 {
        match self {
            Self::Spherical => (std::f64::consts::FRAC_PI_4 - 0.5 * phi).tan(),
            Self::Ellipsoidal { eccentricity, .. } => ancillary::ts(phi.sin_cos(), eccentricity),
        }
    }

    pub fn latitude_from_ts(self, ts: f64) -> f64 {
        let psi = -ts.ln();
        match self {
            Self::Spherical => gudermannian::fwd(psi),
            Self::Ellipsoidal { eccentricity, .. } => {
                ancillary::sinhpsi_to_tanphi(psi.sinh(), eccentricity).atan()
            }
        }
    }
}
