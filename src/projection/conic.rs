use std::f64::consts::FRAC_PI_2;

const STANDARD_PARALLEL_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Copy, Debug)]
pub(crate) enum ConicSetupError {
    PolarStandardParallel,
    OppositeStandardParallels,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum StandardParallels {
    Tangent { phi: f64 },
    Secant { phi1: f64, phi2: f64 },
}

impl StandardParallels {
    pub fn classify(phi1: f64, phi2: f64) -> Self {
        if (phi1 - phi2).abs() < STANDARD_PARALLEL_TOLERANCE {
            Self::Tangent { phi: phi1 }
        } else {
            Self::Secant { phi1, phi2 }
        }
    }

    pub fn validate_conic(self) -> Result<Self, ConicSetupError> {
        match self {
            Self::Tangent { phi } => {
                if phi.abs() >= FRAC_PI_2 || phi.cos().abs() < STANDARD_PARALLEL_TOLERANCE {
                    return Err(ConicSetupError::PolarStandardParallel);
                }
                Ok(self)
            }
            Self::Secant { phi1, phi2 } => {
                if phi1.abs() >= FRAC_PI_2
                    || phi2.abs() >= FRAC_PI_2
                    || phi1.cos().abs() < STANDARD_PARALLEL_TOLERANCE
                    || phi2.cos().abs() < STANDARD_PARALLEL_TOLERANCE
                {
                    return Err(ConicSetupError::PolarStandardParallel);
                }
                if (phi1 + phi2).abs() < STANDARD_PARALLEL_TOLERANCE {
                    return Err(ConicSetupError::OppositeStandardParallels);
                }
                Ok(self)
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum ConicInverse {
    Apex,
    Point { lam: f64, rho: f64 },
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Conic {
    a: f64,
    n: f64,
    rho0: f64,
}

impl Conic {
    pub fn new(a: f64, n: f64, rho0: f64) -> Self {
        Self { a, n, rho0 }
    }

    pub fn n(self) -> f64 {
        self.n
    }

    pub fn polar_point(self) -> (f64, f64) {
        (0.0, self.a * self.rho0)
    }

    pub fn project(self, lam: f64, rho: f64) -> (f64, f64) {
        let theta = self.n * lam;
        let (sin_theta, cos_theta) = theta.sin_cos();
        let x = self.a * rho * sin_theta;
        let y = self.a * (self.rho0 - rho * cos_theta);
        (x, y)
    }

    pub fn inverse(self, x: f64, y: f64) -> ConicInverse {
        let sign = self.n.signum();
        let rho_sin = sign * x / self.a;
        let rho_cos = sign * (self.rho0 - y / self.a);
        let rho = sign * rho_sin.hypot(rho_cos);
        if rho == 0.0 {
            return ConicInverse::Apex;
        }

        let theta = rho_sin.atan2(rho_cos);
        let lam = theta / self.n;
        ConicInverse::Point { lam, rho }
    }
}
