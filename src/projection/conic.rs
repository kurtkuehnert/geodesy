use crate::Error;
use crate::op::ParsedParameters;
use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Conic {
    a: f64,
    n: f64,
    rho0: f64,
}

impl Conic {
    pub(crate) const STANDARD_PARALLEL_TOLERANCE: f64 = 1e-10;

    pub fn cone_constant(
        params: &ParsedParameters,
        phi1: f64,
        phi2: f64,
        strict_polar_standard_parallels: bool,
        secant: impl FnOnce() -> f64,
    ) -> Result<f64, Error> {
        let secant_case = (phi1 - phi2).abs() >= Self::STANDARD_PARALLEL_TOLERANCE;
        let polar_parallel = phi1.abs() >= FRAC_PI_2
            || phi2.abs() >= FRAC_PI_2
            || (strict_polar_standard_parallels
                && (phi1.cos().abs() < Self::STANDARD_PARALLEL_TOLERANCE
                    || phi2.cos().abs() < Self::STANDARD_PARALLEL_TOLERANCE));
        if polar_parallel {
            return Err(Error::BadParam(
                "lat_1/lat_2".to_string(),
                params.name.clone(),
            ));
        }

        if secant_case && (phi1 + phi2).abs() < Self::STANDARD_PARALLEL_TOLERANCE {
            return Err(Error::General(Box::leak(
                format!(
                    "{}: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
                    params.name
                )
                .into_boxed_str(),
            )));
        }

        let n = if secant_case { secant() } else { phi1.sin() };

        if n == 0.0 || !n.is_finite() {
            return Err(Error::General(Box::leak(
                format!(
                    "{}: Invalid value for standard parallels: cone constant is not finite",
                    params.name
                )
                .into_boxed_str(),
            )));
        }

        Ok(n)
    }

    pub fn new(a: f64, n: f64, rho0: f64) -> Self {
        Self { a, n, rho0 }
    }

    pub fn n(self) -> f64 {
        self.n
    }

    pub fn polar_point(self) -> (f64, f64) {
        (0.0, self.a * self.rho0)
    }

    pub fn pole(self) -> (f64, f64) {
        (0.0, self.n.signum() * FRAC_PI_2)
    }

    pub fn project(self, lam: f64, rho: f64) -> (f64, f64) {
        let theta = self.n * lam;
        let (sin_theta, cos_theta) = theta.sin_cos();
        let x = self.a * rho * sin_theta;
        let y = self.a * (self.rho0 - rho * cos_theta);
        (x, y)
    }

    pub fn inverse(self, x: f64, y: f64) -> Option<(f64, f64)> {
        let sign = self.n.signum();
        let rho_sin = sign * x / self.a;
        let rho_cos = sign * (self.rho0 - y / self.a);
        let rho = sign * rho_sin.hypot(rho_cos);
        if rho == 0.0 {
            return None;
        }

        let theta = rho_sin.atan2(rho_cos);
        let lam = theta / self.n;
        Some((lam, rho))
    }
}
