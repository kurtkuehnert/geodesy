use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum AzimuthalAspect {
    Oblique,
    Polar { pole_sign: f64 },
}

impl AzimuthalAspect {
    pub fn classify(lat_0: f64, tolerance: f64) -> Self {
        let abs_lat_0 = lat_0.abs();
        if (abs_lat_0 - FRAC_PI_2).abs() < tolerance {
            Self::Polar {
                pole_sign: if lat_0 < 0.0 { -1.0 } else { 1.0 },
            }
        } else {
            Self::Oblique
        }
    }

    pub fn pole_sign(self) -> Option<f64> {
        match self {
            Self::Polar { pole_sign } => Some(pole_sign),
            Self::Oblique => None,
        }
    }

    pub fn polar_xy(self, lam: f64, rho: f64) -> (f64, f64) {
        let pole_sign = self
            .pole_sign()
            .expect("polar azimuthal helper requires polar aspect");
        let (sin_lam, cos_lam) = lam.sin_cos();
        (rho * sin_lam, -pole_sign * rho * cos_lam)
    }

    pub fn polar_lam(self, x: f64, y: f64) -> f64 {
        let pole_sign = self
            .pole_sign()
            .expect("polar azimuthal helper requires polar aspect");
        x.atan2(-pole_sign * y)
    }
}
