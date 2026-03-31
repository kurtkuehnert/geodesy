use std::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ProjectionAspect {
    Equatorial,
    Oblique,
    NorthPolar,
    SouthPolar,
}

impl ProjectionAspect {
    pub fn classify(lat_0: f64, tolerance: f64) -> Self {
        let abs_lat_0 = lat_0.abs();
        if (abs_lat_0 - FRAC_PI_2).abs() < tolerance {
            if lat_0 < 0.0 {
                Self::SouthPolar
            } else {
                Self::NorthPolar
            }
        } else if abs_lat_0 < tolerance {
            Self::Equatorial
        } else {
            Self::Oblique
        }
    }

    pub fn is_equatorial(self) -> bool {
        matches!(self, Self::Equatorial)
    }

    pub fn is_oblique(self) -> bool {
        matches!(self, Self::Oblique)
    }

    pub fn is_north_polar(self) -> bool {
        matches!(self, Self::NorthPolar)
    }

    pub fn is_south_polar(self) -> bool {
        matches!(self, Self::SouthPolar)
    }

    pub fn is_polar(self) -> bool {
        matches!(self, Self::NorthPolar | Self::SouthPolar)
    }
}

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
}
