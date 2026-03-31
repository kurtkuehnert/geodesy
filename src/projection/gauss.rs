//! Gauss conformal sphere reduction.
//!
//! Attribution:
//! - PROJ 9.8.0 `gauss.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/gauss.cpp>
//! - PROJ 9.8.0 `sterea.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/projections/sterea.cpp>

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const CONVERGENCE_TOLERANCE: f64 = 1e-14;
const K_POLE_TOLERANCE: f64 = 1e-10;
const POLAR_SOURCE_TOLERANCE: f64 = 1e-14;
const MAX_ITER: usize = 20;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Gauss {
    eccentricity: f64,
    c: f64,
    k: f64,
    ratexp: f64,
    pub phic0: f64,
    pub sinc0: f64,
    pub cosc0: f64,
    pub r2: f64,
}

fn srat(esinp: f64, ratexp: f64) -> f64 {
    ((1.0 - esinp) / (1.0 + esinp)).powf(ratexp)
}

impl Gauss {
    pub fn new(eccentricity: f64, phi0: f64) -> Option<Self> {
        let es = eccentricity * eccentricity;
        let sphi = phi0.sin();
        let cphi2 = phi0.cos().powi(2);
        let rc = (1.0 - es).sqrt() / (1.0 - es * sphi * sphi);
        let c = (1.0 + es * cphi2 * cphi2 / (1.0 - es)).sqrt();
        if c == 0.0 {
            return None;
        }

        let phic0 = (sphi / c).asin();
        let ratexp = 0.5 * c * eccentricity;
        let srat_val = srat(eccentricity * sphi, ratexp);
        if srat_val == 0.0 {
            return None;
        }

        let k = if 0.5 * phi0 + FRAC_PI_4 < K_POLE_TOLERANCE {
            1.0 / srat_val
        } else {
            (0.5 * phic0 + FRAC_PI_4).tan() / ((0.5 * phi0 + FRAC_PI_4).tan().powf(c) * srat_val)
        };

        Some(Self {
            eccentricity,
            c,
            k,
            ratexp,
            phic0,
            sinc0: phic0.sin(),
            cosc0: phic0.cos(),
            r2: 2.0 * rc,
        })
    }

    pub fn forward(self, lon: f64, lat: f64) -> (f64, f64) {
        if (lat.abs() - FRAC_PI_2).abs() < POLAR_SOURCE_TOLERANCE {
            return (self.c * lon, lat.signum() * FRAC_PI_2);
        }

        let phi = 2.0
            * (self.k
                * (0.5 * lat + FRAC_PI_4).tan().powf(self.c)
                * srat(self.eccentricity * lat.sin(), self.ratexp))
            .atan()
            - FRAC_PI_2;
        (self.c * lon, phi)
    }

    pub fn inverse(self, lon: f64, lat: f64) -> Option<(f64, f64)> {
        let lam = lon / self.c;
        let num = ((0.5 * lat + FRAC_PI_4).tan() / self.k).powf(1.0 / self.c);
        let mut slp_phi = lat;
        for _ in 0..MAX_ITER {
            let phi = 2.0
                * (num * srat(self.eccentricity * slp_phi.sin(), -0.5 * self.eccentricity)).atan()
                - FRAC_PI_2;
            if (phi - slp_phi).abs() < CONVERGENCE_TOLERANCE {
                return Some((lam, phi));
            }
            slp_phi = phi;
        }
        None
    }
}
