//! Oblique Stereographic Alternative, closely following PROJ's implementation.
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const DEL_TOL: f64 = 1e-14;
const MAX_ITER: usize = 20;

fn srat(esinp: f64, ratexp: f64) -> f64 {
    ((1.0 - esinp) / (1.0 + esinp)).powf(ratexp)
}

fn legacy_w(sinphi: f64, e: f64, n: f64, c: f64) -> f64 {
    let s1 = (1.0 + sinphi) / (1.0 - sinphi);
    let s2 = (1.0 - e * sinphi) / (1.0 + e * sinphi);
    c * (s1 * s2.powf(e)).powf(n)
}

fn gauss_ini(e: f64, phi0: f64) -> Option<(f64, f64, f64, f64, f64)> {
    let es = e * e;
    let sphi = phi0.sin();
    let cphi2 = phi0.cos().powi(2);
    let rc = (1.0 - es).sqrt() / (1.0 - es * sphi * sphi);
    let c = (1.0 + es * cphi2 * cphi2 / (1.0 - es)).sqrt();
    if c == 0.0 {
        return None;
    }
    let chi = (sphi / c).asin();
    let ratexp = 0.5 * c * e;
    let srat_val = srat(e * sphi, ratexp);
    if srat_val == 0.0 {
        return None;
    }
    let k = if 0.5 * phi0 + FRAC_PI_4 < 1e-10 {
        1.0 / srat_val
    } else {
        (0.5 * chi + FRAC_PI_4).tan() / ((0.5 * phi0 + FRAC_PI_4).tan().powf(c) * srat_val)
    };
    Some((c, k, ratexp, chi, rc))
}

fn gauss(lon: f64, lat: f64, c: f64, k: f64, e: f64, ratexp: f64) -> (f64, f64) {
    if (lat.abs() - FRAC_PI_2).abs() < 1e-14 {
        return (c * lon, lat.signum() * FRAC_PI_2);
    }
    let phi = 2.0
        * (k * (0.5 * lat + FRAC_PI_4).tan().powf(c) * srat(e * lat.sin(), ratexp)).atan()
        - FRAC_PI_2;
    let lam = c * lon;
    (lam, phi)
}

fn inv_gauss(lon: f64, lat: f64, c: f64, k: f64, e: f64) -> Option<(f64, f64)> {
    let lam = lon / c;
    let num = ((0.5 * lat + FRAC_PI_4).tan() / k).powf(1.0 / c);
    let mut slp_phi = lat;
    for _ in 0..MAX_ITER {
        let phi = 2.0 * (num * srat(e * slp_phi.sin(), -0.5 * e)).atan() - FRAC_PI_2;
        if (phi - slp_phi).abs() < DEL_TOL {
            return Some((lam, phi));
        }
        slp_phi = phi;
    }
    None
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    if !op.params.boolean("proj_polar") {
        return legacy_fwd(op, operands);
    }

    let Ok(gauss_c) = op.params.real("gauss_c") else {
        return 0;
    };
    let Ok(gauss_k) = op.params.real("gauss_k") else {
        return 0;
    };
    let Ok(ratexp) = op.params.real("gauss_ratexp") else {
        return 0;
    };
    let k0 = op.params.k(0);
    let Ok(r2) = op.params.real("r2") else {
        return 0;
    };
    let Ok(sinc0) = op.params.real("sinc0") else {
        return 0;
    };
    let Ok(cosc0) = op.params.real("cosc0") else {
        return 0;
    };
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let (lam, phi) = gauss(lon - lon_0, lat, gauss_c, gauss_k, e, ratexp);
        let sinc = phi.sin();
        let cosc = phi.cos();
        let cosl = lam.cos();
        let denom = 1.0 + sinc0 * sinc + cosc0 * cosc * cosl;
        if denom == 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let scale = a * k0 * r2 / denom;
        let x = x_0 + scale * cosc * lam.sin();
        let y = y_0 + scale * (cosc0 * sinc - sinc0 * cosc * cosl);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    if !op.params.boolean("proj_polar") {
        return legacy_inv(op, operands);
    }

    let Ok(gauss_c) = op.params.real("gauss_c") else {
        return 0;
    };
    let Ok(gauss_k) = op.params.real("gauss_k") else {
        return 0;
    };
    let k0 = op.params.k(0);
    let Ok(r2) = op.params.real("r2") else {
        return 0;
    };
    let Ok(phic0) = op.params.real("phic0") else {
        return 0;
    };
    let Ok(sinc0) = op.params.real("sinc0") else {
        return 0;
    };
    let Ok(cosc0) = op.params.real("cosc0") else {
        return 0;
    };
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (mut x, mut y) = operands.xy(i);
        x = (x - x_0) / (a * k0);
        y = (y - y_0) / (a * k0);
        let rho = x.hypot(y);

        let (lam, phi) = if rho != 0.0 {
            let c = 2.0 * rho.atan2(r2);
            let sinc = c.sin();
            let cosc = c.cos();
            let phi = (cosc * sinc0 + y * sinc * cosc0 / rho).asin();
            let lam = (x * sinc).atan2(rho * cosc0 * cosc - y * sinc0 * sinc);
            (lam, phi)
        } else {
            (0.0, phic0)
        };

        let Some((lam, phi)) = inv_gauss(lam, phi, gauss_c, gauss_k, e) else {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        };

        operands.set_xy(i, lon_0 + lam, phi);
        successes += 1;
    }
    successes
}

fn legacy_fwd(op: &Op, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(n) = op.params.real("legacy_n") else {
        return 0;
    };
    let Ok(c) = op.params.real("legacy_c") else {
        return 0;
    };
    let k0 = op.params.k(0);
    let Ok(r2) = op.params.real("legacy_r2") else {
        return 0;
    };
    let Ok(sin_chi0) = op.params.real("legacy_sin_chi0") else {
        return 0;
    };
    let Ok(cos_chi0) = op.params.real("legacy_cos_chi0") else {
        return 0;
    };
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    let sinc0 = op.params.real("sinc0").unwrap_or(0.0);
    let cosc0 = op.params.real("cosc0").unwrap_or(0.0);
    let proj_r2 = op.params.real("r2").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        if (lat.abs() - FRAC_PI_2).abs() < 1e-14 {
            let sign = lat.signum();
            let denom = 1.0 + sinc0 * sign;
            if denom == 0.0 {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let scale = a * k0 * proj_r2 / denom;
            let x = x_0;
            let y = y_0 + scale * cosc0 * sign;
            operands.set_xy(i, x, y);
            successes += 1;
            continue;
        }
        let sinphi = lat.sin();
        let ww = legacy_w(sinphi, e, n, c);
        let sin_chi = (ww - 1.0) / (ww + 1.0);
        let chi = sin_chi.asin();
        let cos_chi = chi.cos();
        let d_lambda = n * (lon - lon_0);
        let (sin_dl, cos_dl) = d_lambda.sin_cos();
        let denom = 1.0 + sin_chi * sin_chi0 + cos_chi * cos_chi0 * cos_dl;
        if denom == 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let scale = a * k0 * r2 / denom;
        let x = x_0 + scale * cos_chi * sin_dl;
        let y = y_0 + scale * (sin_chi * cos_chi0 - cos_chi * sin_chi0 * cos_dl);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn legacy_inv(op: &Op, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(n) = op.params.real("legacy_n") else {
        return 0;
    };
    let Ok(c) = op.params.real("legacy_c") else {
        return 0;
    };
    let k0 = op.params.k(0);
    let Ok(r2) = op.params.real("legacy_r2") else {
        return 0;
    };
    let Ok(chi0) = op.params.real("legacy_chi0") else {
        return 0;
    };
    let Ok(g) = op.params.real("legacy_g") else {
        return 0;
    };
    let Ok(h) = op.params.real("legacy_h") else {
        return 0;
    };
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let easting = operands.xy(i).0 - x_0;
        let northing = operands.xy(i).1 - y_0;
        let iang = easting.atan2(h + northing);
        let jang = easting.atan2(g - northing) - iang;
        let chi = chi0 + 2.0 * ((northing - easting * (0.5 * jang).tan()) / (a * k0 * r2)).atan();
        let lon = lon_0 + (jang + 2.0 * iang) / n;

        let sin_chi = chi.sin();
        let psi = 0.5 * ((1.0 + sin_chi) / (c * (1.0 - sin_chi))).ln() / n;
        let mut phi = 2.0 * psi.exp().atan() - FRAC_PI_2;
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let sinphi = phi.sin();
            let psi_i = (0.5 * phi + FRAC_PI_4).tan().ln()
                + (((1.0 - e * sinphi) / (1.0 + e * sinphi)).ln() * e * 0.5);
            let next =
                phi - (psi_i - psi) * phi.cos() * (1.0 - e * e * sinphi * sinphi) / (1.0 - e * e);
            if (next - phi).abs() < DEL_TOL {
                phi = next;
                converged = true;
                break;
            }
            phi = next;
        }
        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        operands.set_xy(i, lon, phi);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let lat_0 = params.lat(0).to_radians();
    let lon_0 = params.lon(0).to_radians();
    let ellps = params.ellps(0);
    let e = ellps.eccentricity();
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let sphi0 = lat_0.sin();
    let cphi0 = lat_0.cos();

    let Some((gauss_c, gauss_k, gauss_ratexp, phic0, rc)) = gauss_ini(e, lat_0) else {
        return Err(Error::Unsupported("sterea gauss setup failed".into()));
    };

    let rho0 = a * (1.0 - es) / (1.0 - es * sphi0 * sphi0).powf(1.5);
    let nu0 = a / (1.0 - es * sphi0 * sphi0).sqrt();
    let r = (rho0 * nu0).sqrt();
    let legacy_n = (1.0 + es * cphi0.powi(4) / (1.0 - es)).sqrt();
    let s1 = (1.0 + sphi0) / (1.0 - sphi0);
    let s2 = (1.0 - e * sphi0) / (1.0 + e * sphi0);
    let w1 = (s1 * s2.powf(e)).powf(legacy_n);
    let sin_chi00 = (w1 - 1.0) / (w1 + 1.0);
    let legacy_c =
        ((legacy_n + sphi0) * (1.0 - sin_chi00)) / ((legacy_n - sphi0) * (1.0 + sin_chi00));
    let w2 = legacy_c * w1;
    let legacy_chi0 = ((w2 - 1.0) / (w2 + 1.0)).asin();
    let legacy_r2 = 2.0 * r / a;
    let k0 = params.k(0);
    let legacy_g = a * k0 * legacy_r2 * (FRAC_PI_4 - 0.5 * legacy_chi0).tan();
    let legacy_h = 2.0 * a * k0 * legacy_r2 * legacy_chi0.tan() + legacy_g;

    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", lon_0);
    if (lat_0.abs() - FRAC_PI_2).abs() < 1e-12 {
        params.boolean.insert("proj_polar");
    }
    params.real.insert("gauss_c", gauss_c);
    params.real.insert("gauss_k", gauss_k);
    params.real.insert("gauss_ratexp", gauss_ratexp);
    params.real.insert("phic0", phic0);
    params.real.insert("sinc0", phic0.sin());
    params.real.insert("cosc0", phic0.cos());
    params.real.insert("r2", 2.0 * rc);
    params.real.insert("legacy_n", legacy_n);
    params.real.insert("legacy_c", legacy_c);
    params.real.insert("legacy_chi0", legacy_chi0);
    params.real.insert("legacy_sin_chi0", legacy_chi0.sin());
    params.real.insert("legacy_cos_chi0", legacy_chi0.cos());
    params.real.insert("legacy_r2", legacy_r2);
    params.real.insert("legacy_g", legacy_g);
    params.real.insert("legacy_h", legacy_h);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sterea_roundtrip_origin() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("sterea lat_0=52.1666666666667 lon_0=19.1666666666667 k_0=0.999714 x_0=500000 y_0=500000 ellps=krass")?;

        let geo = [Coor4D::geo(52.1666666666667, 19.1666666666667, 0., 0.)];
        let projected = [Coor4D::raw(500000.0, 500000.0, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-8);
        Ok(())
    }

    #[test]
    fn sterea_matches_proj_sample() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("sterea lat_0=52.1666666666667 lon_0=19.1666666666667 k_0=0.999714 x_0=500000 y_0=500000 ellps=krass")?;

        let geo = [Coor4D::geo(49.734897261834, 19.151997685019, 0., 0.)];
        let projected = [Coor4D::raw(498942.3343743173, 229504.5906664739, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-11);
        Ok(())
    }
}
