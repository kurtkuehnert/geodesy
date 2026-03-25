//! Oblique Stereographic Alternative.
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const DEL_TOL: f64 = 1e-14;
const MAX_ITER: usize = 20;

fn w(sinphi: f64, e: f64, n: f64, c: f64) -> f64 {
    let s1 = (1.0 + sinphi) / (1.0 - sinphi);
    let s2 = (1.0 - e * sinphi) / (1.0 + e * sinphi);
    c * (s1 * s2.powf(e)).powf(n)
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(n) = op.params.real("n") else { return 0 };
    let Ok(c) = op.params.real("c") else { return 0 };
    let k0 = op.params.k(0);
    let Ok(r2) = op.params.real("r2") else {
        return 0;
    };
    let Ok(sin_chi0) = op.params.real("sin_chi0") else {
        return 0;
    };
    let Ok(cos_chi0) = op.params.real("cos_chi0") else {
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
        let sinphi = lat.sin();
        let ww = w(sinphi, e, n, c);
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

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(n) = op.params.real("n") else { return 0 };
    let Ok(c) = op.params.real("c") else { return 0 };
    let k0 = op.params.k(0);
    let Ok(r2) = op.params.real("r2") else {
        return 0;
    };
    let Ok(chi0) = op.params.real("chi0") else {
        return 0;
    };
    let Ok(g) = op.params.real("g") else {
        return 0;
    };
    let Ok(h) = op.params.real("h") else {
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
    let es = ellps.eccentricity_squared();
    let sphi0 = lat_0.sin();
    let cphi0 = lat_0.cos();
    let a = ellps.semimajor_axis();

    let rho0 = a * (1.0 - es) / (1.0 - es * sphi0 * sphi0).powf(1.5);
    let nu0 = a / (1.0 - es * sphi0 * sphi0).sqrt();
    let r = (rho0 * nu0).sqrt();
    let n = (1.0 + es * cphi0.powi(4) / (1.0 - es)).sqrt();

    let s1 = (1.0 + sphi0) / (1.0 - sphi0);
    let s2 = (1.0 - e * sphi0) / (1.0 + e * sphi0);
    let w1 = (s1 * s2.powf(e)).powf(n);
    let sin_chi00 = (w1 - 1.0) / (w1 + 1.0);
    let c = ((n + sphi0) * (1.0 - sin_chi00)) / ((n - sphi0) * (1.0 + sin_chi00));
    let w2 = c * w1;
    let chi0 = ((w2 - 1.0) / (w2 + 1.0)).asin();
    let sin_chi0 = chi0.sin();
    let cos_chi0 = chi0.cos();
    let r2 = 2.0 * r / a;
    let k0 = params.k(0);
    let g = a * k0 * r2 * (FRAC_PI_4 - 0.5 * chi0).tan();
    let h = 2.0 * a * k0 * r2 * chi0.tan() + g;

    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", lon_0);
    params.real.insert("n", n);
    params.real.insert("c", c);
    params.real.insert("chi0", chi0);
    params.real.insert("sin_chi0", sin_chi0);
    params.real.insert("cos_chi0", cos_chi0);
    params.real.insert("r2", r2);
    params.real.insert("g", g);
    params.real.insert("h", h);

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
