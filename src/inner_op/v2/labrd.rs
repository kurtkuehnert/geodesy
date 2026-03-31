//! Laborde
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const EPS: f64 = 1.0e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let e = ellps.eccentricity();
    let k_rg = op.params.real["k_rg"];
    let p0s = op.params.real["p0s"];
    let a_coeff = op.params.real["a_coeff"];
    let c = op.params.real["c"];
    let ca = op.params.real["ca"];
    let cb = op.params.real["cb"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let v1 = a_coeff * (std::f64::consts::FRAC_PI_4 + 0.5 * lat).tan().ln();
        let t = e * lat.sin();
        let v2 = 0.5 * e * a_coeff * ((1.0 + t) / (1.0 - t)).ln();
        let ps = 2.0 * (v1 - v2 + c).exp().atan() - std::f64::consts::FRAC_PI_2;
        let i1 = ps - p0s;
        let (sinps, cosps) = ps.sin_cos();
        let sinps2 = sinps * sinps;
        let cosps2 = cosps * cosps;
        let i4 = a_coeff * cosps;
        let i2 = 0.5 * a_coeff * i4 * sinps;
        let i3 = i2 * a_coeff * a_coeff * (5.0 * cosps2 - sinps2) / 12.0;
        let i6_base = i4 * a_coeff * a_coeff;
        let i5 = i6_base * (cosps2 - sinps2) / 6.0;
        let i6 = i6_base
            * a_coeff
            * a_coeff
            * (5.0 * cosps2 * cosps2 + sinps2 * (sinps2 - 18.0 * cosps2))
            / 120.0;
        let lam = frame.remove_central_meridian_raw(lon);
        let t2 = lam * lam;
        let mut x = k_rg * lam * (i4 + t2 * (i5 + t2 * i6));
        let mut y = k_rg * (i1 + t2 * (i2 + t2 * i3));
        let x2 = x * x;
        let y2 = y * y;
        let v1 = 3.0 * x * y2 - x * x2;
        let v2 = y * y2 - 3.0 * x2 * y;
        x += ca * v1 + cb * v2;
        y += ca * v2 - cb * v1;
        let (x, y) = frame.apply_false_origin(x * frame.a, y * frame.a);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let e = ellps.eccentricity();
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let k_0 = op.params.k(0);
    let k_rg = op.params.real["k_rg"];
    let p0s = op.params.real["p0s"];
    let a_coeff = op.params.real["a_coeff"];
    let c = op.params.real["c"];
    let ca = op.params.real["ca"];
    let cb = op.params.real["cb"];
    let cc = op.params.real["cc"];
    let cd = op.params.real["cd"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let (x, y) = frame.remove_false_origin(x, y);
        let mut x = x / frame.a;
        let mut y = y / frame.a;
        let x2 = x * x;
        let y2 = y * y;
        let v1 = 3.0 * x * y2 - x * x2;
        let v2 = y * y2 - 3.0 * x2 * y;
        let v3 = x * (5.0 * y2 * y2 + x2 * (-10.0 * y2 + x2));
        let v4 = y * (5.0 * x2 * x2 + y2 * (-10.0 * x2 + y2));
        x += -ca * v1 - cb * v2 + cc * v3 + cd * v4;
        y += cb * v1 - ca * v2 - cd * v3 + cc * v4;
        let ps = p0s + y / k_rg;
        let mut pe = ps + frame.lat_0 - p0s;
        let mut converged = false;
        for _ in 0..20 {
            let v1 = a_coeff * (std::f64::consts::FRAC_PI_4 + 0.5 * pe).tan().ln();
            let tpe = e * pe.sin();
            let v2 = 0.5 * e * a_coeff * ((1.0 + tpe) / (1.0 - tpe)).ln();
            let delta = ps - 2.0 * (v1 - v2 + c).exp().atan() + std::f64::consts::FRAC_PI_2;
            pe += delta;
            if delta.abs() < EPS {
                converged = true;
                break;
            }
        }
        if !converged {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        }

        let t = e * pe.sin();
        let t = 1.0 - t * t;
        let re = one_es / (t * t.sqrt());
        let tan_ps = ps.tan();
        let tan_ps2 = tan_ps * tan_ps;
        let s = k_rg * k_rg;
        let mut d = re * k_0 * k_rg;
        let i7 = tan_ps / (2.0 * d);
        let i8 = tan_ps * (5.0 + 3.0 * tan_ps2) / (24.0 * d * s);
        d = ps.cos() * k_rg * a_coeff;
        let i9 = 1.0 / d;
        d *= s;
        let i10 = (1.0 + 2.0 * tan_ps2) / (6.0 * d);
        let i11 = (5.0 + tan_ps2 * (28.0 + 24.0 * tan_ps2)) / (120.0 * d * s);
        let x2 = x * x;
        let lat = pe + x2 * (-i7 + i8 * x2);
        let lon = frame.lon_0 + x * (i9 + x2 * (-i10 + x2 * i11));
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

#[rustfmt::skip]
 pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "azi", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let ellps = params.ellps(0);
    let es = ellps.eccentricity_squared();
    let lat_0 = params.lat(0);
    if lat_0 == 0.0 {
        return Err(Error::General(
            "Labrd: Invalid value for lat_0: lat_0 should be different from 0",
        ));
    }
    let az = params.real["azi"].to_radians();
    let k_0 = params.real["k_0"];

    let sinp = lat_0.sin();
    let t = 1.0 - es * sinp * sinp;
    let n = 1.0 / t.sqrt();
    let r = (1.0 - es) * n / t;
    let k_rg = k_0 * (n * r).sqrt();
    let p0s = ((r / n).sqrt() * lat_0.tan()).atan();
    let a_coeff = sinp / p0s.sin();
    let t = ellps.eccentricity() * sinp;
    let c = 0.5 * ellps.eccentricity() * a_coeff * ((1.0 + t) / (1.0 - t)).ln()
        - a_coeff * (std::f64::consts::FRAC_PI_4 + 0.5 * lat_0).tan().ln()
        + (std::f64::consts::FRAC_PI_4 + 0.5 * p0s).tan().ln();
    let t2 = az + az;
    let cb = 1.0 / (12.0 * k_rg * k_rg);
    let ca = (1.0 - t2.cos()) * cb;
    let cb = cb * t2.sin();
    let cc = 3.0 * (ca * ca - cb * cb);
    let cd = 6.0 * ca * cb;

    params.real.insert("k_rg", k_rg);
    params.real.insert("p0s", p0s);
    params.real.insert("a_coeff", a_coeff);
    params.real.insert("c", c);
    params.real.insert("ca", ca);
    params.real.insert("cb", cb);
    params.real.insert("cc", cc);
    params.real.insert("cd", cd);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        state: None,
        steps: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labrd_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("labrd ellps=GRS80 lon_0=0.5 lat_0=2")?;

        let geo = [Coor4D::geo(1.0, 2.0, 0.0, 0.0)];
        let expected = [Coor4D::raw(
            166_973.166_090_228,
            -110_536.912_730_266,
            0.0,
            0.0,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-6);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
