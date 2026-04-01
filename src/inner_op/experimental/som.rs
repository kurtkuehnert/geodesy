//! Space Oblique Mercator family (`som`, `misrsom`, `lsat`)
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

const TOL: f64 = 1e-7;

fn parse_angle(value: &str, def: &str, key: &str) -> Result<f64, Error> {
    if let Some(stripped) = value.strip_suffix('r') {
        return stripped
            .parse::<f64>()
            .map_err(|_| Error::BadParam(key.to_string(), def.to_string()));
    }
    value
        .parse::<f64>()
        .map(|v| v.to_radians())
        .map_err(|_| Error::BadParam(key.to_string(), def.to_string()))
}

fn parse_real(value: &str, def: &str, key: &str) -> Result<f64, Error> {
    value
        .parse::<f64>()
        .map_err(|_| Error::BadParam(key.to_string(), def.to_string()))
}

#[allow(clippy::too_many_arguments)]
fn seraz0(
    lam_deg: f64,
    mult: f64,
    p22: f64,
    sa: f64,
    ca: f64,
    t: f64,
    w: f64,
    q: f64,
    xj: f64,
    coeffs: &mut [f64; 5],
) {
    let lam = lam_deg.to_radians();
    let sd = lam.sin();
    let sdsq = sd * sd;
    let s =
        p22 * sa * lam.cos() * ((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq))).sqrt();
    let d1 = 1.0 + q * sdsq;
    let h =
        ((1.0 + q * sdsq) / (1.0 + w * sdsq)).sqrt() * ((1.0 + w * sdsq) / (d1 * d1) - p22 * ca);
    let sq = (xj * xj + s * s).sqrt();
    let fc = mult * (h * xj - s * s) / sq;
    coeffs[2] += fc;
    coeffs[0] += fc * (2.0 * lam).cos();
    coeffs[1] += fc * (4.0 * lam).cos();
    let fc = mult * s * (h + xj) / sq;
    coeffs[3] += fc * lam.cos();
    coeffs[4] += fc * (3.0 * lam).cos();
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let p22 = op.params.real["p22"];
    let sa = op.params.real["sa"];
    let ca = op.params.real["ca"];
    let xj = op.params.real["xj"];
    let q = op.params.real["q"];
    let t = op.params.real["t"];
    let w = op.params.real["w"];
    let rlm = op.params.real["rlm"];
    let rlm2 = op.params.real["rlm2"];
    let a2 = op.params.real["a2"];
    let a4 = op.params.real["a4"];
    let b = op.params.real["b"];
    let c1 = op.params.real["c1"];
    let c3 = op.params.real["c3"];
    let one_es = 1.0 - ellps.eccentricity_squared();
    let es = ellps.eccentricity_squared();
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (mut lon, mut phi) = operands.xy(i);
        lon -= op.params.real["asc_lon"];
        phi = phi.clamp(-FRAC_PI_2, FRAC_PI_2);
        let mut lampp = if phi >= 0.0 {
            FRAC_PI_2
        } else {
            PI + FRAC_PI_2
        };
        let tanphi = phi.tan();
        let mut lamt = 0.0;
        let mut lamdp = 0.0;
        let mut nn = 0;
        let found = loop {
            let sav0 = lampp;
            let lamtp = lon + p22 * lampp;
            let cl = lamtp.cos();
            let fac = if cl < 0.0 {
                lampp + lampp.sin() * FRAC_PI_2
            } else {
                lampp - lampp.sin() * FRAC_PI_2
            };
            let mut sav = sav0;
            let mut l_ok = false;
            for _ in 0..=50 {
                lamt = lon + p22 * sav;
                let c = if lamt.cos().abs() < TOL {
                    lamt.cos().signum() * TOL
                } else {
                    lamt.cos()
                };
                let xlam = (one_es * tanphi * sa + lamt.sin() * ca) / c;
                lamdp = xlam.atan() + fac;
                if (sav.abs() - lamdp.abs()).abs() < TOL {
                    l_ok = true;
                    break;
                }
                sav = lamdp;
            }
            if l_ok || nn >= 2 || (lamdp > rlm && lamdp < rlm2) {
                break l_ok;
            }
            nn += 1;
            lampp = if lamdp <= rlm {
                2.0 * PI + FRAC_PI_2
            } else {
                FRAC_PI_2
            };
        };
        if !found {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let sp = phi.sin();
        let phidp = ((one_es * ca * sp - sa * phi.cos() * lamt.sin())
            / (1.0 - es * sp * sp).sqrt())
        .clamp(-1.0, 1.0)
        .asin();
        let tanph = (FRAC_PI_4 + 0.5 * phidp).tan().ln();
        let sd = lamdp.sin();
        let sdsq = sd * sd;
        let s = p22
            * sa
            * lamdp.cos()
            * ((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq))).sqrt();
        let d = (xj * xj + s * s).sqrt();
        let x = b * lamdp + a2 * (2.0 * lamdp).sin() + a4 * (4.0 * lamdp).sin() - tanph * s / d;
        let y = c1 * sd + c3 * (3.0 * lamdp).sin() + tanph * xj / d;
        operands.set_xy(i, frame.a * x, frame.a * y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let p22 = op.params.real["p22"];
    let sa = op.params.real["sa"];
    let ca = op.params.real["ca"];
    let xj = op.params.real["xj"];
    let q = op.params.real["q"];
    let t = op.params.real["t"];
    let w = op.params.real["w"];
    let a2 = op.params.real["a2"];
    let a4 = op.params.real["a4"];
    let b = op.params.real["b"];
    let c1 = op.params.real["c1"];
    let c3 = op.params.real["c3"];
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let rone_es = 1.0 / one_es;
    let u = op.params.real["u"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let x = operands.xy(i).0 / frame.a;
        let y = operands.xy(i).1 / frame.a;
        let mut lamdp = x / b;
        let mut sd = 0.0;
        let mut s = 0.0;
        for _ in 0..50 {
            let sav = lamdp;
            sd = lamdp.sin();
            let sdsq = sd * sd;
            s = p22
                * sa
                * lamdp.cos()
                * ((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq))).sqrt();
            lamdp = (x + y * s / xj
                - a2 * (2.0 * lamdp).sin()
                - a4 * (4.0 * lamdp).sin()
                - s / xj * (c1 * lamdp.sin() + c3 * (3.0 * lamdp).sin()))
                / b;
            if (lamdp - sav).abs() < TOL {
                break;
            }
        }
        let fac =
            ((1.0 + s * s / (xj * xj)).sqrt() * (y - c1 * sd - c3 * (3.0 * lamdp).sin())).exp();
        let phidp = 2.0 * (fac.atan() - FRAC_PI_4);
        let dd = sd * sd;
        let lamdp_safe = if lamdp.cos().abs() < TOL {
            lamdp - TOL
        } else {
            lamdp
        };
        let spp = phidp.sin();
        let sppsq = spp * spp;
        let denom = 1.0 - sppsq * (1.0 + u);
        if denom == 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let lamt = (((1.0 - sppsq * rone_es) * lamdp_safe.tan() * ca
            - spp * sa * ((1.0 + q * dd) * (1.0 - sppsq) - sppsq * u).sqrt() / lamdp_safe.cos())
            / denom)
            .atan();
        let sl = if lamt >= 0.0 { 1.0 } else { -1.0 };
        let scl = if lamdp_safe.cos() >= 0.0 { 1.0 } else { -1.0 };
        let lamt = lamt - FRAC_PI_2 * (1.0 - scl) * sl;
        let lon = lamt - p22 * lamdp + op.params.real["asc_lon"];
        let lat = if sa.abs() < TOL {
            (spp / (one_es * one_es + es * sppsq).sqrt())
                .clamp(-1.0, 1.0)
                .asin()
        } else {
            ((lamdp.tan() * lamt.cos() - ca * lamt.sin()) / (one_es * sa)).atan()
        };
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

fn setup(
    def: &str,
    params: &mut ParsedParameters,
    alf: f64,
    p22: f64,
    rlm: f64,
) -> Result<Op, Error> {
    if !(0.0..=PI).contains(&alf) {
        return Err(Error::BadParam("inc_angle".to_string(), def.to_string()));
    }
    if p22 < 0.0 {
        return Err(Error::BadParam("ps_rev".to_string(), def.to_string()));
    }
    let ellps = params.ellps(0);
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let rone_es = 1.0 / one_es;
    let sa = alf.sin();
    let mut ca = alf.cos();
    if ca.abs() < 1e-9 {
        ca = 1e-9;
    }
    let esc = es * ca * ca;
    let ess = es * sa * sa;
    let w = (1.0 - esc) * rone_es;
    let w = w * w - 1.0;
    let q = ess * rone_es;
    let t = ess * (2.0 - es) * rone_es * rone_es;
    let u = esc * rone_es;
    let xj = one_es * one_es * one_es;
    let rlm2 = rlm + 2.0 * PI;
    let mut coeffs = [0.0; 5];
    seraz0(0.0, 1.0, p22, sa, ca, t, w, q, xj, &mut coeffs);
    for lam in [9.0, 27.0, 45.0, 63.0, 81.0] {
        seraz0(lam, 4.0, p22, sa, ca, t, w, q, xj, &mut coeffs);
    }
    for lam in [18.0, 36.0, 54.0, 72.0] {
        seraz0(lam, 2.0, p22, sa, ca, t, w, q, xj, &mut coeffs);
    }
    seraz0(90.0, 1.0, p22, sa, ca, t, w, q, xj, &mut coeffs);
    params.real.insert("a2", coeffs[0] / 30.0);
    params.real.insert("a4", coeffs[1] / 60.0);
    params.real.insert("b", coeffs[2] / 30.0);
    params.real.insert("c1", coeffs[3] / 15.0);
    params.real.insert("c3", coeffs[4] / 45.0);
    params.real.insert("p22", p22);
    params.real.insert("sa", sa);
    params.real.insert("ca", ca);
    params.real.insert("q", q);
    params.real.insert("t", t);
    params.real.insert("u", u);
    params.real.insert("w", w);
    params.real.insert("xj", xj);
    params.real.insert("rlm", rlm);
    params.real.insert("rlm2", rlm2);
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params: params.clone(),
        state: None,
        steps: None,
    })
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(
        parameters,
        &[
            OpParameter::Flag { key: "inv" },
            OpParameter::Text {
                key: "ellps",
                default: Some("GRS80"),
            },
        ],
    )?;
    let given = parameters.instantiated_as.split_into_parameters();
    let alf = parse_angle(
        given
            .get("inc_angle")
            .ok_or_else(|| Error::MissingParam("inc_angle".to_string()))?,
        def,
        "inc_angle",
    )?;
    let p22 = parse_real(
        given
            .get("ps_rev")
            .ok_or_else(|| Error::MissingParam("ps_rev".to_string()))?,
        def,
        "ps_rev",
    )?;
    let asc_lon = parse_angle(
        given
            .get("asc_lon")
            .ok_or_else(|| Error::MissingParam("asc_lon".to_string()))?,
        def,
        "asc_lon",
    )?;
    params.real.insert("asc_lon", asc_lon);
    setup(def, &mut params, alf, p22, 0.0)
}

pub fn misr(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(
        parameters,
        &[
            OpParameter::Flag { key: "inv" },
            OpParameter::Text {
                key: "ellps",
                default: Some("GRS80"),
            },
            OpParameter::Natural {
                key: "path",
                default: None,
            },
        ],
    )?;
    let path = params.natural("path")?;
    if path == 0 || path > 233 {
        return Err(Error::BadParam("path".to_string(), def.clone()));
    }
    params.real.insert(
        "asc_lon",
        (129.3056f64).to_radians() - 2.0 * PI / 233.0 * path as f64,
    );
    setup(
        def,
        &mut params,
        (98.30382f64).to_radians(),
        98.88 / 1440.0,
        0.0,
    )
}

pub fn lsat(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(
        parameters,
        &[
            OpParameter::Flag { key: "inv" },
            OpParameter::Text {
                key: "ellps",
                default: Some("GRS80"),
            },
            OpParameter::Natural {
                key: "lsat",
                default: None,
            },
            OpParameter::Natural {
                key: "path",
                default: None,
            },
        ],
    )?;
    let land = params.natural("lsat")?;
    if land == 0 || land > 5 {
        return Err(Error::BadParam("lsat".to_string(), def.clone()));
    }
    let path = params.natural("path")?;
    let max_path = if land <= 3 { 251 } else { 233 };
    if path == 0 || path > max_path {
        return Err(Error::BadParam("path".to_string(), def.clone()));
    }
    let (alf, mut p22, asc_lon, rlm) = if land <= 3 {
        (
            (99.092f64).to_radians(),
            103.2669323 / 1440.0,
            (128.87f64).to_radians() - 2.0 * PI / 251.0 * path as f64,
            PI * (1.0 / 248.0 + 0.5161290322580645),
        )
    } else {
        (
            (98.2f64).to_radians(),
            98.8841202 / 1440.0,
            (129.3f64).to_radians() - 2.0 * PI / 233.0 * path as f64,
            PI * (1.0 / 248.0 + 0.5161290322580645),
        )
    };
    params.real.insert("asc_lon", asc_lon);
    p22 += 0.0;
    setup(def, &mut params, alf, p22, rlm)
}
