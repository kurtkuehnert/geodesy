//! Oblique Mercator
//! Following IOGP Publication 373-7-2 – Geomatics Guidance Note number 7, part 2 – September 2019
//!
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, TAU};

const TOL: f64 = 1e-7;
const EPS: f64 = 1e-10;

fn tsfn(phi: f64, sinphi: f64, e: f64) -> f64 {
    (FRAC_PI_4 - 0.5 * phi).tan() / ((1.0 - e * sinphi) / (1.0 + e * sinphi)).powf(0.5 * e)
}

fn adjlon(mut lon: f64) -> f64 {
    while lon <= -std::f64::consts::PI {
        lon += TAU;
    }
    while lon > std::f64::consts::PI {
        lon -= TAU;
    }
    lon
}

fn aasin(value: f64) -> f64 {
    value.clamp(-1.0, 1.0).asin()
}

fn fwd_central_with_ctx(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    fwd_central(op, operands)
}

fn inv_central_with_ctx(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    inv_central(op, operands)
}

#[allow(non_snake_case)]
fn fwd_central(op: &Op, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let es = ellps.eccentricity_squared();
    let e = es.sqrt();

    let kc = op.params.k(0);
    let fe = op.params.x(0);
    let fn_ = op.params.y(0);
    let ec = fe;
    let nc = fn_;
    let no_rot = op.params.boolean("no_rot");
    let u_0 = op.params.real["u_0"];

    let latc = op.params.real["latc"].to_radians();
    let lonc = op.params.real["lonc"].to_radians();

    let alpha = op.params.real["alpha"];
    let ninety = alpha == 90_f64;
    let alpha = alpha.to_radians();

    let mut gamma_c = op.params.real["gamma_c"];
    let laborde = gamma_c.is_nan();
    gamma_c = gamma_c.to_radians();

    let no_off = op.params.boolean("no_off");
    let mut variant = !no_off;
    if laborde {
        variant = true;
        gamma_c = alpha;
    }

    let (s, c) = latc.sin_cos();
    let b = (1_f64 + c.powi(4) * ellps.second_eccentricity_squared()).sqrt();
    let a = ellps.semimajor_axis() * b * kc * (1_f64 - es).sqrt() / (1.0 - es * s * s);
    let t0 = tsfn(latc, s, e);
    let d = b * (1.0 - es).sqrt() / (c * (1.0 - es * s * s).sqrt());
    let dd = if d < 1.0 { 0.0 } else { (d * d - 1.0).sqrt() };
    let f = d + dd * latc.signum();
    let h = f * t0.powf(b);
    let g = (f - 1.0 / f) / 2.0;
    let gamma_0 = aasin(alpha.sin() / d);
    let lambda_0 = lonc - aasin(g * gamma_0.tan()) / b;
    let _uc = if ninety {
        a * (lonc - lambda_0)
    } else {
        (a / b) * dd.atan2(alpha.cos()) * latc.signum()
    };

    let (s0, c0) = gamma_0.sin_cos();
    let (sc, cc) = gamma_c.sin_cos();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lon = adjlon(lon - lambda_0);
        let slat = lat.sin();

        let t = tsfn(lat, slat, e);
        let q = h / t.powf(b);
        let s = (q - 1.0 / q) / 2.0;
        let t = (q + 1.0 / q) / 2.0;
        let v = (b * lon).sin();
        let u = (s * s0 - v * c0) / t;
        let v_coord = a * ((1.0 - u) / (1.0 + u)).ln() / (2.0 * b);

        let cblon = (b * lon).cos();
        let base_u = a * (s * c0 + v * s0).atan2(cblon) / b;

        if !variant {
            let (x, y) = if no_rot {
                (base_u + fe, v_coord + fn_)
            } else {
                (
                    v_coord * cc + base_u * sc + fe,
                    base_u * cc - v_coord * sc + fn_,
                )
            };
            operands.set_xy(i, x, y);
            successes += 1;
            continue;
        }

        if ninety {
            let base_u = if lon == lambda_0 { 0.0 } else { base_u };
            let (x, y) = if no_rot {
                (base_u + ec, v_coord + nc)
            } else {
                let u = base_u - u_0;
                (v_coord * cc + u * sc + ec, u * cc - v_coord * sc + nc)
            };
            operands.set_xy(i, x, y);
            successes += 1;
            continue;
        }

        let (x, y) = if no_rot {
            (base_u + ec, v_coord + nc)
        } else {
            let u = base_u - u_0;
            (v_coord * cc + u * sc + ec, u * cc - v_coord * sc + nc)
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

#[allow(non_snake_case)]
fn inv_central(op: &Op, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let es = ellps.eccentricity_squared();
    let e = es.sqrt();

    let kc = op.params.k(0);
    let fe = op.params.x(0);
    let fn_ = op.params.y(0);
    let no_rot = op.params.boolean("no_rot");
    let u_0 = op.params.real["u_0"];

    let latc = op.params.real["latc"].to_radians();
    let lonc = op.params.real["lonc"].to_radians();

    let alpha = op.params.real["alpha"];
    let ninety = alpha == 90_f64;
    let alpha = alpha.to_radians();

    let gamma_c = op.params.real["gamma_c"];
    let laborde = gamma_c.is_nan();
    let gamma_c = if laborde { alpha } else { gamma_c.to_radians() };
    let no_off = op.params.boolean("no_off");
    let variant = !no_off || laborde;

    let (s, c) = latc.sin_cos();
    let b = (1_f64 + c.powi(4) * ellps.second_eccentricity_squared()).sqrt();
    let a = ellps.semimajor_axis() * b * kc * (1_f64 - es).sqrt() / (1.0 - es * s * s);
    let t0 = tsfn(latc, s, e);
    let d = b * (1.0 - es).sqrt() / (c * (1.0 - es * s * s).sqrt());
    let dd = if d < 1.0 { 0.0 } else { (d * d - 1.0).sqrt() };
    let f = d + dd * latc.signum();
    let h = f * t0.powf(b);
    let g = (f - 1.0 / f) / 2.0;
    let gamma_0 = aasin(alpha.sin() / d);
    let lambda_0 = lonc - aasin(g * gamma_0.tan()) / b;

    let _uc = if ninety {
        a * (lonc - lambda_0)
    } else {
        (a / b) * dd.atan2(alpha.cos()) * latc.signum()
    };

    let (s0, c0) = gamma_0.sin_cos();
    let (sc, cc) = gamma_c.sin_cos();
    let offset = if variant { u_0 } else { 0.0 };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (easting, northing) = operands.xy(i);

        let (v, u_coord) = if no_rot {
            (northing - fn_, easting - fe)
        } else {
            (
                (easting - fe) * cc - (northing - fn_) * sc,
                (northing - fn_) * cc + (easting - fe) * sc + offset,
            )
        };

        let q = (-b * v / a).exp();
        let s = (q - 1.0 / q) / 2.0;
        let t = (q + 1.0 / q) / 2.0;
        let v = (b * u_coord / a).sin();
        let u = (v * c0 + s * s0) / t;
        let t = (h / ((1.0 + u) / (1.0 - u)).sqrt()).powf(1.0 / b);

        let chi = FRAC_PI_2 - 2.0 * t.atan();
        let f = [
            0.5 + es * (5.0 / 24.0 + es * (1.0 / 12.0 + es * 13.0 / 360.0)),
            es * (7.0 / 48.0 + es * (29.0 / 240.0 + es * 811.0 / 11520.0)),
            es * es * (7.0 / 120.0 + es * 81.0 / 1120.0),
            es * es * es * 4279.0 / 161280.0,
        ];
        let sines = [
            (2.0 * chi).sin(),
            (4.0 * chi).sin(),
            (6.0 * chi).sin(),
            (8.0 * chi).sin(),
        ];

        let lat =
            chi + es * (f[0] * sines[0] + f[1] * sines[1] + f[2] * sines[2] + f[3] * sines[3]);
        let lon = lambda_0 - (s * c0 - v * s0).atan2((b * u_coord / a).cos()) / b;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

// ----- F O R W A R D -----------------------------------------------------------------

#[allow(non_snake_case)]
fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let es = ellps.eccentricity_squared();
    let e = es.sqrt();
    let false_easting = op.params.x(0);
    let false_northing = op.params.y(0);
    let lon_0 = op.params.real["lon_0"];
    let A = op.params.real["A"];
    let B = op.params.real["B"];
    let E = op.params.real["E"];
    let ArB = op.params.real["ArB"];
    let v_pole_n = op.params.real["v_pole_n"];
    let v_pole_s = op.params.real["v_pole_s"];
    let singam = op.params.real["singam"];
    let cosgam = op.params.real["cosgam"];
    let sinrot = op.params.real["sinrot"];
    let cosrot = op.params.real["cosrot"];
    let u_0 = op.params.real["u_0"];
    let no_rot = op.params.boolean("no_rot");

    let mut successes = 0_usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lon = adjlon(lon - lon_0);
        let (u, v) = if (lat.abs() - FRAC_PI_2).abs() > EPS {
            let w = E / tsfn(lat, lat.sin(), e).powf(B);
            let inv_w = 1.0 / w;
            let s = 0.5 * (w - inv_w);
            let t = 0.5 * (w + inv_w);
            let v_term = (B * lon).sin();
            let u = (s * singam - v_term * cosgam) / t;
            if u.abs() >= 1.0 {
                operands.set_xy(i, f64::NAN, f64::NAN);
                continue;
            }
            let v = 0.5 * ArB * ((1.0 - u) / (1.0 + u)).ln();
            let temp = (B * lon).cos();
            let u = if temp.abs() < TOL {
                A * lon
            } else {
                ArB * (s * cosgam + v_term * singam).atan2(temp)
            };
            (u, v)
        } else {
            let v = if lat > 0.0 { v_pole_n } else { v_pole_s };
            let u = ArB * lat;
            (u, v)
        };

        let (x, y) = if no_rot {
            (u + false_easting, v + false_northing)
        } else {
            let u = u - u_0;
            (
                v * cosrot + u * sinrot + false_easting,
                u * cosrot - v * sinrot + false_northing,
            )
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

// ----- I N V E R S E -----------------------------------------------------------------

#[allow(non_snake_case)]
fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let false_easting = op.params.x(0);
    let false_northing = op.params.y(0);
    let lon_0 = op.params.real["lon_0"];
    let B = op.params.real["B"];
    let E = op.params.real["E"];
    let BrA = op.params.real["BrA"];
    let rB = op.params.real["rB"];
    let singam = op.params.real["singam"];
    let cosgam = op.params.real["cosgam"];
    let sinrot = op.params.real["sinrot"];
    let cosrot = op.params.real["cosrot"];
    let u_0 = op.params.real["u_0"];
    let no_rot = op.params.boolean("no_rot");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let x = x - false_easting;
        let y = y - false_northing;
        let (u, v) = if no_rot {
            (x, y)
        } else {
            (y * cosrot + x * sinrot + u_0, x * cosrot - y * sinrot)
        };
        let q = (-BrA * v).exp();
        if q == 0.0 {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        }
        let s = 0.5 * (q - 1.0 / q);
        let t = 0.5 * (q + 1.0 / q);
        let vv = (BrA * u).sin();
        let up = (vv * cosgam + s * singam) / t;
        let (lat, lon) = if (up.abs() - 1.0).abs() < EPS {
            (if up < 0.0 { -FRAC_PI_2 } else { FRAC_PI_2 }, 0.0)
        } else {
            let ts = E / ((1.0 + up) / (1.0 - up)).sqrt();
            let lat = ancillary::pj_phi2(ts.powf(1.0 / B), e);
            let lon = lon_0 - rB * (s * cosgam - vv * singam).atan2((BrA * u).cos());
            (lat, lon)
        };
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

// ----- C O N S T R U C T O R ---------------------------------------------------------

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 18] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "no_rot" },
    OpParameter::Flag { key: "no_off" },
    OpParameter::Flag { key: "variant" },
    OpParameter::Text { key: "ellps",  default: Some("GRS80") },
    OpParameter::Real { key: "latc",  default: Some(0_f64) },
    OpParameter::Real { key: "lonc",  default: Some(0_f64) },
    OpParameter::Real { key: "alpha",  default: Some(f64::NAN) },
    OpParameter::Real { key: "gamma_c",  default: Some(f64::NAN) },
    OpParameter::Real { key: "lat_1", default: Some(f64::NAN) },
    OpParameter::Real { key: "lat_2", default: Some(f64::NAN) },
    OpParameter::Real { key: "lon_1", default: Some(0_f64) },
    OpParameter::Real { key: "lon_2", default: Some(0_f64) },
    OpParameter::Real { key: "x_0",    default: Some(0_f64) },
    OpParameter::Real { key: "y_0",    default: Some(0_f64) },
    OpParameter::Real { key: "k_0",    default: Some(1_f64) },
    OpParameter::Real { key: "lat_0",  default: Some(0_f64) },
    OpParameter::Real { key: "lon_0",  default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    if given.contains_key("a") && !given.contains_key("ellps") {
        let a = given["a"]
            .parse::<f64>()
            .map_err(|_| Error::MissingParam("a".to_string()))?;
        if a <= 0.0 {
            return Err(Error::General(
                "Omerc: Invalid value for a: a must be positive",
            ));
        }

        let rf = if let Some(rf) = given.get("rf") {
            rf.parse::<f64>()
                .map_err(|_| Error::MissingParam("rf".to_string()))?
        } else if let Some(f) = given.get("f") {
            let f = f
                .parse::<f64>()
                .map_err(|_| Error::MissingParam("f".to_string()))?;
            if f == 0.0 { 0.0 } else { 1.0 / f }
        } else if let Some(b) = given.get("b") {
            let b = b
                .parse::<f64>()
                .map_err(|_| Error::MissingParam("b".to_string()))?;
            if b <= 0.0 {
                return Err(Error::General(
                    "Omerc: Invalid value for b: b must be positive",
                ));
            }
            if (a - b).abs() < f64::EPSILON {
                0.0
            } else {
                a / (a - b)
            }
        } else {
            0.0
        };
        params.text.insert("ellps", format!("{a},{rf}"));
    }
    let ellps = params.ellps(0);
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let e = es.sqrt();
    let phi0 = params.real["latc"].to_radians();
    let k0 = params.k(0);
    let mut alpha_c = params.real["alpha"].to_radians();
    let mut gamma = params.real["gamma_c"].to_radians();
    let has_alpha = !params.real["alpha"].is_nan();
    let has_gamma = !params.real["gamma_c"].is_nan();
    let lonc = params.real["lonc"].to_radians();
    let lam1 = params.real["lon_1"].to_radians();
    let phi1 = params.real["lat_1"].to_radians();
    let mut lam2 = params.real["lon_2"].to_radians();
    let phi2 = params.real["lat_2"].to_radians();

    let no_rot = params.boolean("no_rot");
    let no_off = params.boolean("no_off");

    if has_alpha || has_gamma {
        if phi0.abs() >= FRAC_PI_2 - TOL {
            return Err(Error::General(
                "Omerc: Invalid value for lat_0: |lat_0| should be < 90°",
            ));
        }
    } else {
        if phi1.abs() > FRAC_PI_2 - TOL {
            return Err(Error::General(
                "Omerc: Invalid value for lat_1: |lat_1| should be < 90°",
            ));
        }
        if phi2.abs() > FRAC_PI_2 - TOL {
            return Err(Error::General(
                "Omerc: Invalid value for lat_2: |lat_2| should be < 90°",
            ));
        }
        if (phi1 - phi2).abs() <= TOL {
            return Err(Error::General(
                "Omerc: Invalid value for lat_1/lat_2: lat_1 should be different from lat_2",
            ));
        }
        if phi1.abs() <= TOL {
            return Err(Error::General(
                "Omerc: Invalid value for lat_1: lat_1 should be different from 0",
            ));
        }
        if phi0.abs() >= FRAC_PI_2 - TOL {
            return Err(Error::General(
                "Omerc: Invalid value for lat_0: |lat_0| should be < 90°",
            ));
        }
    }

    let com = (1.0 - es).sqrt();
    let (a_scale, b_scale, e_scale, d_scale) = if phi0.abs() > EPS {
        let (sinph0, cosph0) = phi0.sin_cos();
        let con = 1.0 - es * sinph0 * sinph0;
        let b = (1.0 + es * cosph0.powi(4) / (1.0 - es)).sqrt();
        let a_scale = a * b * k0 * com / con;
        let d = b * com / (cosph0 * con.sqrt());
        let mut f = d * d - 1.0;
        if f <= 0.0 {
            f = 0.0;
        } else {
            f = f.sqrt();
            if phi0 < 0.0 {
                f = -f;
            }
        }
        let e0 = (f + d) * tsfn(phi0, sinph0, e).powf(b);
        (a_scale, b, e0, d)
    } else {
        (a * k0, 1.0 / com, 1.0, 1.0)
    };

    let lam0;
    let gamma0;
    if has_alpha || has_gamma {
        if has_alpha {
            gamma0 = aasin(alpha_c.sin() / d_scale);
            if !has_gamma {
                gamma = alpha_c;
            }
        } else {
            gamma0 = gamma;
            let test = d_scale * gamma0.sin();
            if test.abs() > 1.0 {
                return Err(Error::General("Omerc: Invalid value for gamma"));
            }
            alpha_c = aasin(test);
        }
        lam0 = lonc
            - aasin(
                0.5 * ((d_scale * d_scale - 1.0).max(0.0).sqrt() + d_scale
                    - 1.0 / ((d_scale * d_scale - 1.0).max(0.0).sqrt() + d_scale))
                    * gamma0.tan(),
            ) / b_scale;
    } else {
        let h = tsfn(phi1, phi1.sin(), e).powf(b_scale);
        let l = tsfn(phi2, phi2.sin(), e).powf(b_scale);
        let f = e_scale / h;
        let p = (l - h) / (l + h);
        if p == 0.0 {
            return Err(Error::General("Omerc: Invalid value for eccentricity"));
        }
        let mut j = e_scale * e_scale;
        j = (j - l * h) / (j + l * h);
        let con = lam1 - lam2;
        if con < -std::f64::consts::PI {
            lam2 -= TAU;
        } else if con > std::f64::consts::PI {
            lam2 += TAU;
        }
        let lam0_raw =
            0.5 * (lam1 + lam2) - (j * (0.5 * b_scale * (lam1 - lam2)).tan() / p).atan() / b_scale;
        lam0 = adjlon(lam0_raw);
        let denom = f - 1.0 / f;
        if denom == 0.0 {
            return Err(Error::General("Omerc: Invalid value for eccentricity"));
        }
        gamma0 = (2.0 * (b_scale * adjlon(lam1 - lam0)).sin() / denom).atan();
        alpha_c = aasin(d_scale * gamma0.sin());
        gamma = alpha_c;
    }

    let singam = gamma0.sin();
    let cosgam = gamma0.cos();
    let sinrot = gamma.sin();
    let cosrot = gamma.cos();
    let r_b = 1.0 / b_scale;
    let ar_b = a_scale * r_b;
    let br_a = 1.0 / ar_b;
    let ab = a_scale * b_scale;
    let u_0 = if no_off {
        0.0
    } else {
        let mut u0 =
            (ar_b * (((d_scale * d_scale - 1.0).max(0.0).sqrt()) / alpha_c.cos()).atan()).abs();
        if phi0 < 0.0 {
            u0 = -u0;
        }
        u0
    };
    let f = 0.5 * gamma0;
    let v_pole_n = ar_b * (FRAC_PI_4 - f).tan().ln();
    let v_pole_s = ar_b * (FRAC_PI_4 + f).tan().ln();

    params.real.insert("A", a_scale);
    params.real.insert("B", b_scale);
    params.real.insert("E", e_scale);
    params.real.insert("AB", ab);
    params.real.insert("ArB", ar_b);
    params.real.insert("BrA", br_a);
    params.real.insert("rB", r_b);
    params.real.insert("singam", singam);
    params.real.insert("cosgam", cosgam);
    params.real.insert("sinrot", sinrot);
    params.real.insert("cosrot", cosrot);
    params.real.insert("v_pole_n", v_pole_n);
    params.real.insert("v_pole_s", v_pole_s);
    params.real.insert("u_0", u_0);
    params.real.insert("lon_0", lam0);
    if no_rot {
        params.boolean.insert("no_rot");
    }
    let descriptor = if has_alpha || has_gamma {
        OpDescriptor::new(
            def,
            InnerOp(fwd_central_with_ctx),
            Some(InnerOp(inv_central_with_ctx)),
        )
    } else {
        OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)))
    };

    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use float_eq::assert_float_eq;

    #[test]
    fn omerc() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "
            omerc ellps=evrstSS variant
            x_0=590476.87 y_0=442857.65
            latc=4 lonc=115
            k_0=0.99984 alpha=53:18:56.9537 gamma_c=53:07:48.3685
        ";
        // k_0=0.99984 alpha=53.3158204722 gamma_c=53.1301023611
        let op = ctx.op(definition)?;

        // Validation value from EPSG
        let geo = [Coor2D::geo(5.3872535833, 115.8055054444)];
        let projected = [Coor2D::raw(679245.7281740266, 596562.7774687681)];

        // Forward
        let mut operands = geo;

        assert_eq!(1, ctx.apply(op, Fwd, &mut operands)?);
        for i in 0..operands.len() {
            assert_float_eq!(operands[i].0, projected[i].0, abs_all <= 1e-9);
        }

        // Roundtrip
        assert_eq!(1, ctx.apply(op, Inv, &mut operands)?);
        for i in 0..operands.len() {
            assert_float_eq!(operands[i].0, geo[i].0, abs_all <= 1e-9);
        }

        // Forward
        let mut operands = geo;

        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 1e-9);
        }

        // Roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }

        Ok(())
    }

    #[test]
    fn omerc_spherical_alpha_matches_proj_sample() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("omerc a=6400000 latc=45 alpha=35.264383770917604")?;

        let mut operands = [Coor4D::geo(1., 2., 0., 0.)];
        assert_eq!(1, ctx.apply(op, Fwd, &mut operands)?);
        assert_float_eq!(operands[0][0], -3569.825230822232, abs_all <= 1e-3);
        assert_float_eq!(operands[0][1], -5_093_592.310_871_85, abs_all <= 1e-3);

        Ok(())
    }
}
