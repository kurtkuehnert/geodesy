//! Oblique Mercator
//! Following IOGP Publication 373-7-2 – Geomatics Guidance Note number 7, part 2 – September 2019
//!
use crate::authoring::*;

const PARAMETER_TOLERANCE: f64 = 1e-7;
const ANGULAR_TOLERANCE: f64 = 1e-10;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 18] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "no_rot" },
    OpParameter::Flag { key: "no_off" },
    OpParameter::Flag { key: "variant" },
    OpParameter::Text { key: "ellps",   default: Some("GRS80") },
    OpParameter::Real { key: "lat_0",   default: Some(0_f64) },
    OpParameter::Real { key: "lon_0",   default: Some(0_f64) },
    OpParameter::Real { key: "x_0",     default: Some(0_f64) },
    OpParameter::Real { key: "y_0",     default: Some(0_f64) },
    OpParameter::Real { key: "k_0",     default: Some(1_f64) },
    OpParameter::Real { key: "latc",    default: Some(0_f64) },
    OpParameter::Real { key: "lonc",    default: Some(0_f64) },
    OpParameter::Real { key: "alpha",   default: Some(f64::NAN) },
    OpParameter::Real { key: "gamma_c", default: Some(f64::NAN) },
    OpParameter::Real { key: "lat_1",   default: Some(f64::NAN) },
    OpParameter::Real { key: "lon_1",   default: Some(0_f64) },
    OpParameter::Real { key: "lat_2",   default: Some(f64::NAN) },
    OpParameter::Real { key: "lon_2",   default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OmercMode {
    CentralPoint,
    TwoPoint,
}

#[derive(Clone, Copy, Debug)]
struct OmercState {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
    mode: OmercMode,
    latc: f64,
    lonc: f64,
    alpha: f64,
    gamma_c: f64,
    central_laborde: bool,
    no_rot: bool,
    no_off: bool,
    a_scale: f64,
    b_scale: f64,
    e_scale: f64,
    lam0: f64,
    singam: f64,
    cosgam: f64,
    sinrot: f64,
    cosrot: f64,
    ar_b: f64,
    br_a: f64,
    r_b: f64,
    u_0: f64,
    v_pole_n: f64,
    v_pole_s: f64,
}

impl OmercState {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let frame = ProjectionFrame::from_params(params);
        let a = frame.a;
        let es = ellps.eccentricity_squared();
        let e = es.sqrt();
        let latc = if params.given.contains_key("latc") {
            params.real["latc"].to_radians()
        } else {
            frame.lat_0
        };
        let k0 = frame.k_0;
        let mut alpha = params.real["alpha"].to_radians();
        let mut gamma_c = params.real["gamma_c"].to_radians();
        let has_alpha = !params.real["alpha"].is_nan();
        let has_gamma = !params.real["gamma_c"].is_nan();
        let central_laborde = has_alpha && !has_gamma;
        let lonc = if params.given.contains_key("lonc") {
            params.real["lonc"].to_radians()
        } else {
            frame.lon_0
        };
        let lam1 = params.lon(1);
        let phi1 = params.lat(1);
        let mut lam2 = params.lon(2);
        let phi2 = params.lat(2);
        let no_rot = params.boolean("no_rot");
        let no_off = params.boolean("no_off");
        let mode = if has_alpha || has_gamma {
            OmercMode::CentralPoint
        } else {
            OmercMode::TwoPoint
        };

        match mode {
            OmercMode::CentralPoint => {
                if latc.abs() >= FRAC_PI_2 - PARAMETER_TOLERANCE {
                    return Err(Error::General(
                        "Omerc: Invalid value for lat_0: |lat_0| should be < 90°",
                    ));
                }
            }
            OmercMode::TwoPoint => {
                if phi1.abs() > FRAC_PI_2 - PARAMETER_TOLERANCE {
                    return Err(Error::General(
                        "Omerc: Invalid value for lat_1: |lat_1| should be < 90°",
                    ));
                }
                if phi2.abs() > FRAC_PI_2 - PARAMETER_TOLERANCE {
                    return Err(Error::General(
                        "Omerc: Invalid value for lat_2: |lat_2| should be < 90°",
                    ));
                }
                if (phi1 - phi2).abs() <= PARAMETER_TOLERANCE {
                    return Err(Error::General(
                        "Omerc: Invalid value for lat_1/lat_2: lat_1 should be different from lat_2",
                    ));
                }
                if phi1.abs() <= PARAMETER_TOLERANCE {
                    return Err(Error::General(
                        "Omerc: Invalid value for lat_1: lat_1 should be different from 0",
                    ));
                }
                if latc.abs() >= FRAC_PI_2 - PARAMETER_TOLERANCE {
                    return Err(Error::General(
                        "Omerc: Invalid value for lat_0: |lat_0| should be < 90°",
                    ));
                }
            }
        }

        let com = (1.0 - es).sqrt();
        let (a_scale, b_scale, e_scale, d_scale) = if latc.abs() > ANGULAR_TOLERANCE {
            let (sinph0, cosph0) = latc.sin_cos();
            let con = 1.0 - es * sinph0 * sinph0;
            let b = (1.0 + es * cosph0.powi(4) / (1.0 - es)).sqrt();
            let a_scale = a * b * k0 * com / con;
            let d = b * com / (cosph0 * con.sqrt());
            let mut f = d * d - 1.0;
            if f <= 0.0 {
                f = 0.0;
            } else {
                f = f.sqrt();
                if latc < 0.0 {
                    f = -f;
                }
            }
            let e0 = (f + d) * tsfn(latc, sinph0, e).powf(b);
            (a_scale, b, e0, d)
        } else {
            (a * k0, 1.0 / com, 1.0, 1.0)
        };

        let (lam0, gamma0) = match mode {
            OmercMode::CentralPoint => {
                let gamma0 = if has_alpha {
                    let gamma0 = aasin(alpha.sin() / d_scale);
                    if !has_gamma {
                        gamma_c = alpha;
                    }
                    gamma0
                } else {
                    let gamma0 = gamma_c;
                    let test = d_scale * gamma0.sin();
                    if test.abs() > 1.0 {
                        return Err(Error::General("Omerc: Invalid value for gamma"));
                    }
                    alpha = aasin(test);
                    gamma0
                };
                let lam0 = lonc
                    - aasin(
                        0.5 * ((d_scale * d_scale - 1.0).max(0.0).sqrt() + d_scale
                            - 1.0 / ((d_scale * d_scale - 1.0).max(0.0).sqrt() + d_scale))
                            * gamma0.tan(),
                    ) / b_scale;
                (lam0, gamma0)
            }
            OmercMode::TwoPoint => {
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
                let lam0_raw = 0.5 * (lam1 + lam2)
                    - (j * (0.5 * b_scale * (lam1 - lam2)).tan() / p).atan() / b_scale;
                let lam0 = adjlon(lam0_raw);
                let denom = f - 1.0 / f;
                if denom == 0.0 {
                    return Err(Error::General("Omerc: Invalid value for eccentricity"));
                }
                let gamma0 = (2.0 * (b_scale * adjlon(lam1 - lam0)).sin() / denom).atan();
                alpha = aasin(d_scale * gamma0.sin());
                gamma_c = alpha;
                (lam0, gamma0)
            }
        };

        let singam = gamma0.sin();
        let cosgam = gamma0.cos();
        let sinrot = gamma_c.sin();
        let cosrot = gamma_c.cos();
        let r_b = 1.0 / b_scale;
        let ar_b = a_scale * r_b;
        let br_a = 1.0 / ar_b;
        let u_0 = if no_off {
            0.0
        } else {
            let mut value =
                (ar_b * (((d_scale * d_scale - 1.0).max(0.0).sqrt()) / alpha.cos()).atan()).abs();
            if latc < 0.0 {
                value = -value;
            }
            value
        };
        let f = 0.5 * gamma0;
        let v_pole_n = ar_b * (FRAC_PI_4 - f).tan().ln();
        let v_pole_s = ar_b * (FRAC_PI_4 + f).tan().ln();

        Ok(Self {
            frame,
            ellps,
            mode,
            latc,
            lonc,
            alpha,
            gamma_c,
            central_laborde,
            no_rot,
            no_off,
            a_scale,
            b_scale,
            e_scale,
            lam0,
            singam,
            cosgam,
            sinrot,
            cosrot,
            ar_b,
            br_a,
            r_b,
            u_0,
            v_pole_n,
            v_pole_s,
        })
    }
}

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

#[allow(clippy::too_many_arguments)]
fn forward_core(
    lon: f64,
    lat: f64,
    lon_0: f64,
    e: f64,
    b_scale: f64,
    e_scale: f64,
    ar_b: f64,
    a_scale: f64,
    singam: f64,
    cosgam: f64,
    v_pole_n: f64,
    v_pole_s: f64,
) -> Option<(f64, f64)> {
    let lon = adjlon(lon - lon_0);
    if (lat.abs() - FRAC_PI_2).abs() > ANGULAR_TOLERANCE {
        let w = e_scale / tsfn(lat, lat.sin(), e).powf(b_scale);
        let inv_w = 1.0 / w;
        let s = 0.5 * (w - inv_w);
        let t = 0.5 * (w + inv_w);
        let v_term = (b_scale * lon).sin();
        let u = (s * singam - v_term * cosgam) / t;
        if u.abs() >= 1.0 {
            return None;
        }
        let v = 0.5 * ar_b * ((1.0 - u) / (1.0 + u)).ln();
        let temp = (b_scale * lon).cos();
        let base_u = if temp.abs() < PARAMETER_TOLERANCE {
            a_scale * lon
        } else {
            ar_b * (s * cosgam + v_term * singam).atan2(temp)
        };
        Some((base_u, v))
    } else {
        let v = if lat > 0.0 { v_pole_n } else { v_pole_s };
        let base_u = ar_b * lat;
        Some((base_u, v))
    }
}

#[allow(clippy::too_many_arguments)]
fn inverse_core(
    u: f64,
    v: f64,
    e: f64,
    e_scale: f64,
    b_scale: f64,
    br_a: f64,
    r_b: f64,
    singam: f64,
    cosgam: f64,
    lon_0: f64,
) -> Option<(f64, f64)> {
    let q = (-br_a * v).exp();
    if q == 0.0 {
        return None;
    }
    let s = 0.5 * (q - 1.0 / q);
    let t = 0.5 * (q + 1.0 / q);
    let vv = (br_a * u).sin();
    let up = (vv * cosgam + s * singam) / t;
    if (up.abs() - 1.0).abs() < ANGULAR_TOLERANCE {
        let lat = if up < 0.0 { -FRAC_PI_2 } else { FRAC_PI_2 };
        return Some((lon_0, lat));
    }
    let ts = e_scale / ((1.0 + up) / (1.0 - up)).sqrt();
    let lat = ancillary::pj_phi2(ts.powf(1.0 / b_scale), e);
    let lon = lon_0 - r_b * (s * cosgam - vv * singam).atan2((br_a * u).cos());
    Some((lon, lat))
}

#[allow(non_snake_case)]
fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<OmercState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let Some((x, y)) = (match state.mode {
            OmercMode::TwoPoint => {
                match forward_core(
                    lon,
                    lat,
                    state.lam0,
                    state.ellps.eccentricity(),
                    state.b_scale,
                    state.e_scale,
                    state.ar_b,
                    state.a_scale,
                    state.singam,
                    state.cosgam,
                    state.v_pole_n,
                    state.v_pole_s,
                ) {
                    Some((base_u, v)) => {
                        if state.no_rot {
                            Some(state.frame.apply_false_origin(base_u, v))
                        } else {
                            let u = base_u - state.u_0;
                            Some(state.frame.apply_false_origin(
                                v * state.cosrot + u * state.sinrot,
                                u * state.cosrot - v * state.sinrot,
                            ))
                        }
                    }
                    None => None,
                }
            }
            OmercMode::CentralPoint => {
                let es = state.ellps.eccentricity_squared();
                let e = es.sqrt();
                let kc = state.frame.k_0;
                let laborde = state.central_laborde;

                let (s, c) = state.latc.sin_cos();
                let b = (1.0 + c.powi(4) * state.ellps.second_eccentricity_squared()).sqrt();
                let a =
                    state.ellps.semimajor_axis() * b * kc * (1.0 - es).sqrt() / (1.0 - es * s * s);
                let t0 = tsfn(state.latc, s, e);
                let d = b * (1.0 - es).sqrt() / (c * (1.0 - es * s * s).sqrt());
                let dd = if d < 1.0 { 0.0 } else { (d * d - 1.0).sqrt() };
                let f = d + dd * state.latc.signum();
                let h = f * t0.powf(b);
                let g = (f - 1.0 / f) / 2.0;
                let gamma_0 = aasin(state.alpha.sin() / d);
                let lambda_0 = state.lonc - aasin(g * gamma_0.tan()) / b;
                let (s0, c0) = gamma_0.sin_cos();
                let (sc, cc) = (if laborde { state.alpha } else { state.gamma_c }).sin_cos();

                match forward_core(
                    lon,
                    lat,
                    lambda_0,
                    e,
                    b,
                    h,
                    a / b,
                    a,
                    s0,
                    c0,
                    a / b * (FRAC_PI_4 - 0.5 * gamma_0).tan().ln(),
                    a / b * (FRAC_PI_4 + 0.5 * gamma_0).tan().ln(),
                ) {
                    Some((base_u, v)) => {
                        let variant = !state.no_off || state.central_laborde;
                        let base_u = if variant
                            && (state.alpha.abs() - FRAC_PI_2).abs() <= ANGULAR_TOLERANCE
                            && adjlon(lon - lambda_0) == 0.0
                        {
                            0.0
                        } else {
                            base_u
                        };

                        if state.no_rot {
                            Some(state.frame.apply_false_origin(base_u, v))
                        } else {
                            let u = if variant { base_u - state.u_0 } else { base_u };
                            Some(
                                state
                                    .frame
                                    .apply_false_origin(v * cc + u * sc, u * cc - v * sc),
                            )
                        }
                    }
                    None => None,
                }
            }
        }) else {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

#[allow(non_snake_case)]
fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<OmercState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (easting, northing) = operands.xy(i);
        let Some((lon, lat)) = (match state.mode {
            OmercMode::TwoPoint => {
                let (x, y) = state.frame.remove_false_origin(easting, northing);
                let (u, v) = if state.no_rot {
                    (x, y)
                } else {
                    (
                        y * state.cosrot + x * state.sinrot + state.u_0,
                        x * state.cosrot - y * state.sinrot,
                    )
                };
                inverse_core(
                    u,
                    v,
                    state.ellps.eccentricity(),
                    state.e_scale,
                    state.b_scale,
                    state.br_a,
                    state.r_b,
                    state.singam,
                    state.cosgam,
                    state.lam0,
                )
            }
            OmercMode::CentralPoint => {
                let es = state.ellps.eccentricity_squared();
                let e = es.sqrt();
                let kc = state.frame.k_0;

                let (s, c) = state.latc.sin_cos();
                let b = (1.0 + c.powi(4) * state.ellps.second_eccentricity_squared()).sqrt();
                let a =
                    state.ellps.semimajor_axis() * b * kc * (1.0 - es).sqrt() / (1.0 - es * s * s);
                let t0 = tsfn(state.latc, s, e);
                let d = b * (1.0 - es).sqrt() / (c * (1.0 - es * s * s).sqrt());
                let dd = if d < 1.0 { 0.0 } else { (d * d - 1.0).sqrt() };
                let f = d + dd * state.latc.signum();
                let h = f * t0.powf(b);
                let g = (f - 1.0 / f) / 2.0;
                let gamma_0 = aasin(state.alpha.sin() / d);
                let lambda_0 = state.lonc - aasin(g * gamma_0.tan()) / b;
                let (s0, c0) = gamma_0.sin_cos();
                let (x, y) = state.frame.remove_false_origin(easting, northing);
                let (u, v) = if state.no_rot {
                    (x, y)
                } else {
                    let offset = if !state.no_off || state.central_laborde {
                        state.u_0
                    } else {
                        0.0
                    };
                    (
                        y * state.cosrot + x * state.sinrot + offset,
                        x * state.cosrot - y * state.sinrot,
                    )
                };
                inverse_core(u, v, e, h, b, b / a, 1.0 / b, s0, c0, lambda_0)
            }
        }) else {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        };
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    let state = OmercState::new(&params)?;
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op::with_state(descriptor, params, state))
}

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
        let op = ctx.op(definition)?;

        let geo = [Coor2D::geo(5.3872535833, 115.8055054444)];
        let projected = [Coor2D::raw(679245.7281740266, 596562.7774687681)];

        let mut operands = geo;
        assert_eq!(1, ctx.apply(op, Fwd, &mut operands)?);
        for i in 0..operands.len() {
            assert_float_eq!(operands[i].0, projected[i].0, abs_all <= 1e-9);
        }

        assert_eq!(1, ctx.apply(op, Inv, &mut operands)?);
        for i in 0..operands.len() {
            assert_float_eq!(operands[i].0, geo[i].0, abs_all <= 1e-9);
        }

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 1e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }

        Ok(())
    }

    #[test]
    fn omerc_spherical_alpha_matches_proj_sample() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("omerc a=6400000 lat_0=45 alpha=35.264383770917604")?;

        let mut operands = [Coor4D::geo(1., 2., 0., 0.)];
        assert_eq!(1, ctx.apply(op, Fwd, &mut operands)?);
        assert_float_eq!(operands[0][0], -3569.825230822232, abs_all <= 1e-3);
        assert_float_eq!(operands[0][1], -5_093_592.310_871_85, abs_all <= 1e-3);

        Ok(())
    }
}
