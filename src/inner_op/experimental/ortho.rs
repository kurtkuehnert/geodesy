//! Orthographic and Local Orthographic
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const EPS10: f64 = 1e-10;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    NorthPole,
    SouthPole,
    Equatorial,
    Oblique,
}

fn adjlon(lon: f64) -> f64 {
    (lon + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
}

fn ortho_s_inverse_normalized(op: &Op, x: f64, y: f64) -> Option<(f64, f64)> {
    let frame = ProjectionFrame::from_params(&op.params);
    let mode = mode(op);
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(1.0);

    let rh = x.hypot(y);
    let mut sinc = rh;
    if sinc > 1.0 {
        if sinc - 1.0 > EPS10 {
            return None;
        }
        sinc = 1.0;
    }
    let cosc = (1.0 - sinc * sinc).sqrt();
    if rh <= EPS10 {
        return Some((0.0, frame.lat_0));
    }

    let (mut lam, phi) = match mode {
        Mode::NorthPole => {
            let yy = -y;
            (x.atan2(yy), sinc.acos())
        }
        Mode::SouthPole => (x.atan2(y), -sinc.acos()),
        Mode::Equatorial => {
            let mut arg = y * sinc / rh;
            arg = arg.clamp(-1.0, 1.0);
            let phi = arg.asin();
            let yy = cosc * rh;
            let xx = x * sinc;
            let lam = if yy == 0.0 {
                if xx == 0.0 {
                    0.0
                } else if xx < 0.0 {
                    -std::f64::consts::FRAC_PI_2
                } else {
                    std::f64::consts::FRAC_PI_2
                }
            } else {
                xx.atan2(yy)
            };
            (lam, phi)
        }
        Mode::Oblique => {
            let mut arg = cosc * sinph0 + y * sinc * cosph0 / rh;
            arg = arg.clamp(-1.0, 1.0);
            let phi = arg.asin();
            let yy = (cosc - sinph0 * arg) * rh;
            let xx = x * sinc * cosph0;
            let lam = if yy == 0.0 {
                if xx == 0.0 {
                    0.0
                } else if xx < 0.0 {
                    -std::f64::consts::FRAC_PI_2
                } else {
                    std::f64::consts::FRAC_PI_2
                }
            } else {
                xx.atan2(yy)
            };
            (lam, phi)
        }
    };
    lam = adjlon(lam);
    Some((lam, phi))
}

fn mode(op: &Op) -> Mode {
    match op.params.real("mode").unwrap_or(2.0) as i32 {
        0 => Mode::NorthPole,
        1 => Mode::SouthPole,
        2 => Mode::Equatorial,
        _ => Mode::Oblique,
    }
}

fn unrotate(op: &Op, x: f64, y: f64) -> (f64, f64) {
    let alpha = op.params.real("alpha").unwrap_or(0.0);
    let k_0 = op.params.k(0);
    let cosalpha = alpha.cos();
    let sinalpha = alpha.sin();
    (
        (cosalpha * x + sinalpha * y) / k_0,
        (-sinalpha * x + cosalpha * y) / k_0,
    )
}

fn rotate(op: &Op, x: f64, y: f64) -> (f64, f64) {
    let alpha = op.params.real("alpha").unwrap_or(0.0);
    let k_0 = op.params.k(0);
    let cosalpha = alpha.cos();
    let sinalpha = alpha.sin();
    (
        (x * cosalpha - y * sinalpha) * k_0,
        (x * sinalpha + y * cosalpha) * k_0,
    )
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let es = ellps.eccentricity_squared();
    let frame = ProjectionFrame::from_params(&op.params);
    let phi0 = frame.lat_0;
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(1.0);
    let nu0 = op.params.real("nu0").unwrap_or(1.0);
    let mode = mode(op);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = adjlon(frame.remove_central_meridian_raw(lon));
        let cosphi = lat.cos();
        let sinphi = lat.sin();
        let coslam = lam.cos();
        let sinlam = lam.sin();

        let (xp, yp) = if spherical {
            let yp = match mode {
                Mode::Equatorial => {
                    if cosphi * coslam < -EPS10 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    sinphi
                }
                Mode::Oblique => {
                    if sinph0 * sinphi + cosph0 * cosphi * coslam < -EPS10 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    cosph0 * sinphi - sinph0 * cosphi * coslam
                }
                Mode::NorthPole => {
                    if (lat - phi0).abs() - std::f64::consts::FRAC_PI_2 > EPS10 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    -cosphi * coslam
                }
                Mode::SouthPole => {
                    if (lat - phi0).abs() - std::f64::consts::FRAC_PI_2 > EPS10 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    cosphi * coslam
                }
            };
            (cosphi * sinlam, yp)
        } else {
            if sinph0 * sinphi + cosph0 * cosphi * coslam < -EPS10 {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let nu = 1.0 / (1.0 - es * sinphi * sinphi).sqrt();
            let xp = nu * cosphi * sinlam;
            let yp = nu * (sinphi * cosph0 - cosphi * sinph0 * coslam)
                + es * (nu0 * sinph0 - nu * sinphi) * cosph0;
            (xp, yp)
        };

        let (xr, yr) = rotate(op, xp, yp);
        let (x, y) = frame.apply_false_origin(frame.a * xr, frame.a * yr);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let b = ellps.semiminor_axis();
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let frame = ProjectionFrame::from_params(&op.params);
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(1.0);
    let y_shift = op.params.real("y_shift").unwrap_or(0.0);
    let y_scale = op.params.real("y_scale").unwrap_or(1.0);
    let mode = mode(op);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x_raw, y_raw);
        let (x, y) = unrotate(op, x_local / frame.a, y_local / frame.a);

        let (lam, phi) = if spherical {
            match ortho_s_inverse_normalized(op, x, y) {
                Some(v) => v,
                None => {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
            }
        } else {
            match mode {
                Mode::NorthPole | Mode::SouthPole => {
                    let rh2 = x * x + y * y;
                    if rh2 >= 1.0 - 1e-15 {
                        if rh2 - 1.0 > EPS10 {
                            operands.set_coord(i, &Coor4D::nan());
                            continue;
                        }
                        let lam = x.atan2(y * if mode == Mode::NorthPole { -1.0 } else { 1.0 });
                        (lam, 0.0)
                    } else {
                        let phi = (rh2 * one_es / (1.0 - es * rh2)).sqrt().acos()
                            * if mode == Mode::NorthPole { 1.0 } else { -1.0 };
                        let lam = x.atan2(y * if mode == Mode::NorthPole { -1.0 } else { 1.0 });
                        (lam, phi)
                    }
                }
                Mode::Equatorial => {
                    if x * x + (y * (frame.a / b)).powi(2) > 1.0 + 1e-11 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    let sinphi2 = if y == 0.0 {
                        0.0
                    } else {
                        1.0 / (((1.0 - es) / y).powi(2) + es)
                    };
                    if sinphi2 > 1.0 - 1e-11 {
                        (
                            0.0,
                            std::f64::consts::FRAC_PI_2 * if y > 0.0 { 1.0 } else { -1.0 },
                        )
                    } else {
                        let phi = sinphi2.sqrt().asin() * if y > 0.0 { 1.0 } else { -1.0 };
                        let sinlam = x * ((1.0 - es * sinphi2) / (1.0 - sinphi2)).sqrt();
                        let lam = if sinlam.abs() > 1.0 - 1e-15 {
                            std::f64::consts::FRAC_PI_2 * if x > 0.0 { 1.0 } else { -1.0 }
                        } else {
                            sinlam.asin()
                        };
                        (lam, phi)
                    }
                }
                Mode::Oblique => {
                    let yr = (y - y_shift) / y_scale;
                    if x * x + yr * yr > 1.0 + 1e-11 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    let mut lp = match ortho_s_inverse_normalized(op, x, yr) {
                        Some(v) => v,
                        None => {
                            operands.set_coord(i, &Coor4D::nan());
                            continue;
                        }
                    };

                    let mut converged = false;
                    for _ in 0..20 {
                        let cosphi = lp.1.cos();
                        let sinphi = lp.1.sin();
                        let coslam = lp.0.cos();
                        let sinlam = lp.0.sin();
                        let one_minus_es_sinphi2 = 1.0 - es * sinphi * sinphi;
                        let nu = 1.0 / one_minus_es_sinphi2.sqrt();
                        let x_new = nu * cosphi * sinlam;
                        let y_new = nu * (sinphi * cosph0 - cosphi * sinph0 * coslam)
                            + es * (op.params.real("nu0").unwrap_or(1.0) * sinph0 - nu * sinphi)
                                * cosph0;
                        let rho = one_es * nu / one_minus_es_sinphi2;
                        let j11 = -rho * sinphi * sinlam;
                        let j12 = nu * cosphi * coslam;
                        let j21 = rho * (cosphi * cosph0 + sinphi * sinph0 * coslam);
                        let j22 = nu * sinph0 * cosphi * sinlam;
                        let det = j11 * j22 - j12 * j21;
                        let dx = x - x_new;
                        let dy = y - y_new;
                        let dphi = (j22 * dx - j12 * dy) / det;
                        let dlam = (-j21 * dx + j11 * dy) / det;
                        lp.1 += dphi;
                        if lp.1 > std::f64::consts::FRAC_PI_2 {
                            lp.1 =
                                std::f64::consts::FRAC_PI_2 - (lp.1 - std::f64::consts::FRAC_PI_2);
                            lp.0 = adjlon(lp.0 + std::f64::consts::PI);
                        } else if lp.1 < -std::f64::consts::FRAC_PI_2 {
                            lp.1 = -std::f64::consts::FRAC_PI_2
                                + (-std::f64::consts::FRAC_PI_2 - lp.1);
                            lp.0 = adjlon(lp.0 + std::f64::consts::PI);
                        }
                        lp.0 += dlam;
                        if dphi.abs() < 1e-12 && dlam.abs() < 1e-12 {
                            converged = true;
                            break;
                        }
                    }
                    if !converged {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    lp
                }
            }
        };

        operands.set_xy(i, adjlon(frame.lon_0 + lam), phi);
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
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "alpha", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let lat_0 = params.lat(0);
    if lat_0.abs() > std::f64::consts::FRAC_PI_2 + EPS10 {
        return Err(Error::BadParam("lat_0".to_string(), def.clone()));
    }
    params
        .real
        .insert("alpha", params.real("alpha")?.to_radians());

    let sinph0 = lat_0.sin();
    let cosph0 = lat_0.cos();
    params.real.insert("sinph0", sinph0);
    params.real.insert("cosph0", cosph0);
    let mode = if (lat_0.abs() - std::f64::consts::FRAC_PI_2).abs() <= EPS10 {
        if lat_0 < 0.0 { 1.0 } else { 0.0 }
    } else if lat_0.abs() > EPS10 {
        3.0
    } else {
        2.0
    };
    params.real.insert("mode", mode);

    let ellps = params.ellps(0);
    if ellps.flattening() == 0.0 {
        params.boolean.insert("spherical");
    } else {
        let es = ellps.eccentricity_squared();
        let nu0 = 1.0 / (1.0 - es * sinph0 * sinph0).sqrt();
        params.real.insert("nu0", nu0);
        params.real.insert("y_shift", es * nu0 * sinph0 * cosph0);
        params
            .real
            .insert("y_scale", 1.0 / (1.0 - es * cosph0 * cosph0).sqrt());
    }

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
    fn ortho_spherical_equatorial_matches_proj_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("ortho R=1 lat_0=0 lon_0=0")?;
        let geo = [Coor4D::geo(50.0, 10.0, 0., 0.)];
        let projected = [Coor4D::raw(0.1116, 0.7660, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0].x() - projected[0].x()).abs() < 1e-4);
        assert!((operands[0].y() - projected[0].y()).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn ortho_ellipsoidal_matches_proj_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("ortho ellps=WGS84 lat_0=55 lon_0=5")?;
        let geo = [Coor4D::geo(53.80939444444444, 2.12955, 0., 0.)];
        let projected = [Coor4D::raw(-189_011.711, -128_640.567, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-3);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn ortho_local_matches_proj_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("ortho lat_0=37.628969166666664 lon_0=-122.39394166666668 k_0=0.9999968 alpha=27.7927777777777 x_0=0 y_0=0 ellps=GRS80")?;
        let geo = [Coor4D::geo(37.62607694444444, -122.3846388888889, 0., 0.)];
        let projected = [Coor4D::raw(876.13676, 98.97406, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-4);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
