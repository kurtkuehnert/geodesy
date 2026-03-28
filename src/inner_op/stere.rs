//! Stereographic projection, following PROJ's mode-specific implementation.
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const EPS10: f64 = 1e-10;
const TOL: f64 = 1e-8;
const CONV: f64 = 1e-10;
const NITER: usize = 8;

fn ssfn(phi: f64, sinphi: f64, eccen: f64) -> f64 {
    (0.5 * (FRAC_PI_2 + phi)).tan()
        * ((1.0 - eccen * sinphi) / (1.0 + eccen * sinphi)).powf(0.5 * eccen)
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(akm1) = op.params.real("akm1") else {
        return 0;
    };

    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let spherical = e == 0.0;
    let north_polar = op.params.boolean("north_polar");
    let south_polar = op.params.boolean("south_polar");
    let equatorial = op.params.boolean("equatorial");
    let oblique = op.params.boolean("oblique");
    let sin_x1 = op.params.real("sin_x1").unwrap_or(0.0);
    let cos_x1 = op.params.real("cos_x1").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let sinlam = lam.sin();
        let mut coslam = lam.cos();
        let (x, y) = if spherical {
            let sinphi = lat.sin();
            let cosphi = lat.cos();
            if equatorial {
                let denom = 1.0 + cosphi * coslam;
                if denom <= EPS10 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let a = akm1 / denom;
                (a * cosphi * sinlam, a * sinphi)
            } else if oblique {
                let denom = 1.0 + sin_x1 * sinphi + cos_x1 * cosphi * coslam;
                if denom <= EPS10 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let a = akm1 / denom;
                (
                    a * cosphi * sinlam,
                    a * (cos_x1 * sinphi - sin_x1 * cosphi * coslam),
                )
            } else {
                let mut phi = lat;
                if north_polar {
                    coslam = -coslam;
                    phi = -phi;
                }
                if (phi - FRAC_PI_2).abs() < TOL {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let y = akm1 * (FRAC_PI_4 + 0.5 * phi).tan();
                (sinlam * y, y * coslam)
            }
        } else {
            let mut sin_x = 0.0;
            let mut cos_x = 0.0;
            let mut a = 0.0;
            let mut sinphi = lat.sin();

            if oblique || equatorial {
                let xang = 2.0 * ssfn(lat, sinphi, e).atan() - FRAC_PI_2;
                sin_x = xang.sin();
                cos_x = xang.cos();
            }

            let (mut x, y) = if oblique {
                let denom = cos_x1 * (1.0 + sin_x1 * sin_x + cos_x1 * cos_x * coslam);
                if denom == 0.0 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                a = akm1 / denom;
                (a * cos_x, a * (cos_x1 * sin_x - sin_x1 * cos_x * coslam))
            } else if equatorial {
                let denom = 1.0 + cos_x * coslam;
                if denom == 0.0 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                a = akm1 / denom;
                (a * cos_x, a * sin_x)
            } else {
                let mut phi = lat;
                if south_polar {
                    phi = -phi;
                    coslam = -coslam;
                    sinphi = -sinphi;
                }
                let x = if (phi - FRAC_PI_2).abs() < 1e-15 {
                    0.0
                } else {
                    akm1 * ancillary::ts((sinphi, phi.cos()), e)
                };
                (x, -x * coslam)
            };
            x *= sinlam;
            (x, y)
        };

        operands.set_xy(i, x_0 + x, y_0 + y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(akm1) = op.params.real("akm1") else {
        return 0;
    };

    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let spherical = e == 0.0;
    let north_polar = op.params.boolean("north_polar");
    let south_polar = op.params.boolean("south_polar");
    let equatorial = op.params.boolean("equatorial");
    let oblique = op.params.boolean("oblique");
    let sin_x1 = op.params.real("sin_x1").unwrap_or(0.0);
    let cos_x1 = op.params.real("cos_x1").unwrap_or(0.0);
    let lat_0 = op.params.real("lat_0").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (mut x, mut y) = operands.xy(i);
        x -= x_0;
        y -= y_0;
        let rho = x.hypot(y);

        let (lon, lat) = if spherical {
            let c = 2.0 * (rho / akm1).atan();
            let sinc = c.sin();
            let cosc = c.cos();
            if equatorial {
                let lat = if rho.abs() <= EPS10 {
                    0.0
                } else {
                    (y * sinc / rho).asin()
                };
                let lon = if cosc != 0.0 || x != 0.0 {
                    (x * sinc).atan2(cosc * rho)
                } else {
                    0.0
                };
                (lon, lat)
            } else if oblique {
                let lat = if rho.abs() <= EPS10 {
                    lat_0
                } else {
                    (cosc * sin_x1 + y * sinc * cos_x1 / rho).asin()
                };
                let c = cosc - sin_x1 * lat.sin();
                let lon = if c != 0.0 || x != 0.0 {
                    (x * sinc * cos_x1).atan2(c * rho)
                } else {
                    0.0
                };
                (lon, lat)
            } else {
                if north_polar {
                    y = -y;
                }
                let lat = if rho.abs() <= EPS10 {
                    lat_0
                } else if south_polar {
                    (-cosc).asin()
                } else {
                    cosc.asin()
                };
                let lon = if x == 0.0 && y == 0.0 {
                    0.0
                } else {
                    x.atan2(y)
                };
                (lon, lat)
            }
        } else {
            let (tp, mut phi_l, halfpi, halfe) = if oblique || equatorial {
                let tp = 2.0 * (rho * cos_x1).atan2(akm1);
                let cosphi = tp.cos();
                let sinphi = tp.sin();
                let phi_l = if rho == 0.0 {
                    (cosphi * sin_x1).asin()
                } else {
                    (cosphi * sin_x1 + y * sinphi * cos_x1 / rho).asin()
                };
                let tp = (0.5 * (FRAC_PI_2 + phi_l)).tan();
                x *= sinphi;
                y = rho * cos_x1 * cosphi - y * sin_x1 * sinphi;
                (tp, phi_l, FRAC_PI_2, 0.5 * e)
            } else {
                if north_polar {
                    y = -y;
                }
                let tp = -rho / akm1;
                let phi_l = FRAC_PI_2 - 2.0 * tp.atan();
                (tp, phi_l, -FRAC_PI_2, -0.5 * e)
            };

            let mut lat = f64::NAN;
            let mut converged = false;
            for _ in 0..NITER {
                let sinphi = e * phi_l.sin();
                lat = 2.0 * (tp * ((1.0 + sinphi) / (1.0 - sinphi)).powf(halfe)).atan() - halfpi;
                if (phi_l - lat).abs() < CONV {
                    converged = true;
                    break;
                }
                phi_l = lat;
            }
            if !converged {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let lat = if south_polar { -lat } else { lat };
            let lon = if x == 0.0 && y == 0.0 {
                0.0
            } else {
                x.atan2(y)
            };
            (lon, lat)
        };

        operands.set_xy(i, lon_0 + lon, lat);
        successes += 1;
    }

    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lat_ts", default: Some(90_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let lat_0 = params.lat(0).to_radians();
    let mut lat_ts = params.real("lat_ts").unwrap_or(90.0).to_radians();
    if lat_ts.abs() > FRAC_PI_2 + EPS10 {
        return Err(Error::BadParam("lat_ts".to_string(), def.clone()));
    }

    if (lat_0.abs() - FRAC_PI_2).abs() < EPS10 {
        if lat_0.is_sign_negative() {
            params.boolean.insert("south_polar");
            lat_ts = -lat_ts.abs();
        } else {
            params.boolean.insert("north_polar");
            lat_ts = lat_ts.abs();
        }
    } else if lat_0.abs() > EPS10 {
        params.boolean.insert("oblique");
    } else {
        params.boolean.insert("equatorial");
    }

    params.real.insert("lat_0", lat_0);
    params.real.insert("lat_ts", lat_ts);
    params.real.insert("lon_0", params.lon(0).to_radians());

    let ellps = params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    let k_0 = params.k(0);

    let akm1 = if e != 0.0 {
        if params.boolean("north_polar") || params.boolean("south_polar") {
            let lat_ts_abs = lat_ts.abs();
            if (lat_ts_abs - FRAC_PI_2).abs() < EPS10 {
                let num = 2.0 * k_0;
                let den = ((1.0 + e).powf(1.0 + e) * (1.0 - e).powf(1.0 - e)).sqrt();
                a * num / den
            } else {
                let sin_ts = lat_ts_abs.sin();
                let factor = lat_ts_abs.cos() / ancillary::ts((sin_ts, lat_ts_abs.cos()), e);
                a * k_0 * factor / (1.0 - (e * sin_ts).powi(2)).sqrt()
            }
        } else {
            let t = lat_0.sin();
            let xang = 2.0 * ssfn(lat_0, t, e).atan() - FRAC_PI_2;
            let te = e * t;
            params.real.insert("sin_x1", xang.sin());
            params.real.insert("cos_x1", xang.cos());
            2.0 * a * k_0 * lat_0.cos() / (1.0 - te * te).sqrt()
        }
    } else if params.boolean("oblique") {
        params.real.insert("sin_x1", lat_0.sin());
        params.real.insert("cos_x1", lat_0.cos());
        2.0 * a * k_0
    } else if params.boolean("equatorial") {
        2.0 * a * k_0
    } else if (lat_ts.abs() - FRAC_PI_2).abs() >= EPS10 {
        a * lat_ts.abs().cos() / (FRAC_PI_4 - 0.5 * lat_ts.abs()).tan()
    } else {
        2.0 * a * k_0
    };

    params.real.insert("akm1", akm1);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}

pub fn ups(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    let tokens: Vec<&str> = parameters.instantiated_as.split_whitespace().collect();
    let south = tokens.iter().any(|t| *t == "south");
    let has_named = tokens.iter().any(|t| t.starts_with("ellps="));
    let has_a = tokens.iter().any(|t| t.starts_with("a="));
    let has_b = tokens.iter().any(|t| t.starts_with("b="));
    let has_rf = tokens.iter().any(|t| t.starts_with("rf="));
    let has_f = tokens.iter().any(|t| t.starts_with("f="));
    let has_r = tokens.iter().any(|t| t.starts_with("R="));

    if has_r || (has_a && !has_b && !has_rf && !has_f && !has_named) {
        return Err(Error::Unsupported(
            "ups requires an ellipsoidal definition; spherical UPS is not supported".into(),
        ));
    }

    let mut rewritten = vec![
        "stere".to_string(),
        format!("lat_0={}", if south { -90 } else { 90 }),
        "lon_0=0".to_string(),
        "k_0=0.994".to_string(),
        "x_0=2000000".to_string(),
        "y_0=2000000".to_string(),
    ];

    for token in tokens.into_iter().skip(1) {
        if token == "south" {
            continue;
        }
        rewritten.push(token.to_string());
    }
    if !has_named && !has_a {
        rewritten.push("ellps=GRS80".to_string());
    }

    let def = rewritten.join(" ");
    let raw = RawParameters {
        invoked_as: def.clone(),
        instantiated_as: def,
        globals: parameters.globals.clone(),
        recursion_level: parameters.recursion_level,
    };
    new(&raw, ctx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ups_north_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere lat_0=90 lon_0=0 k_0=0.994 x_0=2000000 y_0=2000000 ellps=WGS84")?;

        let geo = [Coor4D::geo(85., 0., 0., 0.)];
        let projected = [Coor4D::raw(2_000_000.0, 1_444_542.608_617_322_5, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-8);
        Ok(())
    }

    #[test]
    fn south_polar_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere lat_0=-90 lat_ts=-71 lon_0=0 ellps=WGS84")?;

        let geo = [Coor4D::geo(-75., 30., 0., 0.)];
        let projected = [Coor4D::raw(
            819_391.619_181_387,
            1_419_227.915_756_798,
            0.,
            0.,
        )];

        let mut operands = projected;
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-6);
        Ok(())
    }

    #[test]
    fn custom_ellipsoid_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere lat_0=90 lat_ts=70 lon_0=-45 ellps=6378273,298.279411123064")?;

        let geo = [Coor4D::geo(80., 45., 0., 0.)];
        let projected = [Coor4D::raw(1_085_943.187_924_962, 0.0, 0., 0.)];

        let mut operands = projected;
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-11);
        Ok(())
    }

    #[test]
    fn stere_equatorial_matches_proj_sample() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere ellps=GRS80")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(
            222_644.854_550_117,
            110_610.883_474_174,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
