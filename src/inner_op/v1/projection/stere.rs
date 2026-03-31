//! Stereographic projection, following PROJ's mode-specific implementation.
use crate::authoring::*;
use crate::projection::{
    ConformalLatitude, ProjectionAspect, ProjectionFrame, spherical_inverse_equatorial,
    spherical_inverse_oblique,
};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const DENOMINATOR_TOLERANCE: f64 = 1e-10;
const POLAR_TOLERANCE: f64 = 1e-8;
const POLAR_SOURCE_TOLERANCE: f64 = 1e-15;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 9] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "variant_c" },
    OpParameter::Text { key: "ellps",  default: Some("GRS80") },
    OpParameter::Real { key: "lat_0",  default: Some(0_f64) },
    OpParameter::Real { key: "lat_ts", default: Some(90_f64) },
    OpParameter::Real { key: "lon_0",  default: Some(0_f64) },
    OpParameter::Real { key: "k_0",    default: Some(1_f64) },
    OpParameter::Real { key: "x_0",    default: Some(0_f64) },
    OpParameter::Real { key: "y_0",    default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug)]
struct StereState {
    frame: ProjectionFrame,
    conformal: ConformalLatitude,
    aspect: ProjectionAspect,
    oblique: Option<(f64, f64)>,
    akm1: f64,
}

impl StereState {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let def = &params.name;
        let lat_0 = params.lat(0);
        let mut lat_ts = params.real("lat_ts").unwrap_or(FRAC_PI_2);
        if lat_ts.abs() > FRAC_PI_2 + DENOMINATOR_TOLERANCE {
            return Err(Error::BadParam("lat_ts".to_string(), def.clone()));
        }

        let aspect = ProjectionAspect::classify(lat_0, POLAR_TOLERANCE);
        if aspect.is_polar() {
            lat_ts = lat_ts.copysign(lat_0);
        }

        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let conformal = ConformalLatitude::new(ellps);
        let e = conformal.eccentricity();
        let k_0 = params.k(0);
        let variant_c = params.boolean("variant_c");

        if variant_c && !aspect.is_polar() {
            return Err(Error::Unsupported(
                "stere variant_c is only supported for polar stereographic aspects".into(),
            ));
        }

        let oblique = if aspect.is_oblique() {
            if e != 0.0 {
                let xang = conformal.reduced(lat_0);
                Some((xang.sin(), xang.cos()))
            } else {
                Some((lat_0.sin(), lat_0.cos()))
            }
        } else {
            None
        };

        let akm1 = if e != 0.0 {
            if aspect.is_polar() {
                let lat_ts_abs = lat_ts.abs();
                if (lat_ts_abs - FRAC_PI_2).abs() < POLAR_TOLERANCE {
                    let num = 2.0 * k_0;
                    let den = ((1.0 + e).powf(1.0 + e) * (1.0 - e).powf(1.0 - e)).sqrt();
                    a * num / den
                } else {
                    let sin_ts = lat_ts_abs.sin();
                    let factor = lat_ts_abs.cos() / conformal.ts(lat_ts_abs);
                    a * k_0 * factor / (1.0 - (e * sin_ts).powi(2)).sqrt()
                }
            } else {
                let te = e * lat_0.sin();
                2.0 * a * k_0 * lat_0.cos() / (1.0 - te * te).sqrt()
            }
        } else {
            if !aspect.is_polar() {
                2.0 * a * k_0
            } else if (lat_ts.abs() - FRAC_PI_2).abs() >= DENOMINATOR_TOLERANCE {
                a * lat_ts.abs().cos() / (FRAC_PI_4 - 0.5 * lat_ts.abs()).tan()
            } else {
                2.0 * a * k_0
            }
        };

        let mut frame = ProjectionFrame::from_params(params);
        if variant_c {
            let lat_ts_abs = lat_ts.abs();
            let sin_ts = lat_ts_abs.sin();
            let cos_ts = lat_ts_abs.cos();
            let rho_f = a * cos_ts / (1.0 - (e * sin_ts).powi(2)).sqrt();
            frame.y_0 += if aspect.is_south_polar() {
                -rho_f
            } else {
                rho_f
            };
        }

        Ok(Self {
            frame,
            conformal,
            aspect,
            oblique,
            akm1,
        })
    }
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<StereState>();
    let spherical = state.conformal.spherical();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = state.frame.remove_central_meridian_raw(lon);
        let sinlam = lam.sin();
        let mut coslam = lam.cos();
        let (x, y) = if spherical {
            let sinphi = lat.sin();
            let cosphi = lat.cos();
            if state.aspect.is_equatorial() {
                let denom = 1.0 + cosphi * coslam;
                if denom <= DENOMINATOR_TOLERANCE {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let a = state.akm1 / denom;
                (a * cosphi * sinlam, a * sinphi)
            } else if let Some((sin_x1, cos_x1)) = state.oblique {
                let denom = 1.0 + sin_x1 * sinphi + cos_x1 * cosphi * coslam;
                if denom <= DENOMINATOR_TOLERANCE {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let a = state.akm1 / denom;
                (
                    a * cosphi * sinlam,
                    a * (cos_x1 * sinphi - sin_x1 * cosphi * coslam),
                )
            } else {
                let mut phi = lat;
                if state.aspect.is_north_polar() {
                    coslam = -coslam;
                    phi = -phi;
                }
                if (phi - FRAC_PI_2).abs() < POLAR_TOLERANCE {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let y = state.akm1 * (FRAC_PI_4 + 0.5 * phi).tan();
                (sinlam * y, y * coslam)
            }
        } else {
            let mut sin_x = 0.0;
            let mut cos_x = 0.0;
            let a;
            let mut phi = lat;
            let mut sinphi = phi.sin();

            if !state.aspect.is_polar() {
                let xang = state.conformal.reduced(phi);
                sin_x = xang.sin();
                cos_x = xang.cos();
            }

            let (mut x, y) = if let Some((sin_x1, cos_x1)) = state.oblique {
                let denom = cos_x1 * (1.0 + sin_x1 * sin_x + cos_x1 * cos_x * coslam);
                if denom == 0.0 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                a = state.akm1 / denom;
                (a * cos_x, a * (cos_x1 * sin_x - sin_x1 * cos_x * coslam))
            } else if state.aspect.is_equatorial() {
                let denom = 1.0 + cos_x * coslam;
                if denom == 0.0 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                a = state.akm1 / denom;
                (a * cos_x, a * sin_x)
            } else {
                if state.aspect.is_south_polar() {
                    phi = -phi;
                    coslam = -coslam;
                    sinphi = -sinphi;
                }
                let x = if (phi - FRAC_PI_2).abs() < POLAR_SOURCE_TOLERANCE {
                    0.0
                } else {
                    state.akm1 * state.conformal.ts(phi)
                };
                (x, -x * coslam)
            };
            x *= sinlam;
            (x, y)
        };

        let (x, y) = state.frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<StereState>();
    let spherical = state.conformal.spherical();
    let lat_0 = state.frame.lat_0;

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (mut x, mut y) = operands.xy(i);
        let (xf, yf) = state.frame.remove_false_origin(x, y);
        x = xf;
        y = yf;
        let rho = x.hypot(y);

        let (lon, lat) = if spherical {
            let c = 2.0 * (rho / state.akm1).atan();
            if state.aspect.is_equatorial() {
                spherical_inverse_equatorial(x, y, rho, c, DENOMINATOR_TOLERANCE, 0.0)
            } else if let Some((sin_x1, cos_x1)) = state.oblique {
                spherical_inverse_oblique(
                    x,
                    y,
                    rho,
                    c,
                    DENOMINATOR_TOLERANCE,
                    lat_0,
                    sin_x1,
                    cos_x1,
                )
            } else {
                let cosc = c.cos();
                if state.aspect.is_north_polar() {
                    y = -y;
                }
                let lat = if rho.abs() <= DENOMINATOR_TOLERANCE {
                    lat_0
                } else if state.aspect.is_south_polar() {
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
            let (chi, x, y) = if let Some((sin_x1, cos_x1)) = state.oblique {
                let tp = 2.0 * (rho * cos_x1).atan2(state.akm1);
                let cosphi = tp.cos();
                let sinphi = tp.sin();
                let chi = if rho == 0.0 {
                    (cosphi * sin_x1).asin()
                } else {
                    (cosphi * sin_x1 + y * sinphi * cos_x1 / rho).asin()
                };
                (chi, x * sinphi, rho * cos_x1 * cosphi - y * sin_x1 * sinphi)
            } else if state.aspect.is_equatorial() {
                let tp = 2.0 * rho.atan2(state.akm1);
                let cosphi = tp.cos();
                let sinphi = tp.sin();
                let chi = if rho == 0.0 {
                    0.0
                } else {
                    (y * sinphi / rho).asin()
                };
                (chi, x * sinphi, rho * cosphi)
            } else {
                let mut y = y;
                if state.aspect.is_north_polar() {
                    y = -y;
                }
                (FRAC_PI_2 - 2.0 * (-rho / state.akm1).atan(), x, y)
            };

            let lat = if state.aspect.is_polar() {
                let lat = state.conformal.geographic(std::f64::consts::PI - chi);
                if state.aspect.is_south_polar() {
                    -lat
                } else {
                    lat
                }
            } else {
                state.conformal.geographic(chi)
            };
            let lon = if x == 0.0 && y == 0.0 {
                0.0
            } else {
                x.atan2(y)
            };
            (lon, lat)
        };

        operands.set_xy(i, state.frame.lon_0 + lon, lat);
        successes += 1;
    }

    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    let state = StereState::new(&params)?;
    Ok(Op::with_state(descriptor, params, state))
}

pub fn ups(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    let tokens: Vec<&str> = parameters.instantiated_as.split_whitespace().collect();
    let south = tokens.contains(&"south");
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
    fn variant_c_matches_epsg_example() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "stere variant_c lat_0=-90 lat_ts=-67 lon_0=140 x_0=300000 y_0=200000 ellps=intl",
        )?;

        let geo = [Coor4D::geo(-66.60522777777778, 140.0714, 0., 0.)];
        let projected = [Coor4D::raw(303_169.522, 244_055.721, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(
            operands[0].hypot2(&projected[0]) < 0.02,
            "got {:?}, expected {:?}",
            operands[0],
            projected[0]
        );
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-8);
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
