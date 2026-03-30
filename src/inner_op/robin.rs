//! Robinson
use crate::authoring::*;

#[derive(Clone, Copy)]
struct Coefs {
    c0: f32,
    c1: f32,
    c2: f32,
    c3: f32,
}

const X: [Coefs; 19] = [
    Coefs {
        c0: 1.0,
        c1: 2.2199e-17,
        c2: -7.15515e-05,
        c3: 3.1103e-06,
    },
    Coefs {
        c0: 0.9986,
        c1: -0.000482243,
        c2: -2.4897e-05,
        c3: -1.3309e-06,
    },
    Coefs {
        c0: 0.9954,
        c1: -0.00083103,
        c2: -4.48605e-05,
        c3: -9.86701e-07,
    },
    Coefs {
        c0: 0.99,
        c1: -0.00135364,
        c2: -5.9661e-05,
        c3: 3.6777e-06,
    },
    Coefs {
        c0: 0.9822,
        c1: -0.00167442,
        c2: -4.49547e-06,
        c3: -5.72411e-06,
    },
    Coefs {
        c0: 0.973,
        c1: -0.00214868,
        c2: -9.03571e-05,
        c3: 1.8736e-08,
    },
    Coefs {
        c0: 0.96,
        c1: -0.00305085,
        c2: -9.00761e-05,
        c3: 1.64917e-06,
    },
    Coefs {
        c0: 0.9427,
        c1: -0.00382792,
        c2: -6.53386e-05,
        c3: -2.6154e-06,
    },
    Coefs {
        c0: 0.9216,
        c1: -0.00467746,
        c2: -0.00010457,
        c3: 4.81243e-06,
    },
    Coefs {
        c0: 0.8962,
        c1: -0.00536223,
        c2: -3.23831e-05,
        c3: -5.43432e-06,
    },
    Coefs {
        c0: 0.8679,
        c1: -0.00609363,
        c2: -0.000113898,
        c3: 3.32484e-06,
    },
    Coefs {
        c0: 0.835,
        c1: -0.00698325,
        c2: -6.40253e-05,
        c3: 9.34959e-07,
    },
    Coefs {
        c0: 0.7986,
        c1: -0.00755338,
        c2: -5.00009e-05,
        c3: 9.35324e-07,
    },
    Coefs {
        c0: 0.7597,
        c1: -0.00798324,
        c2: -3.5971e-05,
        c3: -2.27626e-06,
    },
    Coefs {
        c0: 0.7186,
        c1: -0.00851367,
        c2: -7.01149e-05,
        c3: -8.6303e-06,
    },
    Coefs {
        c0: 0.6732,
        c1: -0.00986209,
        c2: -0.000199569,
        c3: 1.91974e-05,
    },
    Coefs {
        c0: 0.6213,
        c1: -0.010418,
        c2: 8.83923e-05,
        c3: 6.24051e-06,
    },
    Coefs {
        c0: 0.5722,
        c1: -0.00906601,
        c2: 0.000182,
        c3: 6.24051e-06,
    },
    Coefs {
        c0: 0.5322,
        c1: -0.00677797,
        c2: 0.000275608,
        c3: 6.24051e-06,
    },
];

const Y: [Coefs; 19] = [
    Coefs {
        c0: -5.20417e-18,
        c1: 0.0124,
        c2: 1.21431e-18,
        c3: -8.45284e-11,
    },
    Coefs {
        c0: 0.062,
        c1: 0.0124,
        c2: -1.26793e-09,
        c3: 4.22642e-10,
    },
    Coefs {
        c0: 0.124,
        c1: 0.0124,
        c2: 5.07171e-09,
        c3: -1.60604e-09,
    },
    Coefs {
        c0: 0.186,
        c1: 0.0123999,
        c2: -1.90189e-08,
        c3: 6.00152e-09,
    },
    Coefs {
        c0: 0.248,
        c1: 0.0124002,
        c2: 7.10039e-08,
        c3: -2.24e-08,
    },
    Coefs {
        c0: 0.31,
        c1: 0.0123992,
        c2: -2.64997e-07,
        c3: 8.35986e-08,
    },
    Coefs {
        c0: 0.372,
        c1: 0.0124029,
        c2: 9.88983e-07,
        c3: -3.11994e-07,
    },
    Coefs {
        c0: 0.434,
        c1: 0.0123893,
        c2: -3.69093e-06,
        c3: -4.35621e-07,
    },
    Coefs {
        c0: 0.4958,
        c1: 0.0123198,
        c2: -1.02252e-05,
        c3: -3.45523e-07,
    },
    Coefs {
        c0: 0.5571,
        c1: 0.0121916,
        c2: -1.54081e-05,
        c3: -5.82288e-07,
    },
    Coefs {
        c0: 0.6176,
        c1: 0.0119938,
        c2: -2.41424e-05,
        c3: -5.25327e-07,
    },
    Coefs {
        c0: 0.6769,
        c1: 0.011713,
        c2: -3.20223e-05,
        c3: -5.16405e-07,
    },
    Coefs {
        c0: 0.7346,
        c1: 0.0113541,
        c2: -3.97684e-05,
        c3: -6.09052e-07,
    },
    Coefs {
        c0: 0.7903,
        c1: 0.0109107,
        c2: -4.89042e-05,
        c3: -1.04739e-06,
    },
    Coefs {
        c0: 0.8435,
        c1: 0.0103431,
        c2: -6.4615e-05,
        c3: -1.40374e-09,
    },
    Coefs {
        c0: 0.8936,
        c1: 0.00969686,
        c2: -6.4636e-05,
        c3: -8.547e-06,
    },
    Coefs {
        c0: 0.9394,
        c1: 0.00840947,
        c2: -0.000192841,
        c3: -4.2106e-06,
    },
    Coefs {
        c0: 0.9761,
        c1: 0.00616527,
        c2: -0.000256,
        c3: -4.2106e-06,
    },
    Coefs {
        c0: 1.0,
        c1: 0.00328947,
        c2: -0.000319159,
        c3: -4.2106e-06,
    },
];

const FXC: f64 = 0.8487;
const FYC: f64 = 1.3523;
const C1: f64 = 11.459_155_902_616_464;
const RC1: f64 = 0.087_266_462_599_716_47;
const NODES: usize = 18;
const ONEEPS: f64 = 1.000001;
const EPS: f64 = 1e-10;
const MAX_ITER: usize = 100;

fn v(c: Coefs, z: f64) -> f64 {
    let c0 = c.c0 as f64;
    let c1 = c.c1 as f64;
    let c2 = c.c2 as f64;
    let c3 = c.c3 as f64;
    c0 + z * (c1 + z * (c2 + z * c3))
}

fn dv(c: Coefs, z: f64) -> f64 {
    let c1 = c.c1 as f64;
    let c2 = c.c2 as f64;
    let c3 = c.c3 as f64;
    c1 + 2.0 * z * c2 + z * z * 3.0 * c3
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let mut dphi = lat.abs();
        if dphi.is_nan() {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let mut idx = (dphi * C1 + 1e-15).floor() as isize;
        if idx < 0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        if idx as usize >= NODES {
            idx = NODES as isize;
        }
        dphi = (dphi - RC1 * idx as f64).to_degrees();
        let x = x_0 + a * v(X[idx as usize], dphi) * FXC * lam;
        let mut y = a * v(Y[idx as usize], dphi) * FYC;
        if lat < 0.0 {
            y = -y;
        }
        operands.set_xy(i, x, y_0 + y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let x = (x_raw - x_0) / a;
        let y = (y_raw - y_0) / a;

        let mut lam = x / FXC;
        let mut phi = (y / FYC).abs();
        if phi >= 1.0 {
            if phi > ONEEPS {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            phi = if y < 0.0 {
                -std::f64::consts::FRAC_PI_2
            } else {
                std::f64::consts::FRAC_PI_2
            };
            lam /= X[NODES].c0 as f64;
            operands.set_xy(i, lon_0 + lam, phi);
            successes += 1;
            continue;
        }

        let mut idx = (phi * NODES as f64).floor() as usize;
        if idx >= NODES {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        loop {
            if (Y[idx].c0 as f64) > phi {
                idx -= 1;
            } else if (Y[idx + 1].c0 as f64) <= phi {
                idx += 1;
            } else {
                break;
            }
        }
        let tcoefs = Y[idx];
        let mut t = 5.0 * (phi - tcoefs.c0 as f64) / ((Y[idx + 1].c0 - tcoefs.c0) as f64);
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let dt = (v(tcoefs, t) - phi) / dv(tcoefs, t);
            t -= dt;
            if dt.abs() < EPS {
                converged = true;
                break;
            }
        }
        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        phi = (5.0 * idx as f64 + t).to_radians();
        if y < 0.0 {
            phi = -phi;
        }
        lam /= v(X[idx], t);
        if lam.abs() > std::f64::consts::PI {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        operands.set_xy(i, lon_0 + lam, phi);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    params.real.insert("lon_0", params.lon(0).to_radians());

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
    fn robin_matches_proj_gie_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("robin a=6400000")?;
        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(
            189_588.423_282_508,
            107_318.530_350_703,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-9);
        Ok(())
    }
}
