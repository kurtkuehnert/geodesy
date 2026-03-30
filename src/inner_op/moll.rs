//! Mollweide
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const MAX_ITER: usize = 30;
const LOOP_TOL: f64 = 1e-7;

fn moll_constants() -> (f64, f64, f64) {
    let p = std::f64::consts::FRAC_PI_2;
    let p2 = p + p;
    let sp = p.sin();
    let r = (std::f64::consts::TAU * sp / (p2 + p2.sin())).sqrt();
    (2.0 * r / std::f64::consts::PI, r / sp, p2 + p2.sin())
}

fn aasin(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).asin()
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let (c_x, c_y, c_p) = moll_constants();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = frame.lon_delta(lon);
        let mut theta = lat;
        let k = c_p * lat.sin();
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let v = (theta + theta.sin() - k) / (1.0 + theta.cos());
            theta -= v;
            if v.abs() < LOOP_TOL {
                converged = true;
                break;
            }
        }
        if !converged {
            theta = if theta < 0.0 {
                -std::f64::consts::FRAC_PI_2
            } else {
                std::f64::consts::FRAC_PI_2
            };
        } else {
            theta *= 0.5;
        }

        let x = frame.a * c_x * lam * theta.cos();
        let y = frame.a * c_y * theta.sin();
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let (c_x, c_y, c_p) = moll_constants();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x_raw, y_raw);
        let x = x_local / frame.a;
        let y = y_local / frame.a;
        let theta = aasin(y / c_y);
        let cos_theta = theta.cos();
        if cos_theta.abs() < 1e-12 {
            operands.set_xy(i, frame.lon_0, aasin((theta + theta).sin() / c_p));
            successes += 1;
            continue;
        }
        let lam = x / (c_x * cos_theta);
        if lam.abs() > std::f64::consts::PI + 1e-10 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let phi = aasin(((theta + theta) + (theta + theta).sin()) / c_p);
        operands.set_xy(i, frame.lon_0 + lam, phi);
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
    fn moll_matches_proj_gie_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("moll a=6400000")?;
        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(
            201_113.698_641_813,
            124_066.283_433_860,
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
