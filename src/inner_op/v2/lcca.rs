//! Lambert Conformal Conic Alternative
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const DEL_TOL: f64 = 1e-12;
const MAX_ITER: usize = 10;

fn fs(s: f64, c: f64) -> f64 {
    s * (1.0 + s * s * c)
}

fn fsp(s: f64, c: f64) -> f64 {
    1.0 + 3.0 * s * s * c
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let r0 = op.params.real["r0"];
    let l = op.params.real["l"];
    let m0 = op.params.real["m0"];
    let c = op.params.real["c"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let s = ellps.meridian_latitude_to_distance(lat) / frame.a - m0;
        let dr = fs(s, c);
        let r = r0 - dr;
        let theta = frame.remove_central_meridian_raw(lon) * l;
        let (x, y) = frame.apply_false_origin(
            frame.a * frame.k_0 * r * theta.sin(),
            frame.a * frame.k_0 * (r0 - r * theta.cos()),
        );
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let r0 = op.params.real["r0"];
    let l = op.params.real["l"];
    let m0 = op.params.real["m0"];
    let c = op.params.real["c"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = frame.remove_false_origin(operands.xy(i).0, operands.xy(i).1);
        let x = x / (frame.a * frame.k_0);
        let y = y / (frame.a * frame.k_0);
        let theta = x.atan2(r0 - y);
        let dr = y - x * (0.5 * theta).tan();
        let lon = theta / l + frame.lon_0;
        let mut s = dr;
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let dif = (fs(s, c) - dr) / fsp(s, c);
            s -= dif;
            if dif.abs() < DEL_TOL {
                converged = true;
                break;
            }
        }
        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let lat = ellps.meridian_distance_to_latitude((s + m0) * frame.a);
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 9] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: None },
    OpParameter::Real { key: "lat_1", default: Some(0_f64) },
    OpParameter::Real { key: "lat_2", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let phi0 = params.lat(0);
    if phi0.is_nan() {
        return Err(Error::MissingParam("lat_0".to_string()));
    }
    if phi0 == 0.0 {
        return Err(Error::General(
            "Lcca: Invalid value for lat_0: it should be different from 0.",
        ));
    }

    let ellps = params.ellps(0);
    let es = ellps.eccentricity_squared();
    let sinphi0 = phi0.sin();
    let s2p0 = sinphi0 * sinphi0;
    let r = 1.0 / (1.0 - es * s2p0);
    let n0 = r.sqrt();
    let r0m = r * (1.0 - es) * n0;
    let tan0 = phi0.tan();

    params.real.insert("l", sinphi0);
    params.real.insert(
        "m0",
        ellps.meridian_latitude_to_distance(phi0) / ellps.semimajor_axis(),
    );
    params.real.insert("r0", n0 / tan0);
    params.real.insert("c", 1.0 / (6.0 * r0m * n0));

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
    fn lcca_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("lcca ellps=GRS80 lat_0=1 lat_1=0.5 lat_2=2")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(222_605.285_770_237, 67.806_007_272, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        Ok(())
    }
}
