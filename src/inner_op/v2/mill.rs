//! Miller Cylindrical
use crate::authoring::*;
use crate::projection::ProjectionFrame;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let x = frame.a * frame.lon_delta(lon);
        let y = frame.a * 1.25 * (std::f64::consts::FRAC_PI_4 + 0.4 * lat).tan().ln();
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x_raw, y_raw);
        let x = x_local / frame.a;
        let y = y_local / frame.a;
        let lon = frame.lon_0 + x;
        let lat = 2.5 * ((0.8 * y).exp().atan() - std::f64::consts::FRAC_PI_4);
        operands.set_xy(i, lon, lat);
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
    let params = ParsedParameters::new(parameters, &GAMUT)?;
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
    fn mill_matches_proj_gie_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("mill a=6400000")?;
        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(
            223_402.144_255_274,
            111_704.701_754_394,
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
