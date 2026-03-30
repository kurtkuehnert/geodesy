//! Tunisia Mining Grid
use crate::authoring::*;

const PARIS_PM_GRAD: f64 = 2.596_921_3;
const GRAD_TO_RAD: f64 = std::f64::consts::PI / 200.0;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let phi_0 = op.params.real["lat_0_grad"];
    let lon_0 = op.params.real["lon_0_grad"];
    let x_0 = op.params.real["x_0_km"];
    let y_0 = op.params.real["y_0_km"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lat_grad = lat.to_degrees() / 0.9;
        let lon_grad = lon.to_degrees() / 0.9;
        let x_km = x_0 + (lon_grad - lon_0) / 0.012_185;
        let a = if lat_grad > phi_0 {
            0.010_015
        } else {
            0.010_02
        };
        let y_km = y_0 + (lat_grad - phi_0) / a;
        operands.set_xy(i, x_km * 1000.0, y_km * 1000.0);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let phi_0 = op.params.real["lat_0_grad"];
    let lon_0 = op.params.real["lon_0_grad"];
    let x_0 = op.params.real["x_0_km"];
    let y_0 = op.params.real["y_0_km"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let x_km = x / 1000.0;
        let y_km = y / 1000.0;
        let lat_grad = phi_0 + (y_km - y_0) * if y_km > y_0 { 0.010_015 } else { 0.010_02 };
        let lon_grad = lon_0 + (x_km - x_0) * 0.012_185;
        operands.set_xy(i, lon_grad * GRAD_TO_RAD, lat_grad * GRAD_TO_RAD);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Real { key: "lat_0", default: None },
    OpParameter::Real { key: "lon_0", default: None },
    OpParameter::Real { key: "x_0", default: None },
    OpParameter::Real { key: "y_0", default: None },
    OpParameter::Flag { key: "greenwich" },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let mut lon_0_grad = params.real("lon_0")? / 0.9;
    if params.boolean("greenwich") {
        lon_0_grad -= PARIS_PM_GRAD;
    }
    params
        .real
        .insert("lat_0_grad", params.real("lat_0")? / 0.9);
    params.real.insert("lon_0_grad", lon_0_grad);
    params.real.insert("x_0_km", params.real("x_0")? / 1000.0);
    params.real.insert("y_0_km", params.real("y_0")? / 1000.0);

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
    fn tunmg_matches_epsg_example() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("tunmg lat_0=32.93676 lon_0=7.051005 x_0=270000 y_0=360000")?;

        let original = [Coor4D::raw(302_000.0, 598_000.0, 0.0, 0.0)];
        let mut inverse = original;
        ctx.apply(op, Inv, &mut inverse)?;
        assert!((inverse[0][0].to_degrees() - 7.401_933).abs() < 1e-6);
        assert!((inverse[0][1].to_degrees() - 35.081_973).abs() < 1e-6);

        ctx.apply(op, Fwd, &mut inverse)?;
        assert!(inverse[0].hypot2(&original[0]) < 1e-8);
        Ok(())
    }
}
