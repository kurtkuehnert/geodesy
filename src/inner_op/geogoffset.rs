//! Geographic offset.
use crate::authoring::*;

const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let dlon = op.params.real("dlon").unwrap_or(0.0);
    let dlat = op.params.real("dlat").unwrap_or(0.0);
    let dh = op.params.real("dh").unwrap_or(0.0);

    for i in 0..operands.len() {
        let coord = operands.get_coord(i);
        operands.set_xyz(i, coord[0] + dlon, coord[1] + dlat, coord[2] + dh);
    }
    operands.len()
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let dlon = op.params.real("dlon").unwrap_or(0.0);
    let dlat = op.params.real("dlat").unwrap_or(0.0);
    let dh = op.params.real("dh").unwrap_or(0.0);

    for i in 0..operands.len() {
        let coord = operands.get_coord(i);
        operands.set_xyz(i, coord[0] - dlon, coord[1] - dlat, coord[2] - dh);
    }
    operands.len()
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 4] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Real { key: "dlat", default: Some(0_f64) },
    OpParameter::Real { key: "dlon", default: Some(0_f64) },
    OpParameter::Real { key: "dh", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let mut op = Op::basic(parameters, InnerOp(fwd), Some(InnerOp(inv)), &GAMUT)?;
    op.params
        .real
        .insert("dlat", op.params.real("dlat")? * ARCSEC_TO_RAD);
    op.params
        .real
        .insert("dlon", op.params.real("dlon")? * ARCSEC_TO_RAD);
    Ok(op)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geogoffset_matches_proj_gie_examples() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("geogoffset dlon=3600 dlat=-3600 dh=3")?;

        let mut operands = [
            Coor4D::geo(20.0, 10.0, 30.0, 0.0),
            Coor4D::geo(20.0, 10.0, 0.0, 0.0),
        ];
        ctx.apply(op, Fwd, &mut operands)?;

        assert!((operands[0][0].to_degrees() - 11.0).abs() < 1e-12);
        assert!((operands[0][1].to_degrees() - 19.0).abs() < 1e-12);
        assert!((operands[0][2] - 33.0).abs() < 1e-12);
        assert!((operands[1][0].to_degrees() - 11.0).abs() < 1e-12);
        assert!((operands[1][1].to_degrees() - 19.0).abs() < 1e-12);

        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][0].to_degrees() - 10.0).abs() < 1e-12);
        assert!((operands[0][1].to_degrees() - 20.0).abs() < 1e-12);
        assert!((operands[0][2] - 30.0).abs() < 1e-12);
        Ok(())
    }
}
