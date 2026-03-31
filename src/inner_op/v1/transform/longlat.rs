//! Geographic longitude/latitude with optional prime meridian offset.
use crate::authoring::*;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 3] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
];

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let lon_0 = op.params.lon(0);
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        operands.set_xy(i, lon - lon_0, lat);
    }
    operands.len()
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let lon_0 = op.params.lon(0);
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        operands.set_xy(i, lon + lon_0, lat);
    }
    operands.len()
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let mut op = Op::basic(parameters, InnerOp(fwd), Some(InnerOp(inv)), &GAMUT)?;
    let lon_0 = op.params.lon(0);
    op.params.real.insert("lon_0", lon_0);
    Ok(op)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn longlat_applies_prime_meridian_offset() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("longlat lon_0=0.00289027777777778")?;

        let mut operands = [Coor4D::geo(38., 125., 0., 0.)];
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][0].to_degrees() - 125.00289027777778).abs() < 1e-12);

        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][0].to_degrees() - 125.0).abs() < 1e-12);
        Ok(())
    }
}
