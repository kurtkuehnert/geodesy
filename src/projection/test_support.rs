use crate::authoring::*;

const DEFAULT_PROJ_FORWARD_TOLERANCE: f64 = 1e-6;
const DEFAULT_PROJ_ROUNDTRIP_TOLERANCE: f64 = 1e-10;

pub(crate) fn assert_proj_match(
    definition: &str,
    input: Coor4D,
    expected: Coor4D,
) -> Result<(), Error> {
    let mut ctx = Minimal::default();
    let op = ctx.op(definition)?;
    let mut operands = [input];

    ctx.apply(op, Fwd, &mut operands)?;
    assert!(operands[0].hypot2(&expected) < DEFAULT_PROJ_FORWARD_TOLERANCE);

    ctx.apply(op, Inv, &mut operands)?;
    assert!(operands[0].hypot2(&input) < DEFAULT_PROJ_ROUNDTRIP_TOLERANCE);
    Ok(())
}

pub(crate) fn assert_op_err(definition: &str) {
    let mut ctx = Minimal::default();
    assert!(ctx.op(definition).is_err());
}
