use crate::authoring::*;

pub(crate) fn assert_forward(
    definition: &str,
    input: Coor4D,
    expected: Coor4D,
    tolerance: f64,
) -> Result<(), Error> {
    let mut ctx = Minimal::default();
    let op = ctx.op(definition)?;
    let mut operands = [input];

    ctx.apply(op, Fwd, &mut operands)?;
    assert!(operands[0].hypot2(&expected) < tolerance);
    Ok(())
}

pub(crate) fn assert_forward_and_roundtrip(
    definition: &str,
    input: Coor4D,
    expected: Coor4D,
    forward_tolerance: f64,
    roundtrip_tolerance: f64,
) -> Result<(), Error> {
    let mut ctx = Minimal::default();
    let op = ctx.op(definition)?;
    let mut operands = [input];

    ctx.apply(op, Fwd, &mut operands)?;
    assert!(operands[0].hypot2(&expected) < forward_tolerance);

    ctx.apply(op, Inv, &mut operands)?;
    assert!(operands[0].hypot2(&input) < roundtrip_tolerance);
    Ok(())
}
