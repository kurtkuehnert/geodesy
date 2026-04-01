use crate::authoring::*;

const DEFAULT_PROJ_FORWARD_TOLERANCE: f64 = 1e-6;
const DEFAULT_PROJ_ROUNDTRIP_TOLERANCE: f64 = 1e-10;

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

pub(crate) fn assert_proj_match_with_tol(
    definition: &str,
    input: Coor4D,
    expected: Coor4D,
    forward_tolerance: f64,
    roundtrip_tolerance: f64,
) -> Result<(), Error> {
    assert_forward_and_roundtrip(
        definition,
        input,
        expected,
        forward_tolerance,
        roundtrip_tolerance,
    )
}

pub(crate) fn assert_proj_match(
    definition: &str,
    input: Coor4D,
    expected: Coor4D,
) -> Result<(), Error> {
    assert_proj_match_with_tol(
        definition,
        input,
        expected,
        DEFAULT_PROJ_FORWARD_TOLERANCE,
        DEFAULT_PROJ_ROUNDTRIP_TOLERANCE,
    )
}

pub(crate) fn assert_roundtrip(
    definition: &str,
    input: Coor4D,
    tolerance: f64,
) -> Result<(), Error> {
    let mut ctx = Minimal::default();
    let op = ctx.op(definition)?;
    let mut operands = [input];

    ctx.apply(op, Fwd, &mut operands)?;
    ctx.apply(op, Inv, &mut operands)?;
    assert!(operands[0].hypot2(&input) < tolerance);
    Ok(())
}

pub(crate) fn assert_inverse(
    definition: &str,
    input: Coor4D,
    expected: Coor4D,
    tolerance: f64,
) -> Result<(), Error> {
    let mut ctx = Minimal::default();
    let op = ctx.op(definition)?;
    let mut operands = [input];

    ctx.apply(op, Inv, &mut operands)?;
    assert!(operands[0].hypot2(&expected) < tolerance);
    Ok(())
}

pub(crate) fn assert_inverse_rejects(
    definition: &str,
    input: Coor4D,
) -> Result<(), Error> {
    let mut ctx = Minimal::default();
    let op = ctx.op(definition)?;
    let mut operands = [input];

    ctx.apply(op, Inv, &mut operands)?;
    assert!(operands[0][0].is_nan());
    Ok(())
}

pub(crate) fn assert_op_err(definition: &str) {
    let mut ctx = Minimal::default();
    assert!(ctx.op(definition).is_err());
}
