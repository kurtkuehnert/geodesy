/// The unit conversion operator.
/// It has a subset of the conversions supported by PROJ.
/// NB: If no units are specified, the default is meters.
/// ...
/// Conversions are performed by means of a pivot unit.
/// For horizontal conversions, the pivot unit is meters for linear units and radians for angular units.
/// Vertical units always pivot around meters.
/// Unit_A => (meters || radians) => Unit_B
use super::units::UnitParam;
use crate::authoring::*;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 5] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "xy_in", default: Some("m") },
    OpParameter::Text { key: "xy_out", default: Some("m") },
    OpParameter::Text { key: "z_in", default: Some("m") },
    OpParameter::Text { key: "z_out", default: Some("m") },
];

// ----- F O R W A R D -----------------------------------------------------------------

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let mut successes = 0_usize;

    let xy_in_to_pivot = op.params.real("xy_in_to_pivot").unwrap();
    let pivot_to_xy_out = op.params.real("pivot_to_xy_out").unwrap();
    let xy = xy_in_to_pivot * pivot_to_xy_out;

    let z_in_to_pivot = op.params.real("z_in_to_pivot").unwrap();
    let pivot_to_z_out = op.params.real("pivot_to_z_out").unwrap();
    let z = z_in_to_pivot * pivot_to_z_out;

    for i in 0..operands.len() {
        let mut coord = operands.get_coord(i);
        coord[0] *= xy;
        coord[1] *= xy;
        coord[2] *= z;
        operands.set_coord(i, &coord);
        successes += 1;
    }

    successes
}

// ----- I N V E R S E -----------------------------------------------------------------

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let mut successes = 0_usize;

    let xy_in_to_pivot = op.params.real("xy_in_to_pivot").unwrap();
    let pivot_to_xy_out = op.params.real("pivot_to_xy_out").unwrap();
    let xy = xy_in_to_pivot * pivot_to_xy_out;

    let z_in_to_pivot = op.params.real("z_in_to_pivot").unwrap();
    let pivot_to_z_out = op.params.real("pivot_to_z_out").unwrap();
    let z = z_in_to_pivot * pivot_to_z_out;

    for i in 0..operands.len() {
        let mut coord = operands.get_coord(i);
        coord[0] /= xy;
        coord[1] /= xy;
        coord[2] /= z;
        operands.set_coord(i, &coord);
        successes += 1;
    }

    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let xy_in = params.text("xy_in").unwrap();
    let xy_out = params.text("xy_out").unwrap();
    let z_in = params.text("z_in").unwrap();
    let z_out = params.text("z_out").unwrap();

    let Some(xy_in) = UnitParam::lookup(xy_in.as_str()) else {
        return Err(Error::BadParam("xy_in".to_string(), xy_in));
    };
    let Some(xy_out) = UnitParam::lookup(xy_out.as_str()) else {
        return Err(Error::BadParam("xy_out".to_string(), xy_out));
    };
    let Some(z_in) = UnitParam::lookup(z_in.as_str()) else {
        return Err(Error::BadParam("z_in".to_string(), z_in));
    };
    let Some(z_out) = UnitParam::lookup(z_out.as_str()) else {
        return Err(Error::BadParam("z_out".to_string(), z_out));
    };

    if !xy_in.is_compatible_with(xy_out) {
        return Err(Error::BadParam(
            "xy_in/xy_out".to_string(),
            format!(
                "{}->{}",
                params.text("xy_in").unwrap(),
                params.text("xy_out").unwrap()
            ),
        ));
    }
    if !z_in.is_compatible_with(z_out) {
        return Err(Error::BadParam(
            "z_in/z_out".to_string(),
            format!(
                "{}->{}",
                params.text("z_in").unwrap(),
                params.text("z_out").unwrap()
            ),
        ));
    }

    params.real.insert("xy_in_to_pivot", xy_in.multiplier());
    params
        .real
        .insert("pivot_to_xy_out", 1. / xy_out.multiplier());
    params.real.insert("z_in_to_pivot", z_in.multiplier());
    params
        .real
        .insert("pivot_to_z_out", 1. / z_out.multiplier());

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));

    Ok(Op {
        descriptor,
        params,
        state: None,
        steps: None,
    })
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use float_eq::assert_float_eq;

    #[test]
    fn xyz_default_units() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        let op = ctx.op("unitconvert xy_in=us-ft z_in=us-ft")?;

        let mut operands = [Coor4D::raw(5., 5., 5., 1.)];

        // Forward
        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);

        assert_eq!(successes, 1);

        // Inverse + roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][0], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        Ok(())
    }

    #[test]
    fn xy_us_ft_round_trip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        let op = ctx.op("unitconvert xy_in=us-ft xy_out=us-ft")?;

        let mut operands = [Coor4D::raw(5., 5., 1., 1.)];

        // Forward
        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);

        assert_eq!(successes, 1);

        // Inverse + roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][0], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        Ok(())
    }

    #[test]
    fn xyz_us_ft_to_m() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        let op = ctx.op("unitconvert xy_in=us-ft xy_out=m z_in=us-ft z_out=m")?;

        let mut operands = [Coor4D::raw(5., 5., 5., 1.)];

        // Forward
        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);

        assert_eq!(successes, 1);

        // Inverse + roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][0], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        Ok(())
    }

    #[test]
    fn xy_yd_to_m() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        let op = ctx.op("unitconvert xy_in=us-yd xy_out=m")?;

        let mut operands = [Coor4D::raw(1000., 1000., 500., 1.)];

        // Forward
        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 914.40182880, abs_all <= 1e-5);
        assert_float_eq!(operands[0][1], 914.40182880, abs_all <= 1e-5);
        assert_float_eq!(operands[0][2], 500., abs_all <= 1e-5);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-5);

        assert_eq!(successes, 1);

        // Inverse + roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][0], 1000.0, abs_all <= 1e-5);
        assert_float_eq!(operands[0][1], 1000.0, abs_all <= 1e-5);
        assert_float_eq!(operands[0][2], 500., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        Ok(())
    }
    #[test]
    fn xy_grad_to_deg() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        let op = ctx.op("unitconvert xy_in=grad xy_out=deg")?;

        let mut operands = [Coor4D::raw(135.0, 40., 500., 1.)];

        // Forward
        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 121.5, abs_all <= 1e-5);
        assert_float_eq!(operands[0][1], 36.0, abs_all <= 1e-5);
        assert_float_eq!(operands[0][2], 500., abs_all <= 1e-5);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-5);

        assert_eq!(successes, 1);

        // Inverse + roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][0], 135.0, abs_all <= 1e-5);
        assert_float_eq!(operands[0][1], 40.0, abs_all <= 1e-5);
        assert_float_eq!(operands[0][2], 500., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        Ok(())
    }

    #[test]
    fn unknown_unit() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        assert!(ctx.op("unitconvert xy_in=unknown xy_out=deg").is_err());
        Ok(())
    }

    #[test]
    fn numeric_xy_multiplier_to_m() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        let op = ctx.op("unitconvert xy_in=0.3048 xy_out=m")?;

        let mut operands = [Coor4D::raw(10., 20., 5., 1.)];

        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_eq!(successes, 1);
        assert_float_eq!(operands[0][0], 3.048, abs_all <= 1e-12);
        assert_float_eq!(operands[0][1], 6.096, abs_all <= 1e-12);
        assert_float_eq!(operands[0][2], 5., abs_all <= 1e-12);

        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][0], 10., abs_all <= 1e-12);
        assert_float_eq!(operands[0][1], 20., abs_all <= 1e-12);
        assert_float_eq!(operands[0][2], 5., abs_all <= 1e-12);
        Ok(())
    }

    #[test]
    fn rejects_linear_to_angular_xy_conversion() {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        assert!(ctx.op("unitconvert xy_in=m xy_out=rad").is_err());
    }

    #[test]
    fn rejects_angular_to_linear_z_conversion() {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        assert!(ctx.op("unitconvert z_in=rad z_out=m").is_err());
    }

    #[test]
    fn rejects_zero_and_non_finite_numeric_units() {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(new));
        assert!(ctx.op("unitconvert xy_in=0").is_err());
        assert!(ctx.op("unitconvert xy_out=0").is_err());
        assert!(ctx.op("unitconvert xy_in=1e400").is_err());
        assert!(ctx.op("unitconvert z_out=1e400").is_err());
    }
}
