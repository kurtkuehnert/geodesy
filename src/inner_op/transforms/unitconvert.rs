//! Unit conversion operator.
//!
//! Attribution:
//! - PROJ 9.8.0 `unitconvert.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/conversions/unitconvert.cpp>
//!
//! This supports a subset of PROJ's conversions: horizontal and vertical
//! unit scaling through pivot units (meters for linear units, radians for
//! angular units, and meters for vertical units).

use crate::authoring::*;
use crate::units::UnitParam;

#[derive(Clone, Copy, Debug)]
pub(crate) struct UnitConvert {
    xy: f64,
    z: f64,
}

impl PointOp for UnitConvert {
    const NAME: &'static str = "unitconvert";
    const TITLE: &'static str = "Unit Conversion";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "xy_in", default: Some("m") },
        OpParameter::Text { key: "xy_out", default: Some("m") },
        OpParameter::Text { key: "z_in", default: Some("m") },
        OpParameter::Text { key: "z_out", default: Some("m") },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let xy_in_text = params.text("xy_in").unwrap();
        let xy_out_text = params.text("xy_out").unwrap();
        let z_in_text = params.text("z_in").unwrap();
        let z_out_text = params.text("z_out").unwrap();

        let Some(xy_in) = UnitParam::lookup(xy_in_text.as_str()) else {
            return Err(Error::BadParam("xy_in".to_string(), xy_in_text));
        };
        let Some(xy_out) = UnitParam::lookup(xy_out_text.as_str()) else {
            return Err(Error::BadParam("xy_out".to_string(), xy_out_text));
        };
        let Some(z_in) = UnitParam::lookup(z_in_text.as_str()) else {
            return Err(Error::BadParam("z_in".to_string(), z_in_text));
        };
        let Some(z_out) = UnitParam::lookup(z_out_text.as_str()) else {
            return Err(Error::BadParam("z_out".to_string(), z_out_text));
        };

        if !xy_in.is_compatible_with(xy_out) {
            return Err(Error::BadParam(
                "xy_in/xy_out".to_string(),
                format!("{xy_in_text}->{xy_out_text}"),
            ));
        }
        if !z_in.is_compatible_with(z_out) {
            return Err(Error::BadParam(
                "z_in/z_out".to_string(),
                format!("{z_in_text}->{z_out_text}"),
            ));
        }

        Ok(Self {
            xy: xy_in.multiplier() / xy_out.multiplier(),
            z: z_in.multiplier() / z_out.multiplier(),
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D::raw(
            coord[0] * self.xy,
            coord[1] * self.xy,
            coord[2] * self.z,
            coord[3],
        ))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D::raw(
            coord[0] / self.xy,
            coord[1] / self.xy,
            coord[2] / self.z,
            coord[3],
        ))
    }
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use float_eq::assert_float_eq;

    #[test]
    fn xyz_default_units() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        let op = ctx.op("unitconvert xy_in=us-ft z_in=us-ft")?;

        let mut operands = [Coor4D::raw(5., 5., 5., 1.)];

        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        assert_eq!(successes, 1);

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
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        let op = ctx.op("unitconvert xy_in=us-ft xy_out=us-ft")?;

        let mut operands = [Coor4D::raw(5., 5., 1., 1.)];

        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 5., abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1., abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        assert_eq!(successes, 1);

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
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        let op = ctx.op("unitconvert xy_in=us-ft xy_out=m z_in=us-ft z_out=m")?;

        let mut operands = [Coor4D::raw(5., 5., 5., 1.)];

        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][1], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][2], 1.524003048, abs_all <= 1e-9);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-9);
        assert_eq!(successes, 1);

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
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        let op = ctx.op("unitconvert xy_in=us-yd xy_out=m")?;

        let mut operands = [Coor4D::raw(1000., 1000., 500., 1.)];

        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 914.40182880, abs_all <= 1e-5);
        assert_float_eq!(operands[0][1], 914.40182880, abs_all <= 1e-5);
        assert_float_eq!(operands[0][2], 500., abs_all <= 1e-5);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-5);
        assert_eq!(successes, 1);

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
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        let op = ctx.op("unitconvert xy_in=grad xy_out=deg")?;

        let mut operands = [Coor4D::raw(135.0, 40., 500., 1.)];

        let successes = ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][0], 121.5, abs_all <= 1e-5);
        assert_float_eq!(operands[0][1], 36.0, abs_all <= 1e-5);
        assert_float_eq!(operands[0][2], 500., abs_all <= 1e-5);
        assert_float_eq!(operands[0][3], 1., abs_all <= 1e-5);
        assert_eq!(successes, 1);

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
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        assert!(ctx.op("unitconvert xy_in=unknown xy_out=deg").is_err());
        Ok(())
    }

    #[test]
    fn numeric_xy_multiplier_to_m() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
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
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        assert!(ctx.op("unitconvert xy_in=m xy_out=rad").is_err());
    }

    #[test]
    fn rejects_angular_to_linear_z_conversion() {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        assert!(ctx.op("unitconvert z_in=rad z_out=m").is_err());
    }

    #[test]
    fn rejects_zero_and_non_finite_numeric_units() {
        let mut ctx = Minimal::default();
        ctx.register_op("unitconvert", OpConstructor(Op::point::<UnitConvert>));
        assert!(ctx.op("unitconvert xy_in=0").is_err());
        assert!(ctx.op("unitconvert xy_out=0").is_err());
        assert!(ctx.op("unitconvert xy_in=1e400").is_err());
        assert!(ctx.op("unitconvert z_out=1e400").is_err());
    }
}
