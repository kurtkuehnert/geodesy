use super::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct AxisSwap {
    dimensionality: usize,
    pos: [usize; 4],
    sgn: [f64; 4],
}

impl PointOp for AxisSwap {
    const NAME: &'static str = "axisswap";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Series { key: "order", default: None },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        if params.given.contains_key("axis") && params.given.contains_key("order") {
            return Err(Error::BadParam(
                "axis/order".to_string(),
                "axisswap does not allow both axis and order on the same step".to_string(),
            ));
        }

        let order = params.series("order").map_err(|_| {
            Error::BadParam(
                "order".to_string(),
                "axisswap requires an explicit order or axis parameter".to_string(),
            )
        })?;

        if order.len() > 4 {
            return Err(Error::BadParam(
                "order".to_string(),
                "More than 4 indices given".to_string(),
            ));
        }

        let dimensionality = order.len();
        let mut pos = [0_usize, 1, 2, 3];
        let mut sgn = [1.0, 1.0, 1.0, 1.0];
        let mut seen = [false; 4];

        for (index, value) in order.iter().copied().enumerate() {
            let axis = value as i64;
            if (axis as f64) != value
                || axis == 0
                || (axis.unsigned_abs() as usize) > dimensionality
            {
                return Err(Error::BadParam("order".to_string(), value.to_string()));
            }

            let source = value.abs() as usize - 1;
            if seen[source] {
                return Err(Error::BadParam(
                    "order".to_string(),
                    "duplicate axis specified".to_string(),
                ));
            }
            seen[source] = true;

            pos[index] = source;
            sgn[index] = 1.0_f64.copysign(value);
        }

        Ok(Self {
            dimensionality,
            pos,
            sgn,
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let mut out = coord;
        for index in 0..self.dimensionality {
            out[index] = coord[self.pos[index]] * self.sgn[index];
        }
        Some(out)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let mut out = coord;
        for index in 0..self.dimensionality {
            out[self.pos[index]] = coord[index] * self.sgn[index];
        }
        Some(out)
    }
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn four_dim() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("axisswap order=2,1,-3,-4")?;
        let mut operands = [Coor4D([1., 2., 3., 4.])];

        ctx.apply(op, Fwd, &mut operands)?;
        assert_eq!(operands[0], Coor4D([2., 1., -3., -4.]));

        ctx.apply(op, Inv, &mut operands)?;
        assert_eq!(operands[0], Coor4D([1., 2., 3., 4.]));

        Ok(())
    }

    #[test]
    fn no_args_is_error() {
        let mut ctx = Minimal::default();
        assert!(ctx.op("axisswap").is_err());
    }

    #[test]
    fn two_dim() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("axisswap order=2,-1")?;
        let mut operands = [Coor4D([1., 2., 3., 4.])];

        ctx.apply(op, Fwd, &mut operands)?;
        assert_eq!(operands[0], Coor4D([2., -1., 3., 4.]));

        ctx.apply(op, Inv, &mut operands)?;
        assert_eq!(operands[0], Coor4D([1., 2., 3., 4.]));
        Ok(())
    }

    #[test]
    fn bad_parameters() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        assert!(ctx.op("axisswap order=4,4,4,2,-1").is_err());
        assert!(ctx.op("axisswap order=4,-4,2,-1").is_err());
        assert!(ctx.op("axisswap order=2,3").is_err());
        assert!(ctx.op("axisswap order").is_err());
        assert!(ctx.op("axisswap").is_err());

        Ok(())
    }
}
