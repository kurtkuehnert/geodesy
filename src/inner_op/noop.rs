//! The no-operation. Does nothing, and is good at it
use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct NoOp;

impl PointOp for NoOp {
    const NAME: &'static str = "noop";
    const TITLE: &'static str = "No operation";
    const GAMUT: &'static [OpParameter] = &[];

    fn build(_params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self)
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(coord)
    }
}

// ----- T E S T S ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    const GDA94: Coor4D = Coor4D([-4052051.7643, 4212836.2017, -2545106.0245, 0.0]);

    #[test]
    fn no_change() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("noop")?;

        // EPSG:1134 - 3 parameter, ED50/WGS84, s = sqrt(27) m
        let mut operands = [GDA94];

        // Forward
        ctx.apply(op, Fwd, &mut operands)?;
        assert_eq!(operands[0], GDA94);

        // Inverse + roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        assert_eq!(operands[0], GDA94);
        Ok(())
    }
}
