use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct AddOne;

impl PointOp for AddOne {
    const GAMUT: &'static [OpParameter] = &[];

    fn build(_params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self)
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D([coord[0] + 1., coord[1], coord[2], coord[3]]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D([coord[0] - 1., coord[1], coord[2], coord[3]]))
    }
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addone() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("addone")?;
        let mut data = crate::test_data::coor2d();
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 56.);
        assert_eq!(data[1][0], 60.);
        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);
        Ok(())
    }
}
