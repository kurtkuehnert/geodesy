//! Geographic offset.
use crate::authoring::*;

const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;

#[derive(Clone, Copy, Debug)]
pub(crate) struct GeogOffset {
    dlon: f64,
    dlat: f64,
    dh: f64,
}

impl PointOp for GeogOffset {
    const NAME: &'static str = "geogoffset";
    const TITLE: &'static str = "Geographic offsets";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Real { key: "dlat", default: Some(0_f64) },
        OpParameter::Real { key: "dlon", default: Some(0_f64) },
        OpParameter::Real { key: "dh", default: Some(0_f64) },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self {
            dlat: params.real("dlat")? * ARCSEC_TO_RAD,
            dlon: params.real("dlon")? * ARCSEC_TO_RAD,
            dh: params.real("dh")?,
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D::raw(
            coord[0] + self.dlon,
            coord[1] + self.dlat,
            coord[2] + self.dh,
            coord[3],
        ))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        Some(Coor4D::raw(
            coord[0] - self.dlon,
            coord[1] - self.dlat,
            coord[2] - self.dh,
            coord[3],
        ))
    }
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
