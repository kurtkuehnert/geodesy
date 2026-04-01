/// Permanent tide systems
use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct PermTide {
    ellps: Ellipsoid,
    coefficient: f64,
}

impl PointOp for PermTide {
    const NAME: &'static str = "permtide";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Real { key: "k",     default: Some(0.3) },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
        OpParameter::Text { key: "from",  default: None },
        OpParameter::Text { key: "to",    default: None }
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let k = params.real("k")?;

        let Ok(to) = params.text("to") else {
            return Err(Error::MissingParam(
                "permtide: must specify 'to=' as exactly one of {'mean', 'zero', 'free'}"
                    .to_string(),
            ));
        };
        let Ok(from) = params.text("from") else {
            return Err(Error::MissingParam(
                "permtide: must specify 'from=' as exactly one of {'mean', 'zero', 'free'}"
                    .to_string(),
            ));
        };

        let coefficient = match (to.as_str(), from.as_str()) {
            ("mean", "mean") => 0.0,
            ("mean", "zero") => 1.0,
            ("mean", "free") => 1.0 + k,
            ("zero", "zero") => 0.0,
            ("zero", "mean") => -1.0,
            ("zero", "free") => k,
            ("free", "free") => 0.0,
            ("free", "mean") => -(1.0 + k),
            ("free", "zero") => -k,
            _ => f64::NAN,
        };

        if coefficient.is_nan() {
            return Err(Error::BadParam(
                "'to=' or 'from='".to_string(),
                "must be one of {'mean', 'zero', 'free'}".to_string(),
            ));
        }

        Ok(Self { ellps, coefficient })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let phibar = self.ellps.latitude_geographic_to_geocentric(coord[1]);
        let s = phibar.sin();
        Some(Coor4D::raw(
            coord[0],
            coord[1],
            coord[2] + self.coefficient * (-0.198) * (1.5 * s * s - 0.5),
            coord[3],
        ))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let phibar = self.ellps.latitude_geographic_to_geocentric(coord[1]);
        let s = phibar.sin();
        Some(Coor4D::raw(
            coord[0],
            coord[1],
            coord[2] - self.coefficient * (-0.198) * (1.5 * s * s - 0.5),
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
    fn permtide() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        // A test point near Copenhagen
        let pnt = Coor4D::geo(55., 12., 0., 0.);

        // Mean -> zero
        let op = ctx.op("permtide from=mean to=zero ellps=GRS80")?;
        let mut operands = [pnt];
        ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][2], 0.099407199, abs_all <= 1e-8);
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][2], pnt[2], abs_all <= 1e-12);

        // Mean -> free
        let op = ctx.op("permtide from=mean to=free ellps=GRS80")?;
        let mut operands = [pnt];
        ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0][2], 0.1292293587824579, abs_all <= 1e-8);
        ctx.apply(op, Inv, &mut operands)?;
        assert_float_eq!(operands[0][2], pnt[2], abs_all <= 1e-12);

        // Inversion
        let fwd_op = ctx.op("permtide from=mean to=zero ellps=GRS80")?;
        let inv_op = ctx.op("permtide from=zero to=mean ellps=GRS80 inv")?;
        let mut operands = [pnt];
        ctx.apply(fwd_op, Fwd, &mut operands)?;
        let fwd_h = operands[0][2];

        let mut operands = [pnt];
        ctx.apply(inv_op, Fwd, &mut operands)?;
        let inv_h = operands[0][2];
        assert_float_eq!(fwd_h, inv_h, abs_all <= 1e-20);

        // Bad tide system names
        assert!(matches!(
            ctx.op("permtide from=cheese to=zero ellps=GRS80"),
            Err(Error::BadParam(_, _))
        ));

        // Missing tide system names
        assert!(matches!(ctx.op("permtide"), Err(Error::MissingParam(_))));
        assert!(matches!(
            ctx.op("permtide to=zero ellps=GRS80"),
            Err(Error::MissingParam(_))
        ));
        assert!(matches!(
            ctx.op("permtide from=mean ellps=GRS80"),
            Err(Error::MissingParam(_))
        ));
        Ok(())
    }
}
