/// Auxiliary latitudes
use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
enum LatitudeMode {
    Geocentric,
    Reduced,
    Conformal(FourierCoefficients),
    Rectifying(FourierCoefficients),
    Authalic(FourierCoefficients),
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Latitude {
    ellps: Ellipsoid,
    mode: LatitudeMode,
}

impl PointOp for Latitude {
    const NAME: &'static str = "latitude";
    const TITLE: &'static str = "Latitude Conversion";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Flag { key: "geocentric" },
        OpParameter::Flag { key: "reduced" },
        OpParameter::Flag { key: "parametric" },
        OpParameter::Flag { key: "conformal" },
        OpParameter::Flag { key: "authalic" },
        OpParameter::Flag { key: "rectifying" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") }
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);

        let mut modes = [
            params.boolean("geocentric").then_some(LatitudeMode::Geocentric),
            (params.boolean("reduced") || params.boolean("parametric"))
                .then_some(LatitudeMode::Reduced),
            params.boolean("conformal").then(|| {
                LatitudeMode::Conformal(ellps.coefficients_for_conformal_latitude_computations())
            }),
            params.boolean("authalic").then(|| {
                LatitudeMode::Authalic(ellps.coefficients_for_authalic_latitude_computations())
            }),
            params.boolean("rectifying").then(|| {
                LatitudeMode::Rectifying(ellps.coefficients_for_rectifying_latitude_computations())
            }),
        ]
        .into_iter()
        .flatten();

        let Some(mode) = modes.next() else {
            return Err(Error::MissingParam("latitude: must specify exactly one of flags authalic/conformal/geocentric/rectifying/reduced/parametric".to_string()));
        };
        if modes.next().is_some() {
            return Err(Error::MissingParam("latitude: must specify exactly one of flags authalic/conformal/geocentric/rectifying/reduced/parametric".to_string()));
        }

        Ok(Self { ellps, mode })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let mut coord = coord;
        coord[1] = match self.mode {
            LatitudeMode::Geocentric => self.ellps.latitude_geographic_to_geocentric(coord[1]),
            LatitudeMode::Reduced => self.ellps.latitude_geographic_to_reduced(coord[1]),
            LatitudeMode::Conformal(coefficients) => {
                self.ellps.latitude_geographic_to_conformal(coord[1], &coefficients)
            }
            LatitudeMode::Rectifying(coefficients) => {
                self.ellps.latitude_geographic_to_rectifying(coord[1], &coefficients)
            }
            LatitudeMode::Authalic(coefficients) => {
                self.ellps.latitude_geographic_to_authalic(coord[1], &coefficients)
            }
        };
        Some(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let mut coord = coord;
        coord[1] = match self.mode {
            LatitudeMode::Geocentric => self.ellps.latitude_geocentric_to_geographic(coord[1]),
            LatitudeMode::Reduced => self.ellps.latitude_reduced_to_geographic(coord[1]),
            LatitudeMode::Conformal(coefficients) => {
                self.ellps.latitude_conformal_to_geographic(coord[1], &coefficients)
            }
            LatitudeMode::Rectifying(coefficients) => {
                self.ellps.latitude_rectifying_to_geographic(coord[1], &coefficients)
            }
            LatitudeMode::Authalic(coefficients) => {
                self.ellps.latitude_authalic_to_geographic(coord[1], &coefficients)
            }
        };
        Some(coord)
    }
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn latitude() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        // Geocentric
        let op = ctx.op("latitude geocentric ellps=GRS80")?;
        let mut operands = [Coor4D::geo(55., 12., 0., 0.)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 54.818_973_308_324_5).abs() < 1e-12);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 55.).abs() < 1e-12);

        // Reduced (alias parametric)
        let op = ctx.op("latitude reduced ellps=GRS80")?;
        let mut operands = [Coor4D::geo(55., 12., 0., 0.)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 54.909_538_187_092_245).abs() < 1e-12);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 55.).abs() < 1e-12);

        // And vice versa: Parametric (alias Reduced)
        let op = ctx.op("latitude parametric ellps=GRS80")?;
        let mut operands = [Coor4D::geo(55., 12., 0., 0.)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 54.909_538_187_092_245).abs() < 1e-12);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 55.).abs() < 1e-12);

        // Conformal
        let op = ctx.op("latitude conformal ellps=GRS80")?;
        let mut operands = [Coor4D::geo(55., 12., 0., 0.)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 54.819_109_023_689_02).abs() < 1e-12);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 55.).abs() < 1e-12);

        // Rectifying
        let op = ctx.op("latitude rectifying ellps=GRS80")?;
        let mut operands = [Coor4D::geo(55., 12., 0., 0.)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 54.772_351_809_646_84).abs() < 1e-12);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 55.).abs() < 1e-12);

        // Authalic
        let op = ctx.op("latitude authalic ellps=GRS80")?;
        let mut operands = [Coor4D::geo(55., 12., 0., 0.)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 54.879_361_594_517_796).abs() < 1e-12);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][1].to_degrees() - 55.).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn geocentric_matches_proj_mercury_ographic_to_ocentric() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "axisswap order=-2,1 | unitconvert xy_in=deg xy_out=rad | latitude geocentric ellps=2440530,1075.123348017621 | unitconvert xy_in=rad xy_out=deg | axisswap order=2,1",
        )?;

        let source = [Coor4D::raw(-53.419625812798856, 49.743387956501664, 0., 0.)];
        let expected = [Coor4D::raw(-53.3685811656, -49.7433879565, 0., 0.)];

        let mut operands = source;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-9);
        Ok(())
    }
}
