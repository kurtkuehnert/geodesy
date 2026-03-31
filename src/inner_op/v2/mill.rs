//! Miller Cylindrical
use crate::authoring::*;
use crate::projection::{ProjectionFrame, projection_gamut};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Mill {
    frame: ProjectionFrame,
}

impl PointOp for Mill {
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = projection_gamut!();

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self {
            frame: ProjectionFrame::from_params(params),
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let x_local = self.frame.a * self.frame.remove_central_meridian(lon);
        let y_local = self.frame.a * 1.25 * (std::f64::consts::FRAC_PI_4 + 0.4 * lat).tan().ln();
        let (x, y) = self.frame.apply_false_origin(x_local, y_local);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let (x_local, y_local) = self.frame.remove_false_origin(coord[0], coord[1]);
        let lam = x_local / self.frame.a;
        let y = y_local / self.frame.a;
        let lon = self.frame.apply_central_meridian(lam);
        let lat = 2.5 * ((0.8 * y).exp().atan() - std::f64::consts::FRAC_PI_4);
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mill_matches_proj_gie_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("mill a=6400000")?;
        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(
            223_402.144_255_274,
            111_704.701_754_394,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
