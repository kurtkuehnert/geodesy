//! CalCOFI Line/Station conversions
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_4, PI};

const EPS10: f64 = 1e-10;
const DEG_TO_LINE: f64 = 5.0;
const DEG_TO_STATION: f64 = 15.0;
const LINE_TO_RAD: f64 = 0.003_490_658_503_988_659;
const STATION_TO_RAD: f64 = 0.001_163_552_834_662_886_3;
const PT_O_LINE: f64 = 80.0;
const PT_O_STATION: f64 = 60.0;
const PT_O_LAMBDA: f64 = -2.114_466_388_791_13;
const PT_O_PHI: f64 = 0.596_029_939_556_063_5;
const ROTATION_ANGLE: f64 = PI / 6.0;

fn sph_ts(phi: f64) -> f64 {
    (FRAC_PI_4 + 0.5 * phi).tan().ln()
}

fn sph_inv_ts(psi: f64) -> f64 {
    std::f64::consts::FRAC_PI_2 - 2.0 * (-psi).exp().atan()
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let ellipsoidal = e != 0.0;

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        if (lat.abs() - std::f64::consts::FRAC_PI_2).abs() <= EPS10 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let y = if ellipsoidal {
            -ancillary::ts(lat.sin_cos(), e).ln()
        } else {
            sph_ts(lat)
        };
        let oy = if ellipsoidal {
            -ancillary::ts(PT_O_PHI.sin_cos(), e).ln()
        } else {
            sph_ts(PT_O_PHI)
        };
        let l1 = (y - oy) * ROTATION_ANGLE.tan();
        let l2 = -lon - l1 + PT_O_LAMBDA;
        let ry_merc = l2 * ROTATION_ANGLE.cos() * ROTATION_ANGLE.sin() + y;
        let ry = if ellipsoidal {
            ancillary::pj_phi2((-ry_merc).exp(), e)
        } else {
            sph_inv_ts(ry_merc)
        };
        let x = PT_O_LINE - (ry - PT_O_PHI).to_degrees() * DEG_TO_LINE / ROTATION_ANGLE.cos();
        let y = PT_O_STATION + (ry - lat).to_degrees() * DEG_TO_STATION / ROTATION_ANGLE.sin();
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let ellipsoidal = e != 0.0;

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let ry = PT_O_PHI - LINE_TO_RAD * (x - PT_O_LINE) * ROTATION_ANGLE.cos();
        let lat = ry - STATION_TO_RAD * (y - PT_O_STATION) * ROTATION_ANGLE.sin();
        let oymctr = if ellipsoidal {
            -ancillary::ts(PT_O_PHI.sin_cos(), e).ln()
        } else {
            sph_ts(PT_O_PHI)
        };
        let rymctr = if ellipsoidal {
            -ancillary::ts(ry.sin_cos(), e).ln()
        } else {
            sph_ts(ry)
        };
        let xymctr = if ellipsoidal {
            -ancillary::ts(lat.sin_cos(), e).ln()
        } else {
            sph_ts(lat)
        };
        let l1 = (xymctr - oymctr) * ROTATION_ANGLE.tan();
        let l2 = (rymctr - xymctr) / (ROTATION_ANGLE.cos() * ROTATION_ANGLE.sin());
        let lon = PT_O_LAMBDA - (l1 + l2);
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 2] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        state: None,
        steps: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calcofi_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("calcofi ellps=GRS80")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(508.444_872_150, -1_171.764_860_418, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-9);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
