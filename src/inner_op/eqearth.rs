//! Equal Earth
use crate::authoring::*;

const A1: f64 = 1.340_264;
const A2: f64 = -0.081_106;
const A3: f64 = 0.000_893;
const A4: f64 = 0.003_796;
const M: f64 = 0.866_025_403_784_438_6;
const MAX_Y: f64 = 1.317_362_759_157_4;
const EPS: f64 = 1e-11;
const MAX_ITER: usize = 12;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(rqda) = op.params.real("rqda") else {
        return 0;
    };
    let spherical = op.params.boolean("spherical");
    let Ok(qp) = op.params.real("qp") else {
        return 0;
    };
    let authalic = if spherical {
        None
    } else {
        op.params.fourier_coefficients("authalic").ok()
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let sbeta = if spherical {
            lat.sin()
        } else {
            (ancillary::qs(lat.sin(), ellps.eccentricity()) / qp).clamp(-1.0, 1.0)
        };

        let psi = (M * sbeta).asin();
        let psi2 = psi * psi;
        let psi6 = psi2 * psi2 * psi2;
        let x = x_0
            + a * rqda * (lon - lon_0) * psi.cos()
                / (M * (A1 + 3.0 * A2 * psi2 + psi6 * (7.0 * A3 + 9.0 * A4 * psi2)));
        let y = y_0 + a * rqda * psi * (A1 + A2 * psi2 + psi6 * (A3 + A4 * psi2));

        if !spherical {
            let _ = authalic.as_ref().unwrap();
        }

        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(rqda) = op.params.real("rqda") else {
        return 0;
    };
    let spherical = op.params.boolean("spherical");
    let authalic = if spherical {
        None
    } else {
        op.params.fourier_coefficients("authalic").ok()
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let mut x = operands.xy(i).0 - x_0;
        let mut y = operands.xy(i).1 - y_0;
        x /= a * rqda;
        y = (y / (a * rqda)).clamp(-MAX_Y, MAX_Y);

        let mut psi = y;
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let psi2 = psi * psi;
            let psi6 = psi2 * psi2 * psi2;
            let f = psi * (A1 + A2 * psi2 + psi6 * (A3 + A4 * psi2)) - y;
            let fder = A1 + 3.0 * A2 * psi2 + psi6 * (7.0 * A3 + 9.0 * A4 * psi2);
            let delta = f / fder;
            psi -= delta;
            if delta.abs() < EPS {
                converged = true;
                break;
            }
        }

        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        let psi2 = psi * psi;
        let psi6 = psi2 * psi2 * psi2;
        let lon = lon_0
            + M * x * (A1 + 3.0 * A2 * psi2 + psi6 * (7.0 * A3 + 9.0 * A4 * psi2)) / psi.cos();
        let beta = (psi.sin() / M).clamp(-1.0, 1.0).asin();
        let lat = if spherical {
            beta
        } else {
            ellps.latitude_authalic_to_geographic(beta, authalic.as_ref().unwrap())
        };
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 5] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let ellps = params.ellps(0);
    if ellps.flattening() == 0.0 {
        params.boolean.insert("spherical");
        params.real.insert("qp", 2.0);
        params.real.insert("rqda", 1.0);
    } else {
        let qp = ancillary::qs(1.0, ellps.eccentricity());
        let rqda = (0.5 * qp).sqrt();
        params.real.insert("qp", qp);
        params.real.insert("rqda", rqda);
        let authalic = ellps.coefficients_for_authalic_latitude_computations();
        params.fourier_coefficients.insert("authalic", authalic);
    }

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eqearth_origin_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqearth lon_0=-90 ellps=WGS84")?;

        let geo = [Coor4D::geo(0., -90., 0., 0.)];
        let projected = [Coor4D::raw(0.0, 0.0, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-12);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
