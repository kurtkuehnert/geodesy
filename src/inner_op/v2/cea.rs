//! Equal Area Cylindrical
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const EPS: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let frame = ProjectionFrame::from_params(&op.params);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let x = frame.a * frame.k_0 * frame.lon_delta(lon);
        let y = if spherical {
            frame.a * lat.sin() / frame.k_0
        } else {
            frame.a * 0.5 * ancillary::qs(lat.sin(), e) / frame.k_0
        };
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let spherical = op.params.boolean("spherical");
    let qp = op.params.real("qp").unwrap_or(2.0);
    let authalic = op.params.fourier_coefficients("authalic").ok();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x_raw, y_raw);
        let x = x_local / frame.a;
        let y = y_local / frame.a;
        let lon = frame.lon_0 + x / frame.k_0;

        let lat = if spherical {
            let arg = y * frame.k_0;
            let t = arg.abs();
            if t > 1.0 + EPS {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            arg.clamp(-1.0, 1.0).asin()
        } else {
            let arg = 2.0 * y * frame.k_0 / qp;
            let t = arg.abs();
            if t > 1.0 + EPS {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let beta = arg.clamp(-1.0, 1.0).asin();
            ellps.latitude_authalic_to_geographic(beta, authalic.as_ref().unwrap())
        };

        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let ellps = params.ellps(0);

    let lat_ts = params.real("lat_ts")?;
    if lat_ts.abs() > std::f64::consts::FRAC_PI_2 {
        return Err(Error::General(
            "CEA: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }
    if lat_ts != 0.0 {
        let (s, c) = lat_ts.sin_cos();
        let mut k_0 = c;
        if ellps.eccentricity_squared() != 0.0 {
            k_0 /= (1.0 - ellps.eccentricity_squared() * s * s).sqrt();
        }
        params.real.insert("k_0", k_0);
    }

    super::insert_authalic_setup(&mut params);

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
    fn cea_matches_proj_gie_ellipsoidal_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cea ellps=GRS80")?;

        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(
            222_638.981_586_547,
            110_568.812_396_267,
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
