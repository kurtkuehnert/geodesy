//! Equal Area Cylindrical
use crate::authoring::*;

const EPS: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let k_0 = op.params.k(0);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let x = x_0 + a * k_0 * (lon - lon_0);
        let y = if spherical {
            y_0 + a * lat.sin() / k_0
        } else {
            y_0 + a * 0.5 * ancillary::qs(lat.sin(), e) / k_0
        };
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
    let k_0 = op.params.k(0);
    let spherical = op.params.boolean("spherical");
    let qp = op.params.real("qp").unwrap_or(2.0);
    let authalic = op.params.fourier_coefficients("authalic").ok();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let x = (x_raw - x_0) / a;
        let y = (y_raw - y_0) / a;
        let lon = lon_0 + x / k_0;

        let lat = if spherical {
            let arg = y * k_0;
            let t = arg.abs();
            if t > 1.0 + EPS {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            arg.clamp(-1.0, 1.0).asin()
        } else {
            let arg = 2.0 * y * k_0 / qp;
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
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    let ellps = params.ellps(0);

    let lat_ts = params.real("lat_ts")?;
    if lat_ts.abs() > 90.0 {
        return Err(Error::General(
            "CEA: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }
    if lat_ts != 0.0 {
        let phi_ts = lat_ts.to_radians();
        let (s, c) = phi_ts.sin_cos();
        let mut k_0 = c;
        if ellps.eccentricity_squared() != 0.0 {
            k_0 /= (1.0 - ellps.eccentricity_squared() * s * s).sqrt();
        }
        params.real.insert("k_0", k_0);
    }

    params.real.insert("lon_0", params.lon(0).to_radians());
    if ellps.flattening() == 0.0 {
        params.boolean.insert("spherical");
        params.real.insert("qp", 2.0);
    } else {
        let qp = ancillary::qs(1.0, ellps.eccentricity());
        params.real.insert("qp", qp);
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
    fn cea_matches_proj_gie_ellipsoidal_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cea ellps=GRS80")?;

        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(222_638.981_586_547, 110_568.812_396_267, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
