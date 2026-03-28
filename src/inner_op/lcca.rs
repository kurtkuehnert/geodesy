//! Lambert Conformal Conic Alternative
use crate::authoring::*;

const DEL_TOL: f64 = 1e-12;
const MAX_ITER: usize = 10;

fn fs(s: f64, c: f64) -> f64 {
    s * (1.0 + s * s * c)
}

fn fsp(s: f64, c: f64) -> f64 {
    1.0 + 3.0 * s * s * c
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let k_0 = op.params.k(0);
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let r0 = op.params.real["r0"];
    let l = op.params.real["l"];
    let m0 = op.params.real["m0"];
    let c = op.params.real["c"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let s = ellps.meridian_latitude_to_distance(lat) / a - m0;
        let dr = fs(s, c);
        let r = r0 - dr;
        let theta = (lon - lon_0) * l;
        operands.set_xy(
            i,
            x_0 + a * k_0 * r * theta.sin(),
            y_0 + a * k_0 * (r0 - r * theta.cos()),
        );
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let k_0 = op.params.k(0);
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let r0 = op.params.real["r0"];
    let l = op.params.real["l"];
    let m0 = op.params.real["m0"];
    let c = op.params.real["c"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let x = (operands.xy(i).0 - x_0) / (a * k_0);
        let y = (operands.xy(i).1 - y_0) / (a * k_0);
        let theta = x.atan2(r0 - y);
        let dr = y - x * (0.5 * theta).tan();
        let lon = theta / l + lon_0;
        let mut s = dr;
        let mut converged = false;
        for _ in 0..MAX_ITER {
            let dif = (fs(s, c) - dr) / fsp(s, c);
            s -= dif;
            if dif.abs() < DEL_TOL {
                converged = true;
                break;
            }
        }
        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let lat = ellps.meridian_distance_to_latitude((s + m0) * a);
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 9] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: None },
    OpParameter::Real { key: "lat_1", default: Some(0_f64) },
    OpParameter::Real { key: "lat_2", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    let phi0 = params
        .real("lat_0")
        .map_err(|_| Error::MissingParam("lat_0".to_string()))?
        .to_radians();
    if phi0 == 0.0 {
        return Err(Error::General(
            "Lcca: Invalid value for lat_0: it should be different from 0.",
        ));
    }

    let ellps = params.ellps(0);
    let es = ellps.eccentricity_squared();
    let sinphi0 = phi0.sin();
    let s2p0 = sinphi0 * sinphi0;
    let r = 1.0 / (1.0 - es * s2p0);
    let n0 = r.sqrt();
    let r0m = r * (1.0 - es) * n0;
    let tan0 = phi0.tan();

    params.real.insert("lat_0", phi0);
    params.real.insert("lon_0", params.lon(0).to_radians());
    params.real.insert("l", sinphi0);
    params.real.insert("m0", ellps.meridian_latitude_to_distance(phi0) / ellps.semimajor_axis());
    params.real.insert("r0", n0 / tan0);
    params.real.insert("c", 1.0 / (6.0 * r0m * n0));

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
    fn lcca_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("lcca ellps=GRS80 lat_0=1 lat_1=0.5 lat_2=2")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(222_605.285_770_237, 67.806_007_272, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        Ok(())
    }
}
