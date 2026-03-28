//! Transverse Cylindrical Equal Area
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let k_0 = op.params.k(0);
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let x = lat.cos() * lam.sin() / k_0;
        let y = k_0 * (lat.tan().atan2(lam.cos()) - lat_0);
        operands.set_xy(i, x_0 + a * x, y_0 + a * y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let k_0 = op.params.k(0);
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let mut x = (operands.xy(i).0 - x_0) / a;
        let mut y = (operands.xy(i).1 - y_0) / a;
        y = y / k_0 + lat_0;
        x *= k_0;
        let t = (1.0 - x * x).sqrt();
        let lat = (t * y.sin()).asin();
        let lon = x.atan2(t * y.cos()) + lon_0;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
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
    params.real.insert("lat_0", params.lat(0).to_radians());
    params.real.insert("lon_0", params.lon(0).to_radians());

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
    fn tcea_spherical_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("tcea a=6400000")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(223_322.760_576_727, 111_769.145_040_586, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
