//! Azimuthal Equidistant
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let origin = Coor4D::raw(lon_0, lat_0, 0.0, 0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let target = Coor4D::raw(lon, lat, 0.0, 0.0);
        let inv = ellps.geodesic_inv(&origin, &target);
        let azimuth = inv[0];
        let distance = inv[2];
        let x = x_0 + distance * azimuth.sin();
        let y = y_0 + distance * azimuth.cos();
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let origin = Coor4D::raw(lon_0, lat_0, 0.0, 0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let dx = x - x_0;
        let dy = y - y_0;
        let distance = dx.hypot(dy);
        if distance == 0.0 {
            operands.set_xy(i, lon_0, lat_0);
            successes += 1;
            continue;
        }
        let azimuth = dx.atan2(dy);
        let dest = ellps.geodesic_fwd(&origin, azimuth, distance);
        operands.set_xy(i, dest[0], dest[1]);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
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
    fn aeqd_roundtrip_origin() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aeqd lat_0=52 lon_0=-97.5 x_0=8264722.177 y_0=4867518.353 ellps=WGS84")?;

        let geo = [Coor4D::geo(52., -97.5, 0., 0.)];
        let projected = [Coor4D::raw(8_264_722.177, 4_867_518.353, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn aeqd_forward_reference() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aeqd lat_0=52 lon_0=-97.5 x_0=8264722.177 y_0=4867518.353 ellps=WGS84")?;

        let geo = [Coor4D::geo(61.390407158824, -101.971128034161, 0., 0.)];
        let projected = [Coor4D::raw(8_024_875.974_4, 5_920_866.114_0, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-3);
        Ok(())
    }
}
