//! Equidistant Cylindrical (Plate Carree)
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(rc) = op.params.real("rc") else {
        return 0;
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let x = x_0 + rc * (lon - lon_0);
        let y = y_0 + a * (lat - lat_0);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(rc) = op.params.real("rc") else {
        return 0;
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let x = operands.xy(i).0 - x_0;
        let y = operands.xy(i).1 - y_0;
        let lon = lon_0 + x / rc;
        let lat = lat_0 + y / a;
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
    OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let lat_ts = params.real("lat_ts")?;
    if lat_ts.abs() > 90.0 {
        return Err(Error::General(
            "Eqc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }

    let lat_0 = params.lat(0).to_radians();
    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", params.lon(0).to_radians());

    let ellps = params.ellps(0);
    let lat_ts_rad = lat_ts.to_radians();
    let cos_lat_ts = lat_ts_rad.cos();
    if cos_lat_ts <= 0.0 {
        return Err(Error::General(
            "Eqc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }

    params
        .real
        .insert("rc", ellps.semimajor_axis() * cos_lat_ts);

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
    fn eqc_forward_matches_proj_plate_carree() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqc ellps=WGS84")?;

        let geo = [Coor4D::geo(47., 2., 0., 0.)];
        let projected = [Coor4D::raw(222_638.981_586_547_13, 5_232_016.067_283_858, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        Ok(())
    }

    #[test]
    fn eqc_forward_matches_proj_with_offsets_and_lat_ts() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqc ellps=WGS84 lat_ts=30 lat_0=10 lon_0=1 x_0=100 y_0=200")?;

        let geo = [Coor4D::geo(50., 3., 0., 0.)];
        let projected = [Coor4D::raw(192_911.013_926_645_74, 4_452_979.631_730_943, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn eqc_spherical_forward_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqc ellps=6371000,0 lat_ts=60")?;

        let geo = [Coor4D::geo(30., 10., 0., 0.)];
        let projected = [Coor4D::raw(555_974.633_222_793_7, 3_335_847.799_336_7617, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
