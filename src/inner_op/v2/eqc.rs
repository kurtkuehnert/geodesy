//! Equidistant Cylindrical (Plate Carree)
use crate::authoring::*;
use crate::projection::ProjectionFrame;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let rc = op.params.real["rc"];
    let m0 = op.params.real["m0"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let x = rc * frame.lon_delta_raw(lon);
        let y = ellps.meridian_latitude_to_distance(lat) - m0;
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let rc = op.params.real["rc"];
    let m0 = op.params.real["m0"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x, y) = frame.remove_false_origin(x_raw, y_raw);
        let lon = frame.lon_0 + x / rc;
        let lat = ellps.meridian_distance_to_latitude(y + m0);
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
    if lat_ts.abs() > std::f64::consts::FRAC_PI_2 {
        return Err(Error::General(
            "Eqc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }

    let ellps = params.ellps(0);
    let cos_lat_ts = lat_ts.cos();
    if cos_lat_ts <= 0.0 {
        return Err(Error::General(
            "Eqc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }

    params.real.insert(
        "rc",
        ellps.prime_vertical_radius_of_curvature(lat_ts) * cos_lat_ts,
    );
    params.real.insert(
        "m0",
        ellps.meridian_latitude_to_distance(params.lat(0)),
    );

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
    fn eqc_forward_matches_proj_plate_carree() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqc ellps=WGS84")?;

        let geo = [Coor4D::geo(47., 2., 0., 0.)];
        let projected = [Coor4D::raw(
            222_638.981_586_547_13,
            5_207_247.008_955_783,
            0.,
            0.,
        )];

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
        let projected = [Coor4D::raw(
            193_072.560_501_793,
            4_435_192.208_449_775,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-5);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn eqc_spherical_forward_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqc ellps=6371000,0 lat_ts=60")?;

        let geo = [Coor4D::geo(30., 10., 0., 0.)];
        let projected = [Coor4D::raw(
            555_974.633_222_793_7,
            3_335_847.799_336_761_7,
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

    #[test]
    fn eqc_supports_bare_semimajor_axis_sphere() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqc a=6400000")?;

        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(
            223_402.144_255_274,
            111_701.072_127_637,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        Ok(())
    }

    #[test]
    fn mercury_ocentric_eqc_pipeline_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "axisswap order=2,1 | unitconvert xy_in=deg xy_out=rad | inv latitude geocentric ellps=2440530,1075.123348017621 | eqc ellps=2440530,1075.123348017621",
        )?;

        let geo = [Coor4D::raw(45., 0., 0., 0.)];
        let projected = [Coor4D::raw(0., 1_916_463.956_798_680_6, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        Ok(())
    }
}
