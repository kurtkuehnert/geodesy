//! Equidistant Cylindrical (Plate Carree)
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let rc = op.params.real["rc"];
    let m0 = op.params.real["m0"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let x = x_0 + rc * (lon - lon_0);
        let y = y_0 + ellps.meridian_latitude_to_distance(lat) - m0;
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let rc = op.params.real["rc"];
    let m0 = op.params.real["m0"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let x = operands.xy(i).0 - x_0;
        let y = operands.xy(i).1 - y_0;
        let lon = lon_0 + x / rc;
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
    let given = parameters.instantiated_as.split_into_parameters();
    if given.contains_key("a") && !given.contains_key("ellps") {
        let a = given["a"]
            .parse::<f64>()
            .map_err(|_| Error::MissingParam("a".to_string()))?;
        if a <= 0.0 {
            return Err(Error::General(
                "Eqc: Invalid value for a: a must be positive",
            ));
        }

        let rf = if let Some(rf) = given.get("rf") {
            rf.parse::<f64>()
                .map_err(|_| Error::MissingParam("rf".to_string()))?
        } else if let Some(f) = given.get("f") {
            let f = f
                .parse::<f64>()
                .map_err(|_| Error::MissingParam("f".to_string()))?;
            if f == 0.0 { 0.0 } else { 1.0 / f }
        } else if let Some(b) = given.get("b") {
            let b = b
                .parse::<f64>()
                .map_err(|_| Error::MissingParam("b".to_string()))?;
            if b <= 0.0 {
                return Err(Error::General(
                    "Eqc: Invalid value for b: b must be positive",
                ));
            }
            if (a - b).abs() < f64::EPSILON {
                0.0
            } else {
                a / (a - b)
            }
        } else {
            0.0
        };
        params.text.insert("ellps", format!("{a},{rf}"));
    }

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

    params.real.insert(
        "rc",
        ellps.prime_vertical_radius_of_curvature(lat_ts_rad) * cos_lat_ts,
    );
    params
        .real
        .insert("m0", ellps.meridian_latitude_to_distance(lat_0));

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
}
