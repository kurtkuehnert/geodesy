//! Mercator
use crate::authoring::*;
use crate::projection::ProjectionFrame;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps",  default: Some("GRS80") },
    OpParameter::Real { key: "lat_0",  default: Some(0_f64) },
    OpParameter::Real { key: "lon_0",  default: Some(0_f64) },
    OpParameter::Real { key: "x_0",    default: Some(0_f64) },
    OpParameter::Real { key: "y_0",    default: Some(0_f64) },
    OpParameter::Real { key: "k_0",    default: Some(1_f64) },
    OpParameter::Real { key: "lat_ts", default: Some(0_f64) },
];

struct MercState {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<MercState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let easting = state.frame.lon_delta(lon) * state.frame.k_0 * state.frame.a;
        let isometric = state.ellps.latitude_geographic_to_isometric(lat);
        let northing = state.frame.a * state.frame.k_0 * isometric;

        let (x, y) = state.frame.apply_false_origin(easting, northing);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

// ----- I N V E R S E -----------------------------------------------------------------

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<MercState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x, y) = state.frame.remove_false_origin(x_raw, y_raw);
        let lon = x / (state.frame.a * state.frame.k_0) + state.frame.lon_0;
        let psi = y / (state.frame.a * state.frame.k_0);
        let lat = state.ellps.latitude_isometric_to_geographic(psi);
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

// ----- C O N S T R U C T O R ---------------------------------------------------------

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let ellps = params.ellps(0);

    let lat_ts = params.real("lat_ts")?;
    if lat_ts.abs() > std::f64::consts::FRAC_PI_2 {
        return Err(Error::General(
            "Merc: Invalid value for lat_ts: |lat_ts| should be <= 90°",
        ));
    }

    // lat_ts trumps k_0
    if lat_ts != 0.0 {
        let sc = lat_ts.sin_cos();
        let k_0 = sc.1 / (1. - ellps.eccentricity_squared() * sc.0 * sc.0).sqrt();
        params.real.insert("k_0", k_0);
    }

    let state = MercState {
        frame: ProjectionFrame::from_params(&params),
        ellps,
    };
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op::with_state(descriptor, params, state))
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merc() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "merc";
        let op = ctx.op(definition)?;

        // Validation value from PROJ: echo 12 55 0 0 | cct -d18 +proj=merc
        // followed by quadrant tests from PROJ builtins.gie
        let geo = [
            Coor4D::geo(55., 12., 0., 0.),
            Coor4D::geo(1., 2., 0., 0.),
            Coor4D::geo(-1., 2., 0., 0.),
            Coor4D::geo(1., -2., 0., 0.),
            Coor4D::geo(-1., -2., 0., 0.),
        ];

        let projected = [
            Coor4D::raw(1_335_833.889_519_282_8, 7_326_837.714_873_877, 0., 0.),
            Coor4D::raw(222_638.981_586_547, 110_579.965_218_249, 0., 0.),
            Coor4D::raw(222_638.981_586_547, -110_579.965_218_249, 0., 0.),
            Coor4D::raw(-222_638.981_586_547, 110_579.965_218_249, 0., 0.),
            Coor4D::raw(-222_638.981_586_547, -110_579.965_218_249, 0., 0.),
        ];

        // Forward
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 20e-9);
        }

        // Roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 20e-9);
        }

        Ok(())
    }

    #[test]
    fn merc_lat_ts() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "merc lat_ts=56";
        let op = ctx.op(definition)?;

        let geo = [Coor4D::geo(55., 12., 0., 0.)];

        // Validation value from PROJ: echo 12 55 0 0 | cct -d18 +proj=merc +lat_ts=56
        let projected = [Coor4D::raw(
            748_713.257_925_886_8,
            4_106_573.862_841_270_4,
            0.,
            0.,
        )];

        // Forward
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 20e-9);
        }

        // Roundtrip
        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 20e-9);
        }

        Ok(())
    }

    #[test]
    fn merc_with_false_origin_and_central_meridian() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "merc lat_ts=-41 lon_0=100 x_0=1234 y_0=5678 ellps=WGS84";
        let op = ctx.op(definition)?;

        // Validation value from PROJ:
        // echo 141 -33 0 0 | cct -d12 +proj=merc +lat_ts=-41 +lon_0=100 +x_0=1234 +y_0=5678 +ellps=WGS84
        let geo = [Coor4D::geo(-33., 141., 0., 0.)];
        let projected = [Coor4D::raw(
            3_450_776.589_667_497,
            -2_920_802.103_023_224_5,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(
            operands[0].hypot2(&projected[0]) < 20e-6,
            "expected {:?}, got {:?}",
            projected[0],
            operands[0]
        );

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 20e-9);

        Ok(())
    }

    #[test]
    fn merc_with_parser_normalized_prime_meridian() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "merc lon_0=13.5";
        let op = ctx.op(definition)?;

        let geo = [Coor4D::geo(45., 13.5, 0., 0.)];
        let mut projected = geo;
        ctx.apply(op, Fwd, &mut projected)?;
        assert!(projected[0].x().abs() < 1e-12);

        ctx.apply(op, Inv, &mut projected)?;
        assert!(projected[0].hypot2(&geo[0]) < 20e-9);

        Ok(())
    }

    #[test]
    fn merc_spherical_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "merc ellps=6356751.9,0";
        let op = ctx.op(definition)?;

        let mut projected = [Coor4D::raw(10_000.0, 20_000.0, 0.0, 0.0)];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() - 0.090133735615956).abs() < 1e-12);
        assert!((projected[0][1].to_degrees() - 0.180267173822635).abs() < 1e-12);
        Ok(())
    }
}
