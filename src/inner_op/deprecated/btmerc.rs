//! Transverse Mercator, following [Bowring (1989)](crate::bibliography::Bibliography::Bow89)
use crate::authoring::*;
use crate::projection::ProjectionFrame;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
    OpParameter::Real { key: "k_0",   default: Some(1_f64) },
];

#[rustfmt::skip]
pub const UTM_GAMUT: [OpParameter; 4] = [
    OpParameter::Flag    { key: "inv" },
    OpParameter::Flag    { key: "south" },
    OpParameter::Text    { key: "ellps", default: Some("GRS80") },
    OpParameter::Natural { key: "zone",  default: None },
];

struct BtmercState {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
    eps: f64,
}

// Forward transverse mercator, following Bowring (1989)
fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<BtmercState>();

    let mut successes = 0_usize;
    let n = operands.len();
    for i in 0..n {
        let mut coord = operands.get_coord(i);

        let lat = coord[1] + state.frame.lat_0;
        let (s, c) = lat.sin_cos();
        let cc = c * c;
        let ss = s * s;

        let dlon = coord[0] - state.frame.lon_0;
        let oo = dlon * dlon;

        #[allow(non_snake_case)]
        let N = state.ellps.prime_vertical_radius_of_curvature(lat);
        let z = state.eps * dlon.powi(3) * c.powi(5) / 6.;
        let sd2 = (dlon / 2.).sin();

        let theta_2 = (2. * s * c * sd2 * sd2).atan2(ss + cc * dlon.cos());

        // Easting
        let sd = dlon.sin();
        coord[0] = state.frame.x_0
            + state.frame.k_0 * N * ((c * sd).atanh() + z * (1. + oo * (36. * cc - 29.) / 10.));

        // Northing
        let m = state.ellps.meridian_latitude_to_distance(lat);
        let znos4 = z * N * dlon * s / 4.;
        let ecc = 4. * state.eps * cc;
        coord[1] = state.frame.y_0
            + state.frame.k_0 * (m + N * theta_2 + znos4 * (9. + ecc + oo * (20. * cc - 11.)));
        operands.set_coord(i, &coord);
        successes += 1;
    }

    successes
}

// ----- I N V E R S E -----------------------------------------------------------------

// Inverse transverse mercator, following Bowring (1989)
fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<BtmercState>();

    let mut successes = 0_usize;
    let n = operands.len();
    for i in 0..n {
        let mut coord = operands.get_coord(i);
        // Footpoint latitude, i.e. the latitude of a point on the central meridian
        // having the same northing as the point of interest
        let lat = state
            .ellps
            .meridian_distance_to_latitude((coord[1] - state.frame.y_0) / state.frame.k_0);
        let (s, c) = lat.sin_cos();
        let t = s / c;
        let cc = c * c;
        #[allow(non_snake_case)]
        let N = state.ellps.prime_vertical_radius_of_curvature(lat);
        let x = (coord[0] - state.frame.x_0) / (state.frame.k_0 * N);
        let xx = x * x;
        let theta_4 = x.sinh().atan2(c);
        let theta_5 = (t * theta_4.cos()).atan();

        // Latitude
        let xet = xx * xx * state.eps * t / 24.;
        coord[1] = state.frame.lat_0 + (1. + cc * state.eps) * (theta_5 - xet * (9. - 10. * cc))
            - state.eps * cc * lat;

        // Longitude
        let approx = state.frame.lon_0 + theta_4;
        let coef = state.eps / 60. * xx * x * c;
        coord[0] = approx - coef * (10. - 4. * xx / cc + xx * cc);
        operands.set_coord(i, &coord);

        successes += 1;
    }

    successes
}

// ----- C O N S T R U C T O R ---------------------------------------------------------

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    let ellps = params.ellps(0);
    let state = BtmercState {
        frame: ProjectionFrame::from_params(&params),
        ellps,
        eps: ellps.second_eccentricity_squared(),
    };
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op::with_state(descriptor, params, state))
}

pub fn utm(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &UTM_GAMUT)?;

    // The UTM zone should be an integer between 1 and 60
    let zone = params.natural("zone")?;
    if !(1..61).contains(&zone) {
        return Err(Error::General(
            "UTM: 'zone' must be an integer in the interval 1..60",
        ));
    }

    crate::inner_op::apply_utm_defaults(&mut params, zone);

    let ellps = params.ellps(0);
    let state = BtmercState {
        frame: ProjectionFrame::from_params(&params),
        ellps,
        eps: ellps.second_eccentricity_squared(),
    };
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op::with_state(descriptor, params, state))
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn btmerc() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "btmerc k_0=0.9996 lon_0=9 x_0=500000";
        let op = ctx.op(definition)?;

        // Validation values from PROJ:
        // echo 12 55 0 0 | cct -d18 +proj=utm +zone=32 | clip
        #[rustfmt::skip]
        let geo = [
            Coor4D::geo( 55.,  12., 0., 0.),
            Coor4D::geo(-55.,  12., 0., 0.),
            Coor4D::geo( 55., -6., 0., 0.),
            Coor4D::geo(-55., -6., 0., 0.)
        ];

        #[rustfmt::skip]
        let projected = [
            Coor4D::raw( 691_875.632_139_661, 6_098_907.825_005_012, 0., 0.),
            Coor4D::raw( 691_875.632_139_661,-6_098_907.825_005_012, 0., 0.),
            Coor4D::raw(-455_673.814_189_040, 6_198_246.671_090_279, 0., 0.),
            Coor4D::raw(-455_673.814_189_040,-6_198_246.671_090_279, 0., 0.)
        ];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 5e-3);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 10e-8);
        }
        Ok(())
    }

    #[test]
    fn butm() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "butm zone=32";
        let op = ctx.op(definition)?;

        // Validation values from PROJ:
        // echo 12 55 0 0 | cct -d18 +proj=utm +zone=32 | clip
        #[rustfmt::skip]
        let geo = [
            Coor4D::geo( 55.,  12., 0., 0.),
            Coor4D::geo(-55.,  12., 0., 0.),
            Coor4D::geo( 55., -6., 0., 0.),
            Coor4D::geo(-55., -6., 0., 0.)
        ];

        #[rustfmt::skip]
        let projected = [
            Coor4D::raw( 691_875.632_139_661, 6_098_907.825_005_012, 0., 0.),
            Coor4D::raw( 691_875.632_139_661,-6_098_907.825_005_012, 0., 0.),
            Coor4D::raw(-455_673.814_189_040, 6_198_246.671_090_279, 0., 0.),
            Coor4D::raw(-455_673.814_189_040,-6_198_246.671_090_279, 0., 0.)
        ];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 5e-3);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 10e-8);
        }
        Ok(())
    }

    #[test]
    fn butm_south() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("butm zone=32 south")?;

        #[rustfmt::skip]
        let geo = [
            Coor4D::geo( 55.,  12., 0., 0.),
            Coor4D::geo(-55.,  12., 0., 0.),
            Coor4D::geo( 55., -6., 0., 0.),
            Coor4D::geo(-55., -6., 0., 0.)
        ];

        #[rustfmt::skip]
        let projected = [
            Coor4D::raw( 691_875.632_139_661, 1e7+6_098_907.825_005_012, 0., 0.),
            Coor4D::raw( 691_875.632_139_661, 1e7-6_098_907.825_005_012, 0., 0.),
            Coor4D::raw(-455_673.814_189_040, 1e7+6_198_246.671_090_279, 0., 0.),
            Coor4D::raw(-455_673.814_189_040, 1e7-6_198_246.671_090_279, 0., 0.)
        ];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 5e-3);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 10e-8);
        }

        Ok(())
    }
}
