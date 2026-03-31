//! Cassini-Soldner
use crate::authoring::*;
use crate::math::angular;
use crate::projection::ProjectionFrame;

const C1: f64 = 1.0 / 6.0;
const C2: f64 = 1.0 / 120.0;
const C3: f64 = 1.0 / 24.0;
const C4: f64 = 1.0 / 3.0;
const C5: f64 = 1.0 / 15.0;
const INVERSE_CONVERGENCE_TOLERANCE: f64 = 1e-12;
const NUMERICAL_JACOBIAN_STEP: f64 = 1e-6;
const JACOBIAN_DETERMINANT_TOLERANCE: f64 = 1e-24;
const INV_ITER: usize = 15;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
    OpParameter::Flag { key: "hyperbolic" },
];

struct CassState {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
    m0: f64,
    spherical: bool,
    hyperbolic: bool,
    es: f64,
}

#[allow(clippy::too_many_arguments)]
fn fwd_ellipsoidal(
    ellps: &Ellipsoid,
    lon_0: f64,
    x_0: f64,
    y_0: f64,
    m0: f64,
    hyperbolic: bool,
    lon: f64,
    lat: f64,
) -> (f64, f64) {
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let lam = angular::normalize_symmetric(lon - lon_0);
    let (sinphi, cosphi) = lat.sin_cos();
    let m = ellps.meridian_latitude_to_distance(lat) / a;
    let nu_sq = 1.0 / (1.0 - es * sinphi * sinphi);
    let nu = nu_sq.sqrt();
    let tanphi = lat.tan();
    let t = tanphi * tanphi;
    let a1 = lam * cosphi;
    let c = es * cosphi * cosphi / one_es;
    let a2 = a1 * a1;
    let x_norm = nu * a1 * (1.0 - a2 * t * (C1 + (8.0 - t + 8.0 * c) * a2 * C2));
    let mut y_norm = (m - m0) + nu * tanphi * a2 * (0.5 + (5.0 - t + 6.0 * c) * a2 * C3);
    if hyperbolic {
        let rho = nu_sq * one_es * nu;
        y_norm -= y_norm * y_norm * y_norm / (6.0 * rho * nu);
    }
    (x_0 + a * x_norm, y_0 + a * y_norm)
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<CassState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);

        let (x, y) = if state.spherical {
            let lam = state.frame.remove_central_meridian(lon);
            (
                state.frame.a * (lat.cos() * lam.sin()).asin(),
                state.frame.a * (lat.tan().atan2(lam.cos()) - state.frame.lat_0),
            )
        } else {
            fwd_ellipsoidal(
                &state.ellps,
                state.frame.lon_0,
                state.frame.x_0,
                state.frame.y_0,
                state.m0,
                state.hyperbolic,
                lon,
                lat,
            )
        };

        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<CassState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = state
            .frame
            .remove_false_origin(operands.xy(i).0, operands.xy(i).1);

        let (lon, lat) = if state.spherical {
            let dd = y / state.frame.a + state.frame.lat_0;
            (
                state.frame.lon_0 + (x / state.frame.a).tan().atan2(dd.cos()),
                (dd.sin() * (x / state.frame.a).cos()).asin(),
            )
        } else {
            let phi1 = state
                .ellps
                .meridian_distance_to_latitude(state.frame.a * (state.m0 + y / state.frame.a));
            let tanphi1 = phi1.tan();
            let t1 = tanphi1 * tanphi1;
            let sinphi1 = phi1.sin();
            let nu1_sq = 1.0 / (1.0 - state.es * sinphi1 * sinphi1);
            let nu1 = nu1_sq.sqrt();
            let rho1 = nu1_sq * (1.0 - state.es) * nu1;
            let d = x / (state.frame.a * nu1);
            let d2 = d * d;
            let lon = state.frame.lon_0
                + d * (1.0 + t1 * d2 * (-C4 + (1.0 + 3.0 * t1) * d2 * C5)) / phi1.cos();
            let lat = phi1 - (nu1 * tanphi1 / rho1) * d2 * (0.5 - (1.0 + 3.0 * t1) * d2 * C3);
            let (mut lon, mut lat) = (lon, lat);

            for _ in 0..INV_ITER {
                let (fx, fy) = fwd_ellipsoidal(
                    &state.ellps,
                    state.frame.lon_0,
                    state.frame.x_0,
                    state.frame.y_0,
                    state.m0,
                    state.hyperbolic,
                    lon,
                    lat,
                );
                let delta_x = fx - (x + state.frame.x_0);
                let delta_y = fy - (y + state.frame.y_0);
                if delta_x.abs() < INVERSE_CONVERGENCE_TOLERANCE
                    && delta_y.abs() < INVERSE_CONVERGENCE_TOLERANCE
                {
                    break;
                }

                let h_lon = if lon > 0.0 {
                    -NUMERICAL_JACOBIAN_STEP
                } else {
                    NUMERICAL_JACOBIAN_STEP
                };
                let h_lat = if lat > 0.0 {
                    -NUMERICAL_JACOBIAN_STEP
                } else {
                    NUMERICAL_JACOBIAN_STEP
                };
                let (fx_lon, fy_lon) = fwd_ellipsoidal(
                    &state.ellps,
                    state.frame.lon_0,
                    state.frame.x_0,
                    state.frame.y_0,
                    state.m0,
                    state.hyperbolic,
                    lon + h_lon,
                    lat,
                );
                let (fx_lat, fy_lat) = fwd_ellipsoidal(
                    &state.ellps,
                    state.frame.lon_0,
                    state.frame.x_0,
                    state.frame.y_0,
                    state.m0,
                    state.hyperbolic,
                    lon,
                    lat + h_lat,
                );
                let j11 = (fx_lon - fx) / h_lon;
                let j21 = (fy_lon - fy) / h_lon;
                let j12 = (fx_lat - fx) / h_lat;
                let j22 = (fy_lat - fy) / h_lat;
                let det = j11 * j22 - j12 * j21;
                if det.abs() < JACOBIAN_DETERMINANT_TOLERANCE {
                    break;
                }
                let inv11 = j22 / det;
                let inv12 = -j12 / det;
                let inv21 = -j21 / det;
                let inv22 = j11 / det;
                let delta_lon = (delta_x * inv11 + delta_y * inv12).clamp(-0.3, 0.3);
                let delta_lat = (delta_x * inv21 + delta_y * inv22).clamp(-0.3, 0.3);
                lon = (lon - delta_lon).clamp(-std::f64::consts::PI, std::f64::consts::PI);
                lat = (lat - delta_lat)
                    .clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
            }

            (lon, lat)
        };

        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    let frame = ProjectionFrame::from_params(&params);

    let ellps = params.ellps(0);
    let spherical = ellps.flattening() == 0.0;
    let state = CassState {
        frame,
        ellps,
        m0: if spherical {
            frame.lat_0
        } else {
            ellps.meridian_latitude_to_distance(frame.lat_0) / ellps.semimajor_axis()
        },
        spherical,
        hyperbolic: params.boolean("hyperbolic"),
        es: ellps.eccentricity_squared(),
    };

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op::with_state(descriptor, params, state))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cass_origin_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cass lat_0=31.7340969444444 lon_0=35.2120805555556 x_0=170251.555 y_0=126867.909 ellps=6378300.789,293.4663155389811")?;

        let geo = [Coor4D::geo(31.7340969444444, 35.2120805555556, 0., 0.)];
        let projected = [Coor4D::raw(170251.555, 126867.909, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn hyperbolic_cass_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "cass hyperbolic lat_0=-16.25 lon_0=179.333333333333 x_0=251727.9155424 y_0=334519.953768 ellps=6378306.3696,293.4663076567816",
        )?;

        let mut projected = [
            Coor4D::raw(251727.9155424, 334519.953768, 0.0, 0.0),
            Coor4D::raw(261727.9155424, 354519.953768, 0.0, 0.0),
        ];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() - 179.33333333333303).abs() < 1e-10);
        assert!((projected[0][1].to_degrees() + 16.25).abs() < 1e-10);
        assert!((projected[1][0].to_degrees() - 179.42679063431376).abs() < 1e-10);
        assert!((projected[1][1].to_degrees() + 16.06923331216869).abs() < 1e-10);
        Ok(())
    }

    #[test]
    fn wraps_longitude_difference_across_dateline() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "cass hyperbolic lat_0=-16.25 lon_0=179.333333333333 x_0=251727.9155424 y_0=334519.953768 ellps=6378306.3696,293.4663076567816",
        )?;

        let geo = [Coor4D::geo(-16.1, -179.5, 0., 0.)];
        let projected = [Coor4D::raw(
            376_542.242_742_735_4,
            350_765.385_231_374_24,
            0.,
            0.,
        )];
        let ellps = Ellipsoid::new(6_378_306.369_6, 1.0 / 293.466_307_656_781_6);

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 5e-5);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(ellps.distance(&operands[0], &geo[0]) < 1e-8);
        Ok(())
    }

    #[test]
    fn cass_world_inverse_matches_proj_far_from_central_meridian() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cass lat_0=0 lon_0=0 ellps=WGS84")?;

        let mut projected = [Coor4D::raw(
            -58_044.692_087_790_16,
            10667870.723512903,
            0.0,
            0.0,
        )];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() - 179.0).abs() < 1e-8);
        assert!((projected[0][1].to_degrees() - 80.0).abs() < 1e-8);
        Ok(())
    }

    #[test]
    fn cass_world_inverse_matches_proj_far_southern_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cass lat_0=0 lon_0=0 ellps=WGS84")?;

        let mut projected = [Coor4D::raw(
            251_999.407_927_484,
            -11_743_302.693_378_212,
            0.0,
            0.0,
        )];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() + 167.737_751_661_178_2).abs() < 1e-8);
        assert!((projected[0][1].to_degrees() + 76.449_895_966_577_9).abs() < 1e-8);
        Ok(())
    }
}
