//! Albers Equal Area
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::FRAC_PI_2;

const STANDARD_PARALLEL_TOLERANCE: f64 = 1e-10;
const AUTHALIC_LIMIT_TOLERANCE: f64 = 1e-7;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "lat_1", default: None },
    OpParameter::Real { key: "lat_2", default: None },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug)]
struct AeaState {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
    authalic: Option<FourierCoefficients>,
    n: f64,
    c: f64,
    dd: f64,
    rho0: f64,
    qp: f64,
    ec: f64,
    spherical: bool,
}

impl AeaState {
    fn new(params: &ParsedParameters, phi0: f64, phi1: f64, phi2: f64) -> Result<AeaState, Error> {
        if phi1.abs() > FRAC_PI_2 || phi2.abs() > FRAC_PI_2 {
            return Err(Error::BadParam("lat_1/lat_2".to_string(), params.name.clone()));
        }
        if (phi1 + phi2).abs() < STANDARD_PARALLEL_TOLERANCE {
            return Err(Error::General(
                "Aea: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
            ));
        }

        let ellps = params.ellps(0);
        let frame = ProjectionFrame::from_params(params);
        let spherical = ellps.flattening() == 0.0;

        let (sinphi1, cosphi1) = phi1.sin_cos();
        let secant = (phi1 - phi2).abs() >= STANDARD_PARALLEL_TOLERANCE;
        let mut n = sinphi1;

        if spherical {
            if secant {
                n = 0.5 * (n + phi2.sin());
            }
            let n2 = n + n;
            let c = cosphi1 * cosphi1 + n2 * sinphi1;
            let dd = 1.0 / n;
            let rho0 = dd * (c - n2 * phi0.sin()).sqrt();
            return Ok(Self {
                frame,
                ellps,
                authalic: None,
                n,
                c,
                dd,
                rho0,
                qp: 2.0,
                ec: 2.0,
                spherical,
            });
        }

        let e = ellps.eccentricity();
        let es = ellps.eccentricity_squared();
        let m1 = ancillary::pj_msfn((sinphi1, cosphi1), es);
        let q1 = ancillary::qs(sinphi1, e);
        if secant {
            let (sinphi2, cosphi2) = phi2.sin_cos();
            let m2 = ancillary::pj_msfn((sinphi2, cosphi2), es);
            let q2 = ancillary::qs(sinphi2, e);
            n = (m1 * m1 - m2 * m2) / (q2 - q1);
        }
        let authalic = ellps.coefficients_for_authalic_latitude_computations();
        let ec = 1.0 - 0.5 * (1.0 - es) * ((1.0 - e) / (1.0 + e)).ln() / e;
        let c = m1 * m1 + n * q1;
        let dd = 1.0 / n;
        let qp = ancillary::qs(1.0, e);
        let rho0 = dd * (c - n * ancillary::qs(phi0.sin(), e)).sqrt();
        Ok(Self {
            frame,
            ellps,
            authalic: Some(authalic),
            n,
            c,
            dd,
            rho0,
            qp,
            ec,
            spherical,
        })
    }
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<AeaState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let rho_term = if state.spherical {
            state.c - 2.0 * state.n * lat.sin()
        } else {
            state.c - state.n * ancillary::qs(lat.sin(), state.ellps.eccentricity())
        };
        if rho_term < 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let rho = state.dd * rho_term.sqrt();
        let theta = state.frame.lon_delta(lon) * state.n;
        let x = state.frame.a * rho * theta.sin();
        let y = state.frame.a * (state.rho0 - rho * theta.cos());
        let (x, y) = state.frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<AeaState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (mut x, y_local) = state.frame.remove_false_origin(x_raw, y_raw);
        x /= state.frame.a;
        let mut y = state.rho0 - y_local / state.frame.a;
        let mut rho = x.hypot(y);
        if rho != 0.0 {
            if state.n < 0.0 {
                rho = -rho;
                x = -x;
                y = -y;
            }

            let lat = if state.spherical {
                let arg = (state.c - (rho / state.dd).powi(2)) / (2.0 * state.n);
                arg.clamp(-1.0, 1.0).asin()
            } else {
                let authalic = state
                    .authalic
                    .expect("ellipsoidal AEA state must carry authalic coefficients");
                let qs = (state.c - (rho / state.dd).powi(2)) / state.n;
                if (state.ec - qs.abs()).abs() > AUTHALIC_LIMIT_TOLERANCE {
                    if qs.abs() > 2.0 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    let xi = (qs / state.qp).clamp(-1.0, 1.0).asin();
                    state.ellps.latitude_authalic_to_geographic(xi, &authalic)
                } else if qs < 0.0 {
                    -FRAC_PI_2
                } else {
                    FRAC_PI_2
                }
            };

            if lat.is_nan() {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let lon = x.atan2(y) / state.n + state.frame.lon_0;
            operands.set_xy(i, lon, lat);
            successes += 1;
        } else {
            let lat = if state.n > 0.0 { FRAC_PI_2 } else { -FRAC_PI_2 };
            operands.set_xy(i, state.frame.lon_0, lat);
            successes += 1;
        }
    }
    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    let phi0 = params.lat(0);
    let phi1 = params.lat(1);
    let phi2 = params.lat(2);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    let state = AeaState::new(&params, phi0, phi1, phi2)?;
    Ok(Op::with_state(descriptor, params, state))
}

pub fn leac(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;

    let phi0 = params.lat(0);
    let phi1 = if params.boolean("south") {
        -FRAC_PI_2
    } else {
        FRAC_PI_2
    };
    let phi2 = params.lat(1);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    let state = AeaState::new(&params, phi0, phi1, phi2)?;
    Ok(Op::with_state(descriptor, params, state))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aea_origin_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aea lat_0=23 lon_0=-96 lat_1=29.5 lat_2=45.5 ellps=GRS80")?;

        let geo = [Coor4D::geo(23., -96., 0., 0.)];
        let projected = [Coor4D::raw(0.0, 0.0, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn aea_inverse_matches_proj_metric_coordinates() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op =
            ctx.op("aea lat_0=0 lon_0=-120 lat_1=34 lat_2=40.5 x_0=0 y_0=-4000000 ellps=GRS80")?;

        let mut projected = [Coor4D::raw(0.0, -112_982.409_1, 0.0, 0.0)];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() + 120.0).abs() < 1e-8);
        assert!((projected[0][1].to_degrees() - 37.0).abs() < 1e-8);
        Ok(())
    }

    #[test]
    fn aea_spherical_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aea lat_0=40 lon_0=0 lat_1=20 lat_2=60 ellps=6378136.6,0")?;

        let mut projected = [
            Coor4D::raw(0.0, 0.0, 0.0, 0.0),
            Coor4D::raw(10_000.0, 20_000.0, 0.0, 0.0),
        ];
        ctx.apply(op, Inv, &mut projected)?;

        assert!(projected[0][0].abs() < 1e-12);
        assert!((projected[0][1].to_degrees() - 40.0).abs() < 1e-10);
        assert!((projected[1][0].to_degrees() - 0.124940293483244).abs() < 1e-10);
        assert!((projected[1][1].to_degrees() - 40.169004441322194).abs() < 1e-10);
        Ok(())
    }

    #[test]
    fn wraps_longitude_difference_across_dateline() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aea lat_0=50 lon_0=-154 lat_1=55 lat_2=65 ellps=clrk66")?;

        let geo = [Coor4D::geo(60., 179., 0., 0.)];
        let projected = [Coor4D::raw(
            -1_459_959.150_054_334_7,
            1_413_239.948_137_74,
            0.,
            0.,
        )];
        let ellps = Ellipsoid::named("clrk66")?;

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(ellps.distance(&operands[0], &geo[0]) < 1e-8);
        Ok(())
    }
}
