//! Lambert azimuthal equal area: EPSG coordinate operation method 9820, implemented
//! following [IOGP, 2019](crate::Bibliography::Iogp19), pp. 78-80
use crate::authoring::*;
use crate::projection::{AuthalicLatitude, ProjectionAspect, ProjectionFrame};

use std::f64::consts::FRAC_PI_2;

const ASPECT_TOLERANCE: f64 = 1e-10;
const POLAR_DOMAIN_TOLERANCE: f64 = 1e-15;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug)]
struct LaeaState {
    frame: ProjectionFrame,
    authalic: AuthalicLatitude,
    aspect: ProjectionAspect,
    qp: f64,
    rq: f64,
    dd: f64,
    xmf: f64,
    ymf: f64,
    sinb1: f64,
    cosb1: f64,
}

impl LaeaState {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let frame = ProjectionFrame::from_params(params);
        let lat_0 = frame.lat_0;
        let abs_lat_0 = lat_0.abs();

        if lat_0.is_nan() || abs_lat_0 > FRAC_PI_2 + ASPECT_TOLERANCE {
            return Err(Error::BadParam("lat_0".to_string(), params.name.clone()));
        }

        let aspect = ProjectionAspect::classify(lat_0, ASPECT_TOLERANCE);

        let ellps = params.ellps(0);
        let authalic = AuthalicLatitude::new(ellps);
        let e = ellps.eccentricity();
        let es = ellps.eccentricity_squared();
        let qp = authalic.qp();
        let rq = (0.5 * qp).sqrt();

        let (dd, xmf, ymf, sinb1, cosb1) = if aspect.is_polar() {
            (1.0, 1.0, 1.0, 0.0, 1.0)
        } else if aspect.is_equatorial() {
            (rq.recip(), 1.0, 0.5 * qp, 0.0, 1.0)
        } else {
            let (sin_phi_0, cos_phi_0) = lat_0.sin_cos();
            let xi_0 = (ancillary::qs(sin_phi_0, e) / qp).asin();
            let (sinb1, cosb1) = xi_0.sin_cos();
            let dd = cos_phi_0 / ((1.0 - es * sin_phi_0 * sin_phi_0).sqrt() * rq * cosb1);
            (dd, rq * dd, rq / dd, sinb1, cosb1)
        };

        Ok(Self {
            frame,
            authalic,
            aspect,
            qp,
            rq,
            dd,
            xmf,
            ymf,
            sinb1,
            cosb1,
        })
    }
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<LaeaState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let (sin_lon, cos_lon) = state.frame.lon_delta(lon).sin_cos();

        let xi = state.authalic.beta_from_phi(lat);
        let (sin_xi, cos_xi) = xi.sin_cos();
        let authalic_q = sin_xi * state.qp;

        let (x, y) = if state.aspect.is_oblique() {
            let denom = 1.0 + state.sinb1 * sin_xi + state.cosb1 * cos_xi * cos_lon;
            if denom.abs() < ASPECT_TOLERANCE {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let scale = (2.0 / denom).sqrt();
            (
                state.xmf * scale * cos_xi * sin_lon,
                state.ymf * scale * (state.cosb1 * sin_xi - state.sinb1 * cos_xi * cos_lon),
            )
        } else if state.aspect.is_equatorial() {
            let denom = 1.0 + cos_xi * cos_lon;
            if denom.abs() < ASPECT_TOLERANCE {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let scale = (2.0 / denom).sqrt();
            (
                state.xmf * scale * cos_xi * sin_lon,
                state.ymf * scale * sin_xi,
            )
        } else if state.aspect.is_north_polar() {
            if (lat + state.frame.lat_0).abs() < ASPECT_TOLERANCE {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let q = state.qp - authalic_q;
            if q < POLAR_DOMAIN_TOLERANCE {
                (0.0, 0.0)
            } else {
                let root_q = q.sqrt();
                (root_q * sin_lon, -root_q * cos_lon)
            }
        } else {
            if (lat + state.frame.lat_0).abs() < ASPECT_TOLERANCE {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let q = state.qp + authalic_q;
            if q < POLAR_DOMAIN_TOLERANCE {
                (0.0, 0.0)
            } else {
                let root_q = q.sqrt();
                (root_q * sin_lon, root_q * cos_lon)
            }
        };

        let (x, y) = state
            .frame
            .apply_false_origin(state.frame.a * x, state.frame.a * y);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<LaeaState>();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (mut x, mut y) = state.frame.remove_false_origin(x_raw, y_raw);
        x /= state.frame.a;
        y /= state.frame.a;

        let (ab, lam) = if state.aspect.is_oblique() || state.aspect.is_equatorial() {
            x /= state.dd;
            y *= state.dd;

            let rho = x.hypot(y);
            if rho < ASPECT_TOLERANCE {
                operands.set_xy(i, state.frame.lon_0, state.frame.lat_0);
                successes += 1;
                continue;
            }

            let asin_argument = 0.5 * rho / state.rq;
            if asin_argument > 1.0 {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }

            let c = 2.0 * asin_argument.asin();
            let (sin_c, cos_c) = c.sin_cos();
            x *= sin_c;

            if state.aspect.is_oblique() {
                let ab = cos_c * state.sinb1 + y * sin_c * state.cosb1 / rho;
                let y = rho * state.cosb1 * cos_c - y * state.sinb1 * sin_c;
                (ab, x.atan2(y))
            } else {
                let ab = y * sin_c / rho;
                let y = rho * cos_c;
                (ab, x.atan2(y))
            }
        } else {
            if state.aspect.is_north_polar() {
                y = -y;
            }
            let q = x * x + y * y;
            if q == 0.0 {
                operands.set_xy(i, state.frame.lon_0, state.frame.lat_0);
                successes += 1;
                continue;
            }
            let mut ab = 1.0 - q / state.qp;
            if state.aspect.is_south_polar() {
                ab = -ab;
            }
            (ab, x.atan2(y))
        };

        if ab.abs() > 1.0 + ASPECT_TOLERANCE {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        let lat = state.authalic.phi_from_beta(ab.clamp(-1.0, 1.0).asin());
        let lon = state.frame.lon_0 + lam;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    let state = LaeaState::new(&params)?;
    Ok(Op::with_state(descriptor, params, state))
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_eq::assert_float_eq;

    #[test]
    fn laea_oblique() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        let op = ctx.op("laea ellps=GRS80 lat_0=52 lon_0=10 x_0=4321000 y_0=3210000")?;

        let geo = [Coor2D::geo(50.0, 5.0)];
        let projected = [Coor2D::raw(3962799.45, 2999718.85)];

        let mut operands = geo;

        ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0].0, projected[0].0, abs_all <= 0.01);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][0].to_degrees() - 5.0).abs() < 1e-12);
        assert!((operands[0][1].to_degrees() - 50.).abs() < 1e-12);

        let mut operands = [Coor4D::raw(1e30, 1e30, 0.0, 0.0)];
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0][0].is_nan());

        Ok(())
    }

    #[test]
    fn origin() {
        let mut ctx = Minimal::new();
        let op = ctx
            .op("laea lon_0=10 lat_0=52 x_0=4321000 y_0=3210000")
            .unwrap();
        let mut data = [Coor2D::geo(52.0, 10.0)];
        let clone = data;
        ctx.apply(op, Fwd, &mut data).unwrap();
        ctx.apply(op, Inv, &mut data).unwrap();
        assert_eq!(data, clone);
    }

    #[test]
    fn laea_polar_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        let cases = [
            (
                "laea ellps=GRS80 lat_0=90 lon_0=0 x_0=0 y_0=0",
                Coor2D::geo(20.0, 80.0),
            ),
            (
                "laea ellps=GRS80 lat_0=-90 lon_0=0 x_0=0 y_0=0",
                Coor2D::geo(-45.0, -70.0),
            ),
        ];

        for (definition, geographic) in cases {
            let op = ctx.op(definition)?;
            let mut operands = [geographic];
            ctx.apply(op, Fwd, &mut operands)?;
            ctx.apply(op, Inv, &mut operands)?;
            assert!((operands[0][0].to_degrees() - geographic[0].to_degrees()).abs() < 1e-11);
            assert!((operands[0][1].to_degrees() - geographic[1].to_degrees()).abs() < 1e-11);
        }

        Ok(())
    }
}
