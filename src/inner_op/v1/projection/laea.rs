//! Lambert azimuthal equal area: EPSG coordinate operation method 9820, implemented
//! following [IOGP, 2019](crate::Bibliography::Iogp19), pp. 78-80
use crate::authoring::*;
use crate::projection::ProjectionFrame;

use std::f64::consts::FRAC_PI_2;
const EPS10: f64 = 1e-10;

#[derive(Clone, Copy)]
struct LaeaCache {
    qp: f64,
    rq: f64,
    dd: f64,
    xmf: f64,
    ymf: f64,
    sinb1: f64,
    cosb1: f64,
    oblique: bool,
    equatorial: bool,
    north_polar: bool,
    south_polar: bool,
}

impl LaeaCache {
    fn from_params(params: &ParsedParameters) -> Option<Self> {
        Some(Self {
            qp: params.real("qp").ok()?,
            rq: params.real("rq").ok()?,
            dd: params.real("dd").ok()?,
            xmf: params.real("xmf").ok()?,
            ymf: params.real("ymf").ok()?,
            sinb1: params.real("sinb1").unwrap_or(0.0),
            cosb1: params.real("cosb1").unwrap_or(1.0),
            oblique: params.boolean("oblique"),
            equatorial: params.boolean("equatorial"),
            north_polar: params.boolean("north_polar"),
            south_polar: params.boolean("south_polar"),
        })
    }
}

// ----- F O R W A R D -----------------------------------------------------------------

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Some(cache) = LaeaCache::from_params(&op.params) else {
        return 0;
    };

    let frame = ProjectionFrame::from_params(&op.params);
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();

    let mut successes = 0_usize;
    let n = operands.len();

    for i in 0..n {
        let (lon, lat) = operands.xy(i);
        let (sin_lon, cos_lon) = frame.lon_delta(lon).sin_cos();

        let xi = (ancillary::qs(lat.sin(), e) / cache.qp).asin();
        let (sin_xi, cos_xi) = xi.sin_cos();

        let mut q = sin_xi * cache.qp;
        let mut b = match () {
            _ if cache.oblique => 1.0 + cache.sinb1 * sin_xi + cache.cosb1 * cos_xi * cos_lon,
            _ if cache.equatorial => 1.0 + cos_xi * cos_lon,
            _ if cache.north_polar => {
                q = cache.qp - q;
                FRAC_PI_2 + lat
            }
            _ => {
                q += cache.qp;
                lat - FRAC_PI_2
            }
        };

        if b.abs() < EPS10 {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        }

        let (x, y) = if cache.oblique {
            b = (2.0 / b).sqrt();
            (
                cache.xmf * b * cos_xi * sin_lon,
                cache.ymf
                    * b
                    * (cache.cosb1 * sin_xi - cache.sinb1 * cos_xi * cos_lon),
            )
        } else if cache.equatorial {
            b = (2.0 / b).sqrt();
            (cache.xmf * b * cos_xi * sin_lon, cache.ymf * b * sin_xi)
        } else if q >= 1.0e-15 {
            let root_q = q.sqrt();
            let y = cos_lon * if cache.south_polar { root_q } else { -root_q };
            (root_q * sin_lon, y)
        } else {
            (0.0, 0.0)
        };

        let (x, y) = frame.apply_false_origin(frame.a * x, frame.a * y);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

// ----- I N V E R S E -----------------------------------------------------------------

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Some(cache) = LaeaCache::from_params(&op.params) else {
        return 0;
    };
    let Ok(authalic) = op.params.fourier_coefficients("authalic") else {
        return 0;
    };

    let frame = ProjectionFrame::from_params(&op.params);
    let lat_0 = op.params.lat(0);

    let ellps = op.params.ellps(0);
    let mut successes = 0_usize;
    let n = operands.len();

    for i in 0..n {
        let (mut x, mut y) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x, y);
        x = x_local / frame.a;
        y = y_local / frame.a;

        let (ab, lam) = if cache.oblique || cache.equatorial {
            x /= cache.dd;
            y *= cache.dd;
            let rho = x.hypot(y);
            if rho < EPS10 {
                operands.set_xy(i, frame.lon_0, lat_0);
                successes += 1;
                continue;
            }

            let asin_argument = 0.5 * rho / cache.rq;
            if asin_argument > 1.0 {
                operands.set_xy(i, f64::NAN, f64::NAN);
                continue;
            }

            let c = 2.0 * asin_argument.asin();
            let (s_ce, c_ce) = c.sin_cos();
            x *= s_ce;
            if cache.oblique {
                let ab = c_ce * cache.sinb1 + y * s_ce * cache.cosb1 / rho;
                y = rho * cache.cosb1 * c_ce - y * cache.sinb1 * s_ce;
                (ab, x.atan2(y))
            } else {
                let ab = y * s_ce / rho;
                y = rho * c_ce;
                (ab, x.atan2(y))
            }
        } else {
            if cache.north_polar {
                y = -y;
            }
            let q = x * x + y * y;
            if q == 0.0 {
                operands.set_xy(i, frame.lon_0, lat_0);
                successes += 1;
                continue;
            }
            let mut ab = 1.0 - q / cache.qp;
            if cache.south_polar {
                ab = -ab;
            }
            (ab, x.atan2(y))
        };

        if ab.abs() > 1.0 + EPS10 {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        }

        let lat = ellps.latitude_authalic_to_geographic(ab.clamp(-1.0, 1.0).asin(), &authalic);
        let lon = frame.lon_0 + lam;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

// ----- C O N S T R U C T O R ---------------------------------------------------------

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },

    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },

    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let lat_0 = params.lat(0);

    if lat_0.is_nan() {
        warn!("LAEA: Bad central latitude!");
        return Err(Error::BadParam("lat_0".to_string(), def.clone()));
    }

    let t = lat_0.abs();
    if t > FRAC_PI_2 + EPS10 {
        warn!("LAEA: Bad central latitude!");
        return Err(Error::BadParam("lat_0".to_string(), def.clone()));
    }

    let polar = (t - FRAC_PI_2).abs() < EPS10;
    let north = polar && (lat_0 > 0.0);
    let equatorial = !polar && t < EPS10;
    match (polar, equatorial, north) {
        (true, _, true) => params.boolean.insert("north_polar"),
        (true, _, false) => params.boolean.insert("south_polar"),
        (_, true, _) => params.boolean.insert("equatorial"),
        _ => params.boolean.insert("oblique"),
    };

    let ellps = params.ellps(0);
    let es = ellps.eccentricity_squared();
    let e = es.sqrt();
    let qp = super::insert_authalic_setup(&mut params).unwrap_or(2.0);
    let rq = (0.5 * qp).sqrt();
    params.real.insert("rq", rq);

    match (polar, equatorial, north) {
        (true, _, _) => {
            params.real.insert("dd", 1.0);
            params.real.insert("xmf", 1.0);
            params.real.insert("ymf", 1.0);
        }
        (_, true, _) => {
            params.real.insert("dd", rq.recip());
            params.real.insert("xmf", 1.0);
            params.real.insert("ymf", 0.5 * qp);
        }
        _ => {
            let (sin_phi_0, cos_phi_0) = lat_0.sin_cos();
            let xi_0 = (ancillary::qs(sin_phi_0, e) / qp).asin();
            let (sinb1, cosb1) = xi_0.sin_cos();
            let dd = cos_phi_0 / ((1.0 - es * sin_phi_0 * sin_phi_0).sqrt() * rq * cosb1);
            params.real.insert("sinb1", sinb1);
            params.real.insert("cosb1", cosb1);
            params.real.insert("dd", dd);
            params.real.insert("xmf", rq * dd);
            params.real.insert("ymf", rq / dd);
        }
    }
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        state: None,
        steps: None,
    })
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use float_eq::assert_float_eq;

    #[test]
    fn laea_oblique() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        // ETRS-LAEA grid definition
        let op = ctx.op("laea ellps=GRS80 lat_0=52 lon_0=10  x_0=4321000 y_0=3210000")?;

        // The test point from IOGP
        let p = Coor2D::geo(50.0, 5.0);
        let geo = [p];
        let p = Coor2D::raw(3962799.45, 2999718.85);
        let projected = [p];

        let mut operands = geo;

        // Forward
        ctx.apply(op, Fwd, &mut operands)?;
        assert_float_eq!(operands[0].0, projected[0].0, abs_all <= 0.01);
        assert!((operands[0][0] - 3962799.45).abs() < 0.01);
        assert!((operands[0][1] - 2999718.85).abs() < 0.01);
        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][0].to_degrees() - 5.0).abs() < 1e-12);
        assert!((operands[0][1].to_degrees() - 50.).abs() < 1e-12);

        let p = Coor4D::raw(1e30, 1e30, 0.0, 0.0);
        let mut operands = [p];
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0][0].is_nan());

        // Missing test points for the polar aspects

        Ok(())
    }

    // Test the "if rho < EPS10" branch in the inverse case.
    // From this issue: https://github.com/busstoptaktik/geodesy/issues/89
    // reported by @maximkaaa
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
            (
                "laea ellps=6371228,0 lat_0=90 lon_0=180 x_0=0 y_0=0",
                Coor2D::geo(75.0, 160.0),
            ),
        ];

        for (definition, coordinate) in cases {
            let op = ctx.op(definition)?;
            let mut operands = [coordinate];
            let clone = operands;
            ctx.apply(op, Fwd, &mut operands)?;
            ctx.apply(op, Inv, &mut operands)?;
            assert!(
                (clone[0][0] - operands[0][0]).abs() < 1e-11,
                "{definition}: expected lon {}, got {}",
                clone[0][0],
                operands[0][0]
            );
            assert!(
                (clone[0][1] - operands[0][1]).abs() < 1e-11,
                "{definition}: expected lat {}, got {}",
                clone[0][1],
                operands[0][1]
            );
        }

        Ok(())
    }
}
