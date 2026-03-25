//! Polar stereographic, implemented from the polar branch of PROJ's `stere`.
use crate::authoring::*;
use std::f64::consts::FRAC_PI_2;

const EPS10: f64 = 1e-10;
const CONV: f64 = 1e-10;
const NITER: usize = 8;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(akm1) = op.params.real("akm1") else {
        return 0;
    };

    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let south_polar = op.params.boolean("south_polar");

    let mut successes = 0_usize;

    for i in 0..operands.len() {
        let (lon, mut lat) = operands.xy(i);
        let lam = lon - lon_0;
        let sinlam = lam.sin();
        let mut coslam = lam.cos();
        let mut sinphi = lat.sin();

        if south_polar {
            lat = -lat;
            coslam = -coslam;
            sinphi = -sinphi;
        }

        let rho = if (lat - FRAC_PI_2).abs() < 1e-15 {
            0.0
        } else {
            akm1 * ancillary::ts((sinphi, lat.cos()), e)
        };

        let x = x_0 + rho * sinlam;
        let y = y_0 - rho * coslam;
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let Ok(akm1) = op.params.real("akm1") else {
        return 0;
    };

    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let south_polar = op.params.boolean("south_polar");

    let halfpi = -FRAC_PI_2;
    let halfe = -0.5 * e;

    let mut successes = 0_usize;

    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let x = x_raw - x_0;
        let mut y = y_raw - y_0;

        if !south_polar {
            y = -y;
        }

        let rho = x.hypot(y);
        let tp = -rho / akm1;
        let mut phi_l = FRAC_PI_2 - 2.0 * tp.atan();
        let mut lat = f64::NAN;

        for _ in 0..NITER {
            let sinphi = e * phi_l.sin();
            lat = 2.0 * (tp * ((1.0 + sinphi) / (1.0 - sinphi)).powf(halfe)).atan() - halfpi;
            if (phi_l - lat).abs() < CONV {
                break;
            }
            phi_l = lat;
        }

        if lat.is_nan() {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        if south_polar {
            lat = -lat;
        }

        let lon = if x == 0.0 && y == 0.0 {
            lon_0
        } else {
            x.atan2(y) + lon_0
        };

        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lat_ts", default: Some(90_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let lat_0 = params.lat(0).to_radians();
    if (lat_0.abs() - FRAC_PI_2).abs() > EPS10 {
        return Err(Error::Unsupported(
            "stere currently supports only polar aspects".to_string(),
        ));
    }

    let mut lat_ts = params.real("lat_ts").unwrap_or(90.0).to_radians();
    if lat_ts.abs() > FRAC_PI_2 + EPS10 {
        return Err(Error::BadParam("lat_ts".to_string(), def.clone()));
    }
    if lat_0.is_sign_negative() {
        params.boolean.insert("south_polar");
        lat_ts = -lat_ts.abs();
    } else {
        params.boolean.insert("north_polar");
        lat_ts = lat_ts.abs();
    }

    params.real.insert("lat_0", lat_0);
    params.real.insert("lat_ts", lat_ts);
    params.real.insert("lon_0", params.lon(0).to_radians());

    let ellps = params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    let k_0 = params.k(0);

    let lat_ts_abs = lat_ts.abs();
    let akm1 = if (lat_ts_abs - FRAC_PI_2).abs() < EPS10 {
        let num = 2.0 * k_0;
        let den = ((1.0 + e).powf(1.0 + e) * (1.0 - e).powf(1.0 - e)).sqrt();
        a * num / den
    } else {
        let sin_ts = lat_ts_abs.sin();
        let factor = lat_ts_abs.cos() / ancillary::ts((sin_ts, lat_ts_abs.cos()), e);
        a * k_0 * factor / (1.0 - (e * sin_ts).powi(2)).sqrt()
    };
    params.real.insert("akm1", akm1);

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
    fn ups_north_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere lat_0=90 lon_0=0 k_0=0.994 x_0=2000000 y_0=2000000 ellps=WGS84")?;

        let geo = [Coor4D::geo(85., 0., 0., 0.)];
        let projected = [Coor4D::raw(2_000_000.0, 1_444_542.608_617_322_5, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-8);
        Ok(())
    }

    #[test]
    fn south_polar_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere lat_0=-90 lat_ts=-71 lon_0=0 ellps=WGS84")?;

        let geo = [Coor4D::geo(-75., 30., 0., 0.)];
        let projected = [Coor4D::raw(
            819_391.619_181_387,
            1_419_227.915_756_798,
            0.,
            0.,
        )];

        let mut operands = projected;
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-6);
        Ok(())
    }

    #[test]
    fn custom_ellipsoid_inverse_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("stere lat_0=90 lat_ts=70 lon_0=-45 ellps=6378273,298.279411123064")?;

        let geo = [Coor4D::geo(80., 45., 0., 0.)];
        let projected = [Coor4D::raw(1_085_943.187_924_962, 0.0, 0., 0.)];

        let mut operands = projected;
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-11);
        Ok(())
    }
}
