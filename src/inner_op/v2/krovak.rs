//! Krovak and Modified Krovak
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const EPS: f64 = 1.0e-15;
const MAX_ITER: usize = 100;
const UQ: f64 = 1.042_168_563_804_74;
const S0: f64 = 1.370_083_462_815_55;

const MOD_X0: f64 = 1_089_000.0;
const MOD_Y0: f64 = 654_000.0;
const MOD_C1: f64 = 2.946_529_277E-02;
const MOD_C2: f64 = 2.515_965_696E-02;
const MOD_C3: f64 = 1.193_845_912E-07;
const MOD_C4: f64 = -4.668_270_147E-07;
const MOD_C5: f64 = 9.233_980_362E-12;
const MOD_C6: f64 = 1.523_735_715E-12;
const MOD_C7: f64 = 1.696_780_024E-18;
const MOD_C8: f64 = 4.408_314_235E-18;
const MOD_C9: f64 = -8.331_083_518E-24;
const MOD_C10: f64 = -3.689_471_323E-24;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let alpha = op.params.real["alpha"];
    let k = op.params.real["k"];
    let n = op.params.real["n"];
    let rho0 = op.params.real["rho0"];
    let ad = op.params.real["ad"];
    let a = op.params.real["a"];
    let e = op.params.real["e"];
    let modified = op.params.boolean("modified");
    let easting_northing =
        !matches!(op.params.text("axis"), Ok(axis) if axis.eq_ignore_ascii_case("swu"));

    let mut successes = 0usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let gfi = ((1.0 + e * lat.sin()) / (1.0 - e * lat.sin())).powf(alpha * e * 0.5);
        let u = 2.0 * (k * (lat * 0.5 + FRAC_PI_4).tan().powf(alpha) / gfi).atan() - FRAC_PI_2;
        let deltav = -frame.remove_central_meridian_raw(lon) * alpha;

        let s = (ad.cos() * u.sin() + ad.sin() * u.cos() * deltav.cos()).asin();
        let cos_s = s.cos();
        if cos_s < 1e-12 {
            operands.set_xy(i, 0.0, 0.0);
            successes += 1;
            continue;
        }
        let d = (u.cos() * deltav.sin() / cos_s).asin();
        let eps = n * d;
        let rho = rho0 * (S0 * 0.5 + FRAC_PI_4).tan().powf(n) / (s * 0.5 + FRAC_PI_4).tan().powf(n);

        let mut southing = a * rho * eps.cos();
        let mut westing = a * rho * eps.sin();

        if modified {
            let xr = southing - MOD_X0;
            let yr = westing - MOD_Y0;
            let (dx, dy) = mod_krovak_dx_dy(xr, yr);
            southing -= dx;
            westing -= dy;
        }

        let (x, y) = if easting_northing {
            (-westing - frame.x_0, -southing - frame.y_0)
        } else {
            (southing + frame.y_0, westing + frame.x_0)
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let alpha = op.params.real["alpha"];
    let k = op.params.real["k"];
    let n = op.params.real["n"];
    let rho0 = op.params.real["rho0"];
    let ad = op.params.real["ad"];
    let a = op.params.real["a"];
    let e = op.params.real["e"];
    let modified = op.params.boolean("modified");
    let easting_northing =
        !matches!(op.params.text("axis"), Ok(axis) if axis.eq_ignore_ascii_case("swu"));

    let mut successes = 0usize;
    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let (mut westing, mut southing) = if easting_northing {
            (-x - frame.x_0, -y - frame.y_0)
        } else {
            (y - frame.x_0, x - frame.y_0)
        };

        if modified {
            let xr = southing - MOD_X0;
            let yr = westing - MOD_Y0;
            let (dx, dy) = mod_krovak_dx_dy(xr, yr);
            southing += dx;
            westing += dy;
        }

        let rho = (southing * southing + westing * westing).sqrt() / a;
        let eps = westing.atan2(southing);
        let d = eps / S0.sin();
        let s = if rho == 0.0 {
            FRAC_PI_2
        } else {
            2.0 * ((rho0 / rho).powf(1.0 / n) * (S0 * 0.5 + FRAC_PI_4).tan()).atan() - FRAC_PI_2
        };

        let u = (ad.cos() * s.sin() - ad.sin() * s.cos() * d.cos()).asin();
        let deltav = (s.cos() * d.sin() / u.cos()).asin();
        let lon = frame.lon_0 - deltav / alpha;

        let mut lat = u;
        let mut prev = u;
        let mut converged = false;
        for _ in 0..MAX_ITER {
            lat = 2.0
                * ((1.0 / k).powf(1.0 / alpha)
                    * (u * 0.5 + FRAC_PI_4).tan().powf(1.0 / alpha)
                    * ((1.0 + e * prev.sin()) / (1.0 - e * prev.sin())).powf(e * 0.5))
                .atan()
                - FRAC_PI_2;
            if (lat - prev).abs() < EPS {
                converged = true;
                break;
            }
            prev = lat;
        }
        if !converged {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        }

        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

fn mod_krovak_dx_dy(xr: f64, yr: f64) -> (f64, f64) {
    let xr2 = xr * xr;
    let yr2 = yr * yr;
    let xr4 = xr2 * xr2;
    let yr4 = yr2 * yr2;
    let dx = MOD_C1 + MOD_C3 * xr - MOD_C4 * yr - 2.0 * MOD_C6 * xr * yr
        + MOD_C5 * (xr2 - yr2)
        + MOD_C7 * xr * (xr2 - 3.0 * yr2)
        - MOD_C8 * yr * (3.0 * xr2 - yr2)
        + 4.0 * MOD_C9 * xr * yr * (xr2 - yr2)
        + MOD_C10 * (xr4 + yr4 - 6.0 * xr2 * yr2);
    let dy = MOD_C2
        + MOD_C3 * yr
        + MOD_C4 * xr
        + 2.0 * MOD_C5 * xr * yr
        + MOD_C6 * (xr2 - yr2)
        + MOD_C8 * xr * (xr2 - 3.0 * yr2)
        + MOD_C7 * yr * (3.0 * xr2 - yr2)
        - 4.0 * MOD_C10 * xr * yr * (xr2 - yr2)
        + MOD_C9 * (xr4 + yr4 - 6.0 * xr2 * yr2);
    (dx, dy)
}

#[rustfmt::skip]
const GAMUT: [OpParameter; 10] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "modified" },
    OpParameter::Text { key: "ellps", default: Some("bessel") },
    OpParameter::Text { key: "axis", default: Some("enu") },
    OpParameter::Real { key: "lat_0", default: Some(49.5) },
    OpParameter::Real { key: "lon_0", default: Some(24.8333333333333) },
    OpParameter::Real { key: "alpha", default: Some(30.2881397527778) },
    OpParameter::Real { key: "k", default: Some(0.9999) },
    OpParameter::Real { key: "x_0", default: Some(0.0) },
    OpParameter::Real { key: "y_0", default: Some(0.0) },
];

fn new_inner(parameters: &RawParameters, modified: bool) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    if modified {
        params.boolean.insert("modified");
    }

    // Match PROJ parity: this projection family is hardwired to the Bessel ellipsoid.
    let ellps = Ellipsoid::named("bessel")?;
    let es = ellps.eccentricity_squared();
    let e = ellps.eccentricity();
    let frame = ProjectionFrame::from_params(&params);
    let phi_0 = frame.lat_0;
    let k_0 = params.real("k")?;

    let alpha = (1.0 + (es * phi_0.cos().powi(4)) / (1.0 - es)).sqrt();
    let u0 = (phi_0.sin() / alpha).asin();
    let g = ((1.0 + e * phi_0.sin()) / (1.0 - e * phi_0.sin())).powf(alpha * e * 0.5);
    let tan_half = (phi_0 * 0.5 + FRAC_PI_4).tan();
    if tan_half == 0.0 {
        return Err(Error::General(
            "Krovak: Invalid value for lat_0: lat_0 + PI/4 should be different from 0",
        ));
    }
    let k = (u0 * 0.5 + FRAC_PI_4).tan() / tan_half.powf(alpha) * g;
    let n0 = (1.0 - es).sqrt() / (1.0 - es * phi_0.sin().powi(2));
    let n = S0.sin();
    let rho0 = k_0 * n0 / S0.tan();
    let ad = FRAC_PI_2 - UQ;

    params.real.insert("alpha", alpha);
    params.real.insert("k", k);
    params.real.insert("n", n);
    params.real.insert("rho0", rho0);
    params.real.insert("ad", ad);
    params.real.insert("lat_0", phi_0);
    params.real.insert("a", ellps.semimajor_axis());
    params.real.insert("e", e);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        state: None,
        steps: None,
    })
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    new_inner(parameters, false)
}

pub fn modified(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    new_inner(parameters, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn krovak_forward_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "krovak lat_0=49.5 lon_0=24.8333333333333 alpha=30.2881397527778 k=0.9999 x_0=0 y_0=0 ellps=bessel",
        )?;

        let geo = [Coor4D::geo(50.0833, 14.417, 0., 0.)];
        let expected = [Coor4D::raw(
            -743_263.655_907_14,
            -1_043_505.836_542_498_3,
            0.,
            0.,
        )];
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-6);
        Ok(())
    }

    #[test]
    fn krovak_axis_swu_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "krovak axis=swu lat_0=49.5 lon_0=24.8333333333333 alpha=30.2881397527778 k=0.9999 x_0=0 y_0=0 ellps=bessel",
        )?;

        let geo = [Coor4D::geo(50.0833, 14.417, 0., 0.)];
        let expected = [Coor4D::raw(
            1_043_505.836_542_498_3,
            743_263.655_907_14,
            0.,
            0.,
        )];
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-6);
        Ok(())
    }

    #[test]
    fn mod_krovak_roundtrip_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "mod_krovak lat_0=49.5 lon_0=24.8333333333333 alpha=30.2881397222222 k=0.9999 x_0=5000000 y_0=5000000 ellps=bessel",
        )?;

        let geo = [Coor4D::geo(50.0833, 14.417, 0., 0.)];
        let expected = [Coor4D::raw(
            -5_743_263.678_112_479,
            -6_043_505.811_232_34,
            0.,
            0.,
        )];
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
