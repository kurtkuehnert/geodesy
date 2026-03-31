//! Bonne projection
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::FRAC_PI_2;

const EPS10: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let es = ellps.eccentricity_squared();
    let frame = ProjectionFrame::from_params(&op.params);
    let phi1 = op.params.real["phi1"];
    let am1 = op.params.real["am1"];
    let m1 = op.params.real["m1"];
    let cphi1 = op.params.real["cphi1"];
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = frame.remove_central_meridian(lon);
        let (x, y) = if spherical {
            let rh = cphi1 + phi1 - lat;
            if rh.abs() > EPS10 {
                let e = lam * lat.cos() / rh;
                (
                    frame.a * rh * e.sin() + frame.x_0,
                    frame.a * (cphi1 - rh * e.cos()) + frame.y_0,
                )
            } else {
                (frame.x_0, frame.y_0)
            }
        } else {
            let (sinphi, cosphi) = lat.sin_cos();
            let rh = am1 + m1 - ellps.meridian_latitude_to_distance(lat) / frame.a;
            if rh.abs() > EPS10 {
                let e = cosphi * lam / (rh * (1.0 - es * sinphi * sinphi).sqrt());
                (
                    frame.a * rh * e.sin() + frame.x_0,
                    frame.a * (am1 - rh * e.cos()) + frame.y_0,
                )
            } else {
                (frame.x_0, frame.y_0)
            }
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let es = ellps.eccentricity_squared();
    let frame = ProjectionFrame::from_params(&op.params);
    let phi1 = op.params.real["phi1"];
    let am1 = op.params.real["am1"];
    let m1 = op.params.real["m1"];
    let cphi1 = op.params.real["cphi1"];
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x_raw, y_raw);
        let x = x_local / frame.a;
        let y = y_local / frame.a;

        if spherical {
            let yy = cphi1 - y;
            let rh = yy.hypot(x).copysign(phi1);
            let lat = cphi1 + phi1 - rh;
            let abs_phi = lat.abs();
            if abs_phi > FRAC_PI_2 {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let lon = if FRAC_PI_2 - abs_phi <= EPS10 {
                frame.lon_0
            } else {
                let lm = rh / lat.cos();
                let theta = if phi1 > 0.0 {
                    x.atan2(yy)
                } else {
                    (-x).atan2(-yy)
                };
                frame.lon_0 + lm * theta
            };
            operands.set_xy(i, lon, lat);
            successes += 1;
            continue;
        }

        let yy = am1 - y;
        let rh = yy.hypot(x).copysign(phi1);
        let lat = ellps.meridian_distance_to_latitude(frame.a * (am1 + m1 - rh));
        let abs_phi = lat.abs();
        if abs_phi < FRAC_PI_2 {
            let sinphi = lat.sin();
            let lm = rh * (1.0 - es * sinphi * sinphi).sqrt() / lat.cos();
            let theta = if phi1 > 0.0 {
                x.atan2(yy)
            } else {
                (-x).atan2(-yy)
            };
            operands.set_xy(i, frame.lon_0 + lm * theta, lat);
            successes += 1;
            continue;
        }
        if abs_phi - FRAC_PI_2 <= EPS10 {
            operands.set_xy(i, frame.lon_0, lat);
            successes += 1;
            continue;
        }
        operands.set_coord(i, &Coor4D::nan());
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT_Y0: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "spherical" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_1", default: None },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT_Y0)?;

    let phi1 = params.lat(1);
    if phi1.abs() < EPS10 {
        return Err(Error::General("Bonne: |lat_1| should be > 0"));
    }

    params.real.insert("phi1", phi1);

    let ellps = params.ellps(0);
    let spherical = params.boolean("spherical") || super::mark_spherical(&mut params);
    if spherical {
        let cphi1 = if phi1.abs() + EPS10 >= FRAC_PI_2 {
            0.0
        } else {
            phi1.tan().recip()
        };
        params.real.insert("cphi1", cphi1);
        params.real.insert("am1", 0.0);
        params.real.insert("m1", 0.0);
    } else {
        let sinphi1 = phi1.sin();
        let cosphi1 = phi1.cos();
        let m1 = ellps.meridian_latitude_to_distance(phi1) / ellps.semimajor_axis();
        let am1 =
            cosphi1 / ((1.0 - ellps.eccentricity_squared() * sinphi1 * sinphi1).sqrt() * sinphi1);
        params.real.insert("am1", am1);
        params.real.insert("m1", m1);
        params.real.insert("cphi1", 0.0);
    }

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
    fn bonne_roundtrips_basic_point() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("bonne lat_1=39.6666666666667 lon_0=-8.13190611111111 ellps=bessel")?;

        let mut operands = [Coor4D::geo(40.0, -8.0, 0.0, 0.0)];
        ctx.apply(op, Fwd, &mut operands)?;
        let projected = operands[0];
        assert!(projected[0].is_finite());
        assert!(projected[1].is_finite());

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&Coor4D::geo(40.0, -8.0, 0.0, 0.0)) < 1e-8);
        Ok(())
    }
}
