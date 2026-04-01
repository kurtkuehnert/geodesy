//! Equidistant Conic
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::FRAC_PI_2;

const EPS10: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let n = op.params.real["n"];
    let c = op.params.real["c"];
    let rho0 = op.params.real["rho0"];
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let rho = if spherical {
            c - lat
        } else {
            c - ellps.meridian_latitude_to_distance(lat) / frame.a
        };
        let theta = n * frame.remove_central_meridian(lon);
        let x = frame.a * rho * theta.sin();
        let y = frame.a * (rho0 - rho * theta.cos());
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let n = op.params.real["n"];
    let c = op.params.real["c"];
    let rho0 = op.params.real["rho0"];
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (mut x, y_local) = frame.remove_false_origin(x_raw, y_raw);
        x /= frame.a;
        let mut y = rho0 - y_local / frame.a;
        let mut rho = x.hypot(y);
        let (lon, lat) = if rho != 0.0 {
            if n < 0.0 {
                rho = -rho;
                x = -x;
                y = -y;
            }
            let lat = if spherical {
                c - rho
            } else {
                ellps.meridian_distance_to_latitude(frame.a * (c - rho))
            };
            (x.atan2(y) / n + frame.lon_0, lat)
        } else {
            (frame.lon_0, if n > 0.0 { FRAC_PI_2 } else { -FRAC_PI_2 })
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
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "lat_1", default: None },
    OpParameter::Real { key: "lat_2", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let phi0 = params.lat(0);
    let phi1 = params.lat(1);
    let phi2 = params.lat(2);

    if phi1.abs() > FRAC_PI_2 {
        return Err(Error::General(
            "Eqdc: Invalid value for lat_1: |lat_1| should be <= 90°",
        ));
    }
    if phi2.abs() > FRAC_PI_2 {
        return Err(Error::General(
            "Eqdc: Invalid value for lat_2: |lat_2| should be <= 90°",
        ));
    }
    if (phi1 + phi2).abs() < EPS10 {
        return Err(Error::General(
            "Eqdc: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
        ));
    }

    params.real.insert("lat_1", phi1);
    params.real.insert("lat_2", phi2);

    let ellps = params.ellps(0);
    let spherical = crate::inner_op::mark_spherical(&mut params);

    let secant = (phi1 - phi2).abs() >= EPS10;
    let (sinphi1, cosphi1) = phi1.sin_cos();
    let mut n = sinphi1;
    let c;
    let rho0;
    if spherical {
        if secant {
            n = (cosphi1 - phi2.cos()) / (phi2 - phi1);
        }
        if n == 0.0 {
            return Err(Error::General(
                "Eqdc: Invalid value for lat_1 and lat_2: lat_1 + lat_2 should be > 0",
            ));
        }
        c = phi1 + cosphi1 / n;
        rho0 = c - phi0;
    } else {
        let es = ellps.eccentricity_squared();
        let m1 = ancillary::pj_msfn((sinphi1, cosphi1), es);
        let ml1 = ellps.meridian_latitude_to_distance(phi1) / ellps.semimajor_axis();
        if secant {
            let (sinphi2, cosphi2) = phi2.sin_cos();
            let ml2 = ellps.meridian_latitude_to_distance(phi2) / ellps.semimajor_axis();
            if ml1 == ml2 {
                return Err(Error::General("Eqdc: Eccentricity too close to 1"));
            }
            n = (m1 - ancillary::pj_msfn((sinphi2, cosphi2), es)) / (ml2 - ml1);
            if n == 0.0 {
                return Err(Error::General("Eqdc: Invalid value for eccentricity"));
            }
        }
        c = ml1 + m1 / n;
        rho0 = c - ellps.meridian_latitude_to_distance(phi0) / ellps.semimajor_axis();
    }

    params.real.insert("n", n);
    params.real.insert("c", c);
    params.real.insert("rho0", rho0);

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
    fn eqdc_roundtrip_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("eqdc ellps=GRS80 lat_1=0.5 lat_2=2")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(
            222_588.440_269_286,
            110_659.134_907_347,
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
}
