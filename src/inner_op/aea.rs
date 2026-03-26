//! Albers Equal Area
use crate::authoring::*;
use std::f64::consts::FRAC_PI_2;

const EPS10: f64 = 1e-10;
const TOL7: f64 = 1e-7;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(n) = op.params.real("n") else { return 0 };
    let Ok(c) = op.params.real("c") else { return 0 };
    let Ok(dd) = op.params.real("dd") else {
        return 0;
    };
    let spherical = op.params.boolean("spherical");
    let Ok(qp) = op.params.real("qp") else {
        return 0;
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let rho_term = if spherical {
            c - 2.0 * n * lat.sin()
        } else {
            c - n * ancillary::qs(lat.sin(), ellps.eccentricity())
        };
        if rho_term < 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let rho = dd * rho_term.sqrt();
        let theta = (lon - lon_0) * n;
        let x = x_0 + a * rho * theta.sin();
        let y = y_0 + a * (op.params.real("rho0").unwrap() - rho * theta.cos());
        let _ = qp; // retained for constructor symmetry/readability
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(n) = op.params.real("n") else { return 0 };
    let Ok(c) = op.params.real("c") else { return 0 };
    let Ok(dd) = op.params.real("dd") else {
        return 0;
    };
    let Ok(rho0) = op.params.real("rho0") else {
        return 0;
    };
    let Ok(qp) = op.params.real("qp") else {
        return 0;
    };
    let Ok(ec) = op.params.real("ec") else {
        return 0;
    };
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let mut x = (operands.xy(i).0 - x_0) / a;
        let mut y = rho0 - (operands.xy(i).1 - y_0) / a;
        let mut rho = x.hypot(y);
        if rho != 0.0 {
            if n < 0.0 {
                rho = -rho;
                x = -x;
                y = -y;
            }

            let lat = if spherical {
                let arg = (c - (rho / dd).powi(2)) / (2.0 * n);
                arg.clamp(-1.0, 1.0).asin()
            } else {
                let Ok(authalic) = op.params.fourier_coefficients("authalic") else {
                    return 0;
                };
                let qs = (c - (rho / dd).powi(2)) / n;
                if (ec - qs.abs()).abs() > TOL7 {
                    if qs.abs() > 2.0 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    let xi = (qs / qp).clamp(-1.0, 1.0).asin();
                    ellps.latitude_authalic_to_geographic(xi, &authalic)
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
            let lon = x.atan2(y) / n + lon_0;
            operands.set_xy(i, lon, lat);
            successes += 1;
        } else {
            let lat = if n > 0.0 { FRAC_PI_2 } else { -FRAC_PI_2 };
            operands.set_xy(i, lon_0, lat);
            successes += 1;
        }
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
    OpParameter::Real { key: "lat_2", default: None },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let phi1 = params.lat(1).to_radians();
    let phi2 = params.lat(2).to_radians();
    let phi0 = params.lat(0).to_radians();
    if phi1.abs() > FRAC_PI_2 || phi2.abs() > FRAC_PI_2 {
        return Err(Error::BadParam("lat_1/lat_2".to_string(), def.clone()));
    }
    if (phi1 + phi2).abs() < EPS10 {
        return Err(Error::General(
            "Aea: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
        ));
    }

    params.real.insert("lat_0", phi0);
    params.real.insert("lat_1", phi1);
    params.real.insert("lat_2", phi2);
    params.real.insert("lon_0", params.lon(0).to_radians());

    let ellps = params.ellps(0);
    let spherical = ellps.flattening() == 0.0;
    if spherical {
        params.boolean.insert("spherical");
    }

    let (sinphi1, cosphi1) = phi1.sin_cos();
    let secant = (phi1 - phi2).abs() >= EPS10;
    let mut n = sinphi1;

    if spherical {
        if secant {
            n = 0.5 * (n + phi2.sin());
        }
        let n2 = n + n;
        let c = cosphi1 * cosphi1 + n2 * sinphi1;
        let dd = 1.0 / n;
        let rho0 = dd * (c - n2 * phi0.sin()).sqrt();
        params.real.insert("n", n);
        params.real.insert("c", c);
        params.real.insert("dd", dd);
        params.real.insert("rho0", rho0);
        params.real.insert("qp", 2.0);
        params.real.insert("ec", 2.0);
    } else {
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
        let qp = ancillary::qs(1.0, e);
        let ec = 1.0 - 0.5 * (1.0 - es) * ((1.0 - e) / (1.0 + e)).ln() / e;
        let c = m1 * m1 + n * q1;
        let dd = 1.0 / n;
        let rho0 = dd * (c - n * ancillary::qs(phi0.sin(), e)).sqrt();
        params.real.insert("n", n);
        params.real.insert("c", c);
        params.real.insert("dd", dd);
        params.real.insert("rho0", rho0);
        params.real.insert("qp", qp);
        params.real.insert("ec", ec);
        let authalic = ellps.coefficients_for_authalic_latitude_computations();
        params.fourier_coefficients.insert("authalic", authalic);
    }

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

        let mut projected = [Coor4D::raw(0.0, -112_982.4091, 0.0, 0.0)];
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
}
