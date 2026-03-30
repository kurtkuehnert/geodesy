//! Lambert Conic Near-Conformal
use crate::authoring::*;
use crate::projection::ProjectionFrame;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let phi_0 = frame.lat_0;
    let a_coef = op.params.real["A"];
    let aprime = op.params.real["aprime"];
    let bprime = op.params.real["bprime"];
    let cprime = op.params.real["cprime"];
    let dprime = op.params.real["dprime"];
    let eprime = op.params.real["eprime"];
    let r_0 = op.params.real["r_0"];
    let s_0 = op.params.real["s_0"];
    let sin_phi_0 = phi_0.sin();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, phi) = operands.xy(i);
        let s = series_s(phi, aprime, bprime, cprime, dprime, eprime);
        let m = s - s_0;
        let m_cap = frame.k_0 * (m + a_coef * m * m * m);
        let r = r_0 - m_cap;
        let theta = frame.lon_delta(lon) * sin_phi_0;
        let x = r * theta.sin();
        let y = m_cap + r * theta.sin() * (theta / 2.0).tan();
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let phi_0 = frame.lat_0;
    let a_coef = op.params.real["A"];
    let aprime = op.params.real["aprime"];
    let bprime = op.params.real["bprime"];
    let cprime = op.params.real["cprime"];
    let dprime = op.params.real["dprime"];
    let eprime = op.params.real["eprime"];
    let r_0 = op.params.real["r_0"];
    let s_0 = op.params.real["s_0"];
    let sin_phi_0 = phi_0.sin();

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (dx, dy) = frame.remove_false_origin(x_raw, y_raw);
        let theta = dx.atan2(r_0 - dy);
        let r_prime = dx.hypot(r_0 - dy).copysign(phi_0);
        let m_prime = r_0 - r_prime;

        let denom = -frame.k_0 - 3.0 * frame.k_0 * a_coef * m_prime * m_prime;
        if denom == 0.0 || sin_phi_0 == 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let m = m_prime
            - (m_prime - frame.k_0 * m_prime - frame.k_0 * a_coef * m_prime.powi(3)) / denom;
        let mut phi_prime = phi_0 + (m / aprime).to_radians();
        let s_prime = series_s(phi_prime, aprime, bprime, cprime, dprime, eprime);
        let ds_prime = series_ds(phi_prime, aprime, bprime, cprime, dprime, eprime);
        if ds_prime == 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        phi_prime -= (m + s_0 - s_prime) / -ds_prime;
        let lon = frame.lon_0 + theta / sin_phi_0;
        operands.set_xy(i, lon, phi_prime);
        successes += 1;
    }
    successes
}

fn series_s(phi: f64, aprime: f64, bprime: f64, cprime: f64, dprime: f64, eprime: f64) -> f64 {
    aprime * phi.to_degrees() - bprime * (2.0 * phi).sin() + cprime * (4.0 * phi).sin()
        - dprime * (6.0 * phi).sin()
        + eprime * (8.0 * phi).sin()
}

fn series_ds(phi: f64, aprime: f64, bprime: f64, cprime: f64, dprime: f64, eprime: f64) -> f64 {
    aprime * (180.0 / std::f64::consts::PI) - 2.0 * bprime * (2.0 * phi).cos()
        + 4.0 * cprime * (4.0 * phi).cos()
        - 6.0 * dprime * (6.0 * phi).cos()
        + 8.0 * eprime * (8.0 * phi).cos()
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: None },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let ellps = params.ellps(0);
    let a = ellps.semimajor_axis();
    let f = ellps.flattening();
    let phi_0 = params.lat(0);
    let rho_0 = ellps.meridian_radius_of_curvature(phi_0);
    let nu_0 = ellps.prime_vertical_radius_of_curvature(phi_0);
    let n = f / (2.0 - f);
    let a_coef = 1.0 / (6.0 * rho_0 * nu_0);
    let aprime = a
        * (1.0 - n + 5.0 * (n * n - n.powi(3)) / 4.0 + 81.0 * (n.powi(4) - n.powi(5)) / 64.0)
        * std::f64::consts::PI
        / 180.0;
    let bprime =
        3.0 * a * (n - n * n + 7.0 * (n.powi(3) - n.powi(4)) / 8.0 + 55.0 * n.powi(5) / 64.0) / 2.0;
    let cprime = 15.0 * a * (n * n - n.powi(3) + 3.0 * (n.powi(4) - n.powi(5)) / 4.0) / 16.0;
    let dprime = 35.0 * a * (n.powi(3) - n.powi(4) + 11.0 * n.powi(5) / 16.0) / 48.0;
    let eprime = 315.0 * a * (n.powi(4) - n.powi(5)) / 512.0;
    let r_0 = params.k(0) * nu_0 / phi_0.tan();
    let s_0 = series_s(phi_0, aprime, bprime, cprime, dprime, eprime);

    params.real.insert("A", a_coef);
    params.real.insert("aprime", aprime);
    params.real.insert("bprime", bprime);
    params.real.insert("cprime", cprime);
    params.real.insert("dprime", dprime);
    params.real.insert("eprime", eprime);
    params.real.insert("r_0", r_0);
    params.real.insert("s_0", s_0);

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
    fn lccnc_roundtrips_origin_and_nearby_point() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "lccnc lat_0=34.65 lon_0=37.35 k_0=0.9996256 x_0=300000 y_0=300000 ellps=clrk80ign",
        )?;

        let original = [
            Coor4D::geo(34.65, 37.35, 0.0, 0.0),
            Coor4D::geo(35.0, 37.5, 0.0, 0.0),
        ];
        let mut operands = original;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(
            operands
                .iter()
                .all(|coord| coord[0].is_finite() && coord[1].is_finite())
        );

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&original[i]) < 5e-9);
        }
        Ok(())
    }
}
