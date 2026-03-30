//! Colombia Urban
use crate::authoring::*;
use crate::projection::ProjectionFrame;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let es = op.params.real["es"];
    let a = op.params.real["a"];
    let a_coeff = op.params.real["A"];
    let b_coeff = op.params.real["B"];
    let rho0 = op.params.real["rho0"];
    let h0_over_a = op.params.real["h0_over_a"];

    let mut successes = 0usize;
    for i in 0..operands.len() {
        let (lam, phi) = operands.xy(i);
        let lam = frame.lon_delta_raw(lam);
        let cosphi = phi.cos();
        let sinphi = phi.sin();
        let nu = 1.0 / (1.0 - es * sinphi * sinphi).sqrt();
        let lam_nu_cosphi = lam * nu * cosphi;

        let sinphi_m = (0.5 * (phi + frame.lat_0)).sin();
        let rho_m = (1.0 - es) / (1.0 - es * sinphi_m * sinphi_m).powf(1.5);
        let g = 1.0 + h0_over_a / rho_m;

        let x = a * a_coeff * lam_nu_cosphi;
        let y = a * g * rho0 * (frame.lat_delta(phi) + b_coeff * lam_nu_cosphi * lam_nu_cosphi);
        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let es = op.params.real["es"];
    let a = op.params.real["a"];
    let b_coeff = op.params.real["B"];
    let c_coeff = op.params.real["C"];
    let d_coeff = op.params.real["D"];

    let mut successes = 0usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x, y) = frame.remove_false_origin(x_raw, y_raw);

        let phi = frame.lat_0 + y / (a * d_coeff) - b_coeff * (x / (a * c_coeff)).powi(2);
        let sinphi = phi.sin();
        let nu = 1.0 / (1.0 - es * sinphi * sinphi).sqrt();
        let lam = frame.lon_0 + x / (a * c_coeff * nu * phi.cos());
        operands.set_xy(i, lam, phi);
        successes += 1;
    }

    successes
}

#[rustfmt::skip]
const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0.0) },
    OpParameter::Real { key: "lon_0", default: Some(0.0) },
    OpParameter::Real { key: "x_0", default: Some(0.0) },
    OpParameter::Real { key: "y_0", default: Some(0.0) },
    OpParameter::Real { key: "h_0", default: Some(0.0) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let ellps = params.ellps(0);
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let phi0 = params.lat(0);
    let h0 = params.real("h_0")?;
    let h0_over_a = h0 / a;

    let sinphi0 = phi0.sin();
    let nu0 = 1.0 / (1.0 - es * sinphi0 * sinphi0).sqrt();
    let rho0 = (1.0 - es) / (1.0 - es * sinphi0 * sinphi0).powf(1.5);

    params.real.insert("a", a);
    params.real.insert("es", es);
    params.real.insert("h0_over_a", h0_over_a);
    params.real.insert("A", 1.0 + h0_over_a / nu0);
    params.real.insert("rho0", rho0);
    params.real.insert("B", phi0.tan() / (2.0 * rho0 * nu0));
    params.real.insert("C", 1.0 + h0_over_a);
    params
        .real
        .insert("D", rho0 * (1.0 + h0_over_a / (1.0 - es)));

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
    fn col_urban_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "col_urban lat_0=4.68048611111111 lon_0=-74.1465916666667 x_0=92334.879 y_0=109320.965 h_0=2550 ellps=GRS80",
        )?;

        let geo = [Coor4D::geo(4.8, -74.25, 0.0, 0.0)];
        let projected = [Coor4D::raw(80859.033, 122543.174, 0.0, 0.0)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-3);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
