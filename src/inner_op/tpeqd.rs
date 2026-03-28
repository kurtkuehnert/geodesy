//! Two Point Equidistant
use crate::authoring::*;
use std::f64::consts::FRAC_PI_2;

fn aasin(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).asin()
}

fn aacos(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).acos()
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let sp1 = op.params.real["sp1"];
    let cp1 = op.params.real["cp1"];
    let sp2 = op.params.real["sp2"];
    let cp2 = op.params.real["cp2"];
    let ccs = op.params.real["ccs"];
    let cs = op.params.real["cs"];
    let sc = op.params.real["sc"];
    let dlam2 = op.params.real["dlam2"];
    let r2z0 = op.params.real["r2z0"];
    let z02 = op.params.real["z02"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let sp = lat.sin();
        let cp = lat.cos();
        let dl1 = lon + dlam2;
        let dl2 = lon - dlam2;
        let mut z1 = aacos(sp1 * sp + cp1 * cp * dl1.cos());
        let mut z2 = aacos(sp2 * sp + cp2 * cp * dl2.cos());
        z1 *= z1;
        z2 *= z2;

        let t = z1 - z2;
        let x = r2z0 * t;
        let mut y = r2z0 * (4.0 * z02 * z2 - (z02 - t) * (z02 - t)).sqrt();
        if ccs * sp - cp * (cs * dl1.sin() - sc * dl2.sin()) < 0.0 {
            y = -y;
        }
        operands.set_xy(i, x_0 + a * x, y_0 + a * y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let hz0 = op.params.real["hz0"];
    let thz0 = op.params.real["thz0"];
    let rhshz0 = op.params.real["rhshz0"];
    let ca = op.params.real["ca"];
    let sa = op.params.real["sa"];
    let lp0 = op.params.real["lp0"];
    let lamc = op.params.real["lamc"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let x = (operands.xy(i).0 - x_0) / a;
        let y = (operands.xy(i).1 - y_0) / a;
        let cz1 = (y.hypot(x + hz0)).cos();
        let cz2 = (y.hypot(x - hz0)).cos();
        let s = cz1 + cz2;
        let d = cz1 - cz2;
        let mut lon = -d.atan2(s * thz0);
        let mut lat = aacos((thz0 * s).hypot(d) * rhshz0);
        if y < 0.0 {
            lat = -lat;
        }
        let sp = lat.sin();
        let cp = lat.cos();
        lon -= lp0;
        let s_lon = lon.cos();
        lat = aasin(sa * sp + ca * cp * s_lon);
        lon = (cp * lon.sin()).atan2(sa * cp * s_lon - ca * sp) + lamc;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_1", default: None },
    OpParameter::Real { key: "lon_1", default: Some(0_f64) },
    OpParameter::Real { key: "lat_2", default: None },
    OpParameter::Real { key: "lon_2", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    let phi1 = params
        .real("lat_1")
        .map_err(|_| Error::MissingParam("lat_1".to_string()))?
        .to_radians();
    let lam1 = params.lon(1).to_radians();
    let phi2 = params
        .real("lat_2")
        .map_err(|_| Error::MissingParam("lat_2".to_string()))?
        .to_radians();
    let lam2 = params.lon(2).to_radians();

    if phi1 == phi2 && lam1 == lam2 {
        return Err(Error::General(
            "Tpeqd: Invalid value for lat_1/lon_1/lat_2/lon_2: the 2 points should be distinct.",
        ));
    }
    if phi1.abs() == FRAC_PI_2 && phi2.abs() == FRAC_PI_2 {
        return Err(Error::General(
            "Tpeqd: Invalid value for lat_1 and lat_2: their absolute value should be < 90°.",
        ));
    }

    let lam0 = angular::normalize_symmetric(0.5 * (lam1 + lam2));
    let mut dlam2 = angular::normalize_symmetric(lam2 - lam1);
    let cp1 = phi1.cos();
    let cp2 = phi2.cos();
    let sp1 = phi1.sin();
    let sp2 = phi2.sin();
    let cs = cp1 * sp2;
    let sc = sp1 * cp2;
    let ccs = cp1 * cp2 * dlam2.sin();
    let cs_minus_sc_cos_dlam = cs - sc * dlam2.cos();
    let mut z02 = (cp2 * dlam2.sin())
        .hypot(cs_minus_sc_cos_dlam)
        .atan2(sp1 * sp2 + cp1 * cp2 * dlam2.cos());
    if z02 == 0.0 {
        return Err(Error::General(
            "Tpeqd: Invalid value for lat_1 and lat_2: their absolute value should be < 90°.",
        ));
    }
    let hz0 = 0.5 * z02;
    let a12 = (cp2 * dlam2.sin()).atan2(cs_minus_sc_cos_dlam);
    let pp = aasin(cp1 * a12.sin());
    let ca = pp.cos();
    let sa = pp.sin();
    let lp0 = angular::normalize_symmetric((cp1 * a12.cos()).atan2(sp1) - hz0);
    dlam2 *= 0.5;
    let lamc = FRAC_PI_2 - (a12.sin() * sp1).atan2(a12.cos()) - dlam2;
    let thz0 = hz0.tan();
    let rhshz0 = 0.5 / hz0.sin();
    let r2z0 = 0.5 / z02;
    z02 *= z02;

    params.real.insert("lon_1", lam1);
    params.real.insert("lat_1", phi1);
    params.real.insert("lon_2", lam2);
    params.real.insert("lat_2", phi2);
    params.real.insert("cp1", cp1);
    params.real.insert("cp2", cp2);
    params.real.insert("sp1", sp1);
    params.real.insert("sp2", sp2);
    params.real.insert("cs", cs);
    params.real.insert("sc", sc);
    params.real.insert("ccs", ccs);
    params.real.insert("dlam2", dlam2);
    params.real.insert("hz0", hz0);
    params.real.insert("thz0", thz0);
    params.real.insert("rhshz0", rhshz0);
    params.real.insert("ca", ca);
    params.real.insert("sa", sa);
    params.real.insert("lp0", lp0);
    params.real.insert("lamc", lamc);
    params.real.insert("r2z0", r2z0);
    params.real.insert("z02", z02);

    params.real.insert("lon_0", lam0);

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
    fn tpeqd_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("tpeqd ellps=GRS80 lat_1=0.5 lat_2=2")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(
            -27_750.758_831_679,
            -222_599.403_691_777,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-3);
        Ok(())
    }
}
