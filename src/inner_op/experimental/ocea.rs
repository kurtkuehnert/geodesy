//! Oblique Cylindrical Equal Area
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_2, PI};

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let k_0 = op.params.k(0);
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let rok = 1.0 / k_0;
    let rtk = k_0;
    let sinphi = op.params.real["sinphi_p"];
    let cosphi = op.params.real["cosphi_p"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let sin_lam = lam.sin();
        let cos_lam = lam.cos();
        let mut x = ((lat.tan() * cosphi + sinphi * sin_lam) / cos_lam).atan();
        if cos_lam < 0.0 {
            x += PI;
        }
        let y = rok * (sinphi * lat.sin() - cosphi * lat.cos() * sin_lam);
        operands.set_xy(i, x_0 + a * rtk * x, y_0 + a * y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let a = op.params.ellps(0).semimajor_axis();
    let k_0 = op.params.k(0);
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let rok = 1.0 / k_0;
    let rtk = k_0;
    let sinphi_p = op.params.real["sinphi_p"];
    let cosphi_p = op.params.real["cosphi_p"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let mut y = (operands.xy(i).1 - y_0) / a;
        let x = (operands.xy(i).0 - x_0) / (a * rtk);
        y /= rok;
        let t = (1.0 - y * y).sqrt();
        let s = x.sin();
        let lat = (y * sinphi_p + t * cosphi_p * s).asin();
        let lon = (t * sinphi_p * s - y * cosphi_p).atan2(t * x.cos()) + lon_0;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    let phi0 = params.lat(0);
    let mut lam_p;
    let phi_p;

    params.real.insert("lat_0", phi0);

    if let Some(alpha) = given.get("alpha") {
        let alpha = PI
            + alpha
                .parse::<f64>()
                .map_err(|_| Error::BadParam("alpha".to_string(), def.clone()))?
                .to_radians();
        let lonc = given
            .get("lonc")
            .map(String::as_str)
            .unwrap_or("0")
            .parse::<f64>()
            .map_err(|_| Error::BadParam("lonc".to_string(), def.clone()))?
            .to_radians();
        lam_p = (-alpha.cos()).atan2(-phi0.sin() * alpha.sin()) + lonc;
        phi_p = (phi0.cos() * alpha.sin()).asin();
    } else {
        let phi1 = given
            .get("lat_1")
            .ok_or_else(|| Error::MissingParam("lat_1".to_string()))?
            .parse::<f64>()
            .map_err(|_| Error::BadParam("lat_1".to_string(), def.clone()))?
            .to_radians();
        let phi2 = given
            .get("lat_2")
            .ok_or_else(|| Error::MissingParam("lat_2".to_string()))?
            .parse::<f64>()
            .map_err(|_| Error::BadParam("lat_2".to_string(), def.clone()))?
            .to_radians();
        let lam1 = given
            .get("lon_1")
            .map(String::as_str)
            .unwrap_or("0")
            .parse::<f64>()
            .map_err(|_| Error::BadParam("lon_1".to_string(), def.clone()))?
            .to_radians();
        let lam2 = given
            .get("lon_2")
            .map(String::as_str)
            .unwrap_or("0")
            .parse::<f64>()
            .map_err(|_| Error::BadParam("lon_2".to_string(), def.clone()))?
            .to_radians();

        lam_p = (phi1.cos() * phi2.sin() * lam1.cos() - phi1.sin() * phi2.cos() * lam2.cos())
            .atan2(phi1.sin() * phi2.cos() * lam2.sin() - phi1.cos() * phi2.sin() * lam1.sin());
        if (lam1 + FRAC_PI_2).abs() < f64::EPSILON {
            lam_p = -lam_p;
        }
        let cos_lamp_m_lam1 = (lam_p - lam1).cos();
        let tan_phi1 = phi1.tan();
        phi_p = if tan_phi1 == 0.0 {
            if cos_lamp_m_lam1 >= 0.0 {
                -FRAC_PI_2
            } else {
                FRAC_PI_2
            }
        } else {
            (-cos_lamp_m_lam1 / tan_phi1).atan()
        };
    }

    params.real.insert("lon_0", lam_p + FRAC_PI_2);
    params.real.insert("sinphi_p", phi_p.sin());
    params.real.insert("cosphi_p", phi_p.cos());

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
    fn ocea_two_point_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("ocea a=6400000 lat_1=0.5 lat_2=2 lon_1=0 lon_2=0")?;
        let geo = [Coor4D::geo(1., 2., 0., 0.)];
        let projected = [Coor4D::raw(
            19_994_423.837_934_088,
            223_322.760_576_728,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-6);
        Ok(())
    }
}
