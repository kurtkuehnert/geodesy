//! Geostationary Satellite View
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let h = op.params.real["radius_g_1"];
    let flip_axis = op.params.boolean("flip_axis");
    let radius_p = op.params.real["radius_p"];
    let radius_p2 = op.params.real["radius_p2"];
    let radius_p_inv2 = op.params.real["radius_p_inv2"];
    let ellipsoidal = ellps.flattening() != 0.0;
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let (sin_lam, cos_lam) = lam.sin_cos();
        let (sin_phi, cos_phi) = lat.sin_cos();
        let (vx, vy, vz) = if ellipsoidal {
            let phi_gc = (radius_p2 * lat.tan()).atan();
            let (sin_gc, cos_gc) = phi_gc.sin_cos();
            let r = radius_p / (radius_p * cos_gc).hypot(sin_gc);
            (r * cos_lam * cos_gc, r * sin_lam * cos_gc, r * sin_gc)
        } else {
            (cos_lam * cos_phi, sin_lam * cos_phi, sin_phi)
        };

        if ellipsoidal && ((h + 1.0 - vx) * vx - vy * vy - vz * vz * radius_p_inv2) < 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        let tmp = h + 1.0 - vx;
        let (x, y) = if flip_axis {
            (h * (vy / vz.hypot(tmp)).atan(), h * (vz / tmp).atan())
        } else {
            (h * (vy / tmp).atan(), h * (vz / vy.hypot(tmp)).atan())
        };
        operands.set_xy(i, x_0 + a * x, y_0 + a * y);
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
    let h = op.params.real["radius_g_1"];
    let flip_axis = op.params.boolean("flip_axis");
    let radius_p = op.params.real["radius_p"];
    let radius_p_inv2 = op.params.real["radius_p_inv2"];
    let c = op.params.real["c"];
    let ellipsoidal = ellps.flattening() != 0.0;
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let x = (operands.xy(i).0 - x_0) / a;
        let y = (operands.xy(i).1 - y_0) / a;
        let vx = -1.0;
        let (vy, vz) = if flip_axis {
            let vz = (y / h).tan();
            let vy = (x / h).tan() * (1.0 + vz * vz).sqrt();
            (vy, vz)
        } else {
            let vy = (x / h).tan();
            let vz = (y / h).tan() * (1.0 + vy * vy).sqrt();
            (vy, vz)
        };

        let qa = if ellipsoidal {
            vy * vy + (vz / radius_p) * (vz / radius_p) + vx * vx
        } else {
            vy * vy + vz * vz + vx * vx
        };
        let qb = 2.0 * (h + 1.0) * vx;
        let det = qb * qb - 4.0 * qa * c;
        if det < 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        let k = (-qb - det.sqrt()) / (2.0 * qa);
        let vx = h + 1.0 + k * vx;
        let vy = vy * k;
        let vz = vz * k;
        let lon = vy.atan2(vx) + lon_0;
        let lat = if ellipsoidal {
            (radius_p_inv2 * vz / vx.hypot(vy)).atan()
        } else {
            (vz / vx.hypot(vy)).atan()
        };
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Text { key: "sweep", default: Some("y") },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
    OpParameter::Real { key: "h", default: None },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    let ellps = params.ellps(0);
    let h = params.real("h")? / ellps.semimajor_axis();
    if !(0.0..=1.0e10).contains(&h) || h == 0.0 {
        return Err(Error::BadParam("h".to_string(), def.clone()));
    }
    match params.text("sweep")?.as_str() {
        "x" => {
            params.boolean.insert("flip_axis");
        }
        "y" => {}
        _ => return Err(Error::BadParam("sweep".to_string(), def.clone())),
    }

    let one_es = 1.0 - ellps.eccentricity_squared();
    let radius_g = 1.0 + h;
    params.real.insert("radius_g_1", h);
    params.real.insert("c", radius_g * radius_g - 1.0);
    if ellps.flattening() == 0.0 {
        params.real.insert("radius_p", 1.0);
        params.real.insert("radius_p2", 1.0);
        params.real.insert("radius_p_inv2", 1.0);
    } else {
        params.real.insert("radius_p", one_es.sqrt());
        params.real.insert("radius_p2", one_es);
        params.real.insert("radius_p_inv2", 1.0 / one_es);
    }

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}
