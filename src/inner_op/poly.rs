//! American Polyconic
use crate::authoring::*;

const TOL: f64 = 1e-10;
const ITOL: f64 = 1e-12;
const I_ITER: usize = 20;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(ml0) = op.params.real("ml0") else {
        return 0;
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;

        let (x, y) = if lat.abs() <= TOL {
            (a * lam + x_0, a * (-ml0) + y_0)
        } else {
            let (sinphi, cosphi) = lat.sin_cos();
            let ms = if cosphi.abs() > TOL {
                ancillary::pj_msfn((sinphi, cosphi), es) / sinphi
            } else {
                0.0
            };
            let el = lam * sinphi;
            (
                a * (ms * el.sin()) + x_0,
                a * ((ellps.meridian_latitude_to_distance(lat) / a - ml0) + ms * (1.0 - el.cos()))
                    + y_0,
            )
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(ml0) = op.params.real("ml0") else {
        return 0;
    };

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let x = (x_raw - x_0) / a;
        let mut y = (y_raw - y_0) / a;
        y += ml0;

        if y.abs() <= TOL {
            operands.set_xy(i, x + lon_0, 0.0);
            successes += 1;
            continue;
        }

        let r = y * y + x * x;
        let mut lat = y;
        let mut converged = false;

        for _ in 0..I_ITER {
            let (sp, cp) = lat.sin_cos();
            if cp.abs() < ITOL || sp.abs() < ITOL {
                break;
            }
            let s2ph = sp * cp;
            let ml = ellps.meridian_latitude_to_distance(lat) / a;
            let root = (1.0 - es * sp * sp).sqrt();
            let c = sp * root / cp;
            let mlb = ml * ml + r;
            let mlp = one_es / (root * root * root);
            let dphi_num = ml + ml + c * mlb - 2.0 * y * (c * ml + 1.0);
            let dphi_den = es * s2ph * (mlb - 2.0 * y * ml) / c
                + 2.0 * (y - ml) * (c * mlp - 1.0 / s2ph)
                - mlp
                - mlp;
            let dphi = dphi_num / dphi_den;
            lat += dphi;
            if dphi.abs() <= ITOL {
                converged = true;
                break;
            }
        }

        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        let sp = lat.sin();
        let lon = if sp.abs() <= TOL {
            x
        } else {
            let arg = x * lat.tan() * (1.0 - es * sp * sp).sqrt() / sp;
            arg.clamp(-1.0, 1.0).asin()
        } + lon_0;

        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let lat_0 = params.lat(0).to_radians();
    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", params.lon(0).to_radians());

    let ellps = params.ellps(0);
    let ml0 = ellps.meridian_latitude_to_distance(lat_0) / ellps.semimajor_axis();
    params.real.insert("ml0", ml0);

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
    fn poly_forward_and_inverse_origin() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("poly lat_0=0 lon_0=-54 x_0=5000000 y_0=10000000 ellps=GRS80")?;

        let geo = [Coor4D::geo(0., -54., 0., 0.)];
        let projected = [Coor4D::raw(5_000_000.0, 10_000_000.0, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn poly_non_origin_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("poly lat_0=0 lon_0=-54 x_0=5000000 y_0=10000000 ellps=GRS80")?;

        let geo = [Coor4D::raw((-50.0_f64).to_radians(), (-10.0_f64).to_radians(), 0., 0.)];
        let projected = [Coor4D::raw(5_438_546.7142, 8_891_486.8987, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-3);

        let mut inverse_probe = projected;
        ctx.apply(op, Inv, &mut inverse_probe)?;
        assert!(inverse_probe[0][0].is_finite());
        assert!(inverse_probe[0][1].is_finite());
        assert!((inverse_probe[0][0].to_degrees() + 50.0).abs() < 1e-2);
        assert!((inverse_probe[0][1].to_degrees() + 10.0).abs() < 1e-6);
        Ok(())
    }
}
