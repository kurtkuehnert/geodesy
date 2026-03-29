//! Sinusoidal (Sanson-Flamsteed)
use crate::authoring::*;

const EPS10: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let (x, y) = if spherical {
            (lam * lat.cos(), lat)
        } else {
            let s = lat.sin();
            let c = lat.cos();
            (
                lam * c / (1.0 - es * s * s).sqrt(),
                ellps.meridian_latitude_to_distance(lat) / a,
            )
        };
        operands.set_xy(i, x_0 + a * x, y_0 + a * y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let x = (x_raw - x_0) / a;
        let y = (y_raw - y_0) / a;

        let (lon, lat) = if spherical {
            let lat = y;
            let s = lat.abs();
            if s > std::f64::consts::FRAC_PI_2 + EPS10 {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let lon = if (std::f64::consts::FRAC_PI_2 - s).abs() <= EPS10 {
                lon_0
            } else {
                lon_0 + x / lat.cos()
            };
            (lon, lat)
        } else {
            let lat = ellps.meridian_distance_to_latitude(y * a);
            let s = lat.abs();
            if s < std::f64::consts::FRAC_PI_2 {
                let sinphi = lat.sin();
                let lon = lon_0 + x * (1.0 - es * sinphi * sinphi).sqrt() / lat.cos();
                (lon, lat)
            } else if (s - EPS10) < std::f64::consts::FRAC_PI_2 {
                (lon_0, lat)
            } else {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
        };

        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    params.real.insert("lon_0", params.lon(0).to_radians());

    if params.ellps(0).flattening() == 0.0 {
        params.boolean.insert("spherical");
    }

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op { descriptor, params, steps: None })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sinu_matches_proj_gie_ellipsoidal_case() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("sinu ellps=GRS80")?;
        let geo = [Coor4D::geo(1.0, 2.0, 0., 0.)];
        let projected = [Coor4D::raw(222_605.299_539_466, 110_574.388_554_153, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-6);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
