//! Cassini-Soldner
use crate::authoring::*;

const C1: f64 = 1.0 / 6.0;
const C2: f64 = 1.0 / 120.0;
const C3: f64 = 1.0 / 24.0;
const C4: f64 = 1.0 / 3.0;
const C5: f64 = 1.0 / 15.0;
const INV_TOL: f64 = 1e-12;
const INV_ITER: usize = 12;

fn fwd_ellipsoidal(
    ellps: &Ellipsoid,
    lon_0: f64,
    x_0: f64,
    y_0: f64,
    m0: f64,
    lon: f64,
    lat: f64,
) -> (f64, f64) {
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let lam = lon - lon_0;
    let (sinphi, cosphi) = lat.sin_cos();
    let m = ellps.meridian_latitude_to_distance(lat);
    let nu_sq = 1.0 / (1.0 - es * sinphi * sinphi);
    let nu = nu_sq.sqrt();
    let tanphi = lat.tan();
    let t = tanphi * tanphi;
    let a1 = lam * cosphi;
    let c = es * cosphi * cosphi / one_es;
    let a2 = a1 * a1;
    (
        x_0 + a * nu * a1 * (1.0 - a2 * t * (C1 + (8.0 - t + 8.0 * c) * a2 * C2)),
        y_0 + (m - m0) + a * nu * tanphi * a2 * (0.5 + (5.0 - t + 6.0 * c) * a2 * C3),
    )
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(m0) = op.params.real("m0") else {
        return 0;
    };
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;

        let (x, y) = if spherical {
            (
                a * (lat.cos() * lam.sin()).asin(),
                a * (lat.tan().atan2(lam.cos()) - op.params.lat(0)),
            )
        } else {
            fwd_ellipsoidal(&ellps, lon_0, x_0, y_0, m0, lon, lat)
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
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let Ok(m0) = op.params.real("m0") else {
        return 0;
    };
    let lat_0 = op.params.lat(0);
    let spherical = op.params.boolean("spherical");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let x = operands.xy(i).0 - x_0;
        let y = operands.xy(i).1 - y_0;

        let (lon, lat) = if spherical {
            let dd = y / a + lat_0;
            (
                lon_0 + (x / a).tan().atan2(dd.cos()),
                (dd.sin() * (x / a).cos()).asin(),
            )
        } else {
            let phi1 = ellps.meridian_distance_to_latitude(m0 + y);
            let tanphi1 = phi1.tan();
            let t1 = tanphi1 * tanphi1;
            let sinphi1 = phi1.sin();
            let nu1_sq = 1.0 / (1.0 - es * sinphi1 * sinphi1);
            let nu1 = nu1_sq.sqrt();
            let rho1 = nu1_sq * (1.0 - es) * nu1;
            let d = x / (a * nu1);
            let d2 = d * d;
            let mut lon =
                lon_0 + d * (1.0 + t1 * d2 * (-C4 + (1.0 + 3.0 * t1) * d2 * C5)) / phi1.cos();
            let mut lat = phi1 - (nu1 * tanphi1 / rho1) * d2 * (0.5 - (1.0 + 3.0 * t1) * d2 * C3);

            for _ in 0..INV_ITER {
                let (fx, fy) = fwd_ellipsoidal(&ellps, lon_0, x_0, y_0, m0, lon, lat);
                let dx = x + x_0 - fx;
                let dy = y + y_0 - fy;
                if dx.hypot(dy) < INV_TOL {
                    break;
                }

                let h_lon = 1e-8;
                let h_lat = 1e-8;
                let (fx_lon, fy_lon) =
                    fwd_ellipsoidal(&ellps, lon_0, x_0, y_0, m0, lon + h_lon, lat);
                let (fx_lat, fy_lat) =
                    fwd_ellipsoidal(&ellps, lon_0, x_0, y_0, m0, lon, lat + h_lat);
                let j11 = (fx_lon - fx) / h_lon;
                let j21 = (fy_lon - fy) / h_lon;
                let j12 = (fx_lat - fx) / h_lat;
                let j22 = (fy_lat - fy) / h_lat;
                let det = j11 * j22 - j12 * j21;
                if det.abs() < 1e-24 {
                    break;
                }
                lon += (dx * j22 - dy * j12) / det;
                lat += (dy * j11 - dx * j21) / det;
            }

            (lon, lat)
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
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
    OpParameter::Flag { key: "hyperbolic" },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let lat_0 = params.lat(0).to_radians();
    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", params.lon(0).to_radians());

    let ellps = params.ellps(0);
    if ellps.flattening() == 0.0 {
        params.boolean.insert("spherical");
        params.real.insert("m0", lat_0);
    } else {
        params
            .real
            .insert("m0", ellps.meridian_latitude_to_distance(lat_0));
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
    fn cass_origin_roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cass lat_0=31.7340969444444 lon_0=35.2120805555556 x_0=170251.555 y_0=126867.909 ellps=6378300.789,293.4663155389811")?;

        let geo = [Coor4D::geo(31.7340969444444, 35.2120805555556, 0., 0.)];
        let projected = [Coor4D::raw(170251.555, 126867.909, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }
}
