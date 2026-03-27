//! Azimuthal Equidistant
use crate::authoring::*;
use std::f64::consts::{FRAC_PI_2, PI};

const EPS10: f64 = 1e-10;
const TOL: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let origin = Coor4D::raw(lon_0, lat_0, 0.0, 0.0);
    let spherical = op.params.boolean("spherical");
    let north_polar = op.params.boolean("north_polar");
    let _south_polar = op.params.boolean("south_polar");
    let equatorial = op.params.boolean("equatorial");
    let oblique = op.params.boolean("oblique");
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = angular::normalize_symmetric(lon - lon_0);

        let (x, y) = if spherical {
            if equatorial {
                let cosphi = lat.cos();
                let sinphi = lat.sin();
                let coslam = lam.cos();
                let sinlam = lam.sin();
                let mut c = cosphi * coslam;
                if (c.abs() - 1.0).abs() < TOL {
                    if c < 0.0 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    (0.0, 0.0)
                } else {
                    c = c.acos();
                    let k = c / c.sin();
                    (a * k * cosphi * sinlam, a * k * sinphi)
                }
            } else if oblique {
                let cosphi = lat.cos();
                let sinphi = lat.sin();
                let coslam = lam.cos();
                let sinlam = lam.sin();
                let cosphi_x_coslam = cosphi * coslam;
                let mut c = sinph0 * sinphi + cosph0 * cosphi_x_coslam;
                if (c.abs() - 1.0).abs() < TOL {
                    if c < 0.0 {
                        operands.set_coord(i, &Coor4D::nan());
                        continue;
                    }
                    (0.0, 0.0)
                } else {
                    c = c.acos();
                    let k = c / c.sin();
                    (
                        a * k * cosphi * sinlam,
                        a * k * (cosph0 * sinphi - sinph0 * cosphi_x_coslam),
                    )
                }
            } else {
                let mut phi = lat;
                let mut coslam = lam.cos();
                let sinlam = lam.sin();
                if north_polar {
                    phi = -phi;
                    coslam = -coslam;
                }
                if (phi - FRAC_PI_2).abs() < EPS10 {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                let rho = a * (FRAC_PI_2 + phi);
                (rho * sinlam, rho * coslam)
            }
        } else if (equatorial || oblique) && lam.abs() < EPS10 && (lat - lat_0).abs() < EPS10 {
            (0.0, 0.0)
        } else if equatorial && lat.abs() < EPS10 && lat_0.abs() < EPS10 {
            (a * lam, 0.0)
        } else {
            let target = Coor4D::raw(lon, lat, 0.0, 0.0);
            let inv = ellps.geodesic_inv(&origin, &target);
            let azimuth = inv[0];
            let distance = inv[2];
            if !azimuth.is_finite() || !distance.is_finite() {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            (distance * azimuth.sin(), distance * azimuth.cos())
        };

        operands.set_xy(i, x_0 + x, y_0 + y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let lon_0 = op.params.lon(0);
    let lat_0 = op.params.lat(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let origin = Coor4D::raw(lon_0, lat_0, 0.0, 0.0);
    let spherical = op.params.boolean("spherical");
    let north_polar = op.params.boolean("north_polar");
    let south_polar = op.params.boolean("south_polar");
    let equatorial = op.params.boolean("equatorial");
    let oblique = op.params.boolean("oblique");
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let dx = operands.xy(i).0 - x_0;
        let dy = operands.xy(i).1 - y_0;

        let (lon, lat) = if spherical {
            let mut c_rh = dx.hypot(dy) / ellps.semimajor_axis();
            if c_rh > PI {
                if c_rh - EPS10 > PI {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                c_rh = PI;
            }
            if c_rh < EPS10 {
                (0.0, lat_0)
            } else if equatorial || oblique {
                let sinc = c_rh.sin();
                let cosc = c_rh.cos();
                if equatorial {
                    let lat = (dy * sinc / (dx.hypot(dy))).asin();
                    let x = dx * sinc;
                    let y = cosc * c_rh;
                    let lon = if y == 0.0 { 0.0 } else { x.atan2(y) };
                    (lon, lat)
                } else {
                    let rho = dx.hypot(dy);
                    let lat = (cosc * sinph0 + dy * sinc * cosph0 / rho).asin();
                    let y = (cosc - sinph0 * lat.sin()) * c_rh;
                    let x = dx * sinc * cosph0;
                    let lon = if y == 0.0 { 0.0 } else { x.atan2(y) };
                    (lon, lat)
                }
            } else if north_polar {
                (dx.atan2(-dy), FRAC_PI_2 - c_rh)
            } else {
                (dx.atan2(dy), c_rh - FRAC_PI_2)
            }
        } else {
            let distance = dx.hypot(dy);
            if distance < EPS10 {
                (lon_0, lat_0)
            } else {
                let azimuth = dx.atan2(dy);
                let dest = ellps.geodesic_fwd(&origin, azimuth, distance);
                if !dest[0].is_finite() || !dest[1].is_finite() {
                    operands.set_coord(i, &Coor4D::nan());
                    continue;
                }
                (dest[0], dest[1])
            }
        };

        let lon = if spherical {
            let lon = angular::normalize_symmetric(lon_0 + lon);
            if (lon + PI).abs() < 1e-12 { PI } else { lon }
        } else {
            lon
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
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let lat_0 = params.lat(0).to_radians();
    let lon_0 = params.lon(0).to_radians();
    if lat_0.abs() > FRAC_PI_2 + EPS10 {
        return Err(Error::BadParam("lat_0".to_string(), def.clone()));
    }
    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", lon_0);

    let ellps = params.ellps(0);
    if ellps.flattening() == 0.0 {
        params.boolean.insert("spherical");
        if (lat_0.abs() - FRAC_PI_2).abs() < EPS10 {
            if lat_0.is_sign_negative() {
                params.boolean.insert("south_polar");
            } else {
                params.boolean.insert("north_polar");
            }
        } else if lat_0.abs() < EPS10 {
            params.boolean.insert("equatorial");
        } else {
            params.boolean.insert("oblique");
            params.real.insert("sinph0", lat_0.sin());
            params.real.insert("cosph0", lat_0.cos());
        }
    } else if lat_0.abs() < EPS10 {
        params.boolean.insert("equatorial");
    } else if (lat_0.abs() - FRAC_PI_2).abs() < EPS10 {
        if lat_0.is_sign_negative() {
            params.boolean.insert("south_polar");
        } else {
            params.boolean.insert("north_polar");
        }
    } else {
        params.boolean.insert("oblique");
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
    fn aeqd_roundtrip_origin() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aeqd lat_0=52 lon_0=-97.5 x_0=8264722.177 y_0=4867518.353 ellps=WGS84")?;

        let geo = [Coor4D::geo(52., -97.5, 0., 0.)];
        let projected = [Coor4D::raw(8_264_722.177, 4_867_518.353, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn aeqd_forward_reference() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("aeqd lat_0=52 lon_0=-97.5 x_0=8264722.177 y_0=4867518.353 ellps=WGS84")?;

        let geo = [Coor4D::geo(61.390407158824, -101.971128034161, 0., 0.)];
        let projected = [Coor4D::raw(8_024_875.974_4, 5_920_866.114_0, 0., 0.)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 1e-3);
        Ok(())
    }

    #[test]
    fn aeqd_rejects_invalid_lat_0() {
        let mut ctx = Minimal::default();
        assert!(matches!(ctx.op("aeqd R=1 lat_0=91"), Err(Error::BadParam(_, _))));
    }
}
