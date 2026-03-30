//! Azimuthal Equidistant
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::{FRAC_PI_2, PI};

const EPS10: f64 = 1e-10;
const TOL: f64 = 1e-10;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let origin = Coor4D::raw(frame.lon_0, frame.lat_0, 0.0, 0.0);
    let spherical = op.params.boolean("spherical");
    let guam = op.params.boolean("guam");
    let north_polar = op.params.boolean("north_polar");
    let _south_polar = op.params.boolean("south_polar");
    let equatorial = op.params.boolean("equatorial");
    let oblique = op.params.boolean("oblique");
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = frame.lon_delta(lon);

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
                    (frame.a * k * cosphi * sinlam, frame.a * k * sinphi)
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
                        frame.a * k * cosphi * sinlam,
                        frame.a * k * (cosph0 * sinphi - sinph0 * cosphi_x_coslam),
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
                let rho = frame.a * (FRAC_PI_2 + phi);
                (rho * sinlam, rho * coslam)
            }
        } else if guam {
            let es = ellps.eccentricity_squared();
            let cosphi = lat.cos();
            let sinphi = lat.sin();
            let t = 1.0 / (1.0 - es * sinphi * sinphi).sqrt();
            (
                frame.a * lam * cosphi * t,
                ellps.meridian_latitude_to_distance(lat)
                    - ellps.meridian_latitude_to_distance(frame.lat_0)
                    + 0.5 * frame.a * lam * lam * cosphi * sinphi * t,
            )
        } else if (equatorial || oblique) && lam.abs() < EPS10 && (lat - frame.lat_0).abs() < EPS10 {
            (0.0, 0.0)
        } else if equatorial && lat.abs() < EPS10 && frame.lat_0.abs() < EPS10 {
            (frame.a * lam, 0.0)
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

        let (x, y) = frame.apply_false_origin(x, y);
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let origin = Coor4D::raw(frame.lon_0, frame.lat_0, 0.0, 0.0);
    let spherical = op.params.boolean("spherical");
    let guam = op.params.boolean("guam");
    let north_polar = op.params.boolean("north_polar");
    let _south_polar = op.params.boolean("south_polar");
    let equatorial = op.params.boolean("equatorial");
    let oblique = op.params.boolean("oblique");
    let sinph0 = op.params.real("sinph0").unwrap_or(0.0);
    let cosph0 = op.params.real("cosph0").unwrap_or(0.0);

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (dx, dy) = frame.remove_false_origin(operands.xy(i).0, operands.xy(i).1);

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
                (0.0, frame.lat_0)
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
        } else if guam {
            let es = ellps.eccentricity_squared();
            let x_norm = dx / frame.a;
            let y_norm = dy / frame.a;
            let x2 = 0.5 * x_norm * x_norm;
            let m0 = ellps.meridian_latitude_to_distance(frame.lat_0);
            let mut lat = frame.lat_0;
            let mut t;
            for _ in 0..3 {
                t = (1.0 - es * lat.sin().powi(2)).sqrt();
                lat = ellps.meridian_distance_to_latitude(
                    m0 + frame.a * (y_norm - x2 * lat.tan() * t),
                );
            }
            t = (1.0 - es * lat.sin().powi(2)).sqrt();
            (frame.lon_0 + x_norm * t / lat.cos(), lat)
        } else {
            let distance = dx.hypot(dy);
            if distance < EPS10 {
                (frame.lon_0, frame.lat_0)
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
            let lon = angular::normalize_symmetric(frame.lon_0 + lon);
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
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "guam" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let lat_0 = params.lat(0);
    if lat_0.abs() > FRAC_PI_2 + EPS10 {
        return Err(Error::BadParam("lat_0".to_string(), def.clone()));
    }

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
        state: None,
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
        assert!(matches!(
            ctx.op("aeqd R=1 lat_0=91"),
            Err(Error::BadParam(_, _))
        ));
    }

    #[test]
    fn aeqd_guam_matches_proj_gie_example() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "aeqd guam ellps=clrk66 x_0=50000 y_0=50000 lon_0=144.74875069444445 lat_0=13.47246633333333",
        )?;

        let mut operands = [Coor4D::geo(13.33903846111111, 144.63533129166666, 0.0, 0.0)];
        ctx.apply(op, Fwd, &mut operands)?;
        assert!((operands[0][0] - 37712.48).abs() < 0.01);
        assert!((operands[0][1] - 35242.00).abs() < 0.01);

        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][0].to_degrees() - 144.63533129166666).abs() < 1e-8);
        assert!((operands[0][1].to_degrees() - 13.33903846111111).abs() < 1e-8);
        Ok(())
    }
}
