//! American Polyconic
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const TOL: f64 = 1e-10;
const CONV: f64 = 1e-10;
const ITOL: f64 = 1e-12;
const N_ITER: usize = 10;
const I_ITER: usize = 20;
const INV_TOL: f64 = 1e-12;
const INV_ITER: usize = 15;

#[allow(clippy::too_many_arguments)]
fn fwd_ellipsoidal(
    ellps: &Ellipsoid,
    lon_0: f64,
    x_0: f64,
    y_0: f64,
    ml0: f64,
    lon: f64,
    lat: f64,
    rectifying: Option<&FourierCoefficients>,
) -> (f64, f64) {
    let a = ellps.semimajor_axis();
    let es = ellps.eccentricity_squared();
    let lam = lon - lon_0;

    if lat.abs() <= TOL {
        return (a * lam + x_0, a * (-ml0) + y_0);
    }

    let (sinphi, cosphi) = lat.sin_cos();
    let ms = if cosphi.abs() > TOL {
        ancillary::pj_msfn((sinphi, cosphi), es) / sinphi
    } else {
        0.0
    };
    let el = lam * sinphi;
    let ml = rectifying.map_or_else(
        || ellps.meridian_latitude_to_distance(lat) / a,
        |coefficients| ellps.latitude_geographic_to_rectifying(lat, coefficients),
    );
    (
        a * (ms * el.sin()) + x_0,
        a * ((ml - ml0) + ms * (1.0 - el.cos())) + y_0,
    )
}

fn spherical_inverse_seed(x: f64, y: f64, lat_0: f64, lon_0: f64) -> Option<(f64, f64)> {
    let yy = lat_0 + y;
    if yy.abs() <= TOL {
        return Some((x + lon_0, 0.0));
    }
    let mut lat = yy;
    let b = x * x + y * y;
    for _ in 0..N_ITER {
        let tp = lat.tan();
        let dphi =
            (yy * (lat * tp + 1.0) - lat - 0.5 * (lat * lat + b) * tp) / ((lat - yy) / tp - 1.0);
        lat -= dphi;
        if dphi.abs() <= CONV {
            let lon = if lat.sin().abs() <= TOL {
                x + lon_0
            } else {
                (x * lat.tan()).clamp(-1.0, 1.0).asin() / lat.sin() + lon_0
            };
            return Some((lon, lat));
        }
    }
    None
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let Ok(ml0) = op.params.real("ml0") else {
        return 0;
    };
    let spherical = op.params.real["spherical"] != 0.0;
    let rectifying = op.params.fourier_coefficients.get("rectifying");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = frame.lon_delta_raw(lon);

        let (x, y) = if spherical && lat.abs() > TOL {
            let cot = lat.tan().recip();
            let e = lam * lat.sin();
            (
                frame.a * (e.sin() * cot) + frame.x_0,
                frame.a * (frame.lat_delta(lat) + cot * (1.0 - e.cos())) + frame.y_0,
            )
        } else if lat.abs() <= TOL {
            (frame.a * lam + frame.x_0, frame.a * (-ml0) + frame.y_0)
        } else {
            fwd_ellipsoidal(
                &ellps,
                frame.lon_0,
                frame.x_0,
                frame.y_0,
                ml0,
                lon,
                lat,
                rectifying,
            )
        };
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let frame = ProjectionFrame::from_params(&op.params);
    let es = ellps.eccentricity_squared();
    let one_es = 1.0 - es;
    let Ok(ml0) = op.params.real("ml0") else {
        return 0;
    };
    let spherical = op.params.real["spherical"] != 0.0;
    let rectifying = op.params.fourier_coefficients.get("rectifying");

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x_raw, y_raw) = operands.xy(i);
        let (x_local, y_local) = frame.remove_false_origin(x_raw, y_raw);
        let x = x_local / frame.a;
        let mut y = y_local / frame.a;
        y += ml0;

        if y.abs() <= TOL {
            operands.set_xy(i, x + frame.lon_0, 0.0);
            successes += 1;
            continue;
        }

        if spherical {
            let mut lat = y;
            let b = x * x + y * y;
            let mut converged = false;
            for _ in 0..N_ITER {
                let tp = lat.tan();
                let dphi = (y * (lat * tp + 1.0) - lat - 0.5 * (lat * lat + b) * tp)
                    / ((lat - y) / tp - 1.0);
                lat -= dphi;
                if dphi.abs() <= CONV {
                    converged = true;
                    break;
                }
            }
            if !converged {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
            let lon = frame.lon_0 + (x * lat.tan()).clamp(-1.0, 1.0).asin() / lat.sin();
            operands.set_xy(i, lon, lat);
            successes += 1;
            continue;
        }

        let r = y * y + x * x;
        let mut lat = y;
        let mut converged = false;
        let mut domain_error = false;

        for _ in 0..I_ITER {
            let (sp, cp) = lat.sin_cos();
            if cp.abs() < ITOL {
                domain_error = true;
                break;
            }
            let s2ph = sp * cp;
            let ml = rectifying.map_or_else(
                || ellps.meridian_latitude_to_distance(lat) / frame.a,
                |coefficients| ellps.latitude_geographic_to_rectifying(lat, coefficients),
            );
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

        let analytic_lat = lat.clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        let analytic_lon = if !domain_error && converged {
            let sp = lat.sin();
            if sp.abs() <= TOL {
                x + frame.lon_0
            } else {
                (x * lat.tan() * (1.0 - es * sp * sp).sqrt()).asin() / sp + frame.lon_0
            }
        } else {
            f64::NAN
        };
        let spherical_seed = spherical_inverse_seed(x, y, frame.lat_0, frame.lon_0);
        let (mut lon, mut lat_candidate) = if let Some((seed_lon, seed_lat)) = spherical_seed {
            if analytic_lon.is_finite() {
                let (fx_a, fy_a) = fwd_ellipsoidal(
                    &ellps,
                    frame.lon_0,
                    frame.x_0,
                    frame.y_0,
                    ml0,
                    analytic_lon,
                    analytic_lat,
                    rectifying,
                );
                let (fx_s, fy_s) = fwd_ellipsoidal(
                    &ellps,
                    frame.lon_0,
                    frame.x_0,
                    frame.y_0,
                    ml0,
                    seed_lon,
                    seed_lat,
                    rectifying,
                );
                if (fx_a - x_raw).hypot(fy_a - y_raw) <= (fx_s - x_raw).hypot(fy_s - y_raw) {
                    (analytic_lon, analytic_lat)
                } else {
                    (seed_lon, seed_lat)
                }
            } else {
                (seed_lon, seed_lat)
            }
        } else if analytic_lon.is_finite() {
            (analytic_lon, analytic_lat)
        } else {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        };

        let mut refined = false;
        for _ in 0..INV_ITER {
            let (fx, fy) = fwd_ellipsoidal(
                &ellps,
                frame.lon_0,
                frame.x_0,
                frame.y_0,
                ml0,
                lon,
                lat_candidate,
                rectifying,
            );
            let delta_x = fx - x_raw;
            let delta_y = fy - y_raw;
            if delta_x.abs() < INV_TOL && delta_y.abs() < INV_TOL {
                refined = true;
                break;
            }

            let h_lon = if lon > 0.0 { -1e-6 } else { 1e-6 };
            let h_lat = if lat_candidate > 0.0 { -1e-6 } else { 1e-6 };
            let (fx_lon, fy_lon) = fwd_ellipsoidal(
                &ellps,
                frame.lon_0,
                frame.x_0,
                frame.y_0,
                ml0,
                lon + h_lon,
                lat_candidate,
                rectifying,
            );
            let (fx_lat, fy_lat) = fwd_ellipsoidal(
                &ellps,
                frame.lon_0,
                frame.x_0,
                frame.y_0,
                ml0,
                lon,
                lat_candidate + h_lat,
                rectifying,
            );
            let j11 = (fx_lon - fx) / h_lon;
            let j21 = (fy_lon - fy) / h_lon;
            let j12 = (fx_lat - fx) / h_lat;
            let j22 = (fy_lat - fy) / h_lat;
            let det = j11 * j22 - j12 * j21;
            if det.abs() < 1e-24 {
                break;
            }
            let inv11 = j22 / det;
            let inv12 = -j12 / det;
            let inv21 = -j21 / det;
            let inv22 = j11 / det;
            let delta_lon = (delta_x * inv11 + delta_y * inv12).clamp(-0.3, 0.3);
            let delta_lat = (delta_x * inv21 + delta_y * inv22).clamp(-0.3, 0.3);
            lon = (lon - delta_lon).clamp(-std::f64::consts::PI, std::f64::consts::PI);
            lat_candidate = (lat_candidate - delta_lat)
                .clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        }
        if !refined {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }

        operands.set_xy(i, lon, lat_candidate);
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
    let ellps = params.ellps(0);
    let spherical = super::insert_rectifying_setup(&mut params);
    params.real.insert("spherical", spherical as i32 as f64);
    let frame = ProjectionFrame::from_params(&params);
    if !spherical {
        let rectifying = ellps.coefficients_for_rectifying_latitude_computations();
        let ml0 = ellps.latitude_geographic_to_rectifying(frame.lat_0, &rectifying);
        params.real.insert("ml0", ml0);
    } else {
        let ml0 = ellps.meridian_latitude_to_distance(frame.lat_0) / ellps.semimajor_axis();
        params.real.insert("ml0", ml0);
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

        let geo = [Coor4D::raw(
            (-50.0_f64).to_radians(),
            (-10.0_f64).to_radians(),
            0.,
            0.,
        )];
        let projected = [Coor4D::raw(5_438_546.714_2, 8_891_486.898_7, 0., 0.)];

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

    #[test]
    fn poly_inverse_matches_proj_world_polyconic() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("poly lat_0=0 lon_0=0 ellps=WGS84")?;

        // Reference point generated by PROJ:
        //   echo '30 20 0 0' | cct +proj=poly +lat_0=0 +lon_0=0 +ellps=WGS84
        let mut projected = [Coor4D::raw(3_122_659.251_8, 2_492_720.888_8, 0.0, 0.0)];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() - 30.0).abs() < 1e-8);
        assert!((projected[0][1].to_degrees() - 20.0).abs() < 1e-8);
        Ok(())
    }

    #[test]
    fn poly_inverse_matches_proj_brazil_polyconic_point() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("poly lat_0=0 lon_0=-54 x_0=5000000 y_0=10000000 ellps=aust_SA")?;

        // Reference point generated by PROJ:
        //   echo '-50 -10 0 0' | cct +proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=aust_SA
        let mut projected = [Coor4D::raw(5_438_546.714_2, 8_891_486.898_7, 0.0, 0.0)];
        ctx.apply(op, Inv, &mut projected)?;

        assert!((projected[0][0].to_degrees() + 50.0).abs() < 4e-5);
        assert!((projected[0][1].to_degrees() + 10.0).abs() < 4e-5);
        Ok(())
    }

    #[test]
    fn poly_spherical_roundtrip_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("poly lat_0=0 lon_0=0 R=6371000")?;

        let geo = [Coor4D::raw(
            (-42.22015752663208_f64).to_radians(),
            3.7614226119335115_f64.to_radians(),
            0.0,
            0.0,
        )];
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        ctx.apply(op, Inv, &mut operands)?;

        assert!((operands[0][0] - geo[0][0]).abs() < 1e-9);
        assert!((operands[0][1] - geo[0][1]).abs() < 1e-9);
        Ok(())
    }
}
