//! Quadrilateralized Spherical Cube
use crate::authoring::*;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, PI};

const EPS10: f64 = 1e-10;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Face {
    Front,
    Right,
    Back,
    Left,
    Top,
    Bottom,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Area {
    Zero,
    One,
    Two,
    Three,
}

fn shift_longitude_origin(longitude: f64, offset: f64) -> f64 {
    let mut slon = longitude + offset;
    if slon < -PI {
        slon += 2.0 * PI;
    } else if slon > PI {
        slon -= 2.0 * PI;
    }
    slon
}

fn equat_face_theta(phi: f64, y: f64, x: f64) -> (f64, Area) {
    if phi < EPS10 {
        return (0.0, Area::Zero);
    }
    let mut theta = y.atan2(x);
    let area = if theta.abs() <= FRAC_PI_4 {
        Area::Zero
    } else if theta > FRAC_PI_4 && theta <= FRAC_PI_2 + FRAC_PI_4 {
        theta -= FRAC_PI_2;
        Area::One
    } else if theta > FRAC_PI_2 + FRAC_PI_4 || theta <= -(FRAC_PI_2 + FRAC_PI_4) {
        theta = if theta >= 0.0 { theta - PI } else { theta + PI };
        Area::Two
    } else {
        theta += FRAC_PI_2;
        Area::Three
    };
    (theta, area)
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let face = match op.params.natural["face"] {
        0 => Face::Front,
        1 => Face::Right,
        2 => Face::Back,
        3 => Face::Left,
        4 => Face::Top,
        _ => Face::Bottom,
    };
    let omf2 = op.params.real["one_minus_f_squared"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, phi_geod) = operands.xy(i);
        let lat = if ellps.flattening() != 0.0 {
            (omf2 * phi_geod.tan()).atan()
        } else {
            phi_geod
        };

        let (theta, phi, area) = match face {
            Face::Top => {
                let phi = FRAC_PI_2 - lat;
                if lon >= FRAC_PI_4 && lon <= FRAC_PI_2 + FRAC_PI_4 {
                    (lon - FRAC_PI_2, phi, Area::Zero)
                } else if lon > FRAC_PI_2 + FRAC_PI_4 || lon <= -(FRAC_PI_2 + FRAC_PI_4) {
                    ((if lon > 0.0 { lon - PI } else { lon + PI }), phi, Area::One)
                } else if lon > -(FRAC_PI_2 + FRAC_PI_4) && lon <= -FRAC_PI_4 {
                    (lon + FRAC_PI_2, phi, Area::Two)
                } else {
                    (lon, phi, Area::Three)
                }
            }
            Face::Bottom => {
                let phi = FRAC_PI_2 + lat;
                if lon >= FRAC_PI_4 && lon <= FRAC_PI_2 + FRAC_PI_4 {
                    (-lon + FRAC_PI_2, phi, Area::Zero)
                } else if lon < FRAC_PI_4 && lon >= -FRAC_PI_4 {
                    (-lon, phi, Area::One)
                } else if lon < -FRAC_PI_4 && lon >= -(FRAC_PI_2 + FRAC_PI_4) {
                    (-lon - FRAC_PI_2, phi, Area::Two)
                } else {
                    ((if lon > 0.0 { -lon + PI } else { -lon - PI }), phi, Area::Three)
                }
            }
            _ => {
                let longitude = match face {
                    Face::Right => shift_longitude_origin(lon, FRAC_PI_2),
                    Face::Back => shift_longitude_origin(lon, PI),
                    Face::Left => shift_longitude_origin(lon, -FRAC_PI_2),
                    _ => lon,
                };
                let (sinlat, coslat) = lat.sin_cos();
                let (sinlon, coslon) = longitude.sin_cos();
                let q = coslat * coslon;
                let r = coslat * sinlon;
                let s = sinlat;
                match face {
                    Face::Front => {
                        let phi = q.acos();
                        let (theta, area) = equat_face_theta(phi, s, r);
                        (theta, phi, area)
                    }
                    Face::Right => {
                        let phi = r.acos();
                        let (theta, area) = equat_face_theta(phi, s, -q);
                        (theta, phi, area)
                    }
                    Face::Back => {
                        let phi = (-q).acos();
                        let (theta, area) = equat_face_theta(phi, s, -r);
                        (theta, phi, area)
                    }
                    Face::Left => {
                        let phi = (-r).acos();
                        let (theta, area) = equat_face_theta(phi, s, q);
                        (theta, phi, area)
                    }
                    _ => unreachable!(),
                }
            }
        };

        let mut mu = ((12.0 / PI) * (theta + (theta.sin() * FRAC_PI_4.cos()).acos() - FRAC_PI_2)).atan();
        let t = ((1.0 - phi.cos()) / (mu.cos() * mu.cos()) / (1.0 - (1.0 / theta.cos()).atan().cos())).sqrt();
        match area {
            Area::One => mu += FRAC_PI_2,
            Area::Two => mu += PI,
            Area::Three => mu += PI + FRAC_PI_2,
            Area::Zero => {}
        }
        operands.set_xy(i, a * t * mu.cos(), a * t * mu.sin());
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let face = match op.params.natural["face"] {
        0 => Face::Front,
        1 => Face::Right,
        2 => Face::Back,
        3 => Face::Left,
        4 => Face::Top,
        _ => Face::Bottom,
    };
    let b = op.params.real["b"];
    let omf = op.params.real["one_minus_f"];
    let omf2 = op.params.real["one_minus_f_squared"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let x = operands.xy(i).0 / a;
        let y = operands.xy(i).1 / a;
        let nu = (x.hypot(y)).atan();
        let mut mu = y.atan2(x);
        let area = if x >= 0.0 && x >= y.abs() {
            Area::Zero
        } else if y >= 0.0 && y >= x.abs() {
            mu -= FRAC_PI_2;
            Area::One
        } else if x < 0.0 && -x >= y.abs() {
            mu = if mu < 0.0 { mu + PI } else { mu - PI };
            Area::Two
        } else {
            mu += FRAC_PI_2;
            Area::Three
        };

        let t = (PI / 12.0) * mu.tan();
        let theta = (t.sin() / (t.cos() - FRAC_1_SQRT_2)).atan();
        let cosmu = mu.cos();
        let tannu = nu.tan();
        let mut cosphi = 1.0 - cosmu * cosmu * tannu * tannu * (1.0 - (1.0 / theta.cos()).atan().cos());
        cosphi = cosphi.clamp(-1.0, 1.0);

        let (mut lon, mut lat) = match face {
            Face::Top => {
                let phi = cosphi.acos();
                let lat = FRAC_PI_2 - phi;
                let lon = match area {
                    Area::Zero => theta + FRAC_PI_2,
                    Area::One => if theta < 0.0 { theta + PI } else { theta - PI },
                    Area::Two => theta - FRAC_PI_2,
                    Area::Three => theta,
                };
                (lon, lat)
            }
            Face::Bottom => {
                let phi = cosphi.acos();
                let lat = phi - FRAC_PI_2;
                let lon = match area {
                    Area::Zero => -theta + FRAC_PI_2,
                    Area::One => -theta,
                    Area::Two => -theta - FRAC_PI_2,
                    Area::Three => if theta < 0.0 { -theta - PI } else { -theta + PI },
                };
                (lon, lat)
            }
            _ => {
                let q = cosphi;
                let mut t = q * q;
                let mut s = if t >= 1.0 { 0.0 } else { (1.0 - t).sqrt() * theta.sin() };
                t += s * s;
                let mut r = if t >= 1.0 { 0.0 } else { (1.0 - t).sqrt() };
                match area {
                    Area::One => {
                        let tmp = r; r = -s; s = tmp;
                    }
                    Area::Two => {
                        r = -r; s = -s;
                    }
                    Area::Three => {
                        let tmp = r; r = s; s = -tmp;
                    }
                    Area::Zero => {}
                }
                match face {
                    Face::Right => {
                        let tmp = q; let q2 = -r; r = tmp; let q = q2; 
                        let lat = (-s).acos() - FRAC_PI_2;
                        let lon = shift_longitude_origin(r.atan2(q), -FRAC_PI_2);
                        (lon, lat)
                    }
                    Face::Back => {
                        let q = -q; let r = -r;
                        let lat = (-s).acos() - FRAC_PI_2;
                        let lon = shift_longitude_origin(r.atan2(q), -PI);
                        (lon, lat)
                    }
                    Face::Left => {
                        let tmp = q; let q2 = r; r = -tmp; let q = q2;
                        let lat = (-s).acos() - FRAC_PI_2;
                        let lon = shift_longitude_origin(r.atan2(q), FRAC_PI_2);
                        (lon, lat)
                    }
                    Face::Front => {
                        let lat = (-s).acos() - FRAC_PI_2;
                        let lon = r.atan2(q);
                        (lon, lat)
                    }
                    _ => unreachable!(),
                }
            }
        };

        if ellps.flattening() != 0.0 {
            let invert = lat < 0.0;
            let tanphi = lat.tan();
            let xa = b / (tanphi * tanphi + omf2).sqrt();
            lat = ((a * a - xa * xa).sqrt() / (omf * xa)).atan();
            if invert {
                lat = -lat;
            }
        }
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 4] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let given = parameters.instantiated_as.split_into_parameters();
    super::override_ellps_from_proj_params(&mut params, def, &given)?;
    let ellps = params.ellps(0);
    let lat_0 = params.lat(0).to_radians();
    let lon_0 = params.lon(0).to_radians();
    let face = if lat_0 >= FRAC_PI_2 - FRAC_PI_4 / 2.0 {
        Face::Top
    } else if lat_0 <= -(FRAC_PI_2 - FRAC_PI_4 / 2.0) {
        Face::Bottom
    } else if lon_0.abs() <= FRAC_PI_4 {
        Face::Front
    } else if lon_0.abs() <= FRAC_PI_2 + FRAC_PI_4 {
        if lon_0 > 0.0 { Face::Right } else { Face::Left }
    } else {
        Face::Back
    };
    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", lon_0);
    params.natural.insert("face", match face {
        Face::Front => 0,
        Face::Right => 1,
        Face::Back => 2,
        Face::Left => 3,
        Face::Top => 4,
        Face::Bottom => 5,
    });
    if ellps.flattening() != 0.0 {
        let b = ellps.semimajor_axis() * (1.0 - ellps.flattening());
        let omf = 1.0 - ellps.flattening();
        params.real.insert("b", b);
        params.real.insert("one_minus_f", omf);
        params.real.insert("one_minus_f_squared", omf * omf);
    } else {
        params.real.insert("b", ellps.semimajor_axis());
        params.real.insert("one_minus_f", 1.0);
        params.real.insert("one_minus_f_squared", 1.0);
    }
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op { descriptor, params, steps: None })
}
