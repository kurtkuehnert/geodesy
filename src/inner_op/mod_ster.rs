//! Modified stereographic variants used by PROJ's `alsk`, `gs48`, and `gs50`.
use crate::authoring::*;
use std::collections::BTreeMap;
use std::f64::consts::FRAC_PI_2;

const EPSLN: f64 = 1e-12;

#[derive(Clone, Copy)]
struct Complex {
    r: f64,
    i: f64,
}

impl Complex {
    fn add(self, rhs: Self) -> Self {
        Self {
            r: self.r + rhs.r,
            i: self.i + rhs.i,
        }
    }
    fn sub(self, rhs: Self) -> Self {
        Self {
            r: self.r - rhs.r,
            i: self.i - rhs.i,
        }
    }
    fn mul(self, rhs: Self) -> Self {
        Self {
            r: self.r * rhs.r - self.i * rhs.i,
            i: self.r * rhs.i + self.i * rhs.r,
        }
    }
    fn norm2(self) -> f64 {
        self.r * self.r + self.i * self.i
    }
}

fn zpoly1(z: Complex, coeffs: &[Complex]) -> Complex {
    let mut a = coeffs[coeffs.len() - 1];
    for coeff in coeffs[..coeffs.len() - 1].iter().rev() {
        let t = a.r;
        a.r = coeff.r + z.r * t - z.i * a.i;
        a.i = coeff.i + z.r * a.i + z.i * t;
    }
    let t = a.r;
    a.r = z.r * t - z.i * a.i;
    a.i = z.r * a.i + z.i * t;
    a
}

fn zpolyd1(z: Complex, coeffs: &[Complex]) -> (Complex, Complex) {
    let mut a = coeffs[coeffs.len() - 1];
    let mut b = a;
    let mut first = true;
    for coeff in coeffs[..coeffs.len() - 1].iter().rev() {
        if first {
            first = false;
        } else {
            let t = b.r;
            b.r = a.r + z.r * t - z.i * b.i;
            b.i = a.i + z.r * b.i + z.i * t;
        }
        let t = a.r;
        a.r = coeff.r + z.r * t - z.i * a.i;
        a.i = coeff.i + z.r * a.i + z.i * t;
    }
    let t = b.r;
    b.r = a.r + z.r * t - z.i * b.i;
    b.i = a.i + z.r * b.i + z.i * t;
    let t = a.r;
    a.r = z.r * t - z.i * a.i;
    a.i = z.r * a.i + z.i * t;
    (a, b)
}

const GS48: [Complex; 5] = [
    Complex { r: 0.98879, i: 0.0 },
    Complex { r: 0.0, i: 0.0 },
    Complex {
        r: -0.050909,
        i: 0.0,
    },
    Complex { r: 0.0, i: 0.0 },
    Complex {
        r: 0.075528,
        i: 0.0,
    },
];
const ALSK_E: [Complex; 6] = [
    Complex {
        r: 0.9945303,
        i: 0.0,
    },
    Complex {
        r: 0.0052083,
        i: -0.0027404,
    },
    Complex {
        r: 0.0072721,
        i: 0.0048181,
    },
    Complex {
        r: -0.0151089,
        i: -0.1932526,
    },
    Complex {
        r: 0.0642675,
        i: -0.1381226,
    },
    Complex {
        r: 0.3582802,
        i: -0.2884586,
    },
];
const ALSK_S: [Complex; 6] = [
    Complex {
        r: 0.9972523,
        i: 0.0,
    },
    Complex {
        r: 0.0052513,
        i: -0.0041175,
    },
    Complex {
        r: 0.0074606,
        i: 0.0048125,
    },
    Complex {
        r: -0.0153783,
        i: -0.1968253,
    },
    Complex {
        r: 0.0636871,
        i: -0.1408027,
    },
    Complex {
        r: 0.3660976,
        i: -0.2937382,
    },
];
const GS50_E: [Complex; 10] = [
    Complex {
        r: 0.9827497,
        i: 0.0,
    },
    Complex {
        r: 0.0210669,
        i: 0.0053804,
    },
    Complex {
        r: -0.1031415,
        i: -0.0571664,
    },
    Complex {
        r: -0.0323337,
        i: -0.0322847,
    },
    Complex {
        r: 0.0502303,
        i: 0.1211983,
    },
    Complex {
        r: 0.0251805,
        i: 0.0895678,
    },
    Complex {
        r: -0.0012315,
        i: -0.1416121,
    },
    Complex {
        r: 0.0072202,
        i: -0.1317091,
    },
    Complex {
        r: -0.0194029,
        i: 0.0759677,
    },
    Complex {
        r: -0.0210072,
        i: 0.0834037,
    },
];
const GS50_S: [Complex; 10] = [
    Complex {
        r: 0.9842990,
        i: 0.0,
    },
    Complex {
        r: 0.0211642,
        i: 0.0037608,
    },
    Complex {
        r: -0.1036018,
        i: -0.0575102,
    },
    Complex {
        r: -0.0329095,
        i: -0.0320119,
    },
    Complex {
        r: 0.0499471,
        i: 0.1223335,
    },
    Complex {
        r: 0.0260460,
        i: 0.0899805,
    },
    Complex {
        r: 0.0007388,
        i: -0.1435792,
    },
    Complex {
        r: 0.0075848,
        i: -0.1334108,
    },
    Complex {
        r: -0.0216473,
        i: 0.0776645,
    },
    Complex {
        r: -0.0225161,
        i: 0.0853673,
    },
];

fn coeffs_for(key: &str) -> &'static [Complex] {
    match key {
        "alsk_e" => &ALSK_E,
        "alsk_s" => &ALSK_S,
        "gs50_e" => &GS50_E,
        "gs50_s" => &GS50_S,
        _ => &GS48,
    }
}

fn chi_from_phi(phi: f64, e: f64) -> f64 {
    if e == 0.0 {
        return phi;
    }
    let esphi = e * phi.sin();
    2.0 * (((FRAC_PI_2 + phi) * 0.5).tan() * ((1.0 - esphi) / (1.0 + esphi)).powf(0.5 * e)).atan()
        - FRAC_PI_2
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    let lon_0 = op.params.lon(0);
    let coeffs = coeffs_for(&op.params.text("coeffs").unwrap_or_default());
    let schio = op.params.real["schio"];
    let cchio = op.params.real["cchio"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let lam = lon - lon_0;
        let chi = chi_from_phi(lat, e);
        let schi = chi.sin();
        let cchi = chi.cos();
        let denom = 1.0 + schio * schi + cchio * cchi * lam.cos();
        if denom == 0.0 {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let s = 2.0 / denom;
        let p = zpoly1(
            Complex {
                r: s * cchi * lam.sin(),
                i: s * (cchio * schi - schio * cchi * lam.cos()),
            },
            coeffs,
        );
        operands.set_xy(i, a * p.r, a * p.i);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    let lon_0 = op.params.lon(0);
    let phi0 = op.params.lat(0);
    let coeffs = coeffs_for(&op.params.text("coeffs").unwrap_or_default());
    let schio = op.params.real["schio"];
    let cchio = op.params.real["cchio"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let mut p = Complex {
            r: operands.xy(i).0 / a,
            i: operands.xy(i).1 / a,
        };
        let mut converged = false;
        for _ in 0..20 {
            let (mut fxy, fpxy) = zpolyd1(p, coeffs);
            fxy = fxy.sub(Complex {
                r: operands.xy(i).0 / a,
                i: operands.xy(i).1 / a,
            });
            let den = fpxy.norm2();
            let dp = Complex {
                r: -(fxy.r * fpxy.r + fxy.i * fpxy.i) / den,
                i: -(fxy.i * fpxy.r - fxy.r * fpxy.i) / den,
            };
            p = p.add(dp);
            if dp.r.abs() + dp.i.abs() <= EPSLN {
                converged = true;
                break;
            }
        }
        if !converged {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let rh = p.r.hypot(p.i);
        if rh <= EPSLN {
            operands.set_xy(i, lon_0, phi0);
            successes += 1;
            continue;
        }
        let z = 2.0 * (0.5 * rh).atan();
        let sinz = z.sin();
        let cosz = z.cos();
        let chi = (cosz * schio + p.i * sinz * cchio / rh)
            .clamp(-1.0, 1.0)
            .asin();
        let mut phi = chi;
        for _ in 0..20 {
            if e == 0.0 {
                break;
            }
            let esphi = e * phi.sin();
            let dphi = 2.0
                * (((FRAC_PI_2 + chi) * 0.5).tan() * ((1.0 + esphi) / (1.0 - esphi)).powf(0.5 * e))
                    .atan()
                - FRAC_PI_2
                - phi;
            phi += dphi;
            if dphi.abs() <= EPSLN {
                break;
            }
        }
        let lon = (p.r * sinz).atan2(rh * cchio * cosz - p.i * schio * sinz) + lon_0;
        operands.set_xy(i, lon, phi);
        successes += 1;
    }
    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 3] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

fn build(def: &str, coeffs: &str, a: f64, es: f64, lon0: f64, phi0: f64) -> Result<Op, Error> {
    let raw = RawParameters::new(def, &BTreeMap::new());
    let mut params = ParsedParameters::new(&raw, &GAMUT)?;
    params.text.insert(
        "ellps",
        if es == 0.0 {
            format!("{a},0")
        } else {
            format!("{a},{}", 1.0 / (1.0 - (1.0 - es).sqrt()))
        },
    );
    params.real.insert("lon_0", lon0.to_radians());
    params.real.insert("lat_0", phi0.to_radians());
    params.text.insert("coeffs", coeffs.to_string());
    let chi0 = chi_from_phi(phi0.to_radians(), es.sqrt());
    params.real.insert("schio", chi0.sin());
    params.real.insert("cchio", chi0.cos());
    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}

pub fn gs48(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let _ = parameters;
    build("gs48", "gs48", 6_370_997.0, 0.0, -96.0, 39.0)
}

pub fn alsk(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let given = parameters.instantiated_as.split_into_parameters();
    let spherical = given.contains_key("R")
        || given
            .get("ellps")
            .and_then(|name| Ellipsoid::named(name).ok())
            .map(|ellps| ellps.flattening() == 0.0)
            .unwrap_or(false);
    if spherical {
        build("alsk", "alsk_s", 6_370_997.0, 0.0, -152.0, 64.0)
    } else {
        build("alsk", "alsk_e", 6_378_206.4, 0.00676866, -152.0, 64.0)
    }
}

pub fn gs50(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let given = parameters.instantiated_as.split_into_parameters();
    let spherical = given.contains_key("R")
        || given
            .get("ellps")
            .and_then(|name| Ellipsoid::named(name).ok())
            .map(|ellps| ellps.flattening() == 0.0)
            .unwrap_or(false);
    if spherical {
        build("gs50", "gs50_s", 6_370_997.0, 0.0, -120.0, 45.0)
    } else {
        build("gs50", "gs50_e", 6_378_206.4, 0.00676866, -120.0, 45.0)
    }
}
