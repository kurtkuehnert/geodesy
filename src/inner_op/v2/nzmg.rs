//! New Zealand Map Grid
use crate::authoring::*;
use crate::projection::ProjectionFrame;

const EPSLN: f64 = 1e-10;
const SEC5_TO_RAD: f64 = 0.484_813_681_109_536;
const RAD_TO_SEC5: f64 = 2.062_648_062_470_963_6;

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
    let mut acc = coeffs[coeffs.len() - 1];
    for coeff in coeffs[..coeffs.len() - 1].iter().rev() {
        acc = coeff.add(z.mul(acc));
    }
    z.mul(acc)
}

fn zpolyd1(z: Complex, coeffs: &[Complex]) -> (Complex, Complex) {
    let mut f = coeffs[coeffs.len() - 1];
    let mut fp = Complex { r: 0.0, i: 0.0 };
    for coeff in coeffs[..coeffs.len() - 1].iter().rev() {
        fp = f.add(z.mul(fp));
        f = coeff.add(z.mul(f));
    }
    (z.mul(f), f.add(z.mul(fp)))
}

const BF: [Complex; 6] = [
    Complex {
        r: 0.755_785_322_8,
        i: 0.0,
    },
    Complex {
        r: 0.249_204_646,
        i: 0.003_371_507,
    },
    Complex {
        r: -0.001_541_739,
        i: 0.041_058_56,
    },
    Complex {
        r: -0.101_629_07,
        i: 0.017_276_09,
    },
    Complex {
        r: -0.266_234_89,
        i: -0.362_492_18,
    },
    Complex {
        r: -0.687_098_3,
        i: -1.165_196_7,
    },
];

const TPSI: [f64; 10] = [
    0.639_917_507_3,
    -0.135_879_761_3,
    0.063_294_409,
    -0.025_268_53,
    0.011_787_9,
    -0.005_516_1,
    0.002_690_6,
    -0.001_333,
    0.000_67,
    -0.000_34,
];

const TPHI: [f64; 9] = [
    1.562_701_424_3,
    0.518_540_639_8,
    -0.033_330_98,
    -0.105_290_6,
    -0.036_859_4,
    0.007_317,
    0.012_20,
    0.003_94,
    -0.001_3,
];

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let a = op.params.real["a"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let dphi = frame.lat_delta(lat) * RAD_TO_SEC5;
        let mut pr = TPSI[TPSI.len() - 1];
        for coeff in TPSI[..TPSI.len() - 1].iter().rev() {
            pr = coeff + dphi * pr;
        }
        let p = Complex {
            r: dphi * pr,
            i: frame.lon_delta_raw(lon),
        };
        let z = zpoly1(p, &BF);
        operands.set_xy(i, frame.x_0 + a * z.i, frame.y_0 + a * z.r);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let frame = ProjectionFrame::from_params(&op.params);
    let a = op.params.real["a"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let mut p = Complex {
            r: (y - frame.y_0) / a,
            i: (x - frame.x_0) / a,
        };
        let mut converged = false;
        for _ in 0..20 {
            let (f, fp) = zpolyd1(p, &BF);
            let f = f.sub(Complex {
                r: (y - frame.y_0) / a,
                i: (x - frame.x_0) / a,
            });
            let den = fp.norm2();
            let dp = Complex {
                r: -(f.r * fp.r + f.i * fp.i) / den,
                i: -(f.i * fp.r - f.r * fp.i) / den,
            };
            p = p.add(dp);
            if dp.r.abs() + dp.i.abs() <= EPSLN {
                converged = true;
                break;
            }
        }
        if !converged {
            operands.set_xy(i, f64::NAN, f64::NAN);
            continue;
        }

        let mut phi = TPHI[TPHI.len() - 1];
        for coeff in TPHI[..TPHI.len() - 1].iter().rev() {
            phi = coeff + p.r * phi;
        }
        let lat = frame.lat_0 + p.r * phi * SEC5_TO_RAD;
        let lon = frame.lon_0 + p.i;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }

    successes
}

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 7] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Real { key: "lat_0", default: Some(-41_f64) },
    OpParameter::Real { key: "lon_0", default: Some(173_f64) },
    OpParameter::Real { key: "x_0", default: Some(2_510_000.0) },
    OpParameter::Real { key: "y_0", default: Some(6_023_150.0) }, 
    OpParameter::Real { key: "a", default: Some(6_378_388.0) },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
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
    fn nzmg_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("nzmg")?;

        let geo = [Coor4D::geo(1.0, 2.0, 0.0, 0.0)];
        let expected = [Coor4D::raw(
            3_352_675_144.747_425,
            -7_043_205_391.100_244,
            0.0,
            0.0,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-3);
        Ok(())
    }
}
