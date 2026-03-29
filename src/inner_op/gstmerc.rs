//! Gauss-Schreiber Transverse Mercator
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let a = ellps.semimajor_axis();
    let n1 = op.params.real["n1"];
    let n2 = op.params.real["n2"];
    let c = op.params.real["c"];
    let xs = op.params.real["xs"];
    let ys = op.params.real["ys"];
    let lon_0 = op.params.real["lon_0"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let l = n1 * (lon - lon_0);
        let ls = c + n1 * ancillary::ts((-lat).sin_cos(), e).ln();
        let sin_ls1 = l.sin() / ls.cosh();
        let ls1 = ancillary::ts((-sin_ls1.asin()).sin_cos(), 0.0).ln();
        let x = (xs + n2 * ls1) / a;
        let y = (ys + n2 * (ls.sinh() / l.cos()).atan()) / a;
        operands.set_xy(i, x * a, y * a);
        successes += 1;
    }

    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let e = ellps.eccentricity();
    let n1 = op.params.real["n1"];
    let n2 = op.params.real["n2"];
    let c = op.params.real["c"];
    let xs = op.params.real["xs"];
    let ys = op.params.real["ys"];
    let lon_0 = op.params.real["lon_0"];
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let l = (((x - xs) / n2).sinh() / ((y - ys) / n2).cos()).atan();
        let sin_c = ((y - ys) / n2).sin() / (((x - xs) / n2).cosh());
        let lc = ancillary::ts((-sin_c.asin()).sin_cos(), 0.0).ln();
        let lon = lon_0 + l / n1;
        let lat = -ancillary::pj_phi2(((lc - c) / n1).exp(), e);
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
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;
    let ellps = params.ellps(0);
    let es = ellps.eccentricity_squared();
    let e = ellps.eccentricity();
    let a = ellps.semimajor_axis();
    let lat_0 = params.real["lat_0"].to_radians();
    let lon_0 = params.real["lon_0"].to_radians();
    let k_0 = params.real["k_0"];
    let x_0 = params.real["x_0"];
    let y_0 = params.real["y_0"];

    let n1 = (1.0 + es * lat_0.cos().powi(4) / (1.0 - es)).sqrt();
    let phic = (lat_0.sin() / n1).asin();
    let c =
        ancillary::ts((-phic).sin_cos(), 0.0).ln() - n1 * ancillary::ts((-lat_0).sin_cos(), e).ln();
    let n2 = k_0 * a * (1.0 - es).sqrt() / (1.0 - es * lat_0.sin().powi(2));
    let xs = x_0;
    let ys = y_0 - n2 * phic;

    params.real.insert("lat_0", lat_0);
    params.real.insert("lon_0", lon_0);
    params.real.insert("n1", n1);
    params.real.insert("n2", n2);
    params.real.insert("c", c);
    params.real.insert("xs", xs);
    params.real.insert("ys", ys);

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
    fn gstmerc_matches_proj() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("gstmerc ellps=6400000,0")?;

        let geo = [Coor4D::geo(1.0, 2.0, 0.0, 0.0)];
        let expected = [Coor4D::raw(
            223_413.466_406_322,
            111_769.145_040_586,
            0.0,
            0.0,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-6);
        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot2(&geo[0]) < 1e-10);
        Ok(())
    }

    #[test]
    fn gstmerc_respects_false_origin() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "gstmerc lat_0=-21.1166666666667 lon_0=55.5333333333333 k_0=1 x_0=160000 y_0=50000 ellps=intl",
        )?;

        let geo = [Coor4D::raw(
            55.54497939138478_f64.to_radians(),
            (-21.059138278165022_f64).to_radians(),
            0.0,
            0.0,
        )];
        let expected = [Coor4D::raw(161_210.417_590_794_74, 56_369.498_060_546_89, 0.0, 0.0)];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&expected[0]) < 1e-6);

        ctx.apply(op, Inv, &mut operands)?;
        assert!((operands[0][0] - geo[0][0]).abs() < 1e-10);
        assert!((operands[0][1] - geo[0][1]).abs() < 1e-10);
        Ok(())
    }
}
