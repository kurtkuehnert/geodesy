//! Roussilhe Stereographic
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let k_0 = op.params.k(0);
    let a1 = op.params.real["a1"];
    let a2 = op.params.real["a2"];
    let a3 = op.params.real["a3"];
    let a4 = op.params.real["a4"];
    let a5 = op.params.real["a5"];
    let a6 = op.params.real["a6"];
    let b1 = op.params.real["b1"];
    let b2 = op.params.real["b2"];
    let b3 = op.params.real["b3"];
    let b4 = op.params.real["b4"];
    let b5 = op.params.real["b5"];
    let b6 = op.params.real["b6"];
    let b7 = op.params.real["b7"];
    let b8 = op.params.real["b8"];
    let s0 = op.params.real["s0"];
    let es = ellps.eccentricity_squared();
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let (sp, cp) = lat.sin_cos();
        let s = ellps.meridian_latitude_to_distance(lat) / a - s0;
        let s2 = s * s;
        let al = (lon - lon_0) * cp / (1.0 - es * sp * sp).sqrt();
        let al2 = al * al;
        let x = k_0 * al * (1.0 + s2 * (a1 + s2 * a4) - al2 * (a2 + s * a3 + s2 * a5 + al2 * a6));
        let y = k_0
            * (al2 * (b1 + al2 * b4)
                + s * (1.0
                    + al2 * (b3 - al2 * b6)
                    + s2 * (b2 + s2 * b8)
                    + s * al2 * (b5 + s * b7)));
        operands.set_xy(i, x_0 + a * x, y_0 + a * y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let a = ellps.semimajor_axis();
    let lon_0 = op.params.lon(0);
    let x_0 = op.params.x(0);
    let y_0 = op.params.y(0);
    let k_0 = op.params.k(0);
    let c1 = op.params.real["c1"];
    let c2 = op.params.real["c2"];
    let c3 = op.params.real["c3"];
    let c4 = op.params.real["c4"];
    let c5 = op.params.real["c5"];
    let c6 = op.params.real["c6"];
    let c7 = op.params.real["c7"];
    let c8 = op.params.real["c8"];
    let d1 = op.params.real["d1"];
    let d2 = op.params.real["d2"];
    let d3 = op.params.real["d3"];
    let d4 = op.params.real["d4"];
    let d5 = op.params.real["d5"];
    let d6 = op.params.real["d6"];
    let d7 = op.params.real["d7"];
    let d8 = op.params.real["d8"];
    let d9 = op.params.real["d9"];
    let d10 = op.params.real["d10"];
    let d11 = op.params.real["d11"];
    let s0 = op.params.real["s0"];
    let es = ellps.eccentricity_squared();
    let mut successes = 0usize;

    for i in 0..operands.len() {
        let x = (operands.xy(i).0 - x_0) / (a * k_0);
        let y = (operands.xy(i).1 - y_0) / (a * k_0);
        let x2 = x * x;
        let y2 = y * y;
        let al = x
            * (1.0 - c1 * y2
                + x2 * (c2 + c3 * y - c4 * x2 + c5 * y2 - c7 * x2 * y)
                + y2 * (c6 * y2 - c8 * x2 * y));
        let s = s0
            + y * (1.0 + y2 * (-d2 + d8 * y2))
            + x2 * (-d1
                + y * (-d3 + y * (-d5 + y * (-d7 + y * d11)))
                + x2 * (d4 + y * (d6 + y * d10) - x2 * d9));
        let lat = ellps.meridian_distance_to_latitude(s * a);
        let sp = lat.sin();
        let lon = lon_0 + al * (1.0 - es * sp * sp).sqrt() / lat.cos();
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
    let phi0 = params.lat(0).to_radians();
    let ellps = params.ellps(0);
    if ellps.flattening() == 0.0 {
        return Err(Error::General("Roussilhe requires an ellipsoid"));
    }
    let es = ellps.eccentricity_squared();
    let sp = phi0.sin();
    let es2 = es * sp * sp;
    let t = 1.0 - es2;
    let n0 = 1.0 / t.sqrt();
    let rr0_2 = t * t / (1.0 - es);
    let rr0_4 = rr0_2 * rr0_2;
    let tan_phi0 = phi0.tan();
    let tan2 = tan_phi0 * tan_phi0;

    params.real.insert("lat_0", phi0);
    params.real.insert("lon_0", params.lon(0).to_radians());
    params.real.insert(
        "s0",
        ellps.meridian_latitude_to_distance(phi0) / ellps.semimajor_axis(),
    );
    params.real.insert("a1", rr0_2 / 4.0);
    params
        .real
        .insert("a2", rr0_2 * (2.0 * tan2 - 1.0 - 2.0 * es2) / 12.0);
    params
        .real
        .insert("a3", rr0_2 * tan_phi0 * (1.0 + 4.0 * tan2) / (12.0 * n0));
    params.real.insert("a4", rr0_4 / 24.0);
    params
        .real
        .insert("a5", rr0_4 * (-1.0 + tan2 * (11.0 + 12.0 * tan2)) / 24.0);
    params
        .real
        .insert("a6", rr0_4 * (-2.0 + tan2 * (11.0 - 2.0 * tan2)) / 240.0);
    params.real.insert("b1", tan_phi0 / (2.0 * n0));
    params.real.insert("b2", rr0_2 / 12.0);
    params
        .real
        .insert("b3", rr0_2 * (1.0 + 2.0 * tan2 - 2.0 * es2) / 4.0);
    params
        .real
        .insert("b4", rr0_2 * tan_phi0 * (2.0 - tan2) / (24.0 * n0));
    params
        .real
        .insert("b5", rr0_2 * tan_phi0 * (5.0 + 4.0 * tan2) / (8.0 * n0));
    params
        .real
        .insert("b6", rr0_4 * (-2.0 + tan2 * (-5.0 + 6.0 * tan2)) / 48.0);
    params
        .real
        .insert("b7", rr0_4 * (5.0 + tan2 * (19.0 + 12.0 * tan2)) / 24.0);
    params.real.insert("b8", rr0_4 / 120.0);
    params.real.insert("c1", rr0_2 / 4.0);
    params
        .real
        .insert("c2", rr0_2 * (2.0 * tan2 - 1.0 - 2.0 * es2) / 12.0);
    params
        .real
        .insert("c3", rr0_2 * tan_phi0 * (1.0 + tan2) / (3.0 * n0));
    params
        .real
        .insert("c4", rr0_4 * (-3.0 + tan2 * (34.0 + 22.0 * tan2)) / 240.0);
    params
        .real
        .insert("c5", rr0_4 * (4.0 + tan2 * (13.0 + 12.0 * tan2)) / 24.0);
    params.real.insert("c6", rr0_4 / 16.0);
    params.real.insert(
        "c7",
        rr0_4 * tan_phi0 * (11.0 + tan2 * (33.0 + 16.0 * tan2)) / (48.0 * n0),
    );
    params
        .real
        .insert("c8", rr0_4 * tan_phi0 * (1.0 + 4.0 * tan2) / (36.0 * n0));
    params.real.insert("d1", tan_phi0 / (2.0 * n0));
    params.real.insert("d2", rr0_2 / 12.0);
    params
        .real
        .insert("d3", rr0_2 * (2.0 * tan2 + 1.0 - 2.0 * es2) / 4.0);
    params
        .real
        .insert("d4", rr0_2 * tan_phi0 * (1.0 + tan2) / (8.0 * n0));
    params
        .real
        .insert("d5", rr0_2 * tan_phi0 * (1.0 + 2.0 * tan2) / (4.0 * n0));
    params
        .real
        .insert("d6", rr0_4 * (1.0 + tan2 * (6.0 + 6.0 * tan2)) / 16.0);
    params
        .real
        .insert("d7", rr0_4 * tan2 * (3.0 + 4.0 * tan2) / 8.0);
    params.real.insert("d8", rr0_4 / 80.0);
    params.real.insert(
        "d9",
        rr0_4 * tan_phi0 * (-21.0 + tan2 * (178.0 - 26.0 * tan2)) / 720.0,
    );
    params.real.insert(
        "d10",
        rr0_4 * tan_phi0 * (29.0 + tan2 * (86.0 + 48.0 * tan2)) / (96.0 * n0),
    );
    params
        .real
        .insert("d11", rr0_4 * tan_phi0 * (37.0 + 44.0 * tan2) / (96.0 * n0));

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}
