//! Molodensky-Badekas
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    transform_common(op, operands, Direction::Fwd)
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    transform_common(op, operands, Direction::Inv)
}

fn transform_common(op: &Op, operands: &mut dyn CoordinateSet, direction: Direction) -> usize {
    let t = op.params.series("T").unwrap();
    let s = op.params.real("S").unwrap();
    let p = op.params.series("P").unwrap();
    let m = op.params.series("ROTFLAT").unwrap();
    let rotated = op.params.boolean("rotated");
    let rot = [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]];

    for i in 0..operands.len() {
        let mut c = operands.get_coord(i);
        let local_x = c[0] - p[0];
        let local_y = c[1] - p[1];
        let local_z = c[2] - p[2];

        if direction == Direction::Fwd {
            let (x, y, z) = if rotated {
                (
                    local_x * rot[0][0] + local_y * rot[0][1] + local_z * rot[0][2],
                    local_x * rot[1][0] + local_y * rot[1][1] + local_z * rot[1][2],
                    local_x * rot[2][0] + local_y * rot[2][1] + local_z * rot[2][2],
                )
            } else {
                (local_x, local_y, local_z)
            };

            c[0] = p[0] + t[0] + s * x;
            c[1] = p[1] + t[1] + s * y;
            c[2] = p[2] + t[2] + s * z;
            operands.set_coord(i, &c);
            continue;
        }

        let x = (c[0] - p[0] - t[0]) / s;
        let y = (c[1] - p[1] - t[1]) / s;
        let z = (c[2] - p[2] - t[2]) / s;

        if rotated {
            c[0] = p[0] + x * rot[0][0] + y * rot[1][0] + z * rot[2][0];
            c[1] = p[1] + x * rot[0][1] + y * rot[1][1] + z * rot[2][1];
            c[2] = p[2] + x * rot[0][2] + y * rot[1][2] + z * rot[2][2];
        } else {
            c[0] = p[0] + x;
            c[1] = p[1] + y;
            c[2] = p[2] + z;
        }
        operands.set_coord(i, &c);
    }

    operands.len()
}

#[rustfmt::skip]
const GAMUT: [OpParameter; 15] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "exact" },
    OpParameter::Flag { key: "mm" },
    OpParameter::Flag { key: "uas" },
    OpParameter::Text { key: "convention", default: Some("") },
    OpParameter::Real { key: "x", default: Some(0.0) },
    OpParameter::Real { key: "y", default: Some(0.0) },
    OpParameter::Real { key: "z", default: Some(0.0) },
    OpParameter::Real { key: "rx", default: Some(0.0) },
    OpParameter::Real { key: "ry", default: Some(0.0) },
    OpParameter::Real { key: "rz", default: Some(0.0) },
    OpParameter::Real { key: "s", default: Some(0.0) },
    OpParameter::Real { key: "px", default: Some(0.0) },
    OpParameter::Real { key: "py", default: Some(0.0) },
    OpParameter::Real { key: "pz", default: Some(0.0) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let linear_scale = if params.boolean("mm") { 0.001 } else { 1.0 };
    let angular_scale = if params.boolean("uas") {
        std::f64::consts::PI / (3.6e9 * 180.0)
    } else {
        std::f64::consts::PI / (3.6e3 * 180.0)
    };

    let t = [
        linear_scale * params.real("x")?,
        linear_scale * params.real("y")?,
        linear_scale * params.real("z")?,
    ];
    let p = [
        linear_scale * params.real("px")?,
        linear_scale * params.real("py")?,
        linear_scale * params.real("pz")?,
    ];
    let r = [
        angular_scale * params.real("rx")?,
        angular_scale * params.real("ry")?,
        angular_scale * params.real("rz")?,
    ];
    let scale = 1.0 + params.real("s")? * 1e-6;

    let has_transform = t != [0.0, 0.0, 0.0] || r != [0.0, 0.0, 0.0] || scale != 1.0;
    if !has_transform {
        return Err(Error::MissingParam(
            "x/y/z/rx/ry/rz/s".to_string(),
        ));
    }

    let convention = params.text("convention")?;
    let rotated = r != [0.0, 0.0, 0.0];
    let mut position_vector = true;
    if rotated {
        if !matches!(convention.as_str(), "position_vector" | "coordinate_frame") {
            return Err(Error::MissingParam("convention".to_string()));
        }
        if convention == "coordinate_frame" {
            position_vector = false;
        }
        params.boolean.insert("rotated");
    }

    params.series.insert("T", Vec::from(t));
    params.series.insert("P", Vec::from(p));
    params.real.insert("S", scale);

    let rot = rotation_matrix(&r, params.boolean("exact"), position_vector);
    let mut flat = Vec::from(rot[0]);
    flat.extend(rot[1]);
    flat.extend(rot[2]);
    params.series.insert("ROTFLAT", flat);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    Ok(Op {
        descriptor,
        params,
        steps: None,
    })
}

fn rotation_matrix(r: &[f64], exact: bool, position_vector: bool) -> [[f64; 3]; 3] {
    let (rx, ry, rz) = (r[0], r[1], r[2]);
    let (mut sx, mut sy, mut sz) = (rx, ry, rz);
    let (mut cx, mut cy, mut cz) = (1.0, 1.0, 1.0);

    if exact {
        (sx, cx) = rx.sin_cos();
        (sy, cy) = ry.sin_cos();
        (sz, cz) = rz.sin_cos();
    }

    let r11 = cy * cz;
    let mut r12 = cx * sz;
    let mut r13 = -cx * sy * cz;

    let r21 = -cy * sz;
    let mut r22 = cx * cz;
    let mut r23 = sx * cz;

    let r31 = sy;
    let r32 = -sx * cy;
    let r33 = cx * cy;

    if exact {
        r12 += sx * sy * cz;
        r13 += sx * sz;

        r22 -= sx * sy * sz;
        r23 += cx * sy * sz;
    }

    if position_vector {
        return [[r11, r21, r31], [r12, r22, r32], [r13, r23, r33]];
    }
    [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molobadekas_matches_proj_fixture() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "molobadekas convention=coordinate_frame x=-270.933 y=115.599 z=-360.226 rx=-5.266 ry=-1.238 rz=2.381 s=-5.109 px=2464351.59 py=-5783466.61 pz=974809.81",
        )?;

        let input = [Coor4D::raw(2550408.96, -5749912.26, 1054891.11, 0.0)];
        let expected = [Coor4D::raw(2550138.45, -5749799.87, 1054530.82, 0.0)];

        let mut operands = input;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot3(&expected[0]) < 0.01);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(operands[0].hypot3(&input[0]) < 0.01);
        Ok(())
    }

    #[test]
    fn molobadekas_requires_transform_parameters() {
        let mut ctx = Minimal::default();
        assert!(ctx.op("molobadekas").is_err());
    }
}
