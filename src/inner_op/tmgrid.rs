//! Transverse Mercator Zoned Grid System
use crate::authoring::*;

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let x_0 = op.params.x(0);
    let lon_i = op.params.real["lon_i"];
    let zone_width = op.params.real["zone_width"];
    let zone_count = op.params.natural["zone_count"] as i64;
    let qs = op.params.real["scaled_radius"];
    let zb = op.params.real["zb"];
    let conformal = op.params.fourier_coefficients["conformal"];
    let tm = op.params.fourier_coefficients["tm"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (lon, lat) = operands.xy(i);
        let zone = forward_zone(lon.to_degrees(), lon_i.to_degrees(), zone_width.to_degrees(), zone_count);
        let lon_0 = zone_central_meridian(lon_i, zone_width, zone as f64);
        let x_offset = x_0 + zone as f64 * 1_000_000.0;

        let lat = ellps.latitude_geographic_to_conformal(lat, &conformal);
        let lon = lon - lon_0;

        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_lon, cos_lon) = lon.sin_cos();
        let cos_lat_lon = cos_lat * cos_lon;
        let mut lat = sin_lat.atan2(cos_lat_lon);

        let inv_denom_tan_lon = sin_lat.hypot(cos_lat_lon).recip();
        let tan_lon = sin_lon * cos_lat * inv_denom_tan_lon;
        let mut lon = tan_lon.asinh();

        let two_inv_denom_tan_lon = 2.0 * inv_denom_tan_lon;
        let two_inv_denom_tan_lon_square = two_inv_denom_tan_lon * inv_denom_tan_lon;
        let tmp_r = cos_lat_lon * two_inv_denom_tan_lon_square;
        let trig = [sin_lat * tmp_r, cos_lat_lon * tmp_r - 1.0];
        let hyp = [
            tan_lon * two_inv_denom_tan_lon,
            two_inv_denom_tan_lon_square - 1.0,
        ];
        let dc = fourier::complex_sin_optimized_for_tmerc(trig, hyp, &tm.fwd);
        lat += dc[0];
        lon += dc[1];

        operands.set_xy(i, qs * lon + x_offset, qs * lat + zb);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let ellps = op.params.ellps(0);
    let x_0 = op.params.x(0);
    let lon_i = op.params.real["lon_i"];
    let zone_width = op.params.real["zone_width"];
    let zone_count = op.params.natural["zone_count"] as i64;
    let qs = op.params.real["scaled_radius"];
    let zb = op.params.real["zb"];
    let conformal = op.params.fourier_coefficients["conformal"];
    let tm = op.params.fourier_coefficients["tm"];

    let mut successes = 0_usize;
    for i in 0..operands.len() {
        let (x, y) = operands.xy(i);
        let zone = inverse_zone(x, zone_count);
        if !(1..=zone_count).contains(&zone) {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let lon_0 = zone_central_meridian(lon_i, zone_width, zone as f64);
        let x_local = x - x_0 - zone as f64 * 1_000_000.0;

        let mut lon = x_local / qs;
        let mut lat = (y - zb) / qs;
        let dc = fourier::complex_sin([2.0 * lat, 2.0 * lon], &tm.inv);
        lat += dc[0];
        lon += dc[1];
        lon = gudermannian::fwd(lon);

        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_lon, cos_lon) = lon.sin_cos();
        let cos_lat_lon = cos_lat * cos_lon;
        lon = sin_lon.atan2(cos_lat_lon);
        lat = (sin_lat * cos_lon).atan2(sin_lon.hypot(cos_lat_lon));

        operands.set_xy(
            i,
            angular::normalize_symmetric(lon + lon_0),
            ellps.latitude_conformal_to_geographic(lat, &conformal),
        );
        successes += 1;
    }
    successes
}

fn forward_zone(lon_deg: f64, lon_i_deg: f64, zone_width_deg: f64, zone_count: i64) -> i64 {
    let zone = (((lon_deg - lon_i_deg) / zone_width_deg) + 1e-12).floor() as i64;
    zone.rem_euclid(zone_count) + 1
}

fn inverse_zone(easting: f64, zone_count: i64) -> i64 {
    let zone = (easting / 1_000_000.0).floor() as i64;
    zone.clamp(1, zone_count)
}

fn zone_central_meridian(lon_i: f64, zone_width: f64, zone: f64) -> f64 {
    lon_i + zone * zone_width - zone_width / 2.0
}

fn precompute(op: &mut Op) {
    let ellps = op.params.ellps(0);
    let lat_0 = op.params.lat(0).to_radians();
    let y_0 = op.params.y(0);
    let qs = op.params.k(0) * ellps.semimajor_axis() * ellps.normalized_meridian_arc_unit();
    op.params.real.insert("scaled_radius", qs);

    let conformal = ellps.coefficients_for_conformal_latitude_computations();
    op.params.fourier_coefficients.insert("conformal", conformal);

    let n = ellps.third_flattening();
    let tm = fourier_coefficients(n, &TRANSVERSE_MERCATOR);
    op.params.fourier_coefficients.insert("tm", tm);

    let z = ellps.latitude_geographic_to_conformal(lat_0, &conformal);
    let zb = y_0 - qs * (z + fourier::sin(2.0 * z, &tm.fwd));
    op.params.real.insert("zb", zb);
}

#[rustfmt::skip]
const TRANSVERSE_MERCATOR: PolynomialCoefficients = PolynomialCoefficients {
    fwd: [
        [1./2.,   -2./3.,   5./16.,   41./180.,   -127./288.0 ,   7891./37800.],
        [0., 13./48.,   -3./5.,   557./1440.,   281./630.,   -1983433./1935360.],
        [0., 0., 61./240.,  -103./140.,   15061./26880.,   167603./181440.],
        [0., 0., 0., 49561./161280.,   -179./168.,   6601661./7257600.],
        [0., 0., 0., 0., 34729./80640.,   -3418889./1995840.],
        [0., 0., 0., 0., 0., 212378941./319334400.]
    ],
    inv: [
        [-1./2.,   2./3.,   -37./96.,   1./360.,   81./512.,   -96199./604800.],
        [0., -1./48.,   -1./15.,   437./1440.,   -46./105.,   1118711./3870720.],
        [0., 0., -17./480.,   37./840.,   209./4480.,   -5569./90720.],
        [0., 0., 0., -4397./161280.,   11./504.,   830251./7257600.],
        [0., 0., 0., 0., -4583./161280.,   108847./3991680.],
        [0., 0., 0., 0., 0., -20648693./638668800.]
    ]
};

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 8] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_i", default: None },
    OpParameter::Real { key: "zone_width", default: None },
    OpParameter::Real { key: "k_0", default: Some(1_f64) },
    OpParameter::Real { key: "x_0", default: Some(0_f64) },
    OpParameter::Real { key: "y_0", default: Some(0_f64) },
];

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let mut params = ParsedParameters::new(parameters, &GAMUT)?;

    let zone_width_deg = params.real("zone_width")?;
    if zone_width_deg <= 0.0 {
        return Err(Error::General("TM Grid: zone_width must be > 0"));
    }
    let zone_count = (360.0 / zone_width_deg).round() as usize;
    params.real.insert("lat_0", params.lat(0).to_radians());
    params
        .real
        .insert("lon_i", params.real("lon_i")?.to_radians());
    params
        .real
        .insert("zone_width", zone_width_deg.to_radians());
    params.natural.insert("zone_count", zone_count);

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    let mut op = Op {
        descriptor,
        params,
        steps: None,
    };
    precompute(&mut op);
    Ok(op)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tmgrid_roundtrips_across_two_utm_zones() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "tmgrid lat_0=0 lon_i=-180 zone_width=6 k_0=0.9996 x_0=500000 y_0=0 ellps=WGS84",
        )?;

        let original = [
            Coor4D::geo(55.0, 12.0, 0.0, 0.0),
            Coor4D::geo(10.0, -77.0, 0.0, 0.0),
        ];
        let mut operands = original;
        ctx.apply(op, Fwd, &mut operands)?;
        assert_eq!((operands[0][0] / 1_000_000.0).floor() as i64, 33);
        assert_eq!((operands[1][0] / 1_000_000.0).floor() as i64, 18);

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&original[i]) < 5e-6);
        }
        Ok(())
    }
}
