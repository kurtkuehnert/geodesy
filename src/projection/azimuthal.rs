pub(crate) fn spherical_inverse_equatorial(
    x: f64,
    y: f64,
    rho: f64,
    c: f64,
    rho_tolerance: f64,
    origin_lat: f64,
) -> (f64, f64) {
    if rho.abs() <= rho_tolerance {
        return (0.0, origin_lat);
    }

    let (sinc, cosc) = c.sin_cos();
    let lat = (y * sinc / rho).asin();
    let lon = if cosc != 0.0 || x != 0.0 {
        (x * sinc).atan2(cosc * rho)
    } else {
        0.0
    };
    (lon, lat)
}

pub(crate) fn spherical_inverse_oblique(
    x: f64,
    y: f64,
    rho: f64,
    c: f64,
    rho_tolerance: f64,
    origin_lat: f64,
    sinph0: f64,
    cosph0: f64,
) -> (f64, f64) {
    if rho.abs() <= rho_tolerance {
        return (0.0, origin_lat);
    }

    let (sinc, cosc) = c.sin_cos();
    let lat = (cosc * sinph0 + y * sinc * cosph0 / rho).asin();
    let denom = cosc - sinph0 * lat.sin();
    let lon = if denom != 0.0 || x != 0.0 {
        (x * sinc * cosph0).atan2(denom * rho)
    } else {
        0.0
    };
    (lon, lat)
}
