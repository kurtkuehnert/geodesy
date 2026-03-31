use crate::authoring::*;
use crate::inner_op::aeqd::AeqdState;

#[rustfmt::skip]
const GAMUT: [OpParameter; 6] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug)]
struct GuamAeqdState {
    base: AeqdState,
    m1: Option<f64>,
}

impl GuamAeqdState {
    fn new(params: &ParsedParameters) -> Result<Self, Error> {
        let base = AeqdState::new(params)?;
        let m1 = if base.spherical {
            None
        } else {
            Some(base.ellps.meridian_latitude_to_distance(base.frame.lat_0))
        };
        Ok(Self { base, m1 })
    }

    fn guam_fwd(&self, lam: f64, phi: f64) -> (f64, f64) {
        let es = self.base.ellps.eccentricity_squared();
        let cosphi = phi.cos();
        let sinphi = phi.sin();
        let t = 1.0 / (1.0 - es * sinphi * sinphi).sqrt();
        (
            self.base.frame.a * lam * cosphi * t,
            self.base.ellps.meridian_latitude_to_distance(phi) - self.m1.unwrap_or(0.0)
                + 0.5 * self.base.frame.a * lam * lam * cosphi * sinphi * t,
        )
    }

    fn guam_inv(&self, x: f64, y: f64) -> (f64, f64) {
        let es = self.base.ellps.eccentricity_squared();
        let x_norm = x / self.base.frame.a;
        let y_norm = y / self.base.frame.a;
        let x2 = 0.5 * x_norm * x_norm;
        let mut lat = self.base.frame.lat_0;
        for _ in 0..3 {
            let t = (1.0 - es * lat.sin().powi(2)).sqrt();
            lat = self.base.ellps.meridian_distance_to_latitude(
                self.m1.unwrap_or(0.0) + self.base.frame.a * (y_norm - x2 * lat.tan() * t),
            );
        }
        let t = (1.0 - es * lat.sin().powi(2)).sqrt();
        (self.base.frame.lon_0 + x_norm * t / lat.cos(), lat)
    }
}

struct GuamAeqd;

impl PointOp for GuamAeqd {
    type State = GuamAeqdState;
    const GAMUT: &'static [OpParameter] = &GAMUT;

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self::State, Error> {
        GuamAeqdState::new(params)
    }

    fn fwd(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = state.base.frame.remove_central_meridian(lon);
        let (x, y) = if state.base.spherical {
            state.base.spherical_fwd(lam, lat)?
        } else {
            state.guam_fwd(lam, lat)
        };
        let (x, y) = state.base.frame.apply_false_origin(x, y);
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(state: &Self::State, coord: Coor4D) -> Option<Coor4D> {
        let (x, y) = state.base.frame.remove_false_origin(coord[0], coord[1]);
        let (lon, lat) = if state.base.spherical {
            state.base.spherical_inv(x, y)?
        } else {
            state.guam_inv(x, y)
        };
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}

pub fn new(parameters: &RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
    Op::point::<GuamAeqd>(parameters, ctx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::projection::assert_forward_and_roundtrip;

    #[test]
    fn guam_aeqd_matches_proj_gie_example() -> Result<(), Error> {
        assert_forward_and_roundtrip(
            "guam_aeqd ellps=clrk66 x_0=50000 y_0=50000 lon_0=144.74875069444445 lat_0=13.47246633333333",
            Coor4D::geo(13.33903846111111, 144.63533129166666, 0.0, 0.0),
            Coor4D::raw(37712.48, 35242.00, 0.0, 0.0),
            0.01,
            1e-10,
        )
    }
}
