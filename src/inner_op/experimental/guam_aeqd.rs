use crate::authoring::*;
use crate::projection::RectifyingLatitude;

#[derive(Clone, Copy, Debug)]
pub(crate) struct GuamAeqdInner {
    a: f64,
    lat_0: f64,
    ellps: Ellipsoid,
    m1: f64,
    rectifying: RectifyingLatitude,
}

pub(crate) type GuamAeqd = Framed<GuamAeqdInner>;

impl GuamAeqdInner {}

impl FramedProjection for GuamAeqdInner {
    const NAME: &'static str = "guam_aeqd";
    const TITLE: &'static str = "Guam Azimuthal Equidistant";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = framed_gamut!(
        OpParameter::Text { key: "ellps", default: Some("clrk66") },
        OpParameter::Real { key: "lat_0", default: Some(0_f64) },
    );

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        let ellps = params.ellps(0);
        let a = ellps.semimajor_axis();
        let lat_0 = params.lat(0);
        let rectifying = ellps.rectifying();
        let m1 = rectifying.distance_from_latitude(lat_0);

        Ok(Self {
            a,
            lat_0,
            ellps,
            m1,
            rectifying,
        })
    }

    fn fwd(&self, lam: f64, lat: f64) -> Option<(f64, f64)> {
        let es = self.ellps.eccentricity_squared();
        let (sinphi, cosphi) = lat.sin_cos();
        let t = 1.0 / (1.0 - es * sinphi * sinphi).sqrt();
        let d = self.rectifying.distance_from_latitude(lat);
        Some((
            self.a * lam * cosphi * t,
            d - self.m1 + 0.5 * self.a * lam * lam * cosphi * sinphi * t,
        ))
    }

    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let es = self.ellps.eccentricity_squared();
        let x_norm = x / self.a;
        let y_norm = y / self.a;
        let x2 = 0.5 * x_norm * x_norm;
        let lat = (0..3).fold(self.lat_0, |lat, _| {
            let t = (1.0 - es * lat.sin().powi(2)).sqrt();
            let d = self.m1 + self.a * (y_norm - x2 * lat.tan() * t);
            self.rectifying.latitude_from_distance(d)
        });
        let t = (1.0 - es * lat.sin().powi(2)).sqrt();
        Some((x_norm * t / lat.cos(), lat))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn guam_aeqd_matches_proj_gie_example() -> Result<(), Error> {
        assert_proj_match(
            "guam_aeqd ellps=clrk66 x_0=50000 y_0=50000 lon_0=144.74875069444445 lat_0=13.47246633333333",
            Coor4D::geo(13.33903846111111, 144.63533129166666, 0.0, 0.0),
            Coor4D::raw(37_712.482_258_147_284, 35_242.003_289_566_49, 0.0, 0.0),
        )
    }
}
