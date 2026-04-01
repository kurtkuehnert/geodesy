use crate::authoring::*;

macro_rules! framed_gamut {
    ( $( $extra:expr ),* $(,)? ) => {
        &[
            OpParameter::Flag { key: "inv" },
            OpParameter::Real { key: "lon_0", default: Some(0.0) },
            OpParameter::Real { key: "x_0", default: Some(0.0) },
            OpParameter::Real { key: "y_0", default: Some(0.0) },
            $( $extra ),*
        ]
    };
}

pub(crate) use framed_gamut;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Framed<T> {
    lon_0: f64,
    x_0: f64,
    y_0: f64,
    inner: T,
}

pub(crate) trait FramedProjection: Sized + Send + Sync + 'static {
    const NAME: &'static str;
    const GAMUT: &'static [OpParameter];

    fn build(params: &ParsedParameters, ctx: &dyn Context) -> Result<Self, Error>;

    fn fwd(&self, lam: f64, phi: f64) -> Option<(f64, f64)>;
    fn inv(&self, x: f64, y: f64) -> Option<(f64, f64)>;
}

impl<T: FramedProjection> PointOp for Framed<T> {
    const NAME: &'static str = T::NAME;
    const GAMUT: &'static [OpParameter] = T::GAMUT;

    fn build(params: &ParsedParameters, ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self {
            lon_0: params.lon(0),
            x_0: params.x(0),
            y_0: params.y(0),
            inner: T::build(params, ctx)?,
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let (lon, lat) = coord.xy();
        let lam = angular::normalize_symmetric(lon - self.lon_0);
        let (x_local, y_local) = self.inner.fwd(lam, lat)?;
        let x = x_local + self.x_0;
        let y = y_local + self.y_0;
        Some(Coor4D::raw(x, y, coord[2], coord[3]))
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let x = coord[0] - self.x_0;
        let y = coord[1] - self.y_0;
        let (lam, lat) = self.inner.inv(x, y)?;
        let lon = self.lon_0 + lam;
        Some(Coor4D::raw(lon, lat, coord[2], coord[3]))
    }
}
