mod azimuthal;
mod conic;
mod frame;
mod framed;
mod gauss;
mod geodesic;
mod latitudes;
#[cfg(test)]
mod test_support;

pub(crate) use azimuthal::AzimuthalAspect;
pub(crate) use conic::Conic;
pub(crate) use frame::{ProjectionFrame, projection_gamut};
pub(crate) use framed::{Framed, FramedProjection, framed_gamut};
pub(crate) use gauss::Gauss;
pub(crate) use geodesic::GeodesicPath;
pub(crate) use latitudes::{AuthalicLatitude, ConformalLatitude, RectifyingLatitude};
#[cfg(test)]
pub(crate) use test_support::{assert_op_err, assert_proj_match};
