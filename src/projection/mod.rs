mod aspect;
mod azimuthal;
mod conic;
mod frame;
mod framed;
mod gauss;
mod geodesic;
mod latitudes;
#[cfg(test)]
mod test_support;

pub(crate) use aspect::{AzimuthalAspect, ProjectionAspect};
pub(crate) use azimuthal::{spherical_inverse_equatorial, spherical_inverse_oblique};
pub(crate) use conic::{Conic, ConicInverse, ConicSetupError, StandardParallels};
pub(crate) use frame::{projection_gamut, ProjectionFrame};
pub(crate) use framed::{framed_gamut, Framed, FramedProjection};
pub(crate) use gauss::Gauss;
pub(crate) use geodesic::GeodesicPath;
pub(crate) use latitudes::{AuthalicLatitude, ConformalLatitude, MeridianLatitude};
#[cfg(test)]
pub(crate) use test_support::{
    assert_forward_and_roundtrip, assert_inverse, assert_inverse_rejects, assert_op_err,
    assert_proj_match, assert_proj_match_with_tol, assert_roundtrip,
};
