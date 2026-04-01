mod aspect;
mod authalic;
mod azimuthal;
mod conformal;
mod frame;
mod framed;
mod gauss;
mod spherical_geodesic;
#[cfg(test)]
mod test_support;

pub(crate) use aspect::{AzimuthalAspect, ProjectionAspect};
pub(crate) use authalic::AuthalicLatitude;
pub(crate) use azimuthal::{spherical_inverse_equatorial, spherical_inverse_oblique};
pub(crate) use conformal::ConformalLatitude;
pub(crate) use frame::{ProjectionFrame, projection_gamut};
pub(crate) use framed::{Framed, FramedProjection, framed_gamut};
pub(crate) use gauss::Gauss;
pub(crate) use spherical_geodesic::SphericalGeodesic;
#[cfg(test)]
pub(crate) use test_support::{
    assert_forward_and_roundtrip, assert_inverse, assert_inverse_rejects, assert_op_err,
    assert_proj_match, assert_roundtrip,
};
