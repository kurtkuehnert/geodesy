use crate::authoring::ParsedParameters;
use crate::ellipsoid::EllipsoidBase;
use crate::math::angular;

/// Shared cached frame values for geographic -> projected operators.
///
/// This is intentionally small and mechanical: the common origin, false-origin,
/// scale, and semimajor-axis values plus a few helpers for the most repeated
/// runtime adjustments. It is not a generic projection framework.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ProjectionFrame {
    pub lon_0: f64,
    pub lat_0: f64,
    pub x_0: f64,
    pub y_0: f64,
    pub k_0: f64,
    pub a: f64,
}

#[allow(dead_code)]
impl ProjectionFrame {
    /// Build a frame from normalized parsed parameters.
    ///
    /// `lon_0`/`lat_0` are expected to already be in radians when present.
    pub fn from_params(params: &ParsedParameters) -> Self {
        Self {
            lon_0: params.lon(0),
            lat_0: params.lat(0),
            x_0: params.x(0),
            y_0: params.y(0),
            k_0: params.k(0),
            a: params.ellps(0).semimajor_axis(),
        }
    }

    /// Central-meridian delta normalized into [-pi, pi].
    pub fn lon_delta(&self, lon: f64) -> f64 {
        angular::normalize_symmetric(lon - self.lon_0)
    }

    /// Absolute longitude from a central-meridian delta.
    pub fn apply_lon_delta(&self, lam: f64) -> f64 {
        self.lon_0 + lam
    }

    /// Central-meridian delta without wrapping.
    pub fn lon_delta_raw(&self, lon: f64) -> f64 {
        lon - self.lon_0
    }

    /// Latitude delta from the latitude of origin.
    pub fn lat_delta(&self, lat: f64) -> f64 {
        lat - self.lat_0
    }

    /// Apply the cached false origin to local projected coordinates.
    pub fn apply_false_origin(&self, x: f64, y: f64) -> (f64, f64) {
        (x + self.x_0, y + self.y_0)
    }

    /// Remove the cached false origin from projected coordinates.
    pub fn remove_false_origin(&self, x: f64, y: f64) -> (f64, f64) {
        (x - self.x_0, y - self.y_0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::authoring::{Minimal, OpParameter, ParsedParameters, RawParameters};

    const GAMUT: [OpParameter; 6] = [
        OpParameter::Text {
            key: "ellps",
            default: Some("GRS80"),
        },
        OpParameter::Real {
            key: "lat_0",
            default: Some(0.0),
        },
        OpParameter::Real {
            key: "lon_0",
            default: Some(0.0),
        },
        OpParameter::Real {
            key: "x_0",
            default: Some(0.0),
        },
        OpParameter::Real {
            key: "y_0",
            default: Some(0.0),
        },
        OpParameter::Real {
            key: "k_0",
            default: Some(1.0),
        },
    ];

    fn normalized_params(definition: &str) -> ParsedParameters {
        let _ctx = Minimal::default();
        let raw = RawParameters::new(definition, &Default::default());
        ParsedParameters::new(&raw, &GAMUT).expect("params")
    }

    #[test]
    fn projection_frame_reads_normalized_common_values() {
        let params = normalized_params("lat_0=12 lon_0=34 x_0=100 y_0=200 k_0=0.9996 ellps=WGS84");
        let frame = ProjectionFrame::from_params(&params);

        assert!((frame.lat_0 - 12_f64.to_radians()).abs() < 1e-15);
        assert!((frame.lon_0 - 34_f64.to_radians()).abs() < 1e-15);
        assert_eq!(frame.x_0, 100.0);
        assert_eq!(frame.y_0, 200.0);
        assert_eq!(frame.k_0, 0.9996);
        assert!((frame.a - 6_378_137.0).abs() < 1e-9);
    }

    #[test]
    fn projection_frame_helpers_apply_and_remove_false_origin() {
        let params = normalized_params("x_0=500000 y_0=1000000");
        let frame = ProjectionFrame::from_params(&params);

        assert_eq!(frame.apply_false_origin(12.0, 34.0), (500012.0, 1000034.0));
        assert_eq!(frame.remove_false_origin(500012.0, 1000034.0), (12.0, 34.0));
    }

    #[test]
    fn projection_frame_deltas_are_relative_to_cached_origin() {
        let params = normalized_params("lat_0=10 lon_0=179");
        let frame = ProjectionFrame::from_params(&params);
        let lon_delta = frame.lon_delta((-179_f64).to_radians());
        let lat_delta = frame.lat_delta(13_f64.to_radians());
        let raw_lon_delta = frame.lon_delta_raw((-179_f64).to_radians());
        let restored_lon = frame.apply_lon_delta(lon_delta);

        assert!((lon_delta - 2_f64.to_radians()).abs() < 1e-15);
        assert!((raw_lon_delta + 358_f64.to_radians()).abs() < 1e-15);
        assert!((restored_lon - 181_f64.to_radians()).abs() < 1e-15);
        assert!((lat_delta - 3_f64.to_radians()).abs() < 1e-15);
    }
}
