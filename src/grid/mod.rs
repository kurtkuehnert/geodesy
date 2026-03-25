//! Grid characteristics and interpolation.

pub mod gravsoft;
pub mod gsa;
#[cfg(feature = "with_gtg")]
pub mod gtg;
pub mod gtx;
pub mod ntv2;
pub mod unigrid;
use crate::prelude::*;
use std::{fmt::Debug, sync::Arc};

pub trait Grid: Debug + Sync + Send {
    fn bands(&self) -> usize;
    /// Returns true if `coord` is contained by `self` or lies within a margin of
    /// `margin` grid cell units. Typically `margin` should be on the order of 1
    /// If `all_inclusive==true`, a point is considered contained if it is on any
    /// of the grid borders. Otherwise only the westmost and southernmost border
    /// is considered to be within.
    fn contains(&self, coord: Coor4D, margin: f64, all_inclusive: bool) -> bool;

    /// Return the name of the subgrid containing `coord` or, `None` if none do.
    /// Mostly intended for debugging purposes
    fn which_subgrid_contains(&self, coord: Coor4D, margin: f64) -> Option<String>;

    /// Returns `None` if neither the grid, nor any of its sub-grids contain the point.
    /// **Contain** is in the sense of the `contains` method, i.e. the point is
    /// considered contained if it is inside a margin of `margin` grid units of
    /// the grid.
    fn at(&self, ctx: Option<&dyn Context>, at: Coor4D, margin: f64) -> Option<Coor4D>;
}

#[derive(Debug, Clone)]
pub enum GridSource {
    /// if external, we ask the context for the value
    External {
        /// The local, user, or global level
        level: usize,
        /// The offset from the top-of-file to the first byte of the grid-proper,
        /// i.e. after the header. The offset given in the `unigrid.index` file is
        /// the offset to the header, not of the grid-proper.
        /// For unigrid files written strictly sequentially, the offset of the grid
        /// will exceed the offset of the header by 64 bytes. This can, however,
        /// not be counted on.
        offset: usize,
    },
    Internal {
        values: Vec<f32>,
    },
}

/// Grid characteristics and interpolation.
///
/// The actual grid may be part of the `GridSource` struct, or
/// provided externally (presumably by a [Context](crate::context::Context)).
///
/// In principle grid format agnostic, but includes a parser for
/// geodetic grids in the Gravsoft format.
#[derive(Debug, Clone)]
pub struct BaseGrid {
    pub name: String,
    pub header: GridHeader,
    // pub lat_n: f64, // Latitude of the first (typically northernmost) row of the grid
    // pub lat_s: f64, // Latitude of the last (typically southernmost) row of the grid
    pub grid: GridSource,
    pub subgrids: Vec<BaseGrid>, // Not optional, because external grids can have subgrids too
}

impl BaseGrid {
    pub fn new(name: &str, header: GridHeader, grid: GridSource) -> Result<Self, Error> {
        let bands = header.bands;
        let rows = header.rows;
        let cols = header.cols;
        let elements = rows * cols * bands;
        if elements == 0 || bands < 1 {
            return Err(Error::General("Malformed grid - A"));
        }

        if let GridSource::Internal { values } = &grid {
            if elements > values.len() {
                return Err(Error::General("Malformed grid - B"));
            }
        }

        let subgrids = Vec::new();

        Ok(BaseGrid {
            name: name.to_string(),
            header,
            grid,
            subgrids,
        })
    }

    /// Read and decode a grid from a filesystem path.
    ///
    /// The decoder is selected from the file extension.
    pub fn read(path: &str) -> Result<Self, Error> {
        let path: std::path::PathBuf = path.into();
        let grid = std::fs::read(&path)?;
        let id = path
            .file_stem()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default();
        Self::read_bytes(id, &grid, path.extension().and_then(|ext| ext.to_str()))
    }

    /// Decode a grid from raw bytes.
    ///
    /// `name` is used as the logical grid identifier in diagnostics and grid
    /// metadata. `format_hint` is typically a file extension such as `"gsb"`
    /// or `"tif"`.
    pub fn read_bytes(name: &str, bytes: &[u8], format_hint: Option<&str>) -> Result<Self, Error> {
        match format_hint
            .unwrap_or_default()
            .to_ascii_lowercase()
            .as_str()
        {
            #[cfg(feature = "with_gtg")]
            "tif" | "tiff" => gtg::gtg(name, bytes),
            #[cfg(not(feature = "with_gtg"))]
            "tif" | "tiff" => Err(Error::Unsupported(
                "GeoTIFF grid support requires the `with_gtg` feature".to_string(),
            )),
            "gsb" => ntv2::ntv2_grid(bytes),
            "gtx" => gtx::gtx(name, bytes),
            _ => {
                if let Ok(grid) = gsa::gsa(name, bytes) {
                    Ok(grid)
                } else {
                    gravsoft::gravsoft(name, bytes)
                }
            }
        }
    }

    pub fn is_projected(&self) -> bool {
        // If any boundary is outside of [-720; 720], the grid must (by a wide margin) be
        // in projected coordinates
        [
            self.header.lat_n,
            self.header.lat_s,
            self.header.lon_w,
            self.header.lon_e,
        ]
        .iter()
        .any(|h| h.abs() > 7.0)
    }
}

/// Grid metadata: bounding box, quantization, size, dimension
#[derive(Default, Clone, PartialEq)]
pub struct GridHeader {
    pub lat_n: f64, // Latitude of the first (typically northernmost) row of the grid
    pub lat_s: f64, // Latitude of the last (typically southernmost) row of the grid
    pub lon_w: f64, // Longitude of the first (typically westernmost) column of each row
    pub lon_e: f64, // Longitude of the last (typically easternmost) column of each row
    pub dlat: f64,  // Signed distance between two consecutive rows
    pub dlon: f64,  // Signed distance between two consecutive columns
    pub rows: usize,
    pub cols: usize,
    pub bands: usize,
}

impl std::fmt::Debug for GridHeader {
    fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let scaled = if self.is_angular() {
            self.to_degrees()
        } else {
            self.clone()
        };
        fmt.debug_struct("GridHeader")
            .field("lat_n", &scaled.lat_n)
            .field("lat_s", &scaled.lat_s)
            .field("lon_w", &scaled.lon_w)
            .field("lon_e", &scaled.lon_e)
            .field("dlon", &scaled.dlon)
            .field("dlat", &scaled.dlat)
            .field("rows", &scaled.rows)
            .field("cols", &scaled.cols)
            .field("bands", &scaled.bands)
            .finish()
    }
}

impl GridHeader {
    pub fn new(
        n: f64,
        s: f64,
        w: f64,
        e: f64,
        dlat: f64,
        dlon: f64,
        bands: usize,
    ) -> Result<Self, Error> {
        let lat_n = n;
        let lat_s = s;
        let lon_w = w;
        let lon_e = e;

        let dlat = dlat.copysign(lat_s - lat_n);
        let dlon = dlon.copysign(lon_e - lon_w);

        let rows = (((lat_s - lat_n) / dlat).abs() + 1.5).floor() as usize;
        let cols = (((lon_e - lon_w) / dlon).abs() + 1.5).floor() as usize;
        let elements = rows * cols * bands;
        if elements == 0 || bands < 1 {
            return Err(Error::General("Malformed grid"));
        }

        Ok(GridHeader {
            lat_n,
            lat_s,
            lon_w,
            lon_e,
            dlat,
            dlon,
            rows,
            cols,
            bands,
        })
    }

    pub fn to_degrees(&self) -> Self {
        let mut h = self.clone();
        h.lat_n = h.lat_n.to_degrees();
        h.lat_s = h.lat_s.to_degrees();
        h.lon_w = h.lon_w.to_degrees();
        h.lon_e = h.lon_e.to_degrees();
        h.dlat = h.dlat.to_degrees();
        h.dlon = h.dlon.to_degrees();
        h
    }

    /// Heuristic determination of whether a header describes an equiangular grid.
    /// If any (potentially) angular element of the header exceeds 10 radians, we
    /// assume we have linear coordinates. Otherwise assume angular.
    /// Mostly for use in the implementation of the Debug trait, where it is useful
    /// for the comprehension, if we can present the angular elements in degrees, rather
    /// than radians.
    pub fn is_angular(&self) -> bool {
        ![
            self.lat_n, self.lat_s, self.lon_e, self.lon_w, self.dlat, self.dlon,
        ]
        .iter()
        .any(|x| (*x).abs() > 10.)
    }
}

impl Grid for BaseGrid {
    fn bands(&self) -> usize {
        self.header.bands
    }

    /// Determine whether a given coordinate falls within the grid boundaries + margin.
    /// "On the boundary" qualifies as within for westernmost and southernmost, or for
    /// all boundaries if `all_inclusive==true`.
    fn contains(&self, position: Coor4D, margin: f64, all_inclusive: bool) -> bool {
        let (mut lon, lat) = position.xy();

        // We start by assuming that the last row (latitude) is the southernmost
        let mut lat_min = self.header.lat_s;
        let mut lat_max = self.header.lat_n;
        // If it's not, we swap
        if self.header.dlat > 0. {
            (lat_min, lat_max) = (lat_max, lat_min);
        }

        let lat_grace = margin * self.header.dlat.abs();
        lat_min -= lat_grace;
        lat_max += lat_grace;
        if lat != lat.clamp(lat_min, lat_max) {
            return false;
        }

        // The default assumption is the other way round for columns (longitudes)
        let mut lon_min = self.header.lon_w;
        let mut lon_max = self.header.lon_e;

        // If it's not, we swap
        if self.header.dlon < 0. {
            (lon_min, lon_max) = (lon_max, lon_min);
        }

        if !self.is_projected() {
            lon = wrap_longitude_to_interval(lon, lon_min, lon_max);
        }
        let lon_grace = margin * self.header.dlon.abs();
        lon_min -= lon_grace;
        lon_max += lon_grace;
        if lon != lon.clamp(lon_min, lon_max) {
            return false;
        }

        // If we fell through all the way down here, we're inside the grid, but we
        // still need to take care of the boundary conventions
        if (!all_inclusive) && ((lon == lon_max) || (lat == lat_max)) {
            return false;
        }
        true
    }

    fn which_subgrid_contains(&self, coord: Coor4D, margin: f64) -> Option<String> {
        if !self.contains(coord, margin.max(1e-12), true) {
            return None;
        }
        for grid in self.subgrids.iter().rev() {
            if grid.contains(coord, margin, false) {
                return Some(grid.name.clone());
            }
        }
        Some(self.name.clone())
    }

    fn at(&self, ctx: Option<&dyn Context>, at: Coor4D, margin: f64) -> Option<Coor4D> {
        // If we're outside of the main grid, there is no chance
        // we're inside of any sub-grid
        if !self.contains(at, margin, true) {
            return None;
        };

        // We cannot handle unigrids without Context support
        if let GridSource::External { .. } = self.grid
            && ctx.is_none()
        {
            return None;
        }

        // Locate the relevant sub-grid (if any)
        for subgrid in self.subgrids.iter().rev() {
            if subgrid.contains(at, margin, false) {
                return interpolate(subgrid, ctx, at, 0.0);
            }
        }

        // If we're not inside any of the sub-grids,
        // we must be inside the main grid
        interpolate(self, ctx, at, margin)
    }
}

// Helper for `at`. Carries out the actual interpolation. Not intended for
// direct calls, and hence non-pub, and placed outside of the Grid trait
// implementation
//
// Since we store the entire grid in a single vector, the interpolation
// routine here looks strongly like a case of "writing Fortran 77 in Rust".
// It is, however, one of the cases where a more extensive use of abstractions
// leads to a significantly larger code base, much harder to maintain and
// comprehend.
fn interpolate(
    base: &BaseGrid,
    ctx: Option<&dyn Context>,
    at: Coor4D,
    margin: f64,
) -> Option<Coor4D> {
    if !base.contains(at, margin, true) {
        return None;
    };
    if let GridSource::External { .. } = base.grid
        && ctx.is_none()
    {
        return None;
    }

    let head = &base.header;

    // For now, we support top-to-bottom, left-to-right scan order only.
    // This is the common case for most non-block grid formats, with
    // NTv2 the odd man out. But since we normalize the NTv2 scan order
    // during parsing, we just cruise along here
    let dlat = head.dlat.abs();
    let dlon = head.dlon.abs();
    let lon = if base.is_projected() {
        at[0]
    } else {
        wrap_longitude_to_interval(at[0], head.lon_w.min(head.lon_e), head.lon_w.max(head.lon_e))
    };

    // The interpolation coordinate relative to the grid origin
    let rlon = lon - head.lon_w;
    let rlat = head.lat_n - at[1];

    // The (row, column) of the lower left node of the grid cell containing
    // the interpolation coordinate - or, in the case of extrapolation:
    // the nearest cell inside the grid.
    let row = (rlat / dlat).ceil() as i64;
    let col = (rlon / dlon).floor() as i64;

    let col = col.clamp(0_i64, (head.cols - 2) as i64) as usize;
    let row = row.clamp(1_i64, (head.rows - 1) as i64) as usize;

    // Linear array index of the first band element of each corner value
    #[rustfmt::skip]
    let (ll, lr, ul, ur) = (
        head.bands * (head.cols *  row      + col    ),
        head.bands * (head.cols *  row      + col + 1),
        head.bands * (head.cols * (row - 1) + col    ),
        head.bands * (head.cols * (row - 1) + col + 1),
    );

    let ll_lon = head.lon_w + col as f64 * dlon;
    let ll_lat = head.lat_n - row as f64 * dlat;

    // Cell relative, cell unit coordinates in a right handed CS
    let rlon = (lon - ll_lon) / dlon;
    let rlat = (at[1] - ll_lat) / dlat;

    // We cannot return more than 4 bands in a Coor4D,
    // so we ignore any exceeding bands
    let maxbands = head.bands.min(4);

    // Collect the grid values for the corners of the grid cell containing
    // the point of interest
    let mut corners = [Coor4D::nan(); 4];
    let corner_indices = [ll, lr, ul, ur];
    const LL: usize = 0;
    const LR: usize = 1;
    const UL: usize = 2;
    const UR: usize = 3;

    let n = match &base.grid {
        GridSource::External { .. } => {
            ctx.unwrap()
                .get_grid_values(base, &corner_indices, &mut corners)
        }
        GridSource::Internal { values } => {
            for i in 0..maxbands {
                corners[LL][i] = values[ll + i] as f64;
                corners[LR][i] = values[lr + i] as f64;
                corners[UL][i] = values[ul + i] as f64;
                corners[UR][i] = values[ur + i] as f64;
            }
            4
        }
    };

    if n != 4 {
        return None;
    }

    // Interpolate (or extrapolate, if we're outside of the physical grid)
    let mut left = Coor4D::origin();
    for i in 0..maxbands {
        let lower = corners[LL][i];
        let upper = corners[UL][i];
        left[i] = (1. - rlat) * lower + rlat * upper;
    }

    let mut right = Coor4D::origin();
    for i in 0..maxbands {
        let lower = corners[LR][i];
        let upper = corners[UR][i];
        right[i] = (1. - rlat) * lower + rlat * upper;
    }

    let mut result = Coor4D::origin();
    for i in 0..maxbands {
        result[i] = (1. - rlon) * left[i] + rlon * right[i];
    }

    // Synthetic parent grids use NaN-filled payloads as a "no value here"
    // marker, so callers fall through to a matching subgrid instead of
    // interpolating across empty space.
    if (0..maxbands).any(|i| !result[i].is_finite()) {
        return None;
    }

    Some(result)
}

fn wrap_longitude_to_interval(mut lon: f64, lon_min: f64, lon_max: f64) -> f64 {
    let center = 0.5 * (lon_min + lon_max);
    let two_pi = std::f64::consts::TAU;

    while lon - center > std::f64::consts::PI {
        lon -= two_pi;
    }
    while lon - center < -std::f64::consts::PI {
        lon += two_pi;
    }

    lon
}

/// Find the most appropriate grid value from a stack (i.e. slice) of grids.
/// Search the grids in slice order and return the first hit.
/// If no hits are found, try once more, this time adding a half grid-cell
/// margin around each grid
pub fn grids_at(
    ctx: Option<&dyn Context>,
    grids: &[Arc<BaseGrid>],
    coord: Coor4D,
    use_null_grid: bool,
) -> Option<Coor4D> {
    for margin in [0.0, 0.5] {
        for grid in grids.iter() {
            let d = grid.at(ctx, coord, margin);
            if d.is_some() {
                return d;
            }
        }
    }

    if use_null_grid {
        return Some(Coor4D::origin());
    }

    None
}

// ----- T E S T S ------------------------------------------------------------------

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;
    use crate::coordinate::AngularUnits;

    const DATUM_HEADER: GridHeader = GridHeader {
        lat_n: 58f64.to_radians(),
        lat_s: 54f64.to_radians(),
        lon_w:  8f64.to_radians(),
        lon_e: 16f64.to_radians(),
        dlat: -1f64.to_radians(),
        dlon:   1f64.to_radians(),
        rows:   5_usize,
        cols:   9_usize,
        bands:  2,
    };

    const GEOID_HEADER: GridHeader = GridHeader {
        lat_n: 58f64.to_radians(),
        lat_s: 54f64.to_radians(),
        lon_w:  8f64.to_radians(),
        lon_e: 16f64.to_radians(),
        dlat: -1f64.to_radians(),
        dlon:  1f64.to_radians(),
        rows:   5_usize,
        cols:   9_usize,
        bands:  1,
    };

    const DATUM: [f32; 5*2*9] = [
        08., 58., 09., 58., 10., 58., 11., 58., 12., 58., 13., 58., 14., 58., 15., 58., 16., 58.,
        08., 57., 09., 57., 10., 57., 11., 57., 12., 57., 13., 57., 14., 57., 15., 57., 16., 57.,
        08., 56., 09., 56., 10., 56., 11., 56., 12., 56., 13., 56., 14., 56., 15., 56., 16., 56.,
        08., 55., 09., 55., 10., 55., 11., 55., 12., 55., 13., 55., 14., 55., 15., 55., 16., 55.,
        08., 54., 09., 54., 10., 54., 11., 54., 12., 54., 13., 54., 14., 54., 15., 54., 16., 54.,
    ];

    const GEOID: [f32; 5*9] = [
        58.08, 58.09, 58.10, 58.11, 58.12, 58.13, 58.14, 58.15, 58.16,
        57.08, 57.09, 57.10, 57.11, 57.12, 57.13, 57.14, 57.15, 57.16,
        56.08, 56.09, 56.10, 56.11, 56.12, 56.13, 56.14, 56.15, 56.16,
        55.08, 55.09, 55.10, 55.11, 55.12, 55.13, 55.14, 55.15, 55.16,
        54.08, 54.09, 54.10, 54.11, 54.12, 54.13, 54.14, 54.15, 54.16,
    ];

    #[test]
    fn grid_header() -> Result<(), Error> {
        // Create a datum correction grid (2 bands)
        let datum_header = DATUM_HEADER;
        let datum_grid = Vec::from(DATUM);
        let datum = BaseGrid::new(
            "hohoho",
            datum_header,
            GridSource::Internal { values: datum_grid },
        )?;

        // Extrapolation
        let c = Coor4D::geo(50., 100., 0., 0.);
        // ...with output in arcsec
        let d = datum.at(None, c, 100.0).unwrap();

        // The grid is constructed to make the position in degrees equal to
        // the extrapolation value in arcsec.
        // Even for this case of extreme extrapolation, we expect the difference
        // to be insignificant
        let ellps = Ellipsoid::named("GRS80")?;
        let dist = ellps.distance(&c, &d.to_radians());
        assert!(dist < 1e-6);

        // Interpolation
        let c = Coor4D::geo(55.06, 12.03, 0., 0.);
        // Check that we're not extrapolating
        assert!(datum.contains(c, 0.0, true));
        // ...with output in arcsec
        let d = datum.at(None, c, 0.0).unwrap();
        // We can do slightly better for interpolation than for extrapolation,
        // but the grid values are f32, so we have only approx 7 significant
        // figures...
        let dist = ellps.distance(&c, &d.to_radians());
        assert!(dist < 1e-6);

        // Create a geoid grid (1 band)
        let geoid_header = GEOID_HEADER;
        let geoid_grid = Vec::from(GEOID);
        let geoid = BaseGrid::new(
            "geoid",
            geoid_header,
            GridSource::Internal { values: geoid_grid },
        )?;

        let c = Coor4D::geo(58.75, 08.25, 0., 0.);
        assert!(!geoid.contains(c, 0.0, true));
        assert!(geoid.contains(c, 1.0, true));

        let n = geoid.at(None, c, 1.0).unwrap();
        assert!((n[0] - (58.75 + 0.0825)).abs() < 0.0001);
        Ok(())
    }

    #[test]
    fn read_base_grid() {
        // trying to read non-existing files fail with Error::Io
        assert!(matches!(BaseGrid::read(""), Err(Error::Io(_))));
        assert!(matches!(BaseGrid::read("not.existing"), Err(Error::Io(_))));

        // BaseGrid::read goes only after the exact filename specified,
        // it does not chase the file down the "known paths"
        assert!(BaseGrid::read("test.datum").is_err());
        assert!(BaseGrid::read("geodesy/datum/test.datum").is_ok());
    }

    #[test]
    fn projected_grids_do_not_wrap_longitude_like_angles() -> Result<(), Error> {
        let projected = GridHeader {
            lat_n: 20.0,
            lat_s: 0.0,
            lon_w: 1000.0,
            lon_e: 1010.0,
            dlat: -10.0,
            dlon: 10.0,
            rows: 3,
            cols: 2,
            bands: 1,
        };
        let grid = BaseGrid::new(
            "projected",
            projected,
            GridSource::Internal {
                values: vec![0.0; 3 * 2],
            },
        )?;

        assert!(!grid.contains(Coor4D::raw(1020.0, 10.0, 0.0, 0.0), 0.0, true));
        Ok(())
    }

    #[test]
    fn at_prefers_most_specific_subgrid() -> Result<(), Error> {
        let parent = BaseGrid::new(
            "parent",
            DATUM_HEADER,
            GridSource::Internal {
                values: vec![0.0; DATUM.len()],
            },
        )?;
        let child_header = GridHeader::new(
            57f64.to_radians(),
            55f64.to_radians(),
            10f64.to_radians(),
            12f64.to_radians(),
            1f64.to_radians(),
            1f64.to_radians(),
            2,
        )?;
        let child = BaseGrid::new(
            "child",
            child_header,
            GridSource::Internal {
                values: vec![42.0; 3 * 3 * 2],
            },
        )?;

        let mut root = parent;
        root.subgrids.push(child);

        let point = Coor4D::geo(56.5, 10.5, 0.0, 0.0);
        assert_eq!(root.which_subgrid_contains(point, 0.0), Some("child".to_string()));
        let value = root.at(None, point, 0.0).unwrap();
        assert_eq!(value[0], 42.0);
        assert_eq!(value[1], 42.0);
        Ok(())
    }
}

// Additional tests for Grid in src/inner_op/gridshift.rs
