//! Read a minimal subset of PROJ Geodetic TIFF Grids (GTG).
//!
//! This reader currently supports the horizontal geographic case used by the
//! first GeoTIFF-backed `hgridshift` corpus pipelines:
//! - `TYPE=HORIZONTAL_OFFSET`
//! - latitude/longitude offset bands described in arc-seconds
//! - regular grids georeferenced by ModelPixelScale + ModelTiepoint

use super::{BaseGrid, GridHeader, GridKind, GridSource};
use crate::{Error, coord::Coor4D};
use std::io::Cursor;
use tiff::decoder::{Decoder, DecodingResult};
use tiff::tags::Tag;

const GDAL_METADATA_TAG: Tag = Tag::Unknown(42112);

pub fn gtg(name: &str, buf: &[u8]) -> Result<BaseGrid, Error> {
    let mut decoder = Decoder::new(Cursor::new(buf)).map_err(tiff_err)?;
    let mut images: Vec<BaseGrid> = Vec::new();
    let mut metas: Vec<GtgMetadata> = Vec::new();

    loop {
        let fallback = metas.last();
        let (image, meta) = decode_gtg_ifd(&mut decoder, name, fallback)?;
        images.push(image);
        metas.push(meta);
        if !decoder.more_images() {
            break;
        }
        decoder.next_image().map_err(tiff_err)?;
    }

    if images.is_empty() {
        return Err(Error::Unsupported(
            "GTG file did not contain any images".to_string(),
        ));
    }

    // Some files encode horizontal and height offsets in separate IFDs rather
    // than as bands in a single IFD (e.g. NADCON5 Alaska format). Attach the
    // height grid as a runtime companion so it is interpolated directly at the
    // query coordinate rather than being pre-resampled onto the horizontal grid.
    if images.len() == 2
        && images[0].header.bands == 2
        && images[1].header.bands == 1
        && metas[1]
            .grid_type
            .eq_ignore_ascii_case("ELLIPSOIDAL_HEIGHT_OFFSET")
    {
        return attach_height_grid(images.remove(0), images.remove(0));
    }

    if images.len() == 1 {
        return Ok(images.remove(0));
    }

    let mut parent = synthetic_parent_grid(name, &images)?;
    parent.subgrids = images;
    Ok(parent)
}

fn decode_gtg_ifd<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
    default_name: &str,
    fallback_meta: Option<&GtgMetadata>,
) -> Result<(BaseGrid, GtgMetadata), Error> {
    let (cols, rows) = decoder.dimensions().map_err(tiff_err)?;
    let cols = cols as usize;
    let rows = rows as usize;

    let pixel_scale = decoder
        .get_tag_f64_vec(Tag::ModelPixelScaleTag)
        .map_err(tiff_err)?;
    let tiepoint = decoder
        .get_tag_f64_vec(Tag::ModelTiepointTag)
        .map_err(tiff_err)?;
    if pixel_scale.len() < 2 || tiepoint.len() < 6 {
        return Err(Error::Unsupported(
            "GTG reader requires ModelPixelScale and ModelTiepoint tags".to_string(),
        ));
    }

    let metadata = decoder
        .get_tag_ascii_string(GDAL_METADATA_TAG)
        .map_err(tiff_err)?;
    let layout = match decoder
        .find_tag_unsigned::<u16>(Tag::PlanarConfiguration)
        .map_err(tiff_err)?
        .unwrap_or(1)
    {
        1 => SampleLayout::Chunky,
        2 => SampleLayout::Separate,
        other => {
            return Err(Error::Unsupported(format!(
                "GTG reader does not support PlanarConfiguration={other}"
            )));
        }
    };

    let samples_per_pixel = decoder
        .find_tag_unsigned::<u16>(Tag::SamplesPerPixel)
        .map_err(tiff_err)?
        .unwrap_or(1) as usize;
    let raster_type = decode_raster_type(decoder)?;
    let crs_type = decode_crs_type(decoder)?;

    let meta = GtgMetadata::from_xml(&metadata, crs_type, layout, fallback_meta)?;
    let n_out = meta.output_bands(crs_type)?;
    let grid_name = meta
        .grid_name
        .clone()
        .unwrap_or_else(|| default_name.to_string());

    // Compute the minimum number of input bands required.
    let required_bands = match n_out {
        1 => meta.height_band.unwrap_or(0) + 1,
        _ => {
            [meta.lat_band, meta.lon_band, meta.height_band.unwrap_or(0)]
                .iter()
                .max()
                .copied()
                .unwrap_or(0)
                + 1
        }
    };
    if samples_per_pixel < required_bands {
        return Err(Error::Unsupported(format!(
            "GTG reader expected at least {required_bands} float32 bands, found {samples_per_pixel}"
        )));
    }

    let values = read_f32_raster(decoder, layout, rows, cols, samples_per_pixel)?;

    let plane_len = rows * cols;
    if values.len() % plane_len != 0 {
        return Err(Error::Unsupported(format!(
            "GTG reader expected whole float32 bands, found {} samples for {} pixels",
            values.len(),
            plane_len
        )));
    }
    let decoded_bands = values.len() / plane_len;
    if decoded_bands < required_bands {
        return Err(Error::Unsupported(format!(
            "GTG reader expected at least {required_bands} float32 bands, found {decoded_bands}"
        )));
    }

    let (dlat, dlon, lon_w, lat_n, lat_s, lon_e) = match crs_type {
        GridCrsType::Geographic => {
            let dlat = -pixel_scale[1].to_radians();
            let dlon = pixel_scale[0].to_radians();
            let (lon_w, lat_n) = match raster_type {
                RasterType::PixelIsArea => (
                    tiepoint[3].to_radians() + 0.5 * dlon,
                    tiepoint[4].to_radians() + 0.5 * dlat,
                ),
                RasterType::PixelIsPoint => (tiepoint[3].to_radians(), tiepoint[4].to_radians()),
            };
            let lat_s =
                (tiepoint[4] - (rows.saturating_sub(1) as f64) * pixel_scale[1]).to_radians();
            let lon_e =
                (tiepoint[3] + (cols.saturating_sub(1) as f64) * pixel_scale[0]).to_radians();
            let (lat_s, lon_e) = match raster_type {
                RasterType::PixelIsArea => (lat_s + 0.5 * dlat, lon_e + 0.5 * dlon),
                RasterType::PixelIsPoint => (lat_s, lon_e),
            };
            (dlat, dlon, lon_w, lat_n, lat_s, lon_e)
        }
        GridCrsType::Projected => {
            let dlat = -pixel_scale[1];
            let dlon = pixel_scale[0];
            let (lon_w, lat_n) = match raster_type {
                RasterType::PixelIsArea => (tiepoint[3] + 0.5 * dlon, tiepoint[4] + 0.5 * dlat),
                RasterType::PixelIsPoint => (tiepoint[3], tiepoint[4]),
            };
            let lat_s = tiepoint[4] - (rows.saturating_sub(1) as f64) * pixel_scale[1];
            let lon_e = tiepoint[3] + (cols.saturating_sub(1) as f64) * pixel_scale[0];
            let (lat_s, lon_e) = match raster_type {
                RasterType::PixelIsArea => (lat_s + 0.5 * dlat, lon_e + 0.5 * dlon),
                RasterType::PixelIsPoint => (lat_s, lon_e),
            };
            (dlat, dlon, lon_w, lat_n, lat_s, lon_e)
        }
    };
    let header = GridHeader {
        lat_n,
        lat_s,
        lon_w,
        lon_e,
        dlat,
        dlon,
        rows,
        cols,
        bands: n_out,
    };

    let mut grid = Vec::with_capacity(plane_len * n_out);
    match (layout, n_out) {
        // ── 1-band vertical offset ────────────────────────────────────────────
        (SampleLayout::Chunky, 1) => {
            let vb = meta.height_band.unwrap_or(0);
            for node in values.chunks_exact(decoded_bands) {
                grid.push(node[vb]);
            }
        }
        (SampleLayout::Separate, 1) => {
            let vb = meta.height_band.unwrap_or(0);
            let plane = &values[vb * plane_len..(vb + 1) * plane_len];
            grid.extend_from_slice(plane);
        }

        // ── 2-band horizontal offset ──────────────────────────────────────────
        (SampleLayout::Chunky, 2) => {
            for node in values.chunks_exact(decoded_bands) {
                let lon = apply_axis_polarity(node[meta.lon_band], meta.lon_positive_east);
                let lat = apply_axis_polarity(node[meta.lat_band], meta.lat_positive_north);
                grid.push(lon);
                grid.push(lat);
            }
        }
        (SampleLayout::Separate, 2) => {
            let lon_plane = &values[meta.lon_band * plane_len..(meta.lon_band + 1) * plane_len];
            let lat_plane = &values[meta.lat_band * plane_len..(meta.lat_band + 1) * plane_len];
            for i in 0..plane_len {
                grid.push(apply_axis_polarity(lon_plane[i], meta.lon_positive_east));
                grid.push(apply_axis_polarity(lat_plane[i], meta.lat_positive_north));
            }
        }

        // ── 3-band geographic-3D offset (horizontal + ellipsoidal height) ─────
        (SampleLayout::Chunky, 3) => {
            let hb = meta.height_band.unwrap();
            for node in values.chunks_exact(decoded_bands) {
                let lon = apply_axis_polarity(node[meta.lon_band], meta.lon_positive_east);
                let lat = apply_axis_polarity(node[meta.lat_band], meta.lat_positive_north);
                grid.push(lon);
                grid.push(lat);
                grid.push(node[hb]); // height offset in metres — no polarity convention
            }
        }
        (SampleLayout::Separate, 3) => {
            let hb = meta.height_band.unwrap();
            let lon_plane = &values[meta.lon_band * plane_len..(meta.lon_band + 1) * plane_len];
            let lat_plane = &values[meta.lat_band * plane_len..(meta.lat_band + 1) * plane_len];
            let h_plane = &values[hb * plane_len..(hb + 1) * plane_len];
            for i in 0..plane_len {
                grid.push(apply_axis_polarity(lon_plane[i], meta.lon_positive_east));
                grid.push(apply_axis_polarity(lat_plane[i], meta.lat_positive_north));
                grid.push(h_plane[i]);
            }
        }

        // n_out is 1, 2, or 3 — other values are unreachable
        _ => unreachable!("output_bands() returned unexpected value {n_out}"),
    }

    let bias = match (crs_type, n_out) {
        (_, 1) => Coor4D::raw(meta.height_constant, 0.0, 0.0, 0.0),
        (GridCrsType::Projected, _) => Coor4D::raw(meta.lon_constant, meta.lat_constant, 0.0, 0.0),
        (GridCrsType::Geographic, _) => Coor4D::origin(),
    };
    let allow_margin_fallback = !(n_out == 1
        && crs_type == GridCrsType::Geographic
        && meta
            .grid_type
            .eq_ignore_ascii_case("VERTICAL_OFFSET_GEOGRAPHIC_TO_VERTICAL"));
    Ok((
        BaseGrid::new(&grid_name, header, GridSource::Internal { values: grid })?
            .with_kind(meta.grid_kind())
            .with_projected(crs_type == GridCrsType::Projected)
            .with_bias(bias)
            .with_margin_fallback(allow_margin_fallback),
        meta,
    ))
}

fn read_f32_raster<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
    layout: SampleLayout,
    rows: usize,
    cols: usize,
    samples_per_pixel: usize,
) -> Result<Vec<f32>, Error> {
    let (chunk_width, chunk_height) = decoder.chunk_dimensions();
    let chunk_width = chunk_width as usize;
    let chunk_height = chunk_height as usize;
    match layout {
        SampleLayout::Chunky => read_chunky_f32_raster(
            decoder,
            rows,
            cols,
            samples_per_pixel,
            chunk_width,
            chunk_height,
        ),
        SampleLayout::Separate => read_planar_f32_raster(
            decoder,
            rows,
            cols,
            samples_per_pixel,
            chunk_width,
            chunk_height,
        ),
    }
}

fn read_chunky_f32_raster<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
    rows: usize,
    cols: usize,
    samples_per_pixel: usize,
    chunk_width: usize,
    chunk_height: usize,
) -> Result<Vec<f32>, Error> {
    let plane_len = rows * cols;
    let chunks_across = cols.div_ceil(chunk_width);
    let chunks_down = rows.div_ceil(chunk_height);
    let mut values = vec![0.0f32; plane_len * samples_per_pixel];

    for chunk_index in 0..(chunks_across * chunks_down) {
        let decoded = decoder.read_chunk(chunk_index as u32).map_err(tiff_err)?;
        let DecodingResult::F32(chunk) = decoded else {
            return Err(Error::Unsupported(
                "GTG reader currently supports only float32 rasters".to_string(),
            ));
        };

        let (chunk_cols, chunk_rows) = decoder.chunk_data_dimensions(chunk_index as u32);
        let chunk_cols = chunk_cols as usize;
        let chunk_rows = chunk_rows as usize;
        let chunk_x = chunk_index % chunks_across;
        let chunk_y = chunk_index / chunks_across;

        for local_row in 0..chunk_rows {
            let dst_row = chunk_y * chunk_height + local_row;
            let dst_col = chunk_x * chunk_width;
            let src_offset = local_row * chunk_cols * samples_per_pixel;
            let dst_offset = (dst_row * cols + dst_col) * samples_per_pixel;
            let row_len = chunk_cols * samples_per_pixel;
            values[dst_offset..dst_offset + row_len]
                .copy_from_slice(&chunk[src_offset..src_offset + row_len]);
        }
    }

    Ok(values)
}

fn read_planar_f32_raster<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
    rows: usize,
    cols: usize,
    samples_per_pixel: usize,
    chunk_width: usize,
    chunk_height: usize,
) -> Result<Vec<f32>, Error> {
    let plane_len = rows * cols;
    let chunks_across = cols.div_ceil(chunk_width);
    let chunks_down = rows.div_ceil(chunk_height);
    let chunks_per_plane = chunks_across * chunks_down;
    let mut values = vec![0.0f32; plane_len * samples_per_pixel];

    for band in 0..samples_per_pixel {
        for chunk_index in 0..chunks_per_plane {
            let decoded = decoder
                .read_chunk((band * chunks_per_plane + chunk_index) as u32)
                .map_err(tiff_err)?;
            let DecodingResult::F32(chunk) = decoded else {
                return Err(Error::Unsupported(
                    "GTG reader currently supports only float32 rasters".to_string(),
                ));
            };

            let (chunk_cols, chunk_rows) = decoder.chunk_data_dimensions(chunk_index as u32);
            let chunk_cols = chunk_cols as usize;
            let chunk_rows = chunk_rows as usize;
            let chunk_x = chunk_index % chunks_across;
            let chunk_y = chunk_index / chunks_across;

            for local_row in 0..chunk_rows {
                let dst_row = chunk_y * chunk_height + local_row;
                let dst_col = chunk_x * chunk_width;
                let dst_offset = band * plane_len + dst_row * cols + dst_col;
                let src_offset = local_row * chunk_cols;
                values[dst_offset..dst_offset + chunk_cols]
                    .copy_from_slice(&chunk[src_offset..src_offset + chunk_cols]);
            }
        }
    }

    Ok(values)
}

fn synthetic_parent_grid(name: &str, images: &[BaseGrid]) -> Result<BaseGrid, Error> {
    let mut lat_n = f64::NEG_INFINITY;
    let mut lat_s = f64::INFINITY;
    let mut lon_w = f64::INFINITY;
    let mut lon_e = f64::NEG_INFINITY;
    let mut dlat = f64::INFINITY;
    let mut dlon = f64::INFINITY;
    let bands = images[0].header.bands;

    for image in images {
        let header = &image.header;
        lat_n = lat_n.max(header.lat_n);
        lat_s = lat_s.min(header.lat_s);
        lon_w = lon_w.min(header.lon_w);
        lon_e = lon_e.max(header.lon_e);
        dlat = dlat.min(header.dlat.abs());
        dlon = dlon.min(header.dlon.abs());
    }

    let header = GridHeader::new(lat_n, lat_s, lon_w, lon_e, -dlat, dlon, bands)?;
    // This synthetic parent only routes lookups into child pages. Its data is
    // intentionally NaN-filled so interpolation does not fabricate values
    // between disjoint subgrids.
    let values = vec![f32::NAN; header.rows * header.cols * header.bands];
    BaseGrid::new(name, header, GridSource::Internal { values })
        .map(|grid| grid.with_kind(images[0].kind))
}

/// Merge a 2-band horizontal grid and a 1-band height grid into a 3-band grid.
///
/// Used for files that encode horizontal and ellipsoidal height offsets in
/// separate IFDs (e.g. NADCON5 Alaska) rather than as three bands in a single IFD.
/// The height grid may be coarser than the horizontal grid, but must be able
/// to provide an interpolated value at every horizontal grid node.
fn attach_height_grid(horiz: BaseGrid, height: BaseGrid) -> Result<BaseGrid, Error> {
    if horiz.kind != GridKind::HorizontalOffset || height.kind != GridKind::EllipsoidalHeightOffset
    {
        return Err(Error::Unsupported(
            "GTG: attach_height_grid requires horizontal and ellipsoidal-height grids".to_string(),
        ));
    }
    Ok(horiz.with_vertical_companion(height))
}

fn apply_axis_polarity(value: f32, positive_is_forward: bool) -> f32 {
    if positive_is_forward { value } else { -value }
}

fn tiff_err(err: tiff::TiffError) -> Error {
    Error::Unsupported(format!("GTG TIFF decode failed: {err}"))
}

fn decode_raster_type<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
) -> Result<RasterType, Error> {
    let keys = read_geo_keys(decoder)?;
    if keys.len() < 4 {
        return Ok(RasterType::PixelIsPoint);
    }

    let number_of_keys = keys[3] as usize;
    for i in 0..number_of_keys {
        let base = 4 + 4 * i;
        if base + 3 >= keys.len() {
            break;
        }
        if keys[base] != 1025 {
            continue;
        }
        return Ok(match keys[base + 3] {
            1 => RasterType::PixelIsArea,
            2 => RasterType::PixelIsPoint,
            other => {
                return Err(Error::Unsupported(format!(
                    "GTG reader does not support RasterTypeGeoKey={other}"
                )));
            }
        });
    }

    Ok(RasterType::PixelIsPoint)
}

fn decode_crs_type<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
) -> Result<GridCrsType, Error> {
    let keys = read_geo_keys(decoder)?;
    if keys.len() < 4 {
        return Ok(GridCrsType::Geographic);
    }

    let number_of_keys = keys[3] as usize;
    let mut model_type = None;
    let mut geographic = false;
    let mut projected = false;
    for i in 0..number_of_keys {
        let base = 4 + 4 * i;
        if base + 3 >= keys.len() {
            break;
        }
        match keys[base] {
            1024 => model_type = Some(keys[base + 3]),
            2048 => geographic = true,
            3072 => projected = true,
            _ => {}
        }
    }

    if let Some(model_type) = model_type {
        return Ok(match model_type {
            1 => GridCrsType::Projected,
            2 | 3 => GridCrsType::Geographic,
            other => {
                return Err(Error::Unsupported(format!(
                    "GTG reader does not support GTModelTypeGeoKey={other}"
                )));
            }
        });
    }
    if projected {
        return Ok(GridCrsType::Projected);
    }
    if geographic {
        return Ok(GridCrsType::Geographic);
    }
    Ok(GridCrsType::Geographic)
}

fn read_geo_keys<R: std::io::Read + std::io::Seek>(
    decoder: &mut Decoder<R>,
) -> Result<Vec<u16>, Error> {
    decoder
        .get_tag_u16_vec(Tag::GeoKeyDirectoryTag)
        .map_err(tiff_err)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SampleLayout {
    Chunky,
    Separate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RasterType {
    PixelIsArea,
    PixelIsPoint,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GridCrsType {
    Geographic,
    Projected,
}

#[derive(Debug, Clone)]
struct GtgMetadata {
    grid_name: Option<String>,
    lat_band: usize,
    lon_band: usize,
    /// Band index for ellipsoidal height offset (GEOGRAPHIC_3D_OFFSET) or the
    /// single value band (VERTICAL_OFFSET_GEOGRAPHIC_TO_VERTICAL).
    height_band: Option<usize>,
    height_description: Option<String>,
    lat_positive_north: bool,
    lon_positive_east: bool,
    lat_unit: String,
    lon_unit: String,
    height_constant: f64,
    lat_constant: f64,
    lon_constant: f64,
    grid_type: String,
}

impl GtgMetadata {
    fn from_xml(
        xml: &str,
        crs_type: GridCrsType,
        _layout: SampleLayout,
        fallback: Option<&Self>,
    ) -> Result<Self, Error> {
        let grid_type = find_item(xml, "TYPE", None, None)
            .or_else(|| fallback.map(|meta| meta.grid_type.clone()))
            .ok_or_else(|| Error::Unsupported("GTG metadata is missing TYPE".to_string()))?;
        let grid_name = find_item(xml, "grid_name", None, None);

        // Vertical grids have a single offset band with a different description.
        // All other recognised types are horizontal (2-band or 3-band).
        if grid_type.eq_ignore_ascii_case("VERTICAL_OFFSET_GEOGRAPHIC_TO_VERTICAL") {
            // Use whichever band the first DESCRIPTION item describes, or 0.
            let (height_band, height_description) = find_first_band_with_description(xml)
                .or_else(|| {
                    fallback.and_then(|m| {
                        m.height_band
                            .map(|band| (band, m.height_description.clone().unwrap_or_default()))
                    })
                })
                .unwrap_or((0, String::new()));
            return Ok(Self {
                grid_name,
                lat_band: 0,
                lon_band: 0,
                height_band: Some(height_band),
                height_description: Some(height_description),
                lat_positive_north: true,
                lon_positive_east: true,
                lat_unit: "metre".to_string(),
                lon_unit: "metre".to_string(),
                height_constant: find_item(xml, "constant_offset", Some(height_band), None)
                    .and_then(|value| value.parse::<f64>().ok())
                    .or_else(|| fallback.map(|m| m.height_constant))
                    .unwrap_or(0.0),
                lat_constant: 0.0,
                lon_constant: 0.0,
                grid_type,
            });
        }

        // Horizontal grids: DESCRIPTION items are optional — fall back to convention.
        let lat_description = match crs_type {
            GridCrsType::Geographic => "latitude_offset",
            GridCrsType::Projected => "northing_offset",
        };
        let lon_description = match crs_type {
            GridCrsType::Geographic => "longitude_offset",
            GridCrsType::Projected => "easting_offset",
        };
        let lat_band = find_band(xml, lat_description)
            .ok()
            .or_else(|| fallback.map(|m| m.lat_band))
            .unwrap_or(match crs_type {
                GridCrsType::Geographic => 0,
                GridCrsType::Projected => 1,
            });
        let lon_band = find_band(xml, lon_description)
            .ok()
            .or_else(|| fallback.map(|m| m.lon_band))
            .unwrap_or(0);

        // GEOGRAPHIC_3D_OFFSET carries an additional ellipsoidal height offset band.
        let height_band = if grid_type.eq_ignore_ascii_case("GEOGRAPHIC_3D_OFFSET") {
            find_band(xml, "ellipsoidal_height_offset")
                .ok()
                .or_else(|| fallback.and_then(|m| m.height_band))
        } else {
            None
        };

        // UNITTYPE is optional; defaults depend on the grid CRS.
        let lat_unit = find_item(xml, "UNITTYPE", Some(lat_band), Some("unittype"))
            .or_else(|| fallback.map(|m| m.lat_unit.clone()))
            .unwrap_or_else(|| match crs_type {
                GridCrsType::Geographic => "arc-second".to_string(),
                GridCrsType::Projected => "metre".to_string(),
            });
        let lon_unit = find_item(xml, "UNITTYPE", Some(lon_band), Some("unittype"))
            .or_else(|| fallback.map(|m| m.lon_unit.clone()))
            .unwrap_or_else(|| match crs_type {
                GridCrsType::Geographic => "arc-second".to_string(),
                GridCrsType::Projected => "metre".to_string(),
            });

        let lat_constant = find_item(xml, "constant_offset", Some(lat_band), None)
            .and_then(|value| value.parse::<f64>().ok())
            .or_else(|| fallback.map(|m| m.lat_constant))
            .unwrap_or(0.0);
        let lon_constant = find_item(xml, "constant_offset", Some(lon_band), None)
            .and_then(|value| value.parse::<f64>().ok())
            .or_else(|| fallback.map(|m| m.lon_constant))
            .unwrap_or(0.0);

        let lat_positive_north = match find_item(xml, "positive_value", Some(lat_band), None) {
            Some(value) if value.eq_ignore_ascii_case("south") => false,
            Some(_) => true,
            None => fallback.map(|meta| meta.lat_positive_north).unwrap_or(true),
        };
        let lon_positive_east = match find_item(xml, "positive_value", Some(lon_band), None) {
            Some(value) if value.eq_ignore_ascii_case("west") => false,
            Some(_) => true,
            None => fallback.map(|meta| meta.lon_positive_east).unwrap_or(true),
        };

        Ok(Self {
            grid_name,
            lat_band,
            lon_band,
            height_band,
            height_description: height_band.map(|_| "ellipsoidal_height_offset".to_string()),
            lat_positive_north,
            lon_positive_east,
            lat_unit,
            lon_unit,
            height_constant: 0.0,
            lat_constant,
            lon_constant,
            grid_type,
        })
    }

    /// Returns the number of output bands this grid will store:
    /// - 1 for vertical-offset grids
    /// - 2 for horizontal-offset grids
    /// - 3 for 3-D geographic-offset grids (horizontal + ellipsoidal height)
    ///
    /// Also validates that horizontal unit conventions are satisfied.
    fn output_bands(&self, crs_type: GridCrsType) -> Result<usize, Error> {
        match self.grid_type.to_ascii_lowercase().as_str() {
            "horizontal_offset" => {
                self.validate_horizontal_units(crs_type)?;
                Ok(2)
            }
            "geographic_3d_offset" => {
                self.validate_horizontal_units(GridCrsType::Geographic)?;
                Ok(if self.height_band.is_some() { 3 } else { 2 })
            }
            "vertical_offset_geographic_to_vertical" | "ellipsoidal_height_offset" => Ok(1),
            other => Err(Error::Unsupported(format!(
                "GTG reader does not support TYPE={other}"
            ))),
        }
    }

    fn grid_kind(&self) -> GridKind {
        match self.grid_type.to_ascii_lowercase().as_str() {
            "vertical_offset_geographic_to_vertical" => {
                GridKind::VerticalOffsetGeographicToVertical
            }
            "ellipsoidal_height_offset" => GridKind::EllipsoidalHeightOffset,
            "geographic_3d_offset" => GridKind::Geographic3DOffset,
            _ => GridKind::HorizontalOffset,
        }
    }

    fn validate_horizontal_units(&self, crs_type: GridCrsType) -> Result<(), Error> {
        let expected = match crs_type {
            GridCrsType::Geographic => "arc-second",
            GridCrsType::Projected => "metre",
        };
        if !self.lat_unit.eq_ignore_ascii_case(expected)
            || !self.lon_unit.eq_ignore_ascii_case(expected)
        {
            return Err(Error::Unsupported(format!(
                "GTG reader requires {expected} horizontal offsets, found lat={} lon={}",
                self.lat_unit, self.lon_unit
            )));
        }
        Ok(())
    }
}

/// Return the band index and value of the first `DESCRIPTION role="description"`
/// item found, regardless of its name. Used for vertical grids where the band
/// name varies.
fn find_first_band_with_description(xml: &str) -> Option<(usize, String)> {
    for line in xml.lines() {
        let line = line.trim();
        if !line.starts_with("<Item ") {
            continue;
        }
        let Some((open, value)) = split_item(line) else {
            continue;
        };
        if attr_value(open, "name") != Some("DESCRIPTION") {
            continue;
        }
        if attr_value(open, "role") != Some("description") {
            continue;
        }
        if let Some(sample) = attr_value(open, "sample") {
            if let Ok(idx) = sample.parse::<usize>() {
                return Some((idx, value.to_string()));
            }
        }
    }
    None
}

fn find_band(xml: &str, description: &str) -> Result<usize, Error> {
    for line in xml.lines() {
        let line = line.trim();
        if !line.starts_with("<Item ") {
            continue;
        }
        let Some((open, value)) = split_item(line) else {
            continue;
        };
        if attr_value(open, "name") != Some("DESCRIPTION") {
            continue;
        }
        if attr_value(open, "role") != Some("description") {
            continue;
        }
        if value != description {
            continue;
        }
        let Some(sample) = attr_value(open, "sample") else {
            continue;
        };
        return sample.parse::<usize>().map_err(|_| {
            Error::Unsupported(format!("GTG metadata sample index is malformed: {sample}"))
        });
    }
    Err(Error::Unsupported(format!(
        "GTG metadata is missing DESCRIPTION for {description}"
    )))
}

fn find_item(xml: &str, name: &str, sample: Option<usize>, role: Option<&str>) -> Option<String> {
    for line in xml.lines() {
        let line = line.trim();
        if !line.starts_with("<Item ") {
            continue;
        }
        let Some((open, value)) = split_item(line) else {
            continue;
        };
        if attr_value(open, "name") != Some(name) {
            continue;
        }
        if sample.map(|v| v.to_string()).as_deref() != attr_value(open, "sample") {
            continue;
        }
        if role != attr_value(open, "role") {
            continue;
        }
        return Some(value.to_string());
    }
    None
}

fn split_item(line: &str) -> Option<(&str, &str)> {
    let (open, rest) = line.split_once('>')?;
    let (value, _) = rest.rsplit_once('<')?;
    Some((open, value.trim()))
}

fn attr_value<'a>(open_tag: &'a str, attr: &str) -> Option<&'a str> {
    let needle = format!(r#"{attr}=""#);
    let start = open_tag.find(&needle)? + needle.len();
    let end = open_tag[start..].find('"')?;
    Some(&open_tag[start..start + end])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_local_geotiff_grid() -> Result<(), Error> {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../cache/grids/us_noaa_conus.tif"
        );
        let grid = BaseGrid::read(path)?;
        assert_eq!(grid.header.bands, 2);
        assert!(grid.header.rows > 0);
        assert!(grid.header.cols > 0);
        Ok(())
    }
}
