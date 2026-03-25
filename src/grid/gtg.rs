//! Read a minimal subset of PROJ Geodetic TIFF Grids (GTG).
//!
//! This reader currently supports the horizontal geographic case used by the
//! first GeoTIFF-backed `hgridshift` corpus pipelines:
//! - `TYPE=HORIZONTAL_OFFSET`
//! - latitude/longitude offset bands described in arc-seconds
//! - regular grids georeferenced by ModelPixelScale + ModelTiepoint

use super::{BaseGrid, GridHeader, GridSource};
use crate::Error;
use std::io::Cursor;
use tiff::decoder::{Decoder, DecodingResult};
use tiff::tags::Tag;

const GDAL_METADATA_TAG: Tag = Tag::Unknown(42112);

pub fn gtg(name: &str, buf: &[u8]) -> Result<BaseGrid, Error> {
    let mut decoder = Decoder::new(Cursor::new(buf)).map_err(tiff_err)?;
    let mut images = Vec::new();
    let mut previous_meta = None;

    loop {
        let (image, meta) = decode_gtg_ifd(&mut decoder, name, previous_meta.as_ref())?;
        previous_meta = Some(meta);
        images.push(image);
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

    let meta = GtgMetadata::from_xml(&metadata, layout, fallback_meta)?;
    meta.validate_horizontal_arcsec()?;
    let grid_name = meta
        .grid_name
        .clone()
        .unwrap_or_else(|| default_name.to_string());
    let required_bands = meta.lat_band.max(meta.lon_band) + 1;
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

    let lat_n = tiepoint[4].to_radians();
    let lon_w = tiepoint[3].to_radians();
    let dlat = -pixel_scale[1].to_radians();
    let dlon = pixel_scale[0].to_radians();
    let lat_s = (tiepoint[4] - (rows.saturating_sub(1) as f64) * pixel_scale[1]).to_radians();
    let lon_e = (tiepoint[3] + (cols.saturating_sub(1) as f64) * pixel_scale[0]).to_radians();
    let header = GridHeader {
        lat_n,
        lat_s,
        lon_w,
        lon_e,
        dlat,
        dlon,
        rows,
        cols,
        bands: 2,
    };

    let mut grid = Vec::with_capacity(values.len());
    match layout {
        SampleLayout::Chunky => {
            for node in values.chunks_exact(decoded_bands) {
                let lon = apply_axis_polarity(node[meta.lon_band], meta.lon_positive_east);
                let lat = apply_axis_polarity(node[meta.lat_band], meta.lat_positive_north);
                grid.push(lon);
                grid.push(lat);
            }
        }
        SampleLayout::Separate => {
            let lon_plane = &values[meta.lon_band * plane_len..(meta.lon_band + 1) * plane_len];
            let lat_plane = &values[meta.lat_band * plane_len..(meta.lat_band + 1) * plane_len];
            for i in 0..plane_len {
                grid.push(apply_axis_polarity(lon_plane[i], meta.lon_positive_east));
                grid.push(apply_axis_polarity(lat_plane[i], meta.lat_positive_north));
            }
        }
    }

    Ok((
        BaseGrid::new(&grid_name, header, GridSource::Internal { values: grid })?,
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
}

fn apply_axis_polarity(value: f32, positive_is_forward: bool) -> f32 {
    if positive_is_forward { value } else { -value }
}

fn tiff_err(err: tiff::TiffError) -> Error {
    Error::Unsupported(format!("GTG TIFF decode failed: {err}"))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SampleLayout {
    Chunky,
    Separate,
}

#[derive(Debug, Clone)]
struct GtgMetadata {
    grid_name: Option<String>,
    lat_band: usize,
    lon_band: usize,
    lat_positive_north: bool,
    lon_positive_east: bool,
    lat_unit: String,
    lon_unit: String,
    grid_type: String,
}

impl GtgMetadata {
    fn from_xml(xml: &str, _layout: SampleLayout, fallback: Option<&Self>) -> Result<Self, Error> {
        let grid_type = find_item(xml, "TYPE", None, None)
            .or_else(|| fallback.map(|meta| meta.grid_type.clone()))
            .ok_or_else(|| Error::Unsupported("GTG metadata is missing TYPE".to_string()))?;
        let grid_name = find_item(xml, "grid_name", None, None);

        let lat_band = find_band(xml, "latitude_offset").or_else(|_| {
            fallback.map(|meta| meta.lat_band).ok_or_else(|| {
                Error::Unsupported(
                    "GTG metadata is missing DESCRIPTION for latitude_offset".to_string(),
                )
            })
        })?;
        let lon_band = find_band(xml, "longitude_offset").or_else(|_| {
            fallback.map(|meta| meta.lon_band).ok_or_else(|| {
                Error::Unsupported(
                    "GTG metadata is missing DESCRIPTION for longitude_offset".to_string(),
                )
            })
        })?;
        let lat_unit = find_item(xml, "UNITTYPE", Some(lat_band), Some("unittype"))
            .or_else(|| fallback.map(|meta| meta.lat_unit.clone()))
            .ok_or_else(|| {
                Error::Unsupported("GTG metadata is missing latitude UNITTYPE".to_string())
            })?;
        let lon_unit = find_item(xml, "UNITTYPE", Some(lon_band), Some("unittype"))
            .or_else(|| fallback.map(|meta| meta.lon_unit.clone()))
            .ok_or_else(|| {
                Error::Unsupported("GTG metadata is missing longitude UNITTYPE".to_string())
            })?;

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
            lat_positive_north,
            lon_positive_east,
            lat_unit,
            lon_unit,
            grid_type,
        })
    }

    fn validate_horizontal_arcsec(&self) -> Result<(), Error> {
        if !self.grid_type.eq_ignore_ascii_case("HORIZONTAL_OFFSET") {
            return Err(Error::Unsupported(format!(
                "GTG reader currently supports only TYPE=HORIZONTAL_OFFSET, found {}",
                self.grid_type
            )));
        }
        if !self.lat_unit.eq_ignore_ascii_case("arc-second")
            || !self.lon_unit.eq_ignore_ascii_case("arc-second")
        {
            return Err(Error::Unsupported(format!(
                "GTG reader currently supports only arc-second horizontal offsets, found lat={} lon={}",
                self.lat_unit, self.lon_unit
            )));
        }
        Ok(())
    }
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
            "/../../data/grids/us_noaa_conus.tif"
        );
        let grid = BaseGrid::read(path)?;
        assert_eq!(grid.header.bands, 2);
        assert!(grid.header.rows > 0);
        assert!(grid.header.cols > 0);
        Ok(())
    }
}
