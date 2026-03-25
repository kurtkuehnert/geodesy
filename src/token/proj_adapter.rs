use crate::Error;
use crate::ellipsoid::{EllipsoidBase, biaxial::Ellipsoid};
use crate::math::angular;

pub fn tidy_proj(elements: &mut Vec<String>) -> Result<(), Error> {
    if elements.first().map(String::as_str) == Some("hgridshift") {
        elements[0] = "gridshift".to_string();
    }

    // Geodesy only supports ellipsoid definitions as named builtins or ellps=a,rf
    // PROJ has richer support which we try navigate here
    // First we find the indices of ellps, a, rf and R elements
    let mut ellps_def: [Option<usize>; 5] = [None; 5];
    for (i, element) in elements.iter().enumerate() {
        if element.starts_with("ellps=") {
            ellps_def[0] = Some(i);
        }
        if element.starts_with("a=") {
            ellps_def[1] = Some(i);
        }
        if element.starts_with("b=") {
            ellps_def[2] = Some(i);
        }
        if element.starts_with("rf=") {
            ellps_def[3] = Some(i);
        }
        if element.starts_with("R=") {
            ellps_def[4] = Some(i);
        }
    }

    if let [None, Some(a_idx), _, Some(rf_idx), _] = ellps_def {
        let a = elements[a_idx][2..].to_string();
        let rf = elements[rf_idx][3..].to_string();
        elements.push(format!("ellps={a},{rf}").to_string());

        if a_idx > rf_idx {
            elements.remove(a_idx);
            elements.remove(rf_idx);
        } else {
            elements.remove(rf_idx);
            elements.remove(a_idx);
        }
    }

    if let [None, Some(a_idx), Some(b_idx), None, _] = ellps_def {
        let a = elements[a_idx][2..].trim().parse::<f64>();
        let b = elements[b_idx][2..].trim().parse::<f64>();
        if let (Ok(a), Ok(b)) = (a, b) {
            if a > 0.0 && b > 0.0 && b <= a {
                let rf = if (a - b).abs() < f64::EPSILON {
                    0.0
                } else {
                    a / (a - b)
                };
                elements.push(format!("ellps={a},{rf}"));
                if a_idx > b_idx {
                    elements.remove(a_idx);
                    elements.remove(b_idx);
                } else {
                    elements.remove(b_idx);
                    elements.remove(a_idx);
                }
            }
        }
    }

    if let [None, None, None, None, Some(r_idx)] = ellps_def {
        let radius = elements[r_idx][2..].to_string();
        elements.push(format!("ellps={radius},0"));
        elements.remove(r_idx);
    }

    normalize_authalic_sphere(elements)?;

    for (i, element) in elements.iter().enumerate() {
        if let Some(stripped) = element.strip_prefix("k=") {
            elements[i] = "k_0=".to_string() + stripped;
            break;
        }
    }

    normalize_prime_meridian(elements)?;
    normalize_omerc(elements);

    Ok(())
}

fn normalize_authalic_sphere(elements: &mut Vec<String>) -> Result<(), Error> {
    let Some(r_a_idx) = elements.iter().position(|e| e == "R_A") else {
        return Ok(());
    };

    let Some(ellps_idx) = elements.iter().position(|e| e.starts_with("ellps=")) else {
        return Err(Error::Unsupported(
            "parse_proj does not support +R_A without a resolvable ellipsoid".to_string(),
        ));
    };

    let definition = elements[ellps_idx][6..].trim();
    let ellps = Ellipsoid::named(definition).map_err(|_| {
        Error::Unsupported(format!(
            "parse_proj cannot derive authalic sphere from ellipsoid definition: {definition}"
        ))
    })?;
    let radius = authalic_radius(&ellps);
    elements[ellps_idx] = format!("ellps={radius},0");
    elements.remove(r_a_idx);
    Ok(())
}

fn authalic_radius(ellps: &impl EllipsoidBase) -> f64 {
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    if e == 0.0 {
        return a;
    }

    let one_minus_es = 1.0 - ellps.eccentricity_squared();
    let qp = one_minus_es
        * (1.0 / one_minus_es - (0.5 / e) * ((1.0 - e) / (1.0 + e)).ln());
    a * (0.5 * qp).sqrt()
}

pub fn extract_axis_step(elements: &mut Vec<String>) -> Result<Option<String>, Error> {
    let Some(axis_idx) = elements.iter().position(|e| e.starts_with("axis=")) else {
        return Ok(None);
    };

    // PROJ rejects specifying both axis= and order= on the same step
    if elements.iter().any(|e| e.starts_with("order=")) {
        return Err(Error::Unsupported(
            "parse_proj does not support both axis= and order= on the same step".to_string(),
        ));
    }

    let axis = elements[axis_idx][5..].to_ascii_lowercase();
    elements.remove(axis_idx);

    let chars: Vec<char> = axis.chars().collect();
    if chars.len() < 2 {
        return Err(Error::Unsupported(format!(
            "parse_proj does not support malformed axis specifier: {axis}"
        )));
    }

    let mapped: Vec<String> = chars
        .iter()
        .take(3)
        .map(|c| match c {
            'e' => Ok("1".to_string()),
            'w' => Ok("-1".to_string()),
            'n' => Ok("2".to_string()),
            's' => Ok("-2".to_string()),
            'u' => Ok("3".to_string()),
            'd' => Ok("-3".to_string()),
            other => Err(Error::Unsupported(format!(
                "parse_proj does not support axis orientation '{other}' in axis={axis}"
            ))),
        })
        .collect::<Result<Vec<_>, _>>()?;

    // When axis= appears on an axisswap step, inject order= directly so the step
    // is fully specified. For other operators, append a separate axisswap step.
    if elements.first().map(String::as_str) == Some("axisswap") {
        elements.push(format!("order={}", mapped.join(",")));
        return Ok(None);
    }

    if mapped == ["1", "2"] || mapped == ["1", "2", "3"] {
        return Ok(None);
    }

    Ok(Some(format!("axisswap order={}", mapped.join(","))))
}

fn normalize_prime_meridian(elements: &mut Vec<String>) -> Result<(), Error> {
    let Some(pm_idx) = elements.iter().position(|e| e.starts_with("pm=")) else {
        return Ok(());
    };
    let pm = elements[pm_idx][3..].to_string();
    let Some(offset_deg) = angular::parse_prime_meridian(&pm) else {
        return Err(Error::Unsupported(format!(
            "parse_proj does not support prime meridian specifier: {pm}"
        )));
    };

    let is_longlat = matches!(
        elements.first().map(String::as_str),
        Some("longlat" | "latlon" | "latlong" | "lonlat")
    );

    for key in ["lon_0=", "lonc=", "lon_1=", "lon_2="] {
        if is_longlat && key == "lon_0=" {
            continue;
        }
        if let Some(idx) = elements.iter().position(|e| e.starts_with(key)) {
            let value = elements[idx][key.len()..].to_string();
            let angle = angular::parse_sexagesimal(&value);
            if angle.is_nan() {
                return Err(Error::Unsupported(format!(
                    "parse_proj cannot adjust longitude parameter {key}{value} for pm={pm}"
                )));
            }
            elements[idx] = format!("{key}{}", angle + offset_deg);
        }
    }

    if is_longlat {
        if let Some(idx) = elements.iter().position(|e| e.starts_with("lon_0=")) {
            let value = elements[idx][6..].to_string();
            let angle = angular::parse_sexagesimal(&value);
            if angle.is_nan() {
                return Err(Error::Unsupported(format!(
                    "parse_proj cannot adjust longitude parameter lon_0={value} for pm={pm}"
                )));
            }
            elements[idx] = format!("lon_0={}", angle + offset_deg);
        } else {
            elements.push(format!("lon_0={offset_deg}"));
        }
    }

    elements.remove(pm_idx);
    Ok(())
}

fn normalize_omerc(elements: &mut Vec<String>) {
    if elements.first().map(String::as_str) != Some("omerc") {
        return;
    }

    rename_parameter(elements, "lat_0=", "latc=");
    rename_parameter(elements, "gamma=", "gamma_c=");

    let no_uoff =
        remove_parameter_flag(elements, "no_uoff") || remove_parameter_flag(elements, "no_off");
    let has_variant = elements.iter().any(|e| e == "variant");
    let has_central_point_form = elements.iter().any(|e| e.starts_with("lonc="));
    if has_central_point_form && !no_uoff && !has_variant {
        elements.push("variant".to_string());
    }
}

fn rename_parameter(elements: &mut [String], from: &str, to: &str) {
    if elements.iter().any(|e| e.starts_with(to)) {
        return;
    }
    if let Some(idx) = elements.iter().position(|e| e.starts_with(from)) {
        let value = elements[idx][from.len()..].to_string();
        elements[idx] = format!("{to}{value}");
    }
}

fn remove_parameter_flag(elements: &mut Vec<String>, flag: &str) -> bool {
    let before = elements.len();
    elements.retain(|e| e != flag);
    before != elements.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rewrites_hgridshift_to_gridshift() -> Result<(), Error> {
        let mut elements = vec![
            "hgridshift".to_string(),
            "grids=us_noaa_conus.tif".to_string(),
        ];
        tidy_proj(&mut elements)?;
        assert_eq!(elements, vec!["gridshift", "grids=us_noaa_conus.tif"]);
        Ok(())
    }
}
