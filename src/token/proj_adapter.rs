use crate::Error;
use crate::ellipsoid::{EllipsoidBase, biaxial::Ellipsoid};
use crate::math::angular;

pub fn tidy_proj(elements: &mut Vec<String>) -> Result<(), Error> {
    if elements.first().map(String::as_str) == Some("hgridshift") {
        elements[0] = "gridshift".to_string();
    }
    if elements.first().map(String::as_str) == Some("etmerc") {
        elements[0] = "tmerc".to_string();
    }
    normalize_geoc_aliases(elements)?;
    normalize_aeqd_variants(elements);
    normalize_sterec(elements);
    normalize_ups(elements)?;

    // Sphere reduction modifiers (R_A, R_V, R_a, …) must run first so that
    // the named ellipsoid is still in resolvable form when the radius is computed.
    normalize_sphere_reductions(elements)?;

    // Collapse all PROJ ellipsoid parameter variants into geodesy's canonical
    // `ellps=<semimajor_axis>,<reciprocal_flattening>` form.
    normalize_ellipsoid_params(elements)?;

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

fn normalize_aeqd_variants(elements: &mut Vec<String>) {
    if elements.first().map(String::as_str) != Some("aeqd") {
        return;
    }

    if remove_parameter_flag(elements, "guam") {
        elements[0] = "guam_aeqd".to_string();
    }
}

fn normalize_sterec(elements: &mut Vec<String>) {
    if elements.first().map(String::as_str) != Some("stere") {
        return;
    }

    if remove_parameter_flag(elements, "variant_c") {
        elements[0] = "sterec".to_string();
    }
}

fn normalize_geoc_aliases(elements: &mut Vec<String>) -> Result<(), Error> {
    let Some(op) = elements.first().cloned() else {
        return Ok(());
    };

    if op == "geocent" {
        elements[0] = "cart".to_string();
        if let Some(idx) = find_prefix(elements, "lon_0=") {
            let lon_0 = parse_f64(&elements[idx][6..], "lon_0")?;
            if lon_0.abs() > f64::EPSILON {
                return Err(Error::Unsupported(format!(
                    "parse_proj does not support geocent with lon_0={lon_0}"
                )));
            }
            elements.remove(idx);
        }
        return Ok(());
    }

    let had_geoc_flag = if op == "geoc" {
        false
    } else {
        remove_parameter_flag(elements, "geoc")
    };
    let is_longlat = matches!(op.as_str(), "longlat" | "latlon" | "latlong" | "lonlat");
    if op == "geoc" || (is_longlat && had_geoc_flag) {
        elements[0] = "latitude".to_string();
        if is_longlat {
            if let Some(idx) = elements.iter().position(|element| element == "inv") {
                elements.remove(idx);
            } else {
                elements.insert(1, "inv".to_string());
            }
        }
        if !elements.iter().any(|element| element == "geocentric") {
            elements.push("geocentric".to_string());
        }
    }

    Ok(())
}

fn normalize_ups(elements: &mut Vec<String>) -> Result<(), Error> {
    if elements.first().map(String::as_str) != Some("ups") {
        return Ok(());
    }

    let ellps_idx = find_prefix(elements, "ellps=");
    let a_idx = find_prefix(elements, "a=");
    let b_idx = find_prefix(elements, "b=");
    let rf_idx = find_prefix(elements, "rf=");
    let f_idx = find_prefix(elements, "f=");
    let r_idx = find_prefix(elements, "R=");
    let south = remove_parameter_flag(elements, "south");

    if r_idx.is_some()
        || (a_idx.is_some() && b_idx.is_none() && rf_idx.is_none() && f_idx.is_none())
    {
        return Err(Error::Unsupported(
            "ups requires an ellipsoidal definition; spherical UPS is not supported".into(),
        ));
    }

    elements[0] = "stere".to_string();
    elements.push(format!("lat_0={}", if south { -90 } else { 90 }));
    elements.push("lon_0=0".to_string());
    elements.push("k_0=0.994".to_string());
    elements.push("x_0=2000000".to_string());
    elements.push("y_0=2000000".to_string());

    // Keep existing ellipsoid params in place for downstream normalization.
    if ellps_idx.is_none()
        && a_idx.is_none()
        && b_idx.is_none()
        && rf_idx.is_none()
        && f_idx.is_none()
    {
        elements.push("ellps=GRS80".to_string());
    }
    Ok(())
}

/// Collapse all PROJ-style ellipsoid parameter combinations into geodesy's single
/// canonical form: `ellps=<semimajor_axis>,<reciprocal_flattening>`.
///
/// PROJ allows many equivalent ways to define an ellipsoid — named builtins
/// (`ellps=GRS80`), numeric combinations (`a` + `b`, `a` + `rf`, `a` + `f`, `R`
/// for a sphere), or named builtins with an overriding shape parameter. Geodesy
/// only understands named builtins (passed through unchanged) or `ellps=a,rf`.
///
/// Also validates shape parameters and rejects degenerate combinations (e.g.
/// `a ≤ 0`, `b ≤ 0`, `f ≥ 1`) that PROJ rejects at construction time.
fn normalize_ellipsoid_params(elements: &mut Vec<String>) -> Result<(), Error> {
    // Convert e= (eccentricity) or es= (eccentricity squared) to f= so the
    // main match below handles them uniformly.
    if let Some(ei) = find_prefix(elements, "e=") {
        let e = parse_f64(elements[ei][2..].trim(), "e")?;
        if !(0.0..1.0).contains(&e) {
            return Err(Error::Unsupported(format!(
                "ellipsoid: e={e} must be in [0, 1)"
            )));
        }
        elements[ei] = format!("f={}", 1.0 - (1.0 - e * e).sqrt());
    } else if let Some(esi) = find_prefix(elements, "es=") {
        let es = parse_f64(elements[esi][3..].trim(), "es")?;
        if !(0.0..1.0).contains(&es) {
            return Err(Error::Unsupported(format!(
                "ellipsoid: es={es} must be in [0, 1)"
            )));
        }
        elements[esi] = format!("f={}", 1.0 - (1.0 - es).sqrt());
    }

    let named = find_prefix(elements, "ellps=");
    let a = find_prefix(elements, "a=");
    let b = find_prefix(elements, "b=");
    let rf = find_prefix(elements, "rf=");
    let sphere = find_prefix(elements, "R=");
    let f = find_prefix(elements, "f=");

    match (named, a, b, rf, sphere, f) {
        // ── anonymous ellipsoid defined by two numeric shape params ──────────

        // a + rf  →  ellps=a,rf
        (None, Some(ai), _, Some(rfi), _, _) => {
            let a_val = parse_f64(&elements[ai][2..], "a")?;
            if a_val <= 0.0 {
                return Err(Error::Unsupported(format!(
                    "ellipsoid: a={a_val} must be positive"
                )));
            }
            let rf_val = elements[rfi][3..].trim().to_string();
            elements.push(format!("ellps={a_val},{rf_val}"));
            remove_two(elements, ai, rfi);
        }

        // a + b  →  ellps=a,rf  (rf derived; only when 0 < b ≤ a)
        (None, Some(ai), Some(bi), None, _, _) => {
            let a_val = parse_f64(&elements[ai][2..], "a")?;
            let b_val = parse_f64(&elements[bi][2..], "b")?;
            if a_val <= 0.0 || b_val <= 0.0 {
                return Err(Error::Unsupported(format!(
                    "ellipsoid: a={a_val} and b={b_val} must both be positive"
                )));
            }
            if b_val <= a_val {
                let rf_val = rf_from_ab(a_val, b_val);
                elements.push(format!("ellps={a_val},{rf_val}"));
                remove_two(elements, ai, bi);
            }
            // b > a: pass through unchanged — Ellipsoid::named will reject at construction time
        }

        // a + f  →  ellps=a,rf  (rf = 1/f, or 0 for a sphere when f = 0)
        (None, Some(ai), None, None, _, Some(fi)) => {
            let a_val = elements[ai][2..].to_string();
            let f_val = parse_f64(&elements[fi][2..], "f")?;
            let rf_val = rf_from_f(f_val);
            elements.push(format!("ellps={a_val},{rf_val}"));
            remove_two(elements, ai, fi);
        }

        // R  →  ellps=R,0  (explicit sphere radius)
        (None, None, None, None, Some(ri), _) => {
            let r_val = parse_f64(&elements[ri][2..], "R")?;
            if r_val <= 0.0 {
                return Err(Error::Unsupported(format!(
                    "ellipsoid: R={r_val} must be positive"
                )));
            }
            elements.push(format!("ellps={r_val},0"));
            elements.remove(ri);
        }

        // ── named ellipsoid with an overriding shape parameter ───────────────

        // ellps=NAME + b  →  ellps=a_from_name, rf derived from a and b
        (Some(ei), None, Some(bi), None, _, _) => {
            let a_val = resolve_named(elements, ei)?;
            let b_val = parse_f64(&elements[bi][2..], "b")?;
            let rf_val = rf_from_ab(a_val, b_val);
            elements.push(format!("ellps={a_val},{rf_val}"));
            remove_two(elements, ei, bi);
        }

        // ellps=NAME + rf  →  ellps=a_from_name, rf
        // Note: rf=0 is rejected here even though `rf=0` is a valid sphere sentinel
        // in EPSG/PROJ database catalogs — as an *explicit user parameter* PROJ
        // treats it as invalid_op_illegal_arg_value.
        (Some(ei), None, None, Some(rfi), _, _) => {
            let a_val = resolve_named(elements, ei)?;
            let rf_str = elements[rfi][3..].trim();
            let rf_num = parse_f64(rf_str, "rf")?;
            if rf_num == 0.0 {
                return Err(Error::Unsupported(
                    "ellipsoid: rf=0 is not valid as an explicit parameter override".into(),
                ));
            }
            elements.push(format!("ellps={a_val},{rf_str}"));
            remove_two(elements, ei, rfi);
        }

        // ellps=NAME + f  →  ellps=a_from_name, rf = 1/f
        (Some(ei), None, None, None, _, Some(fi)) => {
            let a_val = resolve_named(elements, ei)?;
            let f_val = parse_f64(&elements[fi][2..], "f")?;
            let rf_val = rf_from_f(f_val);
            elements.push(format!("ellps={a_val},{rf_val}"));
            remove_two(elements, ei, fi);
        }

        // Everything else: no normalisation needed (named-only, already a,rf, etc.)
        _ => {}
    }

    // Validate any shape params that were not consumed by the match above.
    // These represent lone overrides or invalid combinations that geodesy would
    // silently ignore; validate here to match PROJ's rejection behaviour.
    for (prefix, label) in [("a=", "a"), ("b=", "b")] {
        if let Some(idx) = find_prefix(elements, prefix) {
            let val = parse_f64(&elements[idx][prefix.len()..], label)?;
            if val <= 0.0 {
                return Err(Error::Unsupported(format!(
                    "ellipsoid: {label}={val} must be positive"
                )));
            }
        }
    }
    if let Some(fi) = find_prefix(elements, "f=") {
        let f_val = parse_f64(&elements[fi][2..], "f")?;
        if !(0.0..1.0).contains(&f_val) {
            return Err(Error::Unsupported(format!(
                "ellipsoid: f={f_val} must be in [0, 1)"
            )));
        }
    }

    Ok(())
}

// ── ellipsoid helpers ────────────────────────────────────────────────────────

/// Return the index of the first element that starts with `prefix`, if any.
fn find_prefix(elements: &[String], prefix: &str) -> Option<usize> {
    elements.iter().position(|e| e.starts_with(prefix))
}

/// Remove two elements at distinct indices `i` and `j` (order-independent).
fn remove_two(elements: &mut Vec<String>, i: usize, j: usize) {
    let (hi, lo) = if i > j { (i, j) } else { (j, i) };
    elements.remove(hi);
    elements.remove(lo);
}

/// Resolve a named ellipsoid (`ellps=NAME`) at index `idx` and return its semi-major axis.
fn resolve_named(elements: &[String], idx: usize) -> Result<f64, Error> {
    let name = elements[idx][6..].trim();
    Ellipsoid::named(name)
        .map(|e| e.semimajor_axis())
        .map_err(|_| Error::Unsupported(format!("parse_proj cannot resolve ellipsoid: {name}")))
}

/// Reciprocal flattening from semi-major and semi-minor axes (`0` for a sphere).
fn rf_from_ab(a: f64, b: f64) -> f64 {
    if (a - b).abs() < f64::EPSILON {
        0.0
    } else {
        a / (a - b)
    }
}

/// Reciprocal flattening from flattening `f` (`0` for a sphere when `f = 0`).
fn rf_from_f(f: f64) -> f64 {
    if f.abs() < f64::EPSILON { 0.0 } else { 1.0 / f }
}

/// Parse a float from a string slice, returning a descriptive error on failure.
fn parse_f64(s: &str, param: &str) -> Result<f64, Error> {
    s.trim()
        .parse()
        .map_err(|_| Error::Unsupported(format!("parse_proj cannot parse {param}= value: {s}")))
}

// ── sphere reduction ─────────────────────────────────────────────────────────

/// Replace a PROJ sphere-reduction modifier with `ellps=<radius>,0`.
///
/// PROJ supports several ways to reduce a reference ellipsoid to an equivalent
/// sphere for operators that only support spherical geometry. The supported
/// modifiers are:
///
/// | Modifier       | Sphere radius                              |
/// |----------------|--------------------------------------------|
/// | `R_A`          | Authalic (equal-area) radius               |
/// | `R_V`          | Volumetric radius `(a²b)^(1/3)`            |
/// | `R_C`          | Conformal sphere radius at `lat_0`         |
/// | `R_a`          | Arithmetic mean `(a + b) / 2`              |
/// | `R_g`          | Geometric mean `√(a·b)`                    |
/// | `R_h`          | Harmonic mean `2ab / (a + b)`              |
/// | `R_lat_a=φ`    | Arithmetic mean of N and M at latitude φ   |
/// | `R_lat_g=φ`    | Geometric mean of N and M at latitude φ    |
///
/// Must be called before [`normalize_ellipsoid_params`] so that any named
/// ellipsoid is still in resolvable (`ellps=NAME`) form.
fn normalize_sphere_reductions(elements: &mut Vec<String>) -> Result<(), Error> {
    let mod_idx = elements.iter().position(|e| {
        matches!(e.as_str(), "R_A" | "R_V" | "R_C" | "R_a" | "R_g" | "R_h")
            || e.starts_with("R_lat_a=")
            || e.starts_with("R_lat_g=")
    });

    let Some(mod_idx) = mod_idx else {
        return Ok(());
    };

    let modifier = elements[mod_idx].clone();

    let ellps = resolve_ellipsoid_for_sphere_reduction(elements, &modifier)?;
    let ellps_idx = find_prefix(elements, "ellps=");
    let a_idx = find_prefix(elements, "a=");
    let b_idx = find_prefix(elements, "b=");
    let rf_idx = find_prefix(elements, "rf=");
    let sphere_idx = find_prefix(elements, "R=");
    let f_idx = find_prefix(elements, "f=");

    let a = ellps.semimajor_axis();
    let b = ellps.semiminor_axis();
    let e2 = ellps.eccentricity_squared();

    let radius = match modifier.as_str() {
        "R_A" => authalic_radius(&ellps),
        "R_V" => (a * a * b).cbrt(),
        "R_C" => {
            // PROJ's current Mercator runtime behavior effectively uses the default
            // latitude of origin here, even when a non-zero lat_0 is supplied.
            let phi0 = if elements.first().map(String::as_str) == Some("merc") {
                0.0
            } else {
                find_prefix(elements, "lat_0=")
                    .map(|idx| parse_f64(&elements[idx][6..], "lat_0"))
                    .transpose()?
                    .unwrap_or(0.0)
                    .to_radians()
            };
            let sin_phi0 = phi0.sin();
            let denom = 1.0 - e2 * sin_phi0 * sin_phi0;
            a * (1.0 - e2).sqrt() / denom
        }
        "R_a" => (a + b) / 2.0,
        "R_g" => (a * b).sqrt(),
        "R_h" => 2.0 * a * b / (a + b),
        _ => {
            let (prefix, is_arith) = if let Some(s) = modifier.strip_prefix("R_lat_a=") {
                (s, true)
            } else if let Some(s) = modifier.strip_prefix("R_lat_g=") {
                (s, false)
            } else {
                return Err(Error::Unsupported(format!(
                    "parse_proj does not support sphere reduction modifier: {modifier}"
                )));
            };
            let phi = parse_f64(prefix, "R_lat latitude")?.to_radians();
            let sin_phi = phi.sin();
            let denom = 1.0 - e2 * sin_phi * sin_phi;
            let n = a / denom.sqrt();
            let m = a * (1.0 - e2) / (denom * denom.sqrt());
            if is_arith {
                (n + m) / 2.0
            } else {
                (n * m).sqrt()
            }
        }
    };

    if let Some(idx) = ellps_idx {
        elements[idx] = format!("ellps={radius},0");
    } else {
        elements.push(format!("ellps={radius},0"));
    }
    let mut indices_to_remove = vec![mod_idx];
    indices_to_remove.extend(
        [a_idx, b_idx, rf_idx, sphere_idx, f_idx]
            .into_iter()
            .flatten(),
    );
    indices_to_remove.sort_unstable();
    indices_to_remove.dedup();
    for idx in indices_to_remove.into_iter().rev() {
        elements.remove(idx);
    }
    Ok(())
}

fn resolve_ellipsoid_for_sphere_reduction(
    elements: &[String],
    modifier: &str,
) -> Result<Ellipsoid, Error> {
    let a_idx = find_prefix(elements, "a=");
    let b_idx = find_prefix(elements, "b=");
    let rf_idx = find_prefix(elements, "rf=");
    let sphere_idx = find_prefix(elements, "R=");
    let f_idx = find_prefix(elements, "f=");
    let named = if let Some(ellps_idx) = find_prefix(elements, "ellps=") {
        let definition = elements[ellps_idx][6..].trim().to_string();
        Some(Ellipsoid::named(&definition).map_err(|_| {
            Error::Unsupported(format!(
                "parse_proj cannot derive sphere radius from ellipsoid definition: {definition}"
            ))
        })?)
    } else {
        None
    };

    match (named, a_idx, b_idx, rf_idx, sphere_idx, f_idx) {
        (Some(base), maybe_a, Some(bi), None, None, None) => {
            let a = if let Some(ai) = maybe_a {
                parse_f64(&elements[ai][2..], "a")?
            } else {
                base.semimajor_axis()
            };
            let b = parse_f64(&elements[bi][2..], "b")?;
            Ellipsoid::named(&format!("{a},{}", rf_from_ab(a, b)))
        }
        (Some(base), maybe_a, None, Some(rfi), None, None) => {
            let a = if let Some(ai) = maybe_a {
                parse_f64(&elements[ai][2..], "a")?
            } else {
                base.semimajor_axis()
            };
            let rf = parse_f64(&elements[rfi][3..], "rf")?;
            Ellipsoid::named(&format!("{a},{rf}"))
        }
        (Some(base), maybe_a, None, None, None, Some(fi)) => {
            let a = if let Some(ai) = maybe_a {
                parse_f64(&elements[ai][2..], "a")?
            } else {
                base.semimajor_axis()
            };
            let f = parse_f64(&elements[fi][2..], "f")?;
            Ellipsoid::named(&format!("{a},{}", rf_from_f(f)))
        }
        (Some(base), Some(ai), None, None, None, None) => {
            let a = parse_f64(&elements[ai][2..], "a")?;
            Ellipsoid::named(&format!(
                "{a},{}",
                if base.flattening() == 0.0 {
                    0.0
                } else {
                    1.0 / base.flattening()
                }
            ))
        }
        (Some(base), None, None, None, None, None) => Ok(base),
        (None, Some(ai), Some(bi), None, None, None) => {
            let a = parse_f64(&elements[ai][2..], "a")?;
            let b = parse_f64(&elements[bi][2..], "b")?;
            Ellipsoid::named(&format!("{a},{}", rf_from_ab(a, b)))
        }
        (None, Some(ai), None, Some(rfi), None, None) => {
            let a = parse_f64(&elements[ai][2..], "a")?;
            let rf = parse_f64(&elements[rfi][3..], "rf")?;
            Ellipsoid::named(&format!("{a},{rf}"))
        }
        (None, None, None, None, Some(ri), None) => {
            let r = parse_f64(&elements[ri][2..], "R")?;
            if r <= 0.0 {
                return Err(Error::Unsupported(format!(
                    "ellipsoid: R={r} must be positive"
                )));
            }
            Ellipsoid::named(&format!("{r},0"))
        }
        (None, Some(ai), None, None, None, Some(fi)) => {
            let a = parse_f64(&elements[ai][2..], "a")?;
            let f = parse_f64(&elements[fi][2..], "f")?;
            Ellipsoid::named(&format!("{a},{}", rf_from_f(f)))
        }
        _ => Err(Error::Unsupported(format!(
            "parse_proj does not support +{modifier} without a resolvable ellipsoid"
        ))),
    }
}

fn authalic_radius(ellps: &impl EllipsoidBase) -> f64 {
    let a = ellps.semimajor_axis();
    let e = ellps.eccentricity();
    if e == 0.0 {
        return a;
    }
    let one_minus_es = 1.0 - ellps.eccentricity_squared();
    let qp = one_minus_es * (1.0 / one_minus_es - (0.5 / e) * ((1.0 - e) / (1.0 + e)).ln());
    a * (0.5 * qp).sqrt()
}

// ── axis / parameter normalisation ──────────────────────────────────────────

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
    let is_inverse = elements.iter().any(|element| element == "inv");
    let is_inverse_bonne =
        is_inverse && matches!(elements.first().map(String::as_str), Some("bonne"));

    if is_inverse_bonne {
        elements.remove(pm_idx);
        return Ok(());
    }

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
    if no_uoff {
        elements.push("no_off".to_string());
    }
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
