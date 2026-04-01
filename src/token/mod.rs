use crate::Error;
use std::collections::BTreeMap;

mod proj_adapter;

/// One grid candidate from a `grids=` list.
///
/// Grid candidates are ordered left-to-right within a step, and may be marked
/// optional by the `@grid` syntax.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GridSpec {
    /// The grid identifier as referenced from the pipeline definition.
    pub name: String,
    /// Whether the grid is optional at operator construction time.
    pub optional: bool,
}

/// Grid requirements for one pipeline step.
///
/// The returned `grids` preserve the original left-to-right priority order from
/// the `grids=` parameter. `null_fallback` records whether `@null` was present
/// in the same list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GridDependency {
    /// Ordered grid candidates for this step.
    pub grids: Vec<GridSpec>,
    /// Whether the step's `grids=` list contains `@null`.
    pub null_fallback: bool,
}

/// Convenience methods for lexical analysis of operator definitions.
/// - For splitting a pipeline into steps
/// - For splitting a step into parameters (i.e. key=value-pairs)
/// - For syntactical normalization by desugaring and elimination of non-significant whitespace
/// - For checking whether a given operator is singular or a pipeline
/// - For checking whether a key is a macro name ("resource name"), and
/// - For accessing the name of a given operator.
pub trait Tokenize {
    /// Remove inline and block comments
    fn remove_comments(&self) -> String;
    /// Rotate prefix modifiers (inv, omit_fwd, omit_inv) to the end of parameter list
    fn handle_prefix_modifiers(&self) -> String;

    /// Split a pipeline definition into steps
    fn split_into_steps(&self) -> Vec<String>;

    /// Split a step/an operation into parameters. Give special treatment
    /// to names and flags:
    /// ```txt
    /// 'foo bar=baz bonk=blue flag' -> ('name=foo', 'bar=baz', 'bonk=blue', 'flag=true')
    /// ```
    fn split_into_parameters(&self) -> BTreeMap<String, String>;

    /// Helper function for 'split_into_steps' and 'split_into_parameters':
    /// Glue syntactical elements together, and separate from each other
    /// by a single space:
    ///
    /// 1. Glue key-value pairs together by omitting whitespace around '=':
    ///    ```txt
    ///    key1= value1            key2    =value2  ->  key1=value1 key2=value2
    ///    ```
    /// 2. Trim whitespace on both sides of the macro sigil ':' and leave only
    ///    one space to the left of the dereference sigil '$':
    ///    ```txt
    ///    foo: bar $ baz -> foo:bar $baz
    ///    ```
    /// 3. Trim whitespace around sequence separators ',' and '|':
    ///    ```txt
    ///     foo | bar baz=bonk   ,    bonk  ->  foo|bar baz=bonk,bonk
    ///    ```
    /// 4. Desugar the one-way sequence separators '<' and '>':
    ///    ```txt
    ///     foo > bar < baz  ->  foo|omit_inv bar|omit_fwd baz
    ///    ```
    fn normalize(&self) -> String;

    fn is_pipeline(&self) -> bool;
    fn is_resource_name(&self) -> bool;
    fn operator_name(&self) -> String;
}

/// Tokenize implementation for string-like objects
impl<T> Tokenize for T
where
    T: AsRef<str>,
{
    fn remove_comments(&self) -> String {
        // Impose some line ending sanity
        let all = self
            .as_ref()
            .trim()
            .replace("\r\n", "\n") // The fenestration company
            .replace('\r', "\n") // The fruit company
            .replace("\n:", "\n") // Line continuation markers
            .to_string();
        // Remove comments
        let mut trimmed = String::new();
        for line in all.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            // Remove comments - both inline and separate lines ("block comments")
            let line: Vec<&str> = line.trim().split('#').collect();
            // Full line comment - just skip
            if line[0].starts_with('#') {
                continue;
            }

            // Inline comment, or no comment at all: Collect everything before `#`
            trimmed += " ";
            trimmed += line[0].trim();
        }
        trimmed
    }

    fn handle_prefix_modifiers(&self) -> String {
        if self.is_pipeline() {
            return self.as_ref().to_string();
        }

        let mut elements: Vec<_> = self.as_ref().split_whitespace().collect();
        if elements.is_empty() {
            return "".to_string();
        }

        // Rotate all desugared modifiers to the end of the list
        let modifiers = ["inv", "omit_fwd", "omit_inv"];
        while modifiers.contains(&elements[0]) {
            elements.rotate_left(1);
        }
        elements.join(" ")
    }

    fn split_into_steps(&self) -> Vec<String> {
        self.remove_comments()
            .normalize()
            // split into steps
            .split('|')
            // remove empty steps
            .filter(|x| !x.trim().is_empty())
            // rotate modifiers and convert &str to String
            .map(|x| x.handle_prefix_modifiers())
            // and turn into Vec<String>
            .collect()
    }

    fn split_into_parameters(&self) -> BTreeMap<String, String> {
        // Remove non-significant whitespace, then split by significant
        let step = self.as_ref().normalize().handle_prefix_modifiers();
        let elements: Vec<_> = step.split_whitespace().collect();
        let mut params = BTreeMap::new();
        if step.is_pipeline() {
            params.insert(String::from("_name"), step);
            return params;
        }

        for element in elements {
            // Split a key=value-pair into key and value parts
            let mut parts: Vec<&str> = element.trim().split('=').collect();
            // Add a boolean true part, to make sure we have a value, even for flags
            // (flags are booleans that are true when specified, false when not)
            parts.push("true");
            assert!(parts.len() > 1);

            // If the first arg is a key-without-value, it is the name of the operator
            if params.is_empty() && parts.len() == 2 {
                params.insert(String::from("_name"), String::from(parts[0]));
                continue;
            }

            params.insert(String::from(parts[0]), String::from(parts[1]));
        }

        params
    }

    fn normalize(&self) -> String {
        // Tweak everything into canonical form
        self.as_ref()
            .trim()
            .trim_matches(':')
            .replace("\n:", "\n")
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .replace("= ", "=")
            .replace(": ", ":")
            .replace(", ", ",")
            .replace("| ", "|")
            .replace("> ", ">")
            .replace("< ", "<")
            .replace(" =", "=")
            .replace(" :", ":")
            .replace(" ,", ",")
            .replace(" |", "|")
            .replace(" >", ">")
            .replace(" <", "<")
            .replace('>', "|omit_inv ")
            .replace('<', "|omit_fwd ")
            .replace("₀=", "_0=")
            .replace("₁=", "_1=")
            .replace("₂=", "_2=")
            .replace("₃=", "_3=")
            .replace("₄=", "_4=")
            .replace("₅=", "_5=")
            .replace("₆=", "_6=")
            .replace("₇=", "_7=")
            .replace("₈=", "_8=")
            .replace("₉=", "_9=")
            .replace("$ ", "$") // But keep " $" as is!
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
    }

    fn is_pipeline(&self) -> bool {
        self.as_ref().contains('|') || self.as_ref().contains('<') || self.as_ref().contains('>')
    }

    fn is_resource_name(&self) -> bool {
        self.operator_name().contains(':')
    }

    fn operator_name(&self) -> String {
        if self.is_pipeline() {
            return "".to_string();
        }
        self.split_into_parameters()
            .get("_name")
            .unwrap_or(&"".to_string())
            .to_string()
    }
}

/// Translate a PROJ string into Rust Geodesy format. Since PROJ is syntactically
/// unrestrictive, we do not try to detect any syntax errors: If the input
/// is so cursed as to be intranslatable, this will become clear when trying to
/// instantiate the result as a Geodesy operator. We do, however, check for and
/// report on two *semantically* refusable cases: First, that PROJ does not support
/// nested pipelines (the nesting must be done indirectly through an init-file),
/// second that Rust Geodesy does not support init-files. Hence no support for
/// any kind of nesting here.
///
/// ## Coordinate convention: PROJ vs Geodesy
///
/// **This is the most important difference between PROJ and Geodesy pipelines.**
///
/// PROJ uses **radians internally** but exposes a **degrees boundary** for
/// geographic operators at the API level:
/// - Single-step operations like `+proj=merc` implicitly accept lon/lat in degrees.
/// - Explicit `+proj=pipeline` strings carry their own
///   `+step +proj=unitconvert +xy_in=deg +xy_out=rad` steps at entry/exit to
///   handle the conversion explicitly — degrees at the boundary, radians inside.
///
/// Geodesy uses **radians throughout** with no implicit degree conversion at any
/// boundary. Callers are responsible for the conversion.
///
/// *parse_proj* bridges this gap for single-step PROJ strings by injecting an
/// explicit `unitconvert xy_in=deg xy_out=rad` step before geographic projection
/// operators, and a matching `unitconvert xy_in=rad xy_out=deg` step after
/// geographic identity operators (`lonlat`). This makes a geodesy pipeline
/// produced from a PROJ string behave identically to PROJ at the degree boundary.
///
/// Explicit `+proj=pipeline` strings are **not** wrapped — they already carry
/// their own unitconvert steps.
///
/// **Note on named angular parameters** (`lat_0=`, `lon_0=`, `latc=`, etc.):
/// These are always specified in degrees in both PROJ and Geodesy. Each operator
/// converts them to radians internally. The `unitconvert` wrapping described above
/// only affects coordinate *values* flowing through the pipeline at runtime, not
/// named parameters.
///
/// **Maintenance note:** The list of operators that require degree wrapping is
/// currently a hardcoded allowlist (`GEO_PROJ_OPS`, `GEO_IDENTITY_OPS`) inside
/// this function. It must be kept in sync as new geographic operators are added to
/// Geodesy. A future improvement would be to encode this as per-operator metadata
/// (analogous to how PROJ tracks input/output coordinate types in its operator
/// registry) so the wrapping is driven by the operator itself rather than this
/// centralized list.
///
/// ## Known differences between PROJ and Rust Geodesy definitions:
///
/// ## Ellipsoid definitions
/// - Geodesy only supports a limited set of builtin ellipsoids OR or definition
///   via semi-major and reverse-flattening parameters  `ellps=a,rf`.
/// - PROJ has [richer ellipsoid](https://proj.org/en/9.3/usage/ellipsoids.html#ellipsoid-size-parameters)
///   support which *parse_proj* provides partial support for.
/// - Specifically if an ellipsoid is defined via `a` and `rf` parameters, *parse_proj*
///   will redefine them as `ellps=a,rf` and remove the `a` and `rf` parameters.
/// - A spherical `R` parameter is also rewritten as `ellps=R,0`.
/// - All other cases supported by PROJ are NOT handled by *parse_proj* and will
///   fail when instantiating the operator.
///
/// ## Scaling via `k` parameter
/// - PROJ still supports the deprecated `k` parameter. Most output from `projinfo` will
///   have the scaling defined as `k` instead of `k_0`.
/// - *parse_proj* will replace `k` with `k_0` whenever it is encountered.
///
pub fn parse_proj(definition: &str, proj_degrees: bool) -> Result<String, Error> {
    // If it doesn't look like a PROJ string, we return it unchanged
    if definition.contains('|') | !definition.contains("proj") {
        return Ok(definition.to_string());
    }
    // Impose some line ending sanity and remove the PROJ '+' prefix
    let all = definition
        .replace("\r\n", "\n")
        .replace('\r', "\n")
        .replace(" +", " ")
        .replace("\n+", " ")
        .trim()
        .trim_start_matches('+')
        .to_string();

    // Collect the PROJ string
    let mut trimmed = String::new();
    for line in all.lines() {
        let line = line.trim();

        // Remove block comments
        let line: Vec<&str> = line.trim().split('#').collect();
        // Full line (block) comment - just skip
        if line[0].starts_with('#') {
            continue;
        }

        // Inline comment, or no comment at all: Collect everything before `#`
        trimmed += " ";
        trimmed += line[0].trim();
    }

    // Now split the text into steps. First make sure we do not match
    //"step" as part of a word (stairSTEPping,  poSTEPileptic, STEPwise,
    // quickSTEP), by making it possible to only search for " step "
    trimmed = " ".to_string() + &trimmed.normalize() + " ";

    // Remove empty steps and other non-significant whitespace
    let steps: Vec<String> = trimmed
        // split into steps
        .split(" step ")
        // remove empty steps
        .filter(|x| !x.trim().trim_start_matches("step ").is_empty())
        // remove spurious 'step step' noise and convert &str to String
        .map(|x| x.trim().trim_start_matches("step ").to_string())
        // turn into Vec<String>
        .collect();

    // For accumulating the pipeline steps converted to geodesy syntax
    let mut geodesy_steps = Vec::new();

    // Geodesy does not support pipeline globals, so we must explicitly
    // insert them in the beginning of the argument list of each step
    let mut pipeline_globals = "".to_string();
    let mut pipeline_is_inverted = false;
    // Track whether the input was an explicit +proj=pipeline (as opposed to a single-step
    // PROJ string or the unofficial multi-step-without-pipeline extension).
    let mut is_explicit_pipeline = false;
    // The operator name of the first (and typically only) non-pipeline step.
    let mut primary_op: Option<String> = None;

    for (step_index, step) in steps.iter().enumerate() {
        let mut elements: Vec<_> = step.split_whitespace().map(|x| x.to_string()).collect();

        // Move the "proj=..." element to the front of the collection, stripped for "proj="
        // and handle the pipeline globals, if any
        for (i, element) in elements.iter().enumerate() {
            // Mutating the Vec we are iterating over may seem dangerous but is
            // OK as we break out of the loop immediately after the mutation
            if element.starts_with("init=") {
                return Err(Error::Unsupported(
                    "parse_proj does not support PROJ init clauses: ".to_string() + step,
                ));
            }

            if element.starts_with("proj=") {
                elements.swap(i, 0);
                elements[0] = elements[0][5..].to_string();

                // In the proj=pipeline case, just collect the globals, without
                // introducing a new step into geodesy_steps
                if elements[0] == "pipeline" {
                    is_explicit_pipeline = true;
                    if step_index != 0 {
                        return Err(Error::Unsupported(
                            "PROJ does not support nested pipelines: ".to_string() + &trimmed,
                        ));
                    }
                    elements.remove(0);

                    // The case of 'inv' in globals must be handled separately, since it indicates
                    // the inversion of the entire pipeline, not just an inversion of each step
                    if elements.contains(&"inv".to_string()) {
                        pipeline_is_inverted = true;
                    }

                    // Remove all cases of 'inv' from the global arguments
                    let pipeline_globals_elements: Vec<String> = elements
                        .join(" ")
                        .trim()
                        .to_string()
                        .split_whitespace()
                        .filter(|x| x.trim() != "inv")
                        .map(|x| x.trim().to_string())
                        .collect();
                    pipeline_globals = pipeline_globals_elements.join(" ").trim().to_string();
                    elements.clear();
                }
                break;
            }
        }

        proj_adapter::tidy_proj(&mut elements)?;
        if primary_op.is_none() && !is_explicit_pipeline && !elements.is_empty() {
            primary_op = Some(elements[0].clone());
        }

        // Skip empty steps, insert pipeline globals, handle step and pipeline
        // inversions, and handle directional omissions (omit_fwd, omit_inv)
        let mut geodesy_step = elements.join(" ").trim().to_string();
        if !geodesy_step.is_empty() {
            if !pipeline_globals.is_empty() {
                elements.insert(1, pipeline_globals.clone());
            }

            let step_is_inverted = elements.contains(&"inv".to_string());
            elements = elements
                .iter()
                .filter(|x| x.as_str() != "inv")
                .map(|x| match x.as_str() {
                    "omit_fwd" => "omit_inv",
                    "omit_inv" => "omit_fwd",
                    _ => x,
                })
                .map(|x| x.to_string())
                .collect();

            let effective_inverted = step_is_inverted != pipeline_is_inverted;
            if effective_inverted {
                elements.insert(1, "inv".to_string());
            }

            let axis_step = proj_adapter::extract_axis_step(&mut elements)?;
            geodesy_step = elements.join(" ").trim().to_string();
            if let Some(axis_step) = axis_step {
                geodesy_step = if effective_inverted {
                    format!("{axis_step} | {geodesy_step}")
                } else {
                    format!("{geodesy_step} | {axis_step}")
                };
            }
            if pipeline_is_inverted {
                geodesy_steps.insert(0, geodesy_step);
            } else {
                geodesy_steps.push(geodesy_step);
            }
        }
    }
    let mut result = geodesy_steps.join(" | ").trim().to_string();

    // For a single-step non-pipeline PROJ string that represents a geographic operator,
    // inject degree↔radian conversion so the resulting geodesy pipeline matches PROJ's
    // external coordinate convention (geographic input/output in degrees).
    //
    // Explicit +proj=pipeline strings are excluded: their steps already operate in
    // PROJ's internal radian convention and include explicit unitconvert steps where needed.
    // Multi-step strings (unofficial extension) are excluded by the steps.len() == 1 guard.
    if proj_degrees && steps.len() == 1 && !is_explicit_pipeline && !result.is_empty() {
        let linear_scale = extract_linear_output_scale(&mut result)?;
        if let Some(op) = primary_op.as_deref() {
            let (input, output) = crate::inner_op::op_domains(op);
            let wrap_input = input == Some(crate::inner_op::CoordDomain::Geographic);
            let wrap_output = output == Some(crate::inner_op::CoordDomain::Geographic);
            let mut wrapped = result.clone();
            if wrap_input {
                wrapped = format!("unitconvert xy_in=deg xy_out=rad | {wrapped}");
            }
            if output == Some(crate::inner_op::CoordDomain::Projected)
                && let Some(scale) = linear_scale
            {
                wrapped = format!("{wrapped} | unitconvert xy_in=m xy_out={scale}");
            }
            return Ok(if wrap_output {
                format!("{wrapped} | unitconvert xy_in=rad xy_out=deg")
            } else {
                wrapped
            });
        }
    }

    Ok(result)
}

fn extract_linear_output_scale(definition: &mut String) -> Result<Option<f64>, Error> {
    let mut scale = None;
    let mut kept = Vec::new();
    let tokens: Vec<_> = definition.split_whitespace().collect();
    for token in tokens {
        if let Some(value) = token.strip_prefix("to_meter=") {
            scale = Some(value.parse::<f64>().map_err(|_| {
                Error::Unsupported(format!("parse_proj cannot parse to_meter scale: {value}"))
            })?);
            continue;
        }
        if let Some(value) = token.strip_prefix("units=") {
            scale = Some(match value {
                "km" => 1000.0,
                "m" => 1.0,
                "dm" => 0.1,
                "cm" => 0.01,
                "mm" => 0.001,
                "kmi" => 1852.0,
                "in" => 0.0254,
                "ft" => 0.3048,
                "yd" => 0.9144,
                "mi" => 1609.344,
                "fath" => 1.8288,
                "ch" => 20.1168,
                "link" => 0.201168,
                "us-in" => 100.0 / 3937.0,
                "us-ft" => 1200.0 / 3937.0,
                "us-yd" => 3600.0 / 3937.0,
                "us-ch" => 79200.0 / 3937.0,
                "us-mi" => 6336000.0 / 3937.0,
                "ind-yd" => 0.91439523,
                "ind-ft" => 0.30479841,
                "ind-ch" => 20.11669506,
                _ => {
                    kept.push(token.to_string());
                    continue;
                }
            });
            continue;
        }
        kept.push(token.to_string());
    }
    *definition = kept.join(" ");
    Ok(scale.filter(|value| (*value - 1.0).abs() > f64::EPSILON))
}

/// Inspect a pipeline definition and return one grid dependency chain per step.
///
/// The returned vector preserves pipeline step order. Steps without a `grids=`
/// parameter are omitted. Within each [`GridDependency`], `grids` preserve the
/// original left-to-right priority order.
///
/// This helper accepts both Rust Geodesy syntax and PROJ syntax.
pub fn grid_dependencies(definition: &str) -> Result<Vec<GridDependency>, Error> {
    let definition = parse_proj(definition, false)?;
    let mut dependencies = Vec::new();

    for step in definition.split_into_steps() {
        let params = step.split_into_parameters();
        let Some(grids) = params.get("grids") else {
            continue;
        };

        let mut chain = Vec::new();
        let mut null_fallback = false;
        for item in grids.split(',') {
            let item = item.trim();
            if item.is_empty() {
                continue;
            }

            let optional = item.starts_with('@');
            let name = item.trim_start_matches('@');
            if name.eq_ignore_ascii_case("null") {
                null_fallback = true;
                continue;
            }

            chain.push(GridSpec {
                name: name.to_string(),
                optional,
            });
        }

        dependencies.push(GridDependency {
            grids: chain,
            null_fallback,
        });
    }

    Ok(dependencies)
}

// ----- T E S T S ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    // Test the fundamental tokenization functionality
    #[test]
    fn token() -> Result<(), Error> {
        // Whitespace normalization
        assert_eq!("foo bar $ baz = bonk".normalize(), "foo bar $baz=bonk");
        assert_eq!(
            "foo |  bar baz  =  bonk, bonk , bonk x₁=42  y₁ = 7".normalize(),
            "foo|bar baz=bonk,bonk,bonk x_1=42 y_1=7"
        );

        // Whitespace agnostic desugaring of '<', '>' into '|omit_fwd', '|omit_inv'
        assert_eq!(
            "  : foo>bar <baz  =  bonk,\n: bonk , bonk<zap".normalize(),
            "foo|omit_inv bar|omit_fwd baz=bonk,bonk,bonk|omit_fwd zap"
        );

        // Splitting a pipeline into steps (and handle prefix modifiers)
        let all = "foo>bar |baz  bonk=bonk, bonk , bonk x₁=42 <zap".split_into_steps();
        assert_eq!(all.len(), 4);
        assert_eq!(all[0], "foo");
        assert_eq!(all[1], "bar omit_inv");
        assert_eq!(all[2], "baz bonk=bonk,bonk,bonk x_1=42");
        assert_eq!(all[3], "zap omit_fwd");

        // Parameter splitting
        let args = "foo bar baz=bonk".split_into_parameters();
        assert_eq!(args["_name"], "foo");
        assert_eq!(args["bar"], "true");
        assert_eq!(args["baz"], "bonk");
        assert_eq!("foo bar baz=bonk".operator_name(), "foo");
        assert_eq!("inv foo bar baz=bonk".operator_name(), "foo");

        // Detection of pipelines and resources
        assert!("foo | bar".is_pipeline());
        assert!("foo > bar".is_pipeline());
        assert!("foo < bar".is_pipeline());
        assert!("foo:bar".is_resource_name());

        // Proper handling of subscripts
        let args = "foo x₁=42".split_into_parameters();
        assert_eq!(args["_name"], "foo");
        assert_eq!(args["x_1"], "42");

        // ... and the operator name
        assert_eq!("foo bar baz=  $bonk".operator_name(), "foo");
        Ok(())
    }

    // The PROJ language provides ample opportunity to explore pathological cases
    #[test]
    fn proj() -> Result<(), Error> {
        // Some trivial, but strangely formatted cases
        assert_eq!(
            parse_proj("+a   =   1 +proj =foo    b= 2  ", false)?,
            "foo a=1 b=2"
        );
        assert_eq!(
            parse_proj("+a   =   1 +proj =foo    +   b= 2  ", false)?,
            "foo a=1 b=2"
        );

        // An invalid PROJ string, that parses into an empty pipeline
        assert_eq!(parse_proj("      proj=", false)?, "");

        // A pipeline with a single step and a global argument
        assert_eq!(
            parse_proj("proj=pipeline +foo=bar +step proj=utm zone=32", false)?,
            "utm foo=bar zone=32"
        );

        // A pipeline with 3 steps and 2 global arguments
        assert_eq!(
            parse_proj(
                "proj=pipeline +foo = bar ellps=GRS80 step proj=cart step proj=helmert s=3 step proj=cart ellps=intl",
                false
            )?,
            "cart foo=bar ellps=GRS80 | helmert foo=bar ellps=GRS80 s=3 | cart foo=bar ellps=GRS80 ellps=intl"
        );

        // Although PROJ would choke on this, we accept steps without an initial proj=pipeline
        assert_eq!(
            parse_proj("proj=utm zone=32 step proj=utm inv zone=32", false)?,
            "utm zone=32 | utm inv zone=32"
        );

        // Check for accidental matching of 'step' - even for a hypothetical 'proj=step arg...'
        // and for args called 'step' (which, however, cannot be flags - must come with a value
        // to be recognized as a key=value pair)
        assert_eq!(
            parse_proj(
                "  +step proj = step step=quickstep step step proj=utm inv zone=32 step proj=stepwise step proj=quickstep",
                false
            )?,
            "step step=quickstep | utm inv zone=32 | stepwise | quickstep"
        );

        // Invert the entire pipeline, turning "zone 32-to-zone 33" into "zone 33-to-zone 32"
        // Also throw a few additional spanners in the works, in the form of some ugly, but
        // PROJ-accepted, syntactical abominations
        assert_eq!(
            parse_proj(
                "inv ellps=intl proj=pipeline ugly=syntax +step inv proj=utm zone=32 step proj=utm zone=33",
                false
            )?,
            "utm inv ellps=intl ugly=syntax zone=33 | utm ellps=intl ugly=syntax zone=32"
        );

        // Check for the proper inversion of directional omissions
        assert_eq!(
            parse_proj(
                "proj=pipeline inv   +step   omit_fwd inv proj=utm zone=32   step   omit_inv proj=utm zone=33",
                false
            )?,
            "utm inv omit_fwd zone=33 | utm omit_inv zone=32"
        );

        // Nested pipelines are not supported...

        // Nested pipelines in PROJ requires an `init=` indirection
        assert!(matches!(
            parse_proj("proj=pipeline step proj=pipeline", false),
            Err(Error::Unsupported(_))
        ));
        // ...but `init` is not supported by Rust Geodesy, since that
        // would require a full implementation of PROJ's resolution
        // system - which would be counter to RG's raison d'etre
        assert!(matches!(
            parse_proj("pipeline step init=another_pipeline step proj=noop", false),
            Err(Error::Unsupported(_))
        ));

        // Room here for testing of additional pathological cases...

        // Now check the sanity of the 'pipeline globals' handling
        let mut ctx = Minimal::default();

        // Check that we get the correct argument value when inserting pipeline globals
        // *at the top of the argument list*. Here: x=1 masquerades as the global value,
        // while x=2 is the step local one, which overwrites the global
        let op = ctx.op("helmert x=1 x=2")?;
        let mut operands = crate::test_data::coor2d();
        assert_eq!(2, ctx.apply(op, Fwd, &mut operands)?);
        assert_eq!(operands[0][0], 57.0);
        assert_eq!(operands[1][0], 61.0);

        Ok(())
    }

    #[test]
    fn tidy_proj() -> Result<(), Error> {
        // Ellipsoid defined with `a` and `rf` parameters instead of ellps
        // (explicit pipeline: no degree wrapping even with proj_degrees=true)
        assert_eq!(
            parse_proj(
                "+proj=pipeline +step +inv +proj=tmerc +a=6378249.145 +rf=293.465 +step +proj=step2",
                true
            )?,
            "tmerc inv ellps=6378249.145,293.465 | step2"
        );

        // Ellipsoid is defined with a builtin
        assert_eq!(
            parse_proj("+proj=tmerc +ellps=GRS80", true)?,
            "unitconvert xy_in=deg xy_out=rad | tmerc ellps=GRS80"
        );

        // Ellipsoid is defined with a builtin but is modified by `a` or `rf`
        // Note we don't remove `a` here even though this modification is not supported in RG
        // it's expected to fail in the operator instantiation
        assert_eq!(
            parse_proj("+proj=tmerc +ellps=GRS80 +a=1", true)?,
            "unitconvert xy_in=deg xy_out=rad | tmerc ellps=GRS80 a=1"
        );

        // Replace occurrences of `k=` with `k_0=`
        assert_eq!(
            parse_proj("+proj=tmerc +k=1.5", true)?,
            "unitconvert xy_in=deg xy_out=rad | tmerc k_0=1.5"
        );

        // Convert spherical radius to a zero-flattening ellipsoid
        assert_eq!(
            parse_proj("+proj=laea +lat_0=90 +R=6371228", true)?,
            "unitconvert xy_in=deg xy_out=rad | laea lat_0=90 ellps=6371228,0"
        );

        assert_eq!(
            parse_proj("+proj=laea +R_A +ellps=clrk66 +lat_0=45 +lon_0=-100", true)?,
            "unitconvert xy_in=deg xy_out=rad | laea ellps=6370997.240632998,0 lat_0=45 lon_0=-100"
        );

        assert_eq!(
            parse_proj("+proj=merc +R_C +a=6378136.6 +b=6356751.9 +lat_0=0", true)?,
            "unitconvert xy_in=deg xy_out=rad | merc lat_0=0 ellps=6356751.9,0"
        );

        assert_eq!(
            parse_proj("+proj=merc +ellps=GRS80 +b=6356750 +R_C +lat_0=0", true)?,
            "unitconvert xy_in=deg xy_out=rad | merc ellps=6356750,0 lat_0=0"
        );
        assert_eq!(
            parse_proj("+proj=merc +R_C +R=6378136.6 +lat_0=0", true)?,
            "unitconvert xy_in=deg xy_out=rad | merc lat_0=0 ellps=6378136.6,0"
        );
        assert_eq!(
            parse_proj("+proj=merc +R_C +ellps=WGS84 +lat_0=45", true)?,
            "unitconvert xy_in=deg xy_out=rad | merc ellps=6356752.314245179,0 lat_0=45"
        );

        assert_eq!(
            parse_proj(
                "+proj=cass +lat_0=10.4416666666667 +lon_0=-61.3333333333333 +x_0=86501.46392052 +y_0=65379.0134283 +to_meter=0.201166195164 +a=6378293.64520876 +rf=294.2606763692733",
                true
            )?,
            "unitconvert xy_in=deg xy_out=rad | cass lat_0=10.4416666666667 lon_0=-61.3333333333333 x_0=86501.46392052 y_0=65379.0134283 ellps=6378293.64520876,294.2606763692733 | unitconvert xy_in=m xy_out=0.201166195164"
        );

        assert_eq!(
            parse_proj("+proj=stere +lat_0=90 +a=6378273 +b=6356889.449", true)?,
            "unitconvert xy_in=deg xy_out=rad | stere lat_0=90 ellps=6378273,298.279411123064"
        );
        assert_eq!(
            parse_proj(
                "+proj=stere +variant_c +lat_0=-90 +lat_ts=-67 +lon_0=140 +x_0=300000 +y_0=200000 +ellps=intl",
                true
            )?,
            "unitconvert xy_in=deg xy_out=rad | stere lat_0=-90 lat_ts=-67 lon_0=140 x_0=300000 y_0=200000 ellps=intl"
        );
        assert_eq!(
            parse_proj("+proj=etmerc +ellps=GRS80", true)?,
            "unitconvert xy_in=deg xy_out=rad | tmerc ellps=GRS80"
        );
        assert_eq!(
            parse_proj(
                "+proj=aeqd +guam +ellps=clrk66 +x_0=50000 +y_0=50000 +lon_0=144.74875069444445 +lat_0=13.47246633333333",
                true
            )?,
            "unitconvert xy_in=deg xy_out=rad | guam_aeqd ellps=clrk66 x_0=50000 y_0=50000 lon_0=144.74875069444445 lat_0=13.47246633333333"
        );
        assert_eq!(
            parse_proj("+proj=ups +ellps=GRS80", true)?,
            "unitconvert xy_in=deg xy_out=rad | stere ellps=GRS80 lat_0=90 lon_0=0 k_0=0.994 x_0=2000000 y_0=2000000"
        );
        assert_eq!(
            parse_proj("+proj=ups +south +ellps=GRS80", true)?,
            "unitconvert xy_in=deg xy_out=rad | stere ellps=GRS80 lat_0=-90 lon_0=0 k_0=0.994 x_0=2000000 y_0=2000000"
        );
        assert_eq!(
            parse_proj("+proj=leac +ellps=GRS80 +lat_1=0 +lat_2=2", true)?,
            "unitconvert xy_in=deg xy_out=rad | leac ellps=GRS80 lat_1=0"
        );

        // Generic axis handling should become explicit axisswap steps
        // (explicit pipeline: no degree wrapping)
        assert_eq!(
            parse_proj(
                "+proj=pipeline +step +inv +proj=tmerc +axis=wsu +lon_0=29",
                true
            )?,
            "axisswap order=-1,-2,3 | tmerc inv lon_0=29"
        );

        // Normalize PROJ omerc parameters to the geodesy operator surface
        assert_eq!(
            parse_proj(
                "+proj=omerc +lat_0=4 +lonc=115 +alpha=53 +gamma=52 +k=1",
                true
            )?,
            "unitconvert xy_in=deg xy_out=rad | omerc latc=4 lonc=115 alpha=53 gamma_c=52 k_0=1 variant"
        );
        assert_eq!(
            parse_proj(
                "+proj=omerc +lat_0=4 +lonc=115 +alpha=53 +gamma=52 +no_uoff +pm=paris",
                true
            )?,
            "unitconvert xy_in=deg xy_out=rad | omerc latc=4 lonc=117.33722916666667 alpha=53 gamma_c=52 no_off"
        );

        assert_eq!(
            parse_proj("+proj=longlat +ellps=bessel +pm=0.00289027777777778", true)?,
            "unitconvert xy_in=deg xy_out=rad | lonlat ellps=bessel lon_0=0.00289027777777778 | unitconvert xy_in=rad xy_out=deg"
        );
        assert_eq!(
            parse_proj("+proj=longlat +ellps=bessel +lon_0=10 +pm=paris", true)?,
            "unitconvert xy_in=deg xy_out=rad | lonlat ellps=bessel lon_0=12.337229166666667 | unitconvert xy_in=rad xy_out=deg"
        );
        assert_eq!(
            parse_proj("proj=latlong datum=potsdam ellps=bessel", true)?,
            "unitconvert xy_in=deg xy_out=rad | lonlat datum=potsdam ellps=bessel | unitconvert xy_in=rad xy_out=deg"
        );
        assert_eq!(
            parse_proj("+proj=geoc +ellps=GRS80", true)?,
            "unitconvert xy_in=deg xy_out=rad | latitude ellps=GRS80 geocentric | unitconvert xy_in=rad xy_out=deg"
        );
        assert_eq!(
            parse_proj(
                "+proj=geocent +a=3396190 +b=3376200 +lon_0=0 +units=m",
                true
            )?,
            "unitconvert xy_in=deg xy_out=rad | cart ellps=3396190,169.8944472236118"
        );
        assert_eq!(
            parse_proj("proj=pipeline step proj=longlat ellps=GRS80 geoc inv", true)?,
            "latitude ellps=GRS80 geocentric"
        );
        assert_eq!(
            parse_proj(
                "+proj=pipeline +step +inv +proj=geoc +a=2440530 +b=2438260",
                true
            )?,
            "latitude inv geocentric ellps=2440530,1075.123348017621"
        );
        Ok(())
    }

    #[test]
    fn grid_dependencies_preserve_order_optionality_and_null() -> Result<(), Error> {
        let deps = grid_dependencies(
            "+proj=pipeline \
             +step +proj=hgridshift +grids=@first.gsb,second.tif,@null \
             +step +proj=utm zone=32 \
             +step +proj=gridshift +grids=third.gtx,@fourth.gsb",
        )?;

        assert_eq!(
            deps,
            vec![
                GridDependency {
                    grids: vec![
                        GridSpec {
                            name: "first.gsb".to_string(),
                            optional: true,
                        },
                        GridSpec {
                            name: "second.tif".to_string(),
                            optional: false,
                        },
                    ],
                    null_fallback: true,
                },
                GridDependency {
                    grids: vec![
                        GridSpec {
                            name: "third.gtx".to_string(),
                            optional: false,
                        },
                        GridSpec {
                            name: "fourth.gsb".to_string(),
                            optional: true,
                        },
                    ],
                    null_fallback: false,
                },
            ]
        );
        Ok(())
    }
}
