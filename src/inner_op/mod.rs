use crate::authoring::*;
// ----- B U I L T I N   O P E R A T O R S ---------------------------------------------

// Install new builtin operators by adding them in the `mod` and
// `BUILTIN_OPERATORS` blocks below

/// The coordinate domain an operator works in at its input or output port.
///
/// Used by [`op_domains`] to let callers (e.g. `parse_proj`) decide whether
/// degree↔radian conversion is needed without maintaining a separate list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordDomain {
    /// Angular coordinates — lon/lat, degrees at the PROJ API boundary.
    Geographic,
    /// Linear 2D coordinates — easting/northing or equivalent projected output.
    Projected,
    /// Linear 3D coordinates — geocentric XYZ.
    Cartesian,
}

struct BuiltinOp {
    name: &'static str,
    #[allow(dead_code)]
    title: &'static str,
    constructor: OpConstructor,
    /// `None` means the operator is domain-agnostic (axisswap, unitconvert, …).
    input_domain: Option<CoordDomain>,
    output_domain: Option<CoordDomain>,
}

impl BuiltinOp {
    const fn new(
        name: &'static str,
        ctor: fn(&RawParameters, &dyn Context) -> Result<Op, Error>,
    ) -> Self {
        Self {
            name,
            title: name,
            constructor: OpConstructor(ctor),
            input_domain: None,
            output_domain: None,
        }
    }

    const fn with_domains(
        name: &'static str,
        ctor: fn(&RawParameters, &dyn Context) -> Result<Op, Error>,
        input_domain: CoordDomain,
        output_domain: CoordDomain,
    ) -> Self {
        Self {
            name,
            title: name,
            constructor: OpConstructor(ctor),
            input_domain: Some(input_domain),
            output_domain: Some(output_domain),
        }
    }

    const fn point<T: PointOp>() -> Self {
        Self {
            name: T::NAME,
            title: T::TITLE,
            constructor: OpConstructor(Op::point::<T>),
            input_domain: None,
            output_domain: None,
        }
    }

    const fn point_with_domains<T: PointOp>(
        input_domain: CoordDomain,
        output_domain: CoordDomain,
    ) -> Self {
        Self {
            name: T::NAME,
            title: T::TITLE,
            constructor: OpConstructor(Op::point::<T>),
            input_domain: Some(input_domain),
            output_domain: Some(output_domain),
        }
    }
}

mod addone;
mod curvature;
mod deflection;
mod deformation;
mod geodesic;
mod gravity;
mod gridshift;
mod helmert;
mod iso6709;
mod molobadekas;
mod molodensky;
mod noop;
mod permtide;
mod deprecated;
mod experimental;
mod projections;
mod stack;
mod transforms;

pub(crate) mod pipeline; // Needed by Op for instantiation

use self::{deprecated::*, experimental::*, projections::*, transforms::*};
#[allow(unused_imports)]
pub(crate) use crate::projection::ProjectionFrame;
use CoordDomain::{Cartesian, Geographic, Projected};

#[rustfmt::skip]
const BUILTIN_OPERATORS: [BuiltinOp; 75] = [
    // Geographic projections: main/canonical set, lon/lat degrees in, projected (linear) out
    BuiltinOp::point_with_domains::<aea::Aea>(Geographic, Projected),
    BuiltinOp::point_with_domains::<aeqd::Aeqd>(Geographic, Projected),
    BuiltinOp::with_domains("cass",        cass::new,        Geographic, Projected),
    BuiltinOp::point_with_domains::<cea::Cea>(Geographic, Projected),
    BuiltinOp::with_domains("etmerc",      tmerc::new,       Geographic, Projected),
    BuiltinOp::point_with_domains::<laea::Laea>(Geographic, Projected),
    BuiltinOp::with_domains("lcc",         lcc::new,         Geographic, Projected),
    BuiltinOp::point_with_domains::<leac::Leac>(Geographic, Projected),
    BuiltinOp::point_with_domains::<merc::Merc>(Geographic, Projected),
    BuiltinOp::with_domains("omerc",       omerc::new,       Geographic, Projected),
    BuiltinOp::point_with_domains::<stere::Stere>(Geographic, Projected),
    BuiltinOp::point_with_domains::<sterec::Sterec>(Geographic, Projected),
    BuiltinOp::point_with_domains::<sterea::Sterea>(Geographic, Projected),
    BuiltinOp::with_domains("tmerc",       tmerc::new,       Geographic, Projected),
    BuiltinOp::with_domains("utm",         tmerc::utm,       Geographic, Projected),
    BuiltinOp::point_with_domains::<ups::Ups>(Geographic, Projected),
    BuiltinOp::point_with_domains::<webmerc::WebMerc>(Geographic, Projected),

    // Geographic projections: experimental / secondary set
    BuiltinOp::with_domains("alsk",        mod_ster::alsk,   Geographic, Projected),
    BuiltinOp::with_domains("bonne",       bonne::new,       Geographic, Projected),
    BuiltinOp::with_domains("calcofi",     calcofi::new,     Geographic, Projected),
    BuiltinOp::with_domains("col_urban",   col_urban::new,   Geographic, Projected),
    BuiltinOp::with_domains("eqc",         eqc::new,         Geographic, Projected),
    BuiltinOp::with_domains("eqdc",        eqdc::new,        Geographic, Projected),
    BuiltinOp::with_domains("eqearth",     eqearth::new,     Geographic, Projected),
    BuiltinOp::with_domains("geos",        geos::new,        Geographic, Projected),
    BuiltinOp::point_with_domains::<guam_aeqd::GuamAeqd>(Geographic, Projected),
    BuiltinOp::with_domains("gs48",        mod_ster::gs48,   Geographic, Projected),
    BuiltinOp::with_domains("gs50",        mod_ster::gs50,   Geographic, Projected),
    BuiltinOp::with_domains("gstmerc",     gstmerc::new,     Geographic, Projected),
    BuiltinOp::with_domains("krovak",      krovak::new,      Geographic, Projected),
    BuiltinOp::with_domains("labrd",       labrd::new,       Geographic, Projected),
    BuiltinOp::with_domains("lccnc",       lccnc::new,       Geographic, Projected),
    BuiltinOp::with_domains("lcca",        lcca::new,        Geographic, Projected),
    BuiltinOp::point_with_domains::<mill::Mill>(Geographic, Projected),
    BuiltinOp::with_domains("moll",        moll::new,        Geographic, Projected), 
    BuiltinOp::with_domains("nzmg",        nzmg::new,        Geographic, Projected),
    BuiltinOp::with_domains("ocea",        ocea::new,        Geographic, Projected),
    BuiltinOp::with_domains("ortho",       ortho::new,       Geographic, Projected),
    BuiltinOp::with_domains("poly",        poly::new,        Geographic, Projected),
    BuiltinOp::with_domains("mod_krovak",  krovak::modified, Geographic, Projected),
    BuiltinOp::with_domains("qsc",         qsc::new,         Geographic, Projected),
    BuiltinOp::with_domains("robin",       robin::new,       Geographic, Projected),
    BuiltinOp::with_domains("rouss",       rouss::new,       Geographic, Projected),
    BuiltinOp::with_domains("sinu",        sinu::new,        Geographic, Projected),
    BuiltinOp::with_domains("som",         som::new,         Geographic, Projected),
    BuiltinOp::with_domains("somerc",      somerc::new,      Geographic, Projected),
    BuiltinOp::with_domains("tcea",        tcea::new,        Geographic, Projected),
    BuiltinOp::with_domains("tmgrid",      tmgrid::new,      Geographic, Projected),
    BuiltinOp::with_domains("tpeqd",       tpeqd::new,       Geographic, Projected),
    BuiltinOp::with_domains("tunmg",       tunmg::new,       Geographic, Projected),
    BuiltinOp::with_domains("lsat",        som::lsat,        Geographic, Projected),
    BuiltinOp::with_domains("misrsom",     som::misr,        Geographic, Projected),

    // Geographic transforms / conversions
    BuiltinOp::point::<adapt::Adapt>(),
    BuiltinOp::point::<axisswap::AxisSwap>(),
    BuiltinOp::point_with_domains::<lonlat::LonLat>(Geographic, Geographic),
    BuiltinOp::point_with_domains::<cart::Cart>(Geographic, Cartesian),
    BuiltinOp::point_with_domains::<latitude::Latitude>(Geographic, Geographic),
    BuiltinOp::point::<unitconvert::UnitConvert>(),

    // Domain-agnostic operators
    BuiltinOp::point::<addone::AddOne>(),
    BuiltinOp::new("curvature",   curvature::new),
    BuiltinOp::new("deflection",  deflection::new),
    BuiltinOp::new("deformation", deformation::new),
    BuiltinOp::new("dm",          iso6709::dm),
    BuiltinOp::new("dms",         iso6709::dms),
    BuiltinOp::new("geodesic",    geodesic::new),
    BuiltinOp::point_with_domains::<geogoffset::GeogOffset>(Geographic, Geographic),
    BuiltinOp::new("gravity",     gravity::new),
    BuiltinOp::new("gridshift",   gridshift::new),
    BuiltinOp::new("helmert",     helmert::new),
    BuiltinOp::new("molobadekas", molobadekas::new),
    BuiltinOp::with_domains("molodensky",  molodensky::new, Geographic, Geographic),
    BuiltinOp::point::<permtide::PermTide>(),

    // Pipeline / flow control
    BuiltinOp::new("pipeline",    pipeline::new),
    BuiltinOp::new("stack",       stack::new),
    BuiltinOp::point::<noop::NoOp>(),
];

#[rustfmt::skip]
const DEPRECATED_OPERATORS: [BuiltinOp; 4] = [
    BuiltinOp::with_domains("btmerc",      btmerc::new,      Geographic, Projected),
    BuiltinOp::with_domains("butm",        btmerc::utm,      Geographic, Projected),
    BuiltinOp::new("pop",         pushpop::pop),
    BuiltinOp::new("push",        pushpop::push),
];
// A BTreeMap would have been a better choice for BUILTIN_OPERATORS, except
// for the annoying fact that it cannot be compile-time const-constructed.

/// Handle instantiation of built-in operators, as defined in `BUILTIN_OPERATORS` above.
pub(crate) fn builtin(name: &str) -> Result<OpConstructor, Error> {
    for p in &BUILTIN_OPERATORS {
        if p.name == name {
            return Ok(p.constructor);
        }
    }
    for p in &DEPRECATED_OPERATORS {
        if p.name == name {
            return Ok(p.constructor);
        }
    }
    Err(Error::NotFound(name.to_string(), String::default()))
}

/// Return the input and output [`CoordDomain`] for a built-in operator.
///
/// Returns `(None, None)` for domain-agnostic operators (axisswap, unitconvert, …)
/// and for unknown operator names.
pub(crate) fn op_domains(name: &str) -> (Option<CoordDomain>, Option<CoordDomain>) {
    BUILTIN_OPERATORS
        .iter()
        .find(|e| e.name == name)
        .or_else(|| DEPRECATED_OPERATORS.iter().find(|e| e.name == name))
        .map(|e| (e.input_domain, e.output_domain))
        .unwrap_or((None, None))
}

#[allow(dead_code)]
pub(crate) fn op_title(name: &str) -> Option<&'static str> {
    BUILTIN_OPERATORS
        .iter()
        .find(|e| e.name == name)
        .or_else(|| DEPRECATED_OPERATORS.iter().find(|e| e.name == name))
        .map(|e| e.title)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_titles_expose_human_readable_projection_names() {
        assert_eq!(op_title("cea"), Some("Equal Area Cylindrical"));
        assert_eq!(op_title("merc"), Some("Mercator"));
        assert_eq!(op_title("lonlat"), Some("Lat/long (Geodetic alias)"));
        assert_eq!(op_title("definitely_not_an_op"), None);
    }
}

pub(crate) fn apply_utm_defaults(params: &mut ParsedParameters, zone: usize) {
    params.real.insert("k_0", 0.9996);
    params
        .real
        .insert("lon_0", (-183.0 + 6.0 * zone as f64).to_radians());
    params.real.insert("lat_0", 0.0);
    params.real.insert("x_0", 500_000.0);
    params.real.insert(
        "y_0",
        if params.boolean("south") {
            10_000_000.0
        } else {
            0.0
        },
    );
}

pub(crate) fn mark_spherical(params: &mut ParsedParameters) -> bool {
    let spherical = params.ellps(0).flattening() == 0.0;
    if spherical {
        params.boolean.insert("spherical");
    }
    spherical
}

pub(crate) fn insert_rectifying_setup(params: &mut ParsedParameters) -> bool {
    let spherical = params.ellps(0).flattening() == 0.0;
    if spherical {
        params.boolean.insert("spherical");
        return true;
    }

    let rectifying = params
        .ellps(0)
        .coefficients_for_rectifying_latitude_computations();
    params.fourier_coefficients.insert("rectifying", rectifying);
    false
}

// ----- S T R U C T   O P C O N S T R U C T O R ---------------------------------------

/// Blueprint for the overall instantiation of an operator.
///
/// OpConstructor needs to be a newtype, rather than a type alias,
/// since we must implement the Debug-trait for OpConstructor (to
/// make auto derive of the Debug-trait work for any derived type).
#[derive(Clone, Copy)]
pub struct OpConstructor(pub fn(args: &RawParameters, ctx: &dyn Context) -> Result<Op, Error>);

// Cannot autoderive the Debug trait
impl core::fmt::Debug for OpConstructor {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "OpConstructor")
    }
}

// ----- S T R U C T   I N N E R O P ---------------------------------------------------

/// Blueprint for the functions doing the actual transformation work.
///
/// InnerOp needs to be a newtype, rather than a type alias, since we
/// must implement the Debug-trait for InnerOp (to make auto derive
/// of the Debug-trait work for any derived type).
pub struct InnerOp(pub fn(op: &Op, ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize);

// Cannot autoderive the Debug trait
impl core::fmt::Debug for InnerOp {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "InnerOp")
    }
}

// Defaults to no_op. Used by the pushpop and stack pseudo-operators
impl Default for InnerOp {
    fn default() -> InnerOp {
        InnerOp(noop_placeholder)
    }
}

fn noop_placeholder(_op: &Op, _ctx: &dyn Context, _operands: &mut dyn CoordinateSet) -> usize {
    // Consider whether this should return an Err-value if used as a placeholder for a
    // non-existing or non-implemented inverse operation
    0
}
