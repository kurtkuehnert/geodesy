use crate::authoring::*;
use std::collections::BTreeMap;

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
            constructor: OpConstructor(ctor),
            input_domain: Some(input_domain),
            output_domain: Some(output_domain),
        }
    }
}

mod adapt;
mod addone;
mod aea;
mod aeqd;
mod axisswap;
mod btmerc;
mod bonne;
mod calcofi;
mod cart;
mod cass;
mod cea;
mod col_urban;
mod curvature;
mod deflection;
mod deformation;
mod eqc;
mod eqdc;
mod eqearth;
mod geodesic;
mod geos;
mod gravity;
mod gridshift;
mod gstmerc;
mod helmert;
mod iso6709;
mod krovak;
mod labrd;
mod laea;
mod latitude;
mod lcc;
mod lccnc;
mod lcca;
mod longlat;
mod merc;
mod mill;
mod mod_ster;
mod moll;
mod molobadekas;
mod molodensky;
mod noop;
mod nzmg;
mod ocea;
mod omerc;
mod ortho;
mod permtide;
pub(crate) mod pipeline; // Needed by Op for instantiation
mod poly;
mod pushpop;
mod qsc;
mod robin;
mod rouss;
mod sinu;
mod som;
mod somerc;
mod stack;
mod stere;
mod sterea;
mod tcea;
mod tmerc;
mod tmgrid;
mod tpeqd;
mod tunmg;
mod unitconvert;
mod units;
mod webmerc;

use CoordDomain::{Cartesian, Geographic, Projected};

#[rustfmt::skip]
const BUILTIN_OPERATORS: [BuiltinOp; 79] = [
    // Geographic projections: lon/lat degrees in, projected (linear) out
    BuiltinOp::with_domains("aea",         aea::new,         Geographic, Projected),
    BuiltinOp::with_domains("aeqd",        aeqd::new,        Geographic, Projected),
    BuiltinOp::with_domains("btmerc",      btmerc::new,      Geographic, Projected),
    BuiltinOp::with_domains("butm",        btmerc::utm,      Geographic, Projected),
    BuiltinOp::with_domains("alsk",        mod_ster::alsk,   Geographic, Projected),
    BuiltinOp::with_domains("bonne",       bonne::new,       Geographic, Projected),
    BuiltinOp::with_domains("cass",        cass::new,        Geographic, Projected),
    BuiltinOp::with_domains("calcofi",     calcofi::new,     Geographic, Projected),
    BuiltinOp::with_domains("cea",         cea::new,         Geographic, Projected),
    BuiltinOp::with_domains("col_urban",   col_urban::new,   Geographic, Projected),
    BuiltinOp::with_domains("eqc",         eqc::new,         Geographic, Projected),
    BuiltinOp::with_domains("eqdc",        eqdc::new,        Geographic, Projected),
    BuiltinOp::with_domains("eqearth",     eqearth::new,     Geographic, Projected),
    BuiltinOp::with_domains("etmerc",      tmerc::new,       Geographic, Projected),
    BuiltinOp::with_domains("geos",        geos::new,        Geographic, Projected),
    BuiltinOp::with_domains("gs48",        mod_ster::gs48,   Geographic, Projected),
    BuiltinOp::with_domains("gs50",        mod_ster::gs50,   Geographic, Projected),
    BuiltinOp::with_domains("gstmerc",     gstmerc::new,     Geographic, Projected),
    BuiltinOp::with_domains("krovak",      krovak::new,      Geographic, Projected),
    BuiltinOp::with_domains("laea",        laea::new,        Geographic, Projected),
    BuiltinOp::with_domains("labrd",       labrd::new,       Geographic, Projected),
    BuiltinOp::with_domains("lcc",         lcc::new,         Geographic, Projected),
    BuiltinOp::with_domains("lccnc",       lccnc::new,       Geographic, Projected),
    BuiltinOp::with_domains("lcca",        lcca::new,        Geographic, Projected),
    BuiltinOp::with_domains("leac",        aea::leac,        Geographic, Projected),
    BuiltinOp::with_domains("merc",        merc::new,        Geographic, Projected),
    BuiltinOp::with_domains("mill",        mill::new,        Geographic, Projected),
    BuiltinOp::with_domains("moll",        moll::new,        Geographic, Projected),
    BuiltinOp::with_domains("nzmg",        nzmg::new,        Geographic, Projected),
    BuiltinOp::with_domains("ocea",        ocea::new,        Geographic, Projected),
    BuiltinOp::with_domains("omerc",       omerc::new,       Geographic, Projected),
    BuiltinOp::with_domains("ortho",       ortho::new,       Geographic, Projected),
    BuiltinOp::with_domains("poly",        poly::new,        Geographic, Projected),
    BuiltinOp::with_domains("mod_krovak",  krovak::modified, Geographic, Projected),
    BuiltinOp::with_domains("qsc",         qsc::new,         Geographic, Projected),
    BuiltinOp::with_domains("robin",       robin::new,       Geographic, Projected),
    BuiltinOp::with_domains("rouss",       rouss::new,       Geographic, Projected),
    BuiltinOp::with_domains("sinu",        sinu::new,        Geographic, Projected),
    BuiltinOp::with_domains("som",         som::new,         Geographic, Projected),
    BuiltinOp::with_domains("somerc",      somerc::new,      Geographic, Projected),
    BuiltinOp::with_domains("stere",       stere::new,       Geographic, Projected),
    BuiltinOp::with_domains("sterea",      sterea::new,      Geographic, Projected),
    BuiltinOp::with_domains("tcea",        tcea::new,        Geographic, Projected),
    BuiltinOp::with_domains("tmerc",       tmerc::new,       Geographic, Projected),
    BuiltinOp::with_domains("tmgrid",      tmgrid::new,      Geographic, Projected),
    BuiltinOp::with_domains("tpeqd",       tpeqd::new,       Geographic, Projected),
    BuiltinOp::with_domains("tunmg",       tunmg::new,       Geographic, Projected),
    BuiltinOp::with_domains("utm",         tmerc::utm,       Geographic, Projected),
    BuiltinOp::with_domains("ups",         stere::ups,       Geographic, Projected),
    BuiltinOp::with_domains("webmerc",     webmerc::new,     Geographic, Projected),

    // Geographic identity: lon/lat degrees in and out
    BuiltinOp::with_domains("longlat",     longlat::new,     Geographic, Geographic),
    BuiltinOp::with_domains("latlon",      longlat::new,     Geographic, Geographic),
    BuiltinOp::with_domains("latlong",     longlat::new,     Geographic, Geographic),
    BuiltinOp::with_domains("lonlat",      longlat::new,     Geographic, Geographic),

    // Geographic to Cartesian: lon/lat degrees in, XYZ meters out
    BuiltinOp::with_domains("cart",        cart::new,        Geographic, Cartesian),

    // Domain-agnostic operators
    BuiltinOp::new("adapt",       adapt::new),
    BuiltinOp::new("addone",      addone::new),
    BuiltinOp::new("axisswap",    axisswap::new),
    BuiltinOp::new("curvature",   curvature::new),
    BuiltinOp::new("deflection",  deflection::new),
    BuiltinOp::new("deformation", deformation::new),
    BuiltinOp::new("dm",          iso6709::dm),
    BuiltinOp::new("dms",         iso6709::dms),
    BuiltinOp::new("geodesic",    geodesic::new),
    BuiltinOp::new("gravity",     gravity::new),
    BuiltinOp::new("gridshift",   gridshift::new),
    BuiltinOp::new("helmert",     helmert::new),
    BuiltinOp::new("latitude",    latitude::new),
    BuiltinOp::with_domains("lsat",        som::lsat,        Geographic, Projected),
    BuiltinOp::with_domains("misrsom",     som::misr,        Geographic, Projected),
    BuiltinOp::new("molobadekas", molobadekas::new),
    BuiltinOp::with_domains("molodensky",  molodensky::new, Geographic, Geographic),
    BuiltinOp::new("permtide",    permtide::new),
    BuiltinOp::new("unitconvert", unitconvert::new),

    // Pipeline / flow control
    BuiltinOp::new("pipeline",    pipeline::new),
    BuiltinOp::new("pop",         pushpop::pop),
    BuiltinOp::new("push",        pushpop::push),
    BuiltinOp::new("stack",       stack::new),
    BuiltinOp::new("noop",        noop::new),
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
        .map(|e| (e.input_domain, e.output_domain))
        .unwrap_or((None, None))
}

pub(crate) fn override_ellps_from_proj_params(
    params: &mut ParsedParameters,
    def: &str,
    given: &BTreeMap<String, String>,
) -> Result<(), Error> {
    if given.contains_key("ellps") {
        return Ok(());
    }

    let a = if let Some(a) = given.get("a") {
        a.parse::<f64>()
            .map_err(|_| Error::BadParam("a".to_string(), def.to_string()))?
    } else if let Some(r) = given.get("R") {
        r.parse::<f64>()
            .map_err(|_| Error::BadParam("R".to_string(), def.to_string()))?
    } else {
        return Ok(());
    };

    if a <= 0.0 {
        return Err(Error::BadParam("a".to_string(), def.to_string()));
    }

    let rf = if let Some(rf) = given.get("rf") {
        rf.parse::<f64>()
            .map_err(|_| Error::BadParam("rf".to_string(), def.to_string()))?
    } else if let Some(f) = given.get("f") {
        let f = f
            .parse::<f64>()
            .map_err(|_| Error::BadParam("f".to_string(), def.to_string()))?;
        if f == 0.0 { 0.0 } else { 1.0 / f }
    } else if let Some(b) = given.get("b") {
        let b = b
            .parse::<f64>()
            .map_err(|_| Error::BadParam("b".to_string(), def.to_string()))?;
        if b <= 0.0 {
            return Err(Error::BadParam("b".to_string(), def.to_string()));
        }
        if (a - b).abs() < f64::EPSILON {
            0.0
        } else {
            a / (a - b)
        }
    } else {
        0.0
    };

    params.text.insert("ellps", format!("{a},{rf}"));
    Ok(())
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
