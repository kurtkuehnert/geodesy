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

    const fn point<T: PointOp>() -> Self {
        Self {
            name: T::NAME,
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
            constructor: OpConstructor(Op::point::<T>),
            input_domain: Some(input_domain),
            output_domain: Some(output_domain),
        }
    }
}

mod adapt;
mod addone;
mod aea;
mod aeqd;
mod btmerc;
mod cart;
mod cass;
mod cea;
mod curvature;
mod deflection;
mod deformation;
mod geodesic;
mod gravity;
mod gridshift;
mod helmert;
mod iso6709;
mod laea;
mod lcc;
mod leac;
mod lonlat;
mod merc;
mod molobadekas;
mod molodensky;
mod noop;
mod omerc;
mod permtide;
mod stack;
mod stere;
mod sterea;
mod sterec;
mod tmerc;
mod unitconvert;
mod ups;
mod webmerc;

mod axisswap;
pub(crate) mod pipeline; // Needed by Op for instantiation
#[path = "deprecated/pushpop.rs"]
mod pushpop;

// Experimental operators
#[path = "experimental/bonne.rs"]
mod bonne;
#[path = "experimental/calcofi.rs"]
mod calcofi;
#[path = "experimental/col_urban.rs"]
mod col_urban;
#[path = "experimental/eqc.rs"]
mod eqc;
#[path = "experimental/eqdc.rs"]
mod eqdc;
#[path = "experimental/eqearth.rs"]
mod eqearth;
#[path = "experimental/geogoffset.rs"]
mod geogoffset;
#[path = "experimental/geos.rs"]
mod geos;
#[path = "experimental/gstmerc.rs"]
mod gstmerc;
#[path = "experimental/guam_aeqd.rs"]
mod guam_aeqd;
#[path = "experimental/krovak.rs"]
mod krovak;
#[path = "experimental/labrd.rs"]
mod labrd;
#[path = "experimental/lcca.rs"]
mod lcca;
#[path = "experimental/lccnc.rs"]
mod lccnc;
#[path = "experimental/mill.rs"]
mod mill;
#[path = "experimental/mod_ster.rs"]
mod mod_ster;
#[path = "experimental/moll.rs"]
mod moll;
#[path = "experimental/nzmg.rs"]
mod nzmg;
#[path = "experimental/ocea.rs"]
mod ocea;
#[path = "experimental/ortho.rs"]
mod ortho;
#[path = "experimental/poly.rs"]
mod poly;
#[path = "experimental/qsc.rs"]
mod qsc;
#[path = "experimental/robin.rs"]
mod robin;
#[path = "experimental/rouss.rs"]
mod rouss;
#[path = "experimental/sinu.rs"]
mod sinu;
#[path = "experimental/som.rs"]
mod som;
#[path = "experimental/somerc.rs"]
mod somerc;
#[path = "experimental/tcea.rs"]
mod tcea;
#[path = "experimental/tmgrid.rs"]
mod tmgrid;
#[path = "experimental/tpeqd.rs"]
mod tpeqd;
#[path = "experimental/tunmg.rs"]
mod tunmg;

mod latitude;
#[allow(unused_imports)]
pub(crate) use crate::projection::ProjectionFrame;
use CoordDomain::{Cartesian, Geographic, Projected};

#[rustfmt::skip]
const BUILTIN_OPERATORS: [BuiltinOp; 79] = [
    // Geographic projections: lon/lat degrees in, projected (linear) out
    BuiltinOp::point_with_domains::<aea::Aea>(Geographic, Projected),
    BuiltinOp::point_with_domains::<aeqd::Aeqd>(Geographic, Projected),
    BuiltinOp::with_domains("btmerc",      btmerc::new,      Geographic, Projected),
    BuiltinOp::with_domains("butm",        btmerc::utm,      Geographic, Projected),
    BuiltinOp::with_domains("alsk",        mod_ster::alsk,   Geographic, Projected),
    BuiltinOp::with_domains("bonne",       bonne::new,       Geographic, Projected),
    BuiltinOp::with_domains("cass",        cass::new,        Geographic, Projected),
    BuiltinOp::with_domains("calcofi",     calcofi::new,     Geographic, Projected),
    BuiltinOp::point_with_domains::<cea::Cea>(Geographic, Projected),
    BuiltinOp::with_domains("col_urban",   col_urban::new,   Geographic, Projected),
    BuiltinOp::with_domains("eqc",         eqc::new,         Geographic, Projected),
    BuiltinOp::with_domains("eqdc",        eqdc::new,        Geographic, Projected),
    BuiltinOp::with_domains("eqearth",     eqearth::new,     Geographic, Projected),
    BuiltinOp::with_domains("etmerc",      tmerc::new,       Geographic, Projected),
    BuiltinOp::with_domains("geos",        geos::new,        Geographic, Projected),
    BuiltinOp::with_domains("guam_aeqd",   guam_aeqd::new,   Geographic, Projected),
    BuiltinOp::with_domains("gs48",        mod_ster::gs48,   Geographic, Projected),
    BuiltinOp::with_domains("gs50",        mod_ster::gs50,   Geographic, Projected),
    BuiltinOp::with_domains("gstmerc",     gstmerc::new,     Geographic, Projected),
    BuiltinOp::with_domains("krovak",      krovak::new,      Geographic, Projected),
    BuiltinOp::point_with_domains::<laea::Laea>(Geographic, Projected),
    BuiltinOp::with_domains("labrd",       labrd::new,       Geographic, Projected),
    BuiltinOp::with_domains("lcc",         lcc::new,         Geographic, Projected),
    BuiltinOp::with_domains("lccnc",       lccnc::new,       Geographic, Projected),
    BuiltinOp::with_domains("lcca",        lcca::new,        Geographic, Projected),
    BuiltinOp::point_with_domains::<leac::Leac>(Geographic, Projected),
    BuiltinOp::point_with_domains::<merc::Merc>(Geographic, Projected),
    BuiltinOp::point_with_domains::<mill::Mill>(Geographic, Projected),
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
    BuiltinOp::point_with_domains::<stere::Stere>(Geographic, Projected),
    BuiltinOp::point_with_domains::<sterec::Sterec>(Geographic, Projected), 
    BuiltinOp::point_with_domains::<sterea::Sterea>(Geographic, Projected),
    BuiltinOp::with_domains("tcea",        tcea::new,        Geographic, Projected),
    BuiltinOp::with_domains("tmerc",       tmerc::new,       Geographic, Projected),
    BuiltinOp::with_domains("tmgrid",      tmgrid::new,      Geographic, Projected),
    BuiltinOp::with_domains("tpeqd",       tpeqd::new,       Geographic, Projected),
    BuiltinOp::with_domains("tunmg",       tunmg::new,       Geographic, Projected),
    BuiltinOp::with_domains("utm",         tmerc::utm,       Geographic, Projected),
    BuiltinOp::point_with_domains::<ups::Ups>(Geographic, Projected),
    BuiltinOp::point_with_domains::<webmerc::WebMerc>(Geographic, Projected),

    // Geographic identity: lon/lat degrees in and out
    BuiltinOp::point_with_domains::<lonlat::LonLat>(Geographic, Geographic),

    // Geographic to Cartesian: lon/lat degrees in, XYZ meters out
    BuiltinOp::point_with_domains::<cart::Cart>(Geographic, Cartesian),

    // Domain-agnostic operators
    BuiltinOp::point::<adapt::Adapt>(),
    BuiltinOp::point::<addone::AddOne>(),
    BuiltinOp::point::<axisswap::AxisSwap>(),
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
    BuiltinOp::point_with_domains::<latitude::Latitude>(Geographic, Geographic),
    BuiltinOp::with_domains("lsat",        som::lsat,        Geographic, Projected),
    BuiltinOp::with_domains("misrsom",     som::misr,        Geographic, Projected),
    BuiltinOp::new("molobadekas", molobadekas::new),
    BuiltinOp::with_domains("molodensky",  molodensky::new, Geographic, Geographic),
    BuiltinOp::point::<permtide::PermTide>(),
    BuiltinOp::point::<unitconvert::UnitConvert>(),

    // Pipeline / flow control
    BuiltinOp::new("pipeline",    pipeline::new),
    BuiltinOp::new("pop",         pushpop::pop),
    BuiltinOp::new("push",        pushpop::push),
    BuiltinOp::new("stack",       stack::new),
    BuiltinOp::point::<noop::NoOp>(),
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
