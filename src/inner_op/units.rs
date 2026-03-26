// Units are taken from PROJ https://github.com/OSGeo/PROJ/blob/master/src/units.c

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnitKind {
    Linear,
    Angular,
}

impl UnitKind {
    /// Only named linear↔angular mismatches are rejected.
    pub fn is_compatible_with(self, other: Self) -> bool {
        self == other
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Unit {
    name: &'static str,
    description: &'static str,
    multiplier: f64,
    kind: UnitKind,
}

impl Unit {
    const fn linear(name: &'static str, description: &'static str, multiplier: f64) -> Self {
        Self {
            name,
            description,
            multiplier,
            kind: UnitKind::Linear,
        }
    }

    const fn angular(name: &'static str, description: &'static str, multiplier: f64) -> Self {
        Self {
            name,
            description,
            multiplier,
            kind: UnitKind::Angular,
        }
    }

    pub fn name(self) -> &'static str {
        self.name
    }
    pub fn multiplier(self) -> f64 {
        self.multiplier
    }
    pub fn kind(self) -> UnitKind {
        self.kind
    }
}

/// A resolved unit parameter: either a named unit with known kind, or a raw numeric scale factor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitParam {
    Named(Unit),
    /// Raw scale factor — no semantic kind, compatible with any unit.
    Numeric(f64),
}

impl UnitParam {
    pub fn lookup(name: &str) -> Option<Self> {
        LINEAR_UNITS
            .iter()
            .chain(ANGULAR_UNITS.iter())
            .find(|u| u.name() == name)
            .copied()
            .map(Self::Named)
            .or_else(|| {
                let m = name.parse::<f64>().ok()?;
                if !m.is_finite() || m <= 0.0 {
                    return None;
                }
                Some(Self::Numeric(m))
            })
    }

    pub fn multiplier(self) -> f64 {
        match self {
            Self::Named(u) => u.multiplier(),
            Self::Numeric(m) => m,
        }
    }

    /// Numeric is compatible with any unit. Only named linear↔angular is rejected.
    pub fn is_compatible_with(self, other: Self) -> bool {
        match (self, other) {
            (Self::Named(a), Self::Named(b)) => a.kind().is_compatible_with(b.kind()),
            _ => true,
        }
    }
}

/// Linear units and their conversion factor to meters.
#[rustfmt::skip]
pub const LINEAR_UNITS: [Unit; 21] = [
    Unit::linear("km",     "Kilometer",                    1000.0),
    Unit::linear("m",      "Meter",                        1.0),
    Unit::linear("dm",     "Decimeter",                    0.1),
    Unit::linear("cm",     "Centimeter",                   0.01),
    Unit::linear("mm",     "Millimeter",                   0.001),
    Unit::linear("kmi",    "International Nautical Mile",  1852.0),
    Unit::linear("in",     "International Inch",           0.0254),
    Unit::linear("ft",     "International Foot",           0.3048),
    Unit::linear("yd",     "International Yard",           0.9144),
    Unit::linear("mi",     "International Statute Mile",   1609.344),
    Unit::linear("fath",   "International Fathom",         1.8288),
    Unit::linear("ch",     "International Chain",          20.1168),
    Unit::linear("link",   "International Link",           0.201168),
    Unit::linear("us-in",  "U.S. Surveyor's Inch",         100.0 / 3937.0),
    Unit::linear("us-ft",  "U.S. Surveyor's Foot",         1200.0 / 3937.0),
    Unit::linear("us-yd",  "U.S. Surveyor's Yard",         3600.0 / 3937.0),
    Unit::linear("us-ch",  "U.S. Surveyor's Chain",        79200.0 / 3937.0),
    Unit::linear("us-mi",  "U.S. Surveyor's Statute Mile", 6336000.0 / 3937.0),
    Unit::linear("ind-yd", "Indian Yard",                  0.91439523),
    Unit::linear("ind-ft", "Indian Foot",                  0.30479841),
    Unit::linear("ind-ch", "Indian Chain",                 20.11669506),
];

/// Angular units and their conversion factor to radians.
#[rustfmt::skip]
pub const ANGULAR_UNITS: [Unit; 3] = [
    Unit::angular("rad",  "Radian",  1.0),
    Unit::angular("deg",  "Degree",  std::f64::consts::PI / 180.0),
    Unit::angular("grad", "Grad",    std::f64::consts::PI / 200.0),
];
