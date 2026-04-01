//! Lambert Conformal Conic
use crate::authoring::*;
use crate::projection::ProjectionFrame;
use std::f64::consts::FRAC_PI_2;

const ANGULAR_TOLERANCE: f64 = 1e-10;

#[rustfmt::skip]
pub const GAMUT: [OpParameter; 10] = [
    OpParameter::Flag { key: "inv" },
    OpParameter::Flag { key: "west_oriented" },
    OpParameter::Text { key: "ellps", default: Some("GRS80") },
    OpParameter::Real { key: "lat_1", default: Some(0_f64) },
    OpParameter::Real { key: "lat_2", default: Some(f64::NAN) },
    OpParameter::Real { key: "lat_0", default: Some(f64::NAN) },
    OpParameter::Real { key: "lon_0", default: Some(0_f64) },
    OpParameter::Real { key: "k_0",   default: Some(1_f64) },
    OpParameter::Real { key: "x_0",   default: Some(0_f64) },
    OpParameter::Real { key: "y_0",   default: Some(0_f64) },
];

#[derive(Clone, Copy, Debug)]
struct LccState {
    frame: ProjectionFrame,
    ellps: Ellipsoid,
    n: f64,
    c: f64,
    rho0: f64,
    west_oriented: bool,
}

impl LccState {
    fn new(
        params: &ParsedParameters,
        frame: ProjectionFrame,
        ellps: Ellipsoid,
        phi1: f64,
        phi2: f64,
        lat_0: f64,
    ) -> Result<Self, Error> {
        let sc = phi1.sin_cos();
        let mut n = sc.0;
        let e = ellps.eccentricity();
        let es = ellps.eccentricity_squared();

        if (phi1 + phi2).abs() < ANGULAR_TOLERANCE {
            return Err(Error::General(
                "Lcc: Invalid value for lat_1 and lat_2: |lat_1 + lat_2| should be > 0",
            ));
        }
        if sc.1.abs() < ANGULAR_TOLERANCE || phi1.abs() >= FRAC_PI_2 {
            return Err(Error::General(
                "Lcc: Invalid value for lat_1: |lat_1| should be < 90°",
            ));
        }
        if phi2.cos().abs() < ANGULAR_TOLERANCE || phi2.abs() >= FRAC_PI_2 {
            return Err(Error::General(
                "Lcc: Invalid value for lat_2: |lat_2| should be < 90°",
            ));
        }

        let m1 = crate::math::ancillary::pj_msfn(sc, es);
        let ml1 = crate::math::ancillary::ts(sc, e);

        if (phi1 - phi2).abs() >= ANGULAR_TOLERANCE {
            let sc = phi2.sin_cos();
            n = (m1 / crate::math::ancillary::pj_msfn(sc, es)).ln();
            if n == 0.0 {
                return Err(Error::General("Lcc: Invalid value for eccentricity"));
            }
            let ml2 = crate::math::ancillary::ts(sc, e);
            let denom = (ml1 / ml2).ln();
            if denom == 0.0 {
                return Err(Error::General("Lcc: Invalid value for eccentricity"));
            }
            n /= denom;
        }

        let c = m1 * ml1.powf(-n) / n;
        let rho0 = if (lat_0.abs() - FRAC_PI_2).abs() > ANGULAR_TOLERANCE {
            c * crate::math::ancillary::ts(lat_0.sin_cos(), e).powf(n)
        } else {
            0.0
        };

        Ok(Self {
            frame,
            ellps,
            n,
            c,
            rho0,
            west_oriented: params.boolean("west_oriented"),
        })
    }
}

fn fwd(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<LccState>();
    let e = state.ellps.eccentricity();
    let mut successes = 0_usize;
    let length = operands.len();

    for i in 0..length {
        let (mut lam, phi) = operands.xy(i);
        lam = state.frame.remove_central_meridian(lam);
        let mut rho = 0.;

        if (phi.abs() - FRAC_PI_2).abs() < ANGULAR_TOLERANCE {
            if phi * state.n <= 0.0 {
                operands.set_coord(i, &Coor4D::nan());
                continue;
            }
        } else {
            rho = state.c * crate::math::ancillary::ts(phi.sin_cos(), e).powf(state.n);
        }
        let sc = (lam * state.n).sin_cos();
        let x = if state.west_oriented {
            state.frame.x_0 - state.frame.a * state.frame.k_0 * rho * sc.0
        } else {
            state.frame.a * state.frame.k_0 * rho * sc.0 + state.frame.x_0
        };
        let y = state.frame.a * state.frame.k_0 * (state.rho0 - rho * sc.1) + state.frame.y_0;
        operands.set_xy(i, x, y);
        successes += 1;
    }
    successes
}

fn inv(op: &Op, _ctx: &dyn Context, operands: &mut dyn CoordinateSet) -> usize {
    let state = op.state::<LccState>();
    let e = state.ellps.eccentricity();
    let mut successes = 0_usize;

    for i in 0..operands.len() {
        let (mut x, mut y) = operands.xy(i);
        x = if state.west_oriented {
            (state.frame.x_0 - x) / (state.frame.a * state.frame.k_0)
        } else {
            (x - state.frame.x_0) / (state.frame.a * state.frame.k_0)
        };
        y = state.rho0 - (y - state.frame.y_0) / (state.frame.a * state.frame.k_0);

        let mut rho = x.hypot(y);

        if rho == 0. {
            let lon = 0.0;
            let lat = FRAC_PI_2.copysign(state.n);
            operands.set_xy(i, lon, lat);
            successes += 1;
            continue;
        }

        if state.n < 0.0 {
            rho = -rho;
            x = -x;
            y = -y;
        }

        let ts0 = (rho / state.c).powf(1.0 / state.n);
        let lat = crate::math::ancillary::pj_phi2(ts0, e);
        if lat.is_infinite() || lat.is_nan() {
            operands.set_coord(i, &Coor4D::nan());
            continue;
        }
        let lon = x.atan2(y) / state.n + state.frame.lon_0;
        operands.set_xy(i, lon, lat);
        successes += 1;
    }
    successes
}

pub fn new(parameters: &RawParameters, _ctx: &dyn Context) -> Result<Op, Error> {
    let def = &parameters.instantiated_as;
    let params = ParsedParameters::new(parameters, &GAMUT)?;
    if !params.given.contains_key("lat_1") {
        return Err(Error::MissingParam("lat_1".to_string()));
    }

    let phi1 = params.lat(1);
    let phi2 = if params.lat(2).is_nan() {
        phi1
    } else {
        params.lat(2)
    };

    let mut lat_0 = params.lat(0);
    if lat_0.is_nan() {
        lat_0 = 0.0;
        if (phi1 - phi2).abs() < ANGULAR_TOLERANCE {
            lat_0 = phi1;
        }
    }

    let descriptor = OpDescriptor::new(def, InnerOp(fwd), Some(InnerOp(inv)));
    let frame = ProjectionFrame::from_params(&params);
    let ellps = params.ellps(0);
    let state = LccState::new(&params, frame, ellps, phi1, phi2, lat_0)?;
    Ok(Op::with_state(descriptor, params, state))
}

// ----- T E S T S ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn one_standard_parallel() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "lcc lat_1=57 lon_0=12";
        let op = ctx.op(definition)?;

        // Validation values from PROJ:
        //     echo 12 55 0 0 | cct -d18 proj=lcc lat_1=57 lon_0=12  -- | clip
        //     echo 10 55 0 0 | cct -d18 proj=lcc lat_1=57 lon_0=12  -- | clip
        //     echo 14 59 0 0 | cct -d18 proj=lcc lat_1=57 lon_0=12  -- | clip

        let geo = [
            Coor4D::geo(55., 12., 0., 0.),
            Coor4D::geo(55., 10., 0., 0.),
            Coor4D::geo(59., 14., 0., 0.),
        ];

        let projected = [
            Coor4D::raw(-0.000000000101829246, -222_728.122_307_816_05, 0., 0.),
            Coor4D::raw(-128_046.472_438_652_24, -220_853.700_160_506_4, 0., 0.),
            Coor4D::raw(115_005.414_566_200_68, 224_484.514_376_338_9, 0., 0.),
        ];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 2e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }
        Ok(())
    }

    #[test]
    fn two_standard_parallels() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "lcc lat_1=33 lat_2=45 lon_0=10";
        let op = ctx.op(definition)?;

        // Validation value from PROJ:
        // echo 12 40 0 0 | cct -d12 proj=lcc lat_1=33 lat_2=45 lon_0=10 -- | clip
        let geo = [Coor4D::geo(40., 12., 0., 0.)];
        let projected = [Coor4D::raw(
            169_863.026_093_938_3,
            4_735_925.219_292_451,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 9e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }
        Ok(())
    }

    #[test]
    fn one_standard_parallel_and_latitudinal_offset() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "lcc lat_1=39 lat_0=35 lon_0=10";
        let op = ctx.op(definition)?;

        // Validation value from PROJ:
        // echo 12 40 0 0 | cct -d12 proj=lcc lat_1=39 lat_0=35 lon_0=10 -- | clip
        let geo = [Coor4D::geo(40., 12., 0., 0.)];
        let projected = [Coor4D([
            170_800.011_728_740_65,
            557_172.361_112_929_4,
            0.,
            0.,
        ])];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 2e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 2e-9);
        }
        Ok(())
    }

    #[test]
    fn two_standard_parallels_and_latitudinal_offset() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10";
        let op = ctx.op(definition)?;

        // Validation value from PROJ:
        // echo 12 40 0 0 | cct -d12 proj=lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10 -- | clip
        let geo = [Coor4D::geo(40., 12., 0., 0.)];
        let projected = [Coor4D::raw(
            169_863.026_093_938_36,
            554_155.440_793_916_6,
            0.,
            0.,
        )];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 2e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }
        Ok(())
    }

    #[test]
    fn two_sp_lat_offset_xy_offset() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10 x_0=12345 y_0=67890";
        let op = ctx.op(definition)?;

        // Validation value from PROJ:
        // echo 12 40 0 0 | cct -d12 proj=lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10  x_0=12345 y_0=67890 -- | clip
        let geo = [Coor4D::geo(40., 12., 0., 0.)];
        let projected = [Coor4D([
            182_208.026_093_938_3,
            622_045.440_793_916_6,
            0.,
            0.,
        ])];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 2e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }
        Ok(())
    }

    #[test]
    fn two_sp_lat_offset_xy_offset_scaling() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let definition = "lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10 x_0=12345 y_0=67890 k_0=0.99";
        let op = ctx.op(definition)?;

        // Validation value from PROJ:
        // echo 12 40 0 0 | cct -d12 proj=lcc lat_1=33 lat_2=45 lat_0=35 lon_0=10  x_0=12345 y_0=67890 k_0=0.99 -- | clip
        let geo = [Coor4D::geo(40., 12., 0., 0.)];
        let projected = [Coor4D([
            180_509.395_832_998_9,
            616_503.886_385_977_5,
            0.,
            0.,
        ])];

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&projected[i]) < 2e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..operands.len() {
            assert!(operands[i].hypot2(&geo[i]) < 1e-9);
        }
        Ok(())
    }

    #[test]
    fn wraps_longitude_difference_across_dateline() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op(
            "lcc lat_0=51 lon_0=-176 lat_1=53.8333333333333 lat_2=51.8333333333333 x_0=1000000 y_0=0 ellps=GRS80",
        )?;

        let geo = [Coor4D::geo(52., 179., 0., 0.)];
        let projected = [Coor4D::raw(
            656_902.480_027_398_8,
            123_207.313_340_269_51,
            0.,
            0.,
        )];
        let ellps = Ellipsoid::named("GRS80")?;

        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        assert!(operands[0].hypot2(&projected[0]) < 2e-8);

        ctx.apply(op, Inv, &mut operands)?;
        assert!(ellps.distance(&operands[0], &geo[0]) < 1e-8);
        Ok(())
    }

    #[test]
    fn west_oriented_flips_westing_around_false_easting() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let standard = ctx.op("lcc lat_1=57 lat_0=57 lon_0=12 x_0=500000 y_0=0")?;
        let west = ctx.op("lcc west_oriented lat_1=57 lat_0=57 lon_0=12 x_0=500000 y_0=0")?;
        let geo = [Coor4D::geo(55., 10., 0., 0.)];

        let mut standard_xy = geo;
        ctx.apply(standard, Fwd, &mut standard_xy)?;

        let mut west_xy = geo;
        ctx.apply(west, Fwd, &mut west_xy)?;

        assert!((west_xy[0][0] - (1_000_000.0 - standard_xy[0][0])).abs() < 1e-9);
        assert!((west_xy[0][1] - standard_xy[0][1]).abs() < 1e-9);

        ctx.apply(west, Inv, &mut west_xy)?;
        assert!(west_xy[0].hypot2(&geo[0]) < 1e-9);
        Ok(())
    }
}
