//! Geographic to cartesian (and vice versa) conversion.
//!
//! Attribution:
//! - PROJ 9.8.0 `cart.cpp`:
//!   <https://github.com/OSGeo/PROJ/blob/9.8.0/src/conversions/cart.cpp>
//! - geodesy's shared Bowring-style geographic/cartesian helper:
//!   [geocart.rs](crate::ellipsoid::geocart)

use crate::authoring::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Cart {
    ellps: Ellipsoid,
}

impl PointOp for Cart {
    const NAME: &'static str = "cart";
    const TITLE: &'static str = "Geodetic to cartesian conversion";
    #[rustfmt::skip]
    const GAMUT: &'static [OpParameter] = &[
        OpParameter::Flag { key: "inv" },
        OpParameter::Text { key: "ellps", default: Some("GRS80") },
    ];

    fn build(params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
        Ok(Self {
            ellps: params.ellps(0),
        })
    }

    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        let coord = self.ellps.cartesian(&coord);
        (!coord.0.iter().any(|c| c.is_nan())).then_some(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        let coord = self.ellps.geographic(&coord);
        (!coord.0.iter().any(|c| c.is_nan())).then_some(coord)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        let op = ctx.op("cart")?;

        let geo = [
            Coor4D::geo(85., 0., 100000., 0.),
            Coor4D::geo(55., 10., -100000., 0.),
            Coor4D::geo(25., 20., 0., 0.),
            Coor4D::geo(0., -20., 0., 0.),
            Coor4D::geo(-25., 20., 10., 0.),
            Coor4D::geo(-25., -20., 10., 0.),
            Coor4D::geo(25., -20., 10., 0.),
        ];

        let cart = [
            Coor4D::raw(566_462.633_537_476_8, 0.0, 6_432_020.333_690_127, 0.0),
            Coor4D::raw(
                3_554_403.475_871_930_4,
                626_737.233_120_170_7,
                5_119_468.318_659_256,
                0.,
            ),
            Coor4D::raw(
                5_435_195.382_145_216,
                1_978_249.336_521_975_5,
                2_679_074.462_877_277_8,
                0.,
            ),
            Coor4D::raw(5_993_488.273_261_571, -2_181_451.330_890_750_5, 0., 0.),
            Coor4D::raw(
                5_435_203.898_652_612,
                1_978_252.436_277_167_4,
                -2_679_078.689_059_895,
                0.,
            ),
            Coor4D::raw(
                5_435_203.898_652_612,
                -1_978_252.436_277_167_4,
                -2_679_078.689_059_895,
                0.,
            ),
            Coor4D::raw(
                5_435_203.898_652_612,
                -1_978_252.436_277_167_4,
                2_679_078.689_059_895,
                0.,
            ),
        ];

        let e = Ellipsoid::default();
        let mut operands = geo;
        ctx.apply(op, Fwd, &mut operands)?;
        for i in 0..4 {
            assert!(operands[i].hypot3(&cart[i]) < 20e-9);
        }

        ctx.apply(op, Inv, &mut operands)?;
        for i in 0..5 {
            assert!(e.distance(&operands[i], &geo[i]) < 1e-4);
        }

        Ok(())
    }
}
