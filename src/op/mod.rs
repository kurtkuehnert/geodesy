mod op_descriptor;
mod parameter;
mod parsed_parameters;
mod raw_parameters;

use crate::authoring::*;
use std::any::Any;
use std::collections::BTreeMap;
use std::fmt::{self, Debug};

pub use op_descriptor::OpDescriptor;
pub use parameter::OpParameter;
pub use parsed_parameters::ParsedParameters;
pub use raw_parameters::RawParameters;

pub trait PointOp: Send + Sync + 'static {
    const NAME: &'static str;
    const TITLE: &'static str;
    const GAMUT: &'static [OpParameter];

    fn build(params: &ParsedParameters, ctx: &dyn Context) -> Result<Self, Error>
    where
        Self: Sized;
    fn fwd(&self, coord: Coor4D) -> Option<Coor4D>;
    fn inv(&self, coord: Coor4D) -> Option<Coor4D>;
}

trait DynPointOp: Send + Sync {
    fn fwd(&self, coord: Coor4D) -> Option<Coor4D>;
    fn inv(&self, coord: Coor4D) -> Option<Coor4D>;
}

struct PointOpInstance<T: PointOp>(T);

impl<T: PointOp> DynPointOp for PointOpInstance<T> {
    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.fwd(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.inv(coord)
    }
}

struct PointRuntime(Box<dyn DynPointOp>);

impl PointRuntime {
    fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.fwd(coord)
    }

    fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
        self.0.inv(coord)
    }
}

/// The defining parameters and functions for an operator
pub struct Op {
    pub descriptor: OpDescriptor,
    pub params: ParsedParameters,
    pub state: Option<Box<dyn Any + Send + Sync>>,
    pub steps: Option<Vec<Op>>,
}

impl Debug for Op {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Op")
            .field("descriptor", &self.descriptor)
            .field("params", &self.params)
            .field("has_state", &self.state.is_some())
            .field("has_point_runtime", &self.point_runtime().is_some())
            .field("steps", &self.steps)
            .finish()
    }
}

impl Op {
    fn point_runtime(&self) -> Option<&PointRuntime> {
        self.state
            .as_ref()
            .and_then(|state| state.downcast_ref::<PointRuntime>())
    }

    // operate fwd/inv, taking operator inversion into account.
    pub fn apply(
        &self,
        ctx: &dyn Context,
        operands: &mut dyn CoordinateSet,
        direction: Direction,
    ) -> usize {
        let going_forward = direction == Direction::Fwd;

        if let Some(point_op) = self.point_runtime() {
            let use_fwd = self.descriptor.inverted != going_forward;
            let mut successes = 0_usize;
            for i in 0..operands.len() {
                let coord = operands.get_coord(i);
                let output = if use_fwd {
                    point_op.fwd(coord)
                } else {
                    point_op.inv(coord)
                };
                if let Some(output) = output {
                    operands.set_coord(i, &output);
                    successes += 1;
                } else {
                    operands.set_coord(i, &Coor4D::nan());
                }
            }
            return successes;
        }

        // We use the .fwd method if we're either not inverted and going forward
        // or inverted and not going forward
        if self.descriptor.inverted != going_forward {
            return self.descriptor.fwd.0(self, ctx, operands);
        }

        // Otherwise, we use the .inv method, if it exists
        if let Some(inv) = &self.descriptor.inv {
            return inv.0(self, ctx, operands);
        }

        // If it doesn't exist, we do nothing, and tell it by stomping on the
        // operands, and reporting zero successes
        operands.stomp();
        0
    }

    pub fn new(definition: &str, ctx: &dyn Context) -> Result<Op, Error> {
        let globals = ctx.globals();
        let definition = definition
            .remove_comments()
            .normalize()
            .handle_prefix_modifiers();
        let parameters = RawParameters::new(&definition, &globals);
        Self::op(parameters, ctx)
    }

    // Helper for implementation of `InnerOp`s: Instantiate an `Op` for the simple
    // (and common) case, where the `InnerOp` constructor does not need to set any
    // other parameters than the ones defined by the instantiation parameter
    // arguments.
    pub fn basic(
        parameters: &RawParameters,
        fwd: InnerOp,
        inv: Option<InnerOp>,
        gamut: &[OpParameter],
    ) -> Result<Op, Error> {
        let def = parameters.instantiated_as.as_str();
        let params = ParsedParameters::new(parameters, gamut)?;
        let descriptor = OpDescriptor::new(def, fwd, inv);
        Ok(Op {
            descriptor,
            params,
            state: None,
            steps: None,
        })
    }

    pub fn with_state<T: Any + Send + Sync>(
        descriptor: OpDescriptor,
        params: ParsedParameters,
        state: T,
    ) -> Op {
        Op {
            descriptor,
            params,
            state: Some(Box::new(state)),
            steps: None,
        }
    }

    pub fn point<T: PointOp + 'static>(
        parameters: &RawParameters,
        ctx: &dyn Context,
    ) -> Result<Op, Error> {
        let def = parameters.instantiated_as.as_str();
        let params = ParsedParameters::new(parameters, T::GAMUT)?;
        let state = T::build(&params, ctx)?;
        let descriptor = OpDescriptor::new(def, InnerOp::default(), Some(InnerOp::default()));
        Ok(Op {
            descriptor,
            params,
            state: Some(Box::new(PointRuntime(Box::new(PointOpInstance::<T>(
                state,
            ))))),
            steps: None,
        })
    }

    pub fn state<T: Any + Send + Sync>(&self) -> &T {
        if self.state.is_none() {
            panic!(
                "operator '{}' is missing expected state {}",
                self.descriptor.instantiated_as,
                std::any::type_name::<T>()
            );
        }

        self.state
            .as_ref()
            .and_then(|state| state.downcast_ref::<T>())
            .unwrap_or_else(|| {
                panic!(
                    "operator '{}' carries state, but not expected type {}",
                    self.descriptor.instantiated_as,
                    std::any::type_name::<T>()
                )
            })
    }

    // Instantiate the actual operator, taking into account the relative order
    // of precendence between pipelines, user defined operators, macros, and
    // built-in operators
    #[allow(clippy::self_named_constructors)]
    pub fn op(parameters: RawParameters, ctx: &dyn Context) -> Result<Op, Error> {
        if parameters.recursion_level > 100 {
            return Err(Error::Recursion(
                parameters.invoked_as,
                parameters.instantiated_as,
            ));
        }

        let name = parameters.instantiated_as.operator_name();

        // A pipeline?
        if parameters.instantiated_as.is_pipeline() {
            return super::inner_op::pipeline::new(&parameters, ctx);
        }

        // A user defined operator?
        if !name.is_resource_name() {
            if let Ok(constructor) = ctx.get_op(&name) {
                return constructor.0(&parameters, ctx)?.handle_op_inversion();
            }
        }
        // If the step is a macro (i.e. potentially an embedded pipeline),
        // we take a copy of the arguments from the macro invocation and enter
        // them into the globals, then enter the expanded macro into the
        // definition field of the paramters, and recursively call Op::op().
        else if let Ok(macro_definition) = ctx.get_resource(&name) {
            let macro_definition = macro_definition
                .remove_comments()
                .normalize()
                .handle_prefix_modifiers();
            // Is the macro called in inverse mode? Search for whitespace-delimited
            // "inv" in order to avoid matching INVariant, subINVolution, etc.
            let inverted = parameters.instantiated_as.contains(" inv ")
                || parameters.instantiated_as.ends_with(" inv");

            let mut globals = parameters.globals.clone();
            globals.remove("_name");
            globals.extend(parameters.instantiated_as.split_into_parameters());

            // Inversion of macros is handled by the `handle_inversion()` method:
            // We do not yet know whether the macro is really a pipeline. So at
            // this level, we remove the `inv` from the globals, to avoid poisoning
            // the pipeline at the single step level
            globals.remove("inv");

            let parameters = RawParameters {
                invoked_as: parameters.invoked_as.clone(),
                instantiated_as: macro_definition,
                globals,
                recursion_level: parameters.recursion_level + 1,
            };

            return Op::op(parameters, ctx)?.handle_inversion(inverted);
        }

        // A built in operator?
        if let Ok(constructor) = super::inner_op::builtin(&name) {
            return constructor.0(&parameters, ctx)?.handle_op_inversion();
        }

        Err(Error::NotFound(
            name,
            ": ".to_string() + &parameters.instantiated_as,
        ))
    }

    fn is_invertible(&self) -> bool {
        self.descriptor.inv.is_some()
    }

    fn handle_op_inversion(self) -> Result<Op, Error> {
        let inverted = self.params.boolean("inv");
        self.handle_inversion(inverted)
    }

    fn handle_inversion(mut self, inverted: bool) -> Result<Op, Error> {
        if self.is_invertible() {
            if inverted {
                self.descriptor.inverted = !self.descriptor.inverted;
            }
            return Ok(self);
        }

        if inverted {
            return Err(Error::NonInvertible(self.descriptor.instantiated_as));
        }

        Ok(self)
    }

    pub fn is_pipeline(&self) -> bool {
        self.steps.is_some()
    }

    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        if self.is_pipeline() {
            self.steps.as_ref().unwrap().len()
        } else {
            1
        }
    }
}

// ----- T E S T S ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    struct AddOnePointOp;

    impl PointOp for AddOnePointOp {
        const NAME: &'static str = "addone-point-op-test";
        const TITLE: &'static str = "Add One Point Op Test";
        const GAMUT: &'static [OpParameter] = &[];

        fn build(_params: &ParsedParameters, _ctx: &dyn Context) -> Result<Self, Error> {
            Ok(Self)
        }

        fn fwd(&self, coord: Coor4D) -> Option<Coor4D> {
            Some(Coor4D([coord[0] + 1.0, coord[1], coord[2], coord[3]]))
        }

        fn inv(&self, coord: Coor4D) -> Option<Coor4D> {
            Some(Coor4D([coord[0] - 1.0, coord[1], coord[2], coord[3]]))
        }
    }

    // Test the fundamental Op-functionality: That we can actually instantiate
    // an Op, and invoke its forward and backward operational modes
    #[test]
    fn basic() -> Result<(), Error> {
        let mut ctx = Minimal::default();

        // Try to invoke garbage as a user defined Op
        assert!(matches!(Op::new("_foo", &ctx), Err(Error::NotFound(_, _))));

        // Check forward and inverse operation
        let op = ctx.op("addone")?;
        let mut data = crate::test_data::coor2d();
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 56.);
        assert_eq!(data[1][0], 60.);
        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // Also for an inverted operator: check forward and inverse operation
        let op = ctx.op("addone inv ")?;
        let mut data = crate::test_data::coor2d();
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 54.);
        assert_eq!(data[1][0], 58.);
        // Corner case: " inv " vs. "inv"
        let op = ctx.op("addone inv")?;
        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        Ok(())
    }

    #[test]
    fn point_runtime() -> Result<(), Error> {
        let ctx = Minimal::default();
        let raw = RawParameters::new("addone_point", &ctx.globals());
        let op = Op::point::<AddOnePointOp>(&raw, &ctx)?;

        let mut data = crate::test_data::coor2d();
        assert_eq!(op.apply(&ctx, &mut data, Fwd), 2);
        assert_eq!(data[0][0], 56.);
        assert_eq!(data[1][0], 60.);

        assert_eq!(op.apply(&ctx, &mut data, Inv), 2);
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);
        Ok(())
    }

    // Test that the recursion-breaker works properly, by defining two mutually
    // dependent macros: `foo:bar=foo:baz` and `foo:baz=foo:bar`, and checking
    // that the instantiation fails with an `Error::Recursion(...)`
    #[test]
    fn nesting() -> Result<(), Error> {
        let mut ctx = Minimal::default();
        ctx.register_resource("foo:baz", "foo:bar");
        ctx.register_resource("foo:bar", "foo:baz");

        assert_eq!("foo:baz", ctx.get_resource("foo:bar")?);
        assert_eq!("foo:bar", ctx.get_resource("foo:baz")?);

        assert!(matches!(ctx.op("foo:baz"), Err(Error::Recursion(_, _))));
        Ok(())
    }

    #[test]
    fn pipeline() -> Result<(), Error> {
        let mut data = crate::test_data::coor2d();
        let mut ctx = Minimal::default();
        let op = ctx.op("addone|addone|addone")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 58.);
        assert_eq!(data[1][0], 62.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        Ok(())
    }

    #[test]
    fn macro_expansion() -> Result<(), Error> {
        let mut data = crate::test_data::coor2d();
        let mut ctx = Minimal::default();
        ctx.register_resource("sub:one", "addone inv");
        let op = ctx.op("addone|sub:one|addone")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 56.);
        assert_eq!(data[1][0], 60.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        Ok(())
    }

    #[test]
    fn macro_expansion_inverted() -> Result<(), Error> {
        let mut data = crate::test_data::coor2d();
        let mut ctx = Minimal::default();
        ctx.register_resource("sub:one", "addone inv");
        let op = ctx.op("addone|sub:one inv|addone")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 58.);
        assert_eq!(data[1][0], 62.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        Ok(())
    }

    #[test]
    fn macro_expansion_with_embedded_pipeline() -> Result<(), Error> {
        let mut data = crate::test_data::coor2d();
        let mut ctx = Minimal::default();
        ctx.register_resource("sub:three", "addone inv|addone inv|addone inv");
        let op = ctx.op("addone|sub:three")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 53.);
        assert_eq!(data[1][0], 57.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        let op = ctx.op("addone|sub:three inv")?;
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 59.);
        assert_eq!(data[1][0], 63.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        Ok(())
    }

    #[test]
    fn macro_expansion_with_defaults() -> Result<(), Error> {
        let mut data = crate::test_data::coor2d();
        let mut ctx = Minimal::default();

        // A macro providing a default value of 1 for the x parameter
        ctx.register_resource("helmert:one", "helmert x=(1)");

        // Instantiating the macro without parameters - getting the default
        let op = ctx.op("helmert:one")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 56.);
        assert_eq!(data[1][0], 60.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // Instantiating the macro with parameters - overwriting the default
        // For good measure, check that it also works inside a pipeline
        let op = ctx.op("addone|helmert:one x=2")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 58.);
        assert_eq!(data[1][0], 62.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // Overwrite the default, and provide additional args
        let op = ctx.op("helmert:one x=2 inv")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 53.);
        assert_eq!(data[1][0], 57.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // Overwrite the default, and provide additional args in a pipeline
        let op = ctx.op("addone|helmert:one inv x=2")?;

        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 54.);
        assert_eq!(data[1][0], 58.);

        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // A macro providing a default value of 1 for the x parameter, unless
        // a macro-parameter called eggs is given
        ctx.register_resource("helmert:won", "helmert x=$eggs(1)");

        // Instantiating the macro without parameters - getting the default
        let op = ctx.op("helmert:won")?;
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 56.);
        assert_eq!(data[1][0], 60.);
        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // Instantiating the macro with eggs = 2
        let op = ctx.op("helmert:won eggs=2")?;
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0][0], 57.);
        assert_eq!(data[1][0], 61.);
        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0][0], 55.);
        assert_eq!(data[1][0], 59.);

        // A macro taking an argument, ham, without any default provided
        ctx.register_resource("helmert:ham", "helmert x=$ham");

        // Instantiating the macro without arguments - getting an error
        assert!(matches!(ctx.op("helmert:ham"), Err(Error::Syntax(_))));

        // Now instantiating the macro with ham = 2
        let op = ctx.op("helmert:ham ham=2")?;
        ctx.apply(op, Fwd, &mut data)?;
        assert_eq!(data[0].x(), 57.);
        assert_eq!(data[1].x(), 61.);
        ctx.apply(op, Inv, &mut data)?;
        assert_eq!(data[0].x(), 55.);
        assert_eq!(data[1].x(), 59.);

        Ok(())
    }

    #[test]
    fn steps() -> Result<(), Error> {
        let steps = "  |\n#\n | |foo bar = baz |   bonk : bonk  $ bonk ||| ".split_into_steps();
        assert_eq!(steps.len(), 2);
        assert_eq!(steps[0], "foo bar=baz");
        assert_eq!(steps[1], "bonk:bonk $bonk");
        let steps = "\n\r\r\n    ||| | \n\n\r\n\r  |  \n\r\r \n  ".split_into_steps();
        assert_eq!(steps.len(), 0);

        Ok(())
    }
}
