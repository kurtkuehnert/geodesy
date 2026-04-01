# Frame Wrapper Notes

Shelved idea for later:

- A plain `frame` pipeline step is the wrong abstraction for projection framing.
- The real behavior is wrapper-shaped:
  - frame start
  - core projection math
  - frame end
- This is needed because central meridian handling and false origin handling may both
  apply, but on different sides of the wrapped operator.

Promising syntax ideas to revisit:

- `merc @frame(lon_0=12 x_0=500000 y_0=0)`
- `noop @frame(lon_0=12)` as a cleaner replacement for PROJ `longlat`

Open design questions:

- How should wrapper syntax be represented in the geodesy string grammar?
- Should wrappers be limited to specific operator domain contracts, or be generally composable?
- Can the existing PROJ degree-boundary wrapping concept use the same mechanism?

Current decision:

- Keep `lonlat` as the single canonical geographic identity operator in geodesy.
- Normalize PROJ aliases (`longlat`, `latlon`, `latlong`) to `lonlat` in the PROJ adapter.
- Do not keep `frame` as a standalone builtin operator for now.

## Reusable Refactor Patterns

These are the general patterns that have worked well for the current operator cleanup.
Future agents can usually follow these first and only deviate with a clear reason.

### Point Ops

- Prefer `PointOp` for genuinely pointwise leaf operators.
- Keep the canonical external name on the operator type via `PointOp::NAME`.
- Inline `GAMUT` on the trait impl with `#[rustfmt::skip]` when the gamut is only used once.
- Use `type State = Self`-style thinking even though the trait no longer has a separate state type:
  the operator struct should hold exactly the runtime state needed by `fwd`/`inv`.
- Keep control-flow/meta operators like `pipeline` and `stack` off the point-op path.

### Framed Projections

- Use `FramedProjection` when a projection follows the common shape:
  - remove central meridian
  - do local projection math
  - apply false origin
- The wrapper owns `inv`, `lon_0`, `x_0`, and `y_0`.
- Projection-local params stay on the inner implementation.
- Use `framed_gamut!(...)` for framed projections, and list only the projection-specific params inside it.
- Use a searchable `*Inner` implementation type plus a type alias for the real operator:
  - `pub(crate) struct CeaInner`
  - `pub(crate) type Cea = Framed<CeaInner>`
- If another module needs the inner implementation, make that `*Inner` type `pub(crate)`.

### Projection Helpers

- Keep cached conversion state separate from call-site math where possible.
- `AuthalicLatitude` is the current model for this:
  - `ellps.authalic()`
  - `geographic_to_authalic`
  - `authalic_to_geographic`
  - `q_from_geographic`
  - `geographic_from_q`
- Prefer directional method names over `from_*` instance methods.

### Testing

- Prefer PROJ-backed golden tests over self-referential roundtrip-only tests.
- For invertible operators, the default test shape is:
  - one PROJ forward golden
  - plus inverse roundtrip in the same helper
- Use `assert_proj_match(...)` as the standard parity helper.
- Keep only the minimal set of cases that covers:
  - spherical vs ellipsoidal paths when both exist
  - major aspect/branch differences
  - one representative singularity or wrapping edge case when relevant
  - constructor/parameter rejection
- Use `assert_op_err(...)` for simple invalid-definition tests.

### Compatibility Surface

- Prefer one canonical geodesy name for an operator.
- Keep spelling aliases and PROJ compatibility names in the PROJ adapter where possible.
- Do not preserve duplicate internal names unless they represent meaningfully different algorithms.
