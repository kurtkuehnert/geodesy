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
