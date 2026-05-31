# Problem 1 Solver Status

## Current Solver Boundary

### Route A

Route A is the reliable baseline and fallback solver for Problem 1.

Current Route A pipeline:

1. nearest-node linear Taylor seed
2. target-time q-sheet continuation
3. residual-first Newton refinement

Engineering target:

- `residual_tolerance_seconds = 1e-2`

Current role:

- reliable baseline
- Route A fallback
- correctness reference for Route B comparisons

### Route B

Route B is the fast query candidate path.

Current Route B pipeline:

1. nearest-node projected tangent finite-difference Hessian
2. quadratic prediction
3. target-time q-sheet continuation
4. one-step residual correction

Current default Hessian mode:

- `ProjectedTangentFiniteDifference`

Important boundaries:

- Route B does not call full Newton
- Projected mode does one residual-based projection per stencil
- query-time Route B still may fail with `residual_after_correction_too_large`
- when Route B is invalid, fallback to Route A

### Legacy Tangent Mode

`TangentFiniteDifference` is retained as:

- legacy mode
- diagnostic mode
- ablation mode

It is not the default Route B path.

### Why High-Level Hybrid Is Not Yet Productionized

Current test-only hybrid behavior is:

1. try Route B first
2. if Route B fails, fallback to Route A

This is not yet productionized because:

- online Hessian is not the target path
- online Route B + Route A fallback is slower than pure Route A
- the online projected Hessian build dominates Route B query time

### Engineering Direction

The next Problem 1 engineering direction is:

1. precomputed Hessian cache
2. root-table build-time precomputation
3. query from cached Hessian
4. Route A fallback only when cached Route B is invalid

In other words:

- online Hessian is not the target path
- precomputed Hessian is the intended Route B acceleration path
