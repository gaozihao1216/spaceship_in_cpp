# Route A Final Design

## Goal

Route A is the nearest-node query baseline for `Problem1RootTable`:

1. choose the nearest table node
2. attach first derivatives on the node branch
3. linearly predict `alpha`
4. refine with residual-first Newton
5. dedup and sort candidates by physical branch semantics

This is the intended fast baseline before any future fallback strategy.

## Newton Policy

Route A uses residual-first Newton refinement.

Engineering target:

- `residual_tolerance_seconds = 1e-2`

Rationale:

- Problem1 flight times are typically `1e7` to `1e8` seconds
- microsecond-level residuals are unnecessary for downstream trajectory search
- second-to-centisecond scale residuals are already much tighter than the practical search resolution

Residual-first Newton behavior:

1. evaluate residual first
2. if `abs(residual_seconds) <= residual_tolerance_seconds`, accept immediately
3. only if residual is not yet acceptable, evaluate derivatives and take a Newton step
4. if residual increases, stop early as a safety guard
5. `max_iterations` is only a safety cap
6. `alpha_step_tolerance` is diagnostic, not the main success criterion

## Branch Identity

`q` is **not** branch identity.

`q` is the target planet time-sheet label, and it can change along the same physical branch as parameters move.

Route A branch identity and dedup should therefore be based on:

- same transfer revolution `k`
- close `time_of_flight_seconds`
- small wrapped alpha distance

That means:

- do not use `q` as a required dedup key
- do not use `q` as the main matching key

## Target-Time Q-Sheet Continuation

Because `q` is a sheet label, Route A should choose `q` by **target-time continuity**, not by residual minimization.

For a source branch:

1. compute linear predicted `alpha`
2. define `source_time_reference`:
   - use `source_branch.target_time_seconds` if finite
   - otherwise use `source_branch.time_of_flight_seconds`
3. for each legal `q`:
   - evaluate `evaluate_problem1_root_residual(..., alpha_linear, k, q)`
   - if valid, compute
     `continuity_error = abs(target_time_seconds(q) - source_time_reference)`
4. choose the `q` with the smallest continuity error

Important:

- do **not** choose `q` by smallest residual
- do **not** use greedy `q` rescue
- do **not** use global alpha winding as the primary `q` correction mechanism

The correct interpretation is:

- `alpha` determines geometry
- `q` determines the target-time sheet
- sheet selection should follow target-time continuity

## Route A Source Pipeline

For each source branch from the nearest node:

1. attach derivatives at the node
2. linearly predict
   `alpha_linear = normalize(source_alpha + d_alpha/dnu_A * dnu_A + d_alpha/dnu_B * dnu_B + d_alpha/dtheta_A * dtheta_A)`
3. select `q` by target-time sheet continuity
4. refine with residual-first Newton using that `(k, q_selected, alpha_linear)`
5. if valid, push candidate
6. global dedup by:
   - same `k`
   - close `time_of_flight_seconds`
   - small wrapped alpha distance
7. group by `k`
8. sort by `time_of_flight_seconds`

## Count Mismatch Handling

If candidate count and expected branch count differ inside a `k` group:

- if there are extra candidates:
  keep the physically reasonable sequence after sorting and conservatively discard extras
- if there are missing candidates:
  mark the group as incomplete / non-admissible

Statistics should still evaluate the first `min(exact_count, candidate_count)` per-k pairs.
Count mismatch should be recorded, but it should not suppress already valid matched pairs.

## Current Baseline

Current validated engineering baseline on the virtual `2°` nearest-node experiment:

- `near_node = 90 / 90 = 1.0`
- `mid_cell = 90 / 95 ≈ 0.947368`
- `physical_launch = 85 / 85 = 1.0`

These numbers are for:

- residual-first Newton
- `residual_tolerance_seconds = 1e-2`
- target-time `q` sheet continuation
- per-k time-order engineering matching

Engineering match criterion:

- `residual_seconds <= 1e-2`
- `time_error <= 1 second`
- `wrapped_alpha_error <= 1e-3 rad`

## Remaining Failures

The remaining failures are not a reason to further harden Route A itself.

Most of them are better interpreted as:

- branch count changes inside the cell
- non-admissible interpolation cells
- branch creation / annihilation regions

These should be handled as:

- incomplete / non-admissible cells
- exact solve fallback
- or future topology-aware continuation logic

They should **not** be force-fixed inside the nearest-node Route A baseline.

