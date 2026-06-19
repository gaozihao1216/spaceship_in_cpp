# Problem 1 Solver Status

## Current Solver Boundary

The active implementation is the direct Problem 1 residual solver in
`include/spaceship_cpp/problem1/problem1.hpp` and
`src/problem1/problem1.cpp`.

Current direct-solve pipeline:

1. scan encounter angle candidates
2. detect residual sign changes for each revolution pair
3. refine accepted roots by bisection
4. return physically valid transfer candidates sorted by flight time

Current role:

- correctness baseline for Problem 1
- source of transfer-time candidates used by tests and diagnostics
- replacement direction after removing the root-table precomputation path

The remaining `Problem1Table` module is a smaller endpoint geometry and
transfer-time table experiment. It does not own root precomputation.
