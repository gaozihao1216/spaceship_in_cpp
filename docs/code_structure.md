# Code Structure

This document describes the current source layout after removing the root-table precomputation path.

## Top Level

- `CMakeLists.txt`: Builds the `spaceship_cpp` static library, the two diagnostics apps, and the test suite.
- `README.md`: Short project overview and build commands.
- `apps/problem1_solve_diagnostics.cpp`: Command-line diagnostics for the direct Problem 1 solver.
- `apps/problem1_table_diagnostics.cpp`: Command-line diagnostics for the remaining Problem 1 endpoint table.
- `scripts/build.sh`: Project build helper.
- `scripts/run_tests.sh`: Test runner helper.

## Public Headers

- `include/spaceship_cpp/common/common.hpp`: Shared constants, angle normalization, numeric helpers, and orbit anomaly helpers.
- `include/spaceship_cpp/common/orbit_math.hpp`: Orbital math helpers used by Problem 1 and trajectory code.
- `include/spaceship_cpp/common/time_utils.hpp`: Time conversion helpers.
- `include/spaceship_cpp/config/global_config.hpp`: Global configuration accessors.
- `include/spaceship_cpp/core_types/core_types.hpp`: Shared core type definitions.
- `include/spaceship_cpp/planet_params/planet_params.hpp`: Planet identifiers, physical parameters, and state-at-time queries.
- `include/spaceship_cpp/problem1/orbit_time_integral.hpp`: Inline orbital time integral helpers for Problem 1.
- `include/spaceship_cpp/problem1/problem1.hpp`: Direct Problem 1 residual evaluation and bisection solver API.
- `include/spaceship_cpp/problem1/problem1_diagnostics.hpp`: Diagnostic result types and helpers for Problem 1.
- `include/spaceship_cpp/problem1/problem1_table.hpp`: Endpoint geometry and transfer-time table API retained for experiments.
- `include/spaceship_cpp/problem2/problem2.hpp`: Umbrella header for the remaining Problem 2 API.
- `include/spaceship_cpp/problem2/problem2_slingshot.hpp`: Problem 2 flyby/slingshot residual equations.
- `include/spaceship_cpp/trajectory/flyby_physics.hpp`: Physical flyby feasibility calculations.
- `include/spaceship_cpp/trajectory/orbit_velocity.hpp`: Orbit velocity and planet-relative speed helpers.
- `include/spaceship_cpp/bfs/bfs.hpp`: Placeholder BFS umbrella header.
- `include/spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp`: Conversion between global periapsis angle and Problem 2 local angle.
- `include/spaceship_cpp/bfs/trajectory_search_state.hpp`: Basic trajectory search state and edge data containers.

## Source Layout

- `src/common/common.cpp`: Implementation file for shared common utilities.
- `src/common/orbit_math.cpp`: Implementation file for orbit math utilities.
- `src/common/time_utils.cpp`: Implementation file for time conversion utilities.
- `src/config/global_config.cpp`: Global config implementation.
- `src/core_types/core_types.cpp`: Core type implementation placeholder.
- `src/planet_params/planet_params.cpp`: Planet parameter tables and time-dependent planet state implementation.
- `src/problem1/problem1.cpp`: Direct Problem 1 residual and solve implementation.
- `src/problem1/diagnostics/problem1_diagnostics.cpp`: Problem 1 diagnostics implementation.
- `src/problem1/tables/problem1_table.cpp`: Endpoint geometry and transfer-time table implementation.
- `src/problem2/problem2.cpp`: Problem 2 umbrella implementation.
- `src/problem2/flyby/problem2_slingshot.cpp`: Problem 2 slingshot residual implementation.
- `src/trajectory/flyby_physics.cpp`: Flyby physical feasibility implementation.
- `src/trajectory/orbit_velocity.cpp`: Orbit velocity helper implementation.
- `src/bfs/bfs.cpp`: BFS umbrella implementation placeholder.
- `src/bfs/adapters/problem2_angle_frame_adapter.cpp`: Angle-frame adapter implementation.

## Tests

- `tests/CMakeLists.txt`: Registers only tests that match the current non-root-table code path.
- `tests/test_build.cpp`: Basic build/link smoke test.
- `tests/test_global_config.cpp`: Global config checks.
- `tests/test_orbit_time_integral.cpp`: Orbit time integral checks.
- `tests/test_orbit_velocity_helpers.cpp`: Orbit velocity helper checks.
- `tests/test_planet_params.cpp`: Planet parameter checks.
- `tests/test_planet_position.cpp`: Planet state-at-time checks.
- `tests/test_problem1_residual.cpp`: Problem 1 residual checks.
- `tests/test_problem1_solve.cpp`: Direct Problem 1 solver checks.
- `tests/test_problem1_table.cpp`: Endpoint table checks.
- `tests/test_problem1_table_diagnostics.cpp`: Problem 1 table diagnostics checks.
- `tests/test_problem1_table_theory.cpp`: Problem 1 table theory checks.
- `tests/test_problem2_slingshot_equations.cpp`: Problem 2 slingshot equation checks.
- `tests/test_time_utils.cpp`: Time utility checks.

## Removed Root-Table Path

The root-table generator, 2-degree table loader, nearest-node query layer, and BFS/Problem2 table-backed search modules have been removed from the active code tree and CMake targets. Current Problem 1 work should use the direct solver or the remaining endpoint transfer-time table.
