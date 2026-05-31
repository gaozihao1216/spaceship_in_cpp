#pragma once

#include "spaceship_cpp/bfs/trajectory_search_state.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"

namespace spaceship_cpp::bfs {

TrajectorySearchExpansionResult expand_trajectory_state_by_gravity_assist(
    const problem1::Problem1RootTable2DegLoader& loader,
    const TrajectorySearchState& state,
    planet_params::PlanetId target_planet,
    const problem2::Problem2GravityAssistSolverOptions& options
);

}  // namespace spaceship_cpp::bfs
