#pragma once

#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"

#include <limits>
#include <vector>

namespace spaceship_cpp::bfs {

struct FixedSequenceThetaSearchOptions {
    int beam_width = 20;

    double max_launch_v_inf = 7000.0;

    double time_weight_m_per_s_per_day = 0.0;

    double max_abs_slingshot_residual = 1e-8;
    double max_abs_problem1_residual_seconds = 1e-2;

    double max_arrival_time = std::numeric_limits<double>::infinity();

    trajectory::FlybyPhysicalFeasibilityOptions flyby_physical_options;
};

FixedLaunchThetaSearchResult run_fixed_sequence_theta_search_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    double initial_theta,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceThetaSearchOptions& sequence_options
);

}  // namespace spaceship_cpp::bfs
