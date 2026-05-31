#pragma once

#include "spaceship_cpp/bfs/fixed_launch_time_theta_scan.hpp"
#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"

#include <vector>

namespace spaceship_cpp::bfs {

FixedLaunchTimeThetaScanResult scan_fixed_sequence_over_theta_candidates_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    const std::vector<ThetaCandidate>& theta_candidates,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceThetaSearchOptions& sequence_options,
    const FixedLaunchTimeThetaScanOptions& scan_options
);

}  // namespace spaceship_cpp::bfs
