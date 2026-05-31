#pragma once

#include "spaceship_cpp/bfs/trajectory_search_state.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct InitialLaunchExpansionOptions {
    int max_transfer_revolution = 1;
    int max_target_revolution = 1;

    double problem1_residual_tolerance_seconds = 1e-2;
    int problem1_max_newton_iterations = 80;

    double max_launch_v_inf = std::numeric_limits<double>::infinity();

    double time_weight_m_per_s_per_day = 0.0;
};

struct InitialLaunchCandidate {
    bool valid = false;

    planet_params::PlanetId launch_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId target_planet = planet_params::PlanetId::Mercury;

    double launch_time = 0.0;
    double arrival_time = 0.0;
    double transfer_time_seconds = 0.0;

    double initial_theta = 0.0;
    double alpha = 0.0;

    double transfer_eccentricity = 0.0;
    double transfer_periapsis_angle = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    double launch_v_inf = 0.0;
    double problem1_residual_seconds = 0.0;
    bool periapsis_flipped_from_input_theta = false;
    double transfer_p_resolution_error = 0.0;

    TrajectorySearchEdge edge;
    TrajectorySearchState next_state;
};

struct InitialLaunchExpansionResult {
    bool ok = false;
    std::string error_message;

    std::vector<InitialLaunchCandidate> candidates;
    std::vector<std::string> target_error_messages;

    int target_planet_count = 0;
    int problem1_query_count = 0;
    int raw_branch_count = 0;
    int accepted_candidate_count = 0;
    int launch_v_inf_pruned_count = 0;
    int periapsis_raw_selected_count = 0;
    int periapsis_flipped_selected_count = 0;

    double min_launch_v_inf = 0.0;
    double max_launch_v_inf = 0.0;
    double max_transfer_p_resolution_error = 0.0;
};

InitialLaunchExpansionResult expand_initial_launch_with_problem1_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    double initial_theta,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const InitialLaunchExpansionOptions& options
);

double compute_arrival_v_inf_from_state(const TrajectorySearchState& state);

}  // namespace spaceship_cpp::bfs
