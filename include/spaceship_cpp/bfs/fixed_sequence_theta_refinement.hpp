#pragma once

#include "spaceship_cpp/bfs/fixed_launch_time_theta_scan.hpp"
#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct FixedSequenceThetaRefinementOptions {
    double theta_center = 0.0;
    double initial_half_width = 0.05;

    int samples_per_round = 41;
    int refinement_rounds = 3;

    double shrink_factor = 0.35;

    bool include_center = true;

    int beam_width = 10;
    double max_launch_v_inf = 7000.0;
    double time_weight_m_per_s_per_day = 0.0;

    double max_abs_slingshot_residual = 1e-8;
    double max_abs_problem1_residual_seconds = 1e-2;
    double max_arrival_time = std::numeric_limits<double>::infinity();
};

struct FixedSequenceThetaRefinementRoundSummary {
    int round_index = 0;

    double theta_left = 0.0;
    double theta_right = 0.0;
    double theta_center = 0.0;
    double half_width = 0.0;

    int theta_sample_count = 0;
    int terminal_theta_count = 0;

    bool best_valid = false;
    double best_theta = 0.0;
    double best_score = std::numeric_limits<double>::infinity();
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();

    double wall_time_ms = 0.0;
};

struct FixedSequenceThetaRefinementResult {
    bool ok = false;
    std::string error_message;

    std::vector<planet_params::PlanetId> sequence;
    double launch_time = 0.0;

    std::vector<FixedSequenceThetaRefinementRoundSummary> round_summaries;

    bool best_valid = false;
    double best_theta = 0.0;
    double best_score = std::numeric_limits<double>::infinity();
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();

    FixedLaunchTimeThetaScanSolution best_solution;
};

FixedSequenceThetaRefinementResult refine_fixed_sequence_theta_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceThetaRefinementOptions& refinement_options
);

}  // namespace spaceship_cpp::bfs
