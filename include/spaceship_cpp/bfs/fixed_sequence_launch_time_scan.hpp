#pragma once

#include "spaceship_cpp/bfs/fixed_sequence_theta_refinement.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct FixedSequenceLaunchTimeScanOptions {
    double launch_time_center = 0.0;
    double initial_half_width_seconds = 0.0;

    int launch_time_sample_count = 11;

    bool refine_theta_at_each_time = true;

    double theta_center = 2.579506094977;
    double theta_half_width = 0.03;
    int theta_samples_per_round = 21;
    int theta_refinement_rounds = 2;
    double theta_shrink_factor = 0.4;

    int beam_width = 10;
    double max_launch_v_inf = 7000.0;
    double time_weight_m_per_s_per_day = 0.0;
};

struct FixedSequenceLaunchTimeScanSampleSummary {
    bool valid = false;

    double launch_time = 0.0;
    double launch_time_offset_seconds = 0.0;

    bool best_valid = false;
    double best_theta = 0.0;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();

    double wall_time_ms = 0.0;
};

struct FixedSequenceLaunchTimeScanResult {
    bool ok = false;
    std::string error_message;

    std::vector<planet_params::PlanetId> sequence;

    std::vector<FixedSequenceLaunchTimeScanSampleSummary> samples;

    bool best_valid = false;
    double best_launch_time = 0.0;
    double best_theta = 0.0;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();

    FixedSequenceThetaRefinementResult best_refinement_result;
};

FixedSequenceLaunchTimeScanResult scan_fixed_sequence_launch_time_locally_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceLaunchTimeScanOptions& options
);

}  // namespace spaceship_cpp::bfs
