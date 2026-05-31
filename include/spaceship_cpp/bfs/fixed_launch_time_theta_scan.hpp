#pragma once

#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/bfs/theta_candidate_selector.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

enum class ThetaScanSelectionMode {
    LowestEstimatedVInf,
    StratifiedByTheta,
    Hybrid
};

struct FixedLaunchTimeThetaScanOptions {
    int max_theta_to_scan = 20;

    ThetaScanSelectionMode theta_selection_mode = ThetaScanSelectionMode::Hybrid;
    int hybrid_lowest_vinf_count = 6;
    int hybrid_stratified_count = 14;

    int top_solution_count = 10;
    int top_sequence_count = 10;

    bool deduplicate_nearby_theta = true;
    double duplicate_theta_tolerance = 1e-8;

    bool print_all_theta_summaries = true;
};

struct FixedLaunchTimeThetaScanThetaSummary {
    bool valid = false;

    double theta = 0.0;
    double estimated_min_launch_v_inf = std::numeric_limits<double>::infinity();

    int initial_candidate_count = 0;
    int terminal_solution_count = 0;
    int node_count = 0;

    bool has_terminal = false;
    double best_score = std::numeric_limits<double>::infinity();
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();

    int best_node_index = -1;
    std::vector<planet_params::PlanetId> best_sequence;
};

struct FixedLaunchTimeThetaScanSolution {
    bool valid = false;

    double theta = 0.0;
    double estimated_min_launch_v_inf = std::numeric_limits<double>::infinity();

    FixedLaunchThetaTerminalSolution terminal_solution;
    std::vector<TrajectorySearchEdge> path;
    std::vector<planet_params::PlanetId> sequence;
};

struct FixedLaunchTimeThetaScanSequenceSummary {
    bool valid = false;

    std::vector<planet_params::PlanetId> sequence;

    int occurrence_count = 0;

    double best_theta = 0.0;
    double best_score = std::numeric_limits<double>::infinity();
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();
};

struct FixedLaunchTimeThetaScanResult {
    bool ok = false;
    std::string error_message;

    double launch_time = 0.0;

    std::vector<FixedLaunchTimeThetaScanThetaSummary> theta_summaries;
    std::vector<FixedLaunchTimeThetaScanSolution> top_solutions;
    std::vector<FixedLaunchTimeThetaScanSequenceSummary> top_sequences;

    int theta_candidate_count = 0;
    int theta_scanned_count = 0;
    int terminal_theta_count = 0;
    int total_terminal_solution_count = 0;

    FixedLaunchTimeThetaScanSolution best_solution;
};

const char* theta_scan_selection_mode_name(ThetaScanSelectionMode mode);

std::vector<ThetaCandidate> select_fixed_launch_time_theta_scan_candidates(
    const std::vector<ThetaCandidate>& theta_candidates,
    const FixedLaunchTimeThetaScanOptions& options
);

std::vector<planet_params::PlanetId> sequence_from_path(
    planet_params::PlanetId launch_planet,
    const std::vector<TrajectorySearchEdge>& path
);

FixedLaunchTimeThetaScanResult scan_fixed_launch_time_theta_candidates_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    planet_params::PlanetId terminal_planet,
    double launch_time,
    const std::vector<ThetaCandidate>& theta_candidates,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const std::vector<planet_params::PlanetId>& allowed_transfer_planets,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedLaunchThetaSearchOptions& search_options,
    const FixedLaunchTimeThetaScanOptions& scan_options
);

}  // namespace spaceship_cpp::bfs
