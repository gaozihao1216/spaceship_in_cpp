#pragma once

#include "spaceship_cpp/bfs/initial_launch_expansion.hpp"
#include "spaceship_cpp/bfs/problem2_branch_profile.hpp"
#include "spaceship_cpp/bfs/trajectory_search_state.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct FixedLaunchThetaSearchOptions {
    int max_depth = 4;
    int beam_width = 20;
    bool enable_beam_pruning = true;

    double max_launch_v_inf = 7000.0;

    double time_weight_m_per_s_per_day = 0.0;

    double max_abs_slingshot_residual = 1e-8;
    double max_abs_problem1_residual_seconds = 1e-2;

    double max_arrival_time = std::numeric_limits<double>::infinity();

    trajectory::FlybyPhysicalFeasibilityOptions flyby_physical_options;

    bool continue_after_reaching_terminal = true;

    int profile_max_node_count = 0;
    double profile_max_wall_time_ms = 0.0;
};

struct FixedLaunchThetaSearchNode {
    bool valid = false;

    TrajectorySearchState state;

    int parent_index = -1;
    TrajectorySearchEdge incoming_edge;

    double partial_score = 0.0;
};

struct FixedLaunchThetaTerminalSolution {
    bool valid = false;

    int node_index = -1;

    double launch_v_inf = 0.0;
    double arrival_v_inf = 0.0;
    double total_delta_v = 0.0;
    double total_flight_time_seconds = 0.0;
    double score = 0.0;

    int depth = 0;
};

struct FixedLaunchThetaSearchLayerSummary {
    int depth = 0;

    int input_state_count = 0;
    int attempted_expansion_count = 0;
    int generated_edge_count = 0;
    int accepted_edge_count = 0;
    int terminal_solution_count_after_layer = 0;
    int output_state_count = 0;

    double layer_wall_time_ms = 0.0;
};

struct OpenSearchDepthWidthStats {
    int depth = 0;
    int input_state_count = 0;
    int attempted_expansion_count = 0;
    int raw_edge_count_total = 0;
    int accepted_edge_count_total = 0;
    int reject_flyby_physical_infeasible_count = 0;
    int zero_raw_edge_solve_count = 0;
    int zero_accepted_edge_solve_count = 0;
    int terminal_solution_count_at_depth = 0;
    int output_state_count_before_beam = 0;
    int output_state_count_after_beam = 0;
    int beam_pruned_count = 0;
    double layer_wall_time_ms = 0.0;
};

struct BfsExpansionAttemptProfile {
    int depth = 0;
    int parent_index = -1;
    planet_params::PlanetId from_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId to_planet = planet_params::PlanetId::Mercury;
    double parent_time = 0.0;
    double parent_incoming_e = 0.0;
    double parent_incoming_theta = 0.0;
    double outer_expansion_wall_ms = 0.0;
    int raw_edge_count = 0;
    int accepted_edge_count = 0;
    int flyby_reject_count = 0;
    bool zero_raw = false;
    bool zero_accepted = false;
    problem2::Problem2GravityAssistSolverProfile problem2_solver_profile;
};

struct FixedLaunchThetaSearchResult {
    bool ok = false;
    std::string error_message;

    double launch_time = 0.0;
    double initial_theta = 0.0;

    std::vector<FixedLaunchThetaSearchNode> nodes;
    std::vector<int> frontier_node_indices;
    std::vector<FixedLaunchThetaTerminalSolution> terminal_solutions;
    std::vector<FixedLaunchThetaSearchLayerSummary> layer_summaries;
    std::vector<OpenSearchDepthWidthStats> depth_width_stats;
    std::vector<BfsExpansionAttemptProfile> expansion_attempt_profiles;
    std::vector<Problem2SolveBranchStats> problem2_solve_branch_stats;
    std::vector<Problem2LayerBranchStats> problem2_layer_branch_stats;

    int initial_candidate_count = 0;
    int launch_v_inf_pruned_count = 0;
    std::vector<std::string> initial_target_error_messages;

    FixedLaunchThetaTerminalSolution best_terminal_solution;

    bool stopped_by_node_limit = false;
    bool stopped_by_wall_time_limit = false;
};

FixedLaunchThetaSearchResult run_fixed_launch_theta_beam_search_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    planet_params::PlanetId terminal_planet,
    double launch_time,
    double initial_theta,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const std::vector<planet_params::PlanetId>& allowed_transfer_planets,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedLaunchThetaSearchOptions& search_options
);

std::vector<TrajectorySearchEdge> reconstruct_fixed_launch_theta_path(
    const FixedLaunchThetaSearchResult& result,
    int node_index
);

}  // namespace spaceship_cpp::bfs
