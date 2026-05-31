#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

enum class Problem1FallbackMode {
    DirectSolve,
    CellVertexRouteA,
    ScoutNoDirectFallback
};

struct Problem1NearestNodeQueryOptions {
    double residual_tolerance_seconds = 1e-2;
    int max_newton_iterations = 80;
    bool fallback_direct_solve = true;
    Problem1FallbackMode fallback_mode = Problem1FallbackMode::DirectSolve;
    bool cell_vertex_routeA_enabled = false;
    int cell_vertex_seed_policy = 1;
    int cell_vertex_top_k_vertices = 2;
    int cell_vertex_max_branch_attempts = 8;
    bool scout_no_direct_fallback = false;
};

struct Problem1NearestNodeQueryProfile {
    int table_load_node_count = 0;
    int table_loaded_branch_count = 0;

    int seed_branch_count = 0;
    int refine_attempt_count = 0;
    int refine_success_count = 0;
    int refine_fail_count = 0;

    int derivative_attach_attempt_count = 0;
    int derivative_attach_success_count = 0;
    int derivative_attach_fail_count = 0;

    int dedup_input_count = 0;
    int dedup_output_count = 0;

    int fallback_direct_solve_count = 0;
    int fallback_direct_branch_count = 0;

    int table_offset_build_count = 0;

    int cell_vertex_routeA_attempt_count = 0;
    int cell_vertex_routeA_success_count = 0;
    int cell_vertex_routeA_failure_count = 0;
    int cell_vertex_routeA_no_valid_vertex_count = 0;
    int cell_vertex_routeA_no_valid_branch_count = 0;
    int cell_vertex_routeA_validation_failure_count = 0;

    int cell_vertex_loaded_vertex_count = 0;
    int cell_vertex_selected_vertex_count = 0;
    int cell_vertex_attempted_branch_count = 0;
    int direct_fallback_avoided_count = 0;
    int direct_fallback_skipped_count = 0;

    double table_load_ms = 0.0;
    double table_offset_build_ms = 0.0;
    double table_node_read_ms = 0.0;
    double seed_generation_ms = 0.0;
    double refine_ms = 0.0;
    double derivative_attach_ms = 0.0;
    double dedup_ms = 0.0;
    double fallback_direct_solve_ms = 0.0;
    double cell_vertex_selection_ms = 0.0;
    double cell_vertex_routeA_ms = 0.0;
    double cell_vertex_validation_ms = 0.0;
    double total_ms = 0.0;
};

struct Problem1NearestNodeQueryResult {
    bool valid = false;
    std::string invalid_reason;

    int nearest_nu_A_index = 0;
    int nearest_nu_B_index = 0;
    int nearest_theta_A_index = 0;
    long long nearest_linear_index = 0;

    double delta_nu_A = 0.0;
    double delta_nu_B = 0.0;
    double delta_theta_A = 0.0;

    std::vector<Problem1SolutionBranch> branches;

    int table_seed_count = 0;
    int refined_success_count = 0;
    int refined_fail_count = 0;
    bool used_direct_solve_fallback = false;

    Problem1NearestNodeQueryProfile profile;
};

Problem1NearestNodeQueryResult query_problem1_from_2deg_nearest_node(
    const Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_transfer_revolution,
    int max_target_revolution,
    const Problem1NearestNodeQueryOptions& options = {}
);

}  // namespace spaceship_cpp::problem1
