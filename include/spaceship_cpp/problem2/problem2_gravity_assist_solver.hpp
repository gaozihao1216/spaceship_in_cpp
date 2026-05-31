#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::problem2 {

struct Problem2GravityAssistSolverOptions {
    int theta_sample_count = 64;
    int max_transfer_revolution = 1;
    int max_target_revolution = 1;

    double problem1_residual_tolerance_seconds = 1e-2;
    int problem1_max_newton_iterations = 80;

    double bisection_residual_tolerance = 1e-8;
    double bisection_theta_tolerance = 1e-10;
    int bisection_max_iterations = 80;

    bool allow_roundoff_boundary_relaxed_residual = true;
    double boundary_roundoff_abs_tolerance_m = 1e-3;
    double boundary_roundoff_rel_tolerance = 1e-12;

    // true: high-quality search, recursively subdivides topology-change intervals.
    // false: faster scan, only searches initially stable intervals and may miss roots.
    bool topology_adaptive_enabled = true;
    int topology_max_depth = 10;
    double topology_epsilon = 1e-5;

    bool adaptive_interval_trace_enabled = false;

    bool one_sided_stable_interval_enabled = false;
    int one_sided_validation_sample_count = 1;
    double one_sided_max_residual_tolerance = 1e-8;
    double one_sided_branch_identity_tolerance = 1e-8;

    bool adaptive_midpoint_routeA_from_endpoint_enabled = false;
    int adaptive_midpoint_seed_policy = 0;

    bool use_problem1_solve_for_problem2_seed = false;

    problem1::Problem1FallbackMode problem1_fallback_mode = problem1::Problem1FallbackMode::DirectSolve;
    bool problem1_cell_vertex_routeA_enabled = false;
    int problem1_cell_vertex_seed_policy = 1;
    int problem1_cell_vertex_top_k_vertices = 2;
    int problem1_cell_vertex_max_branch_attempts = 8;
    bool problem1_scout_no_direct_fallback = false;
};

enum class Problem2ResidualSource {
    Strict,
    BoundaryAmbiguousRoundoff
};

struct Problem2GravityAssistSolution {
    bool valid = false;

    double theta_prime = 0.0;
    double alpha = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;

    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;

    double slingshot_residual = 0.0;
    double problem1_residual_seconds = 0.0;

    Problem2ResidualSource residual_source = Problem2ResidualSource::Strict;
    bool boundary_ambiguous = false;

    int bisection_iterations = 0;
    double final_theta_width = 0.0;

    bool origin_was_topology_change = false;
};

struct Problem2GravityAssistSolverSummary {
    int theta_sample_count = 0;

    int initial_stable_interval_count = 0;
    int initial_topology_change_interval_count = 0;

    int adaptive_stable_subinterval_count = 0;
    int adaptive_topology_split_count = 0;
    int topology_transition_core_skipped_count = 0;
    int max_topology_recursion_depth_reached = 0;

    int sign_change_candidate_count = 0;
    int continuation_stable_candidate_count = 0;
    int bisection_attempt_count = 0;
    int bisection_success_count = 0;
    int dedup_success_count = 0;

    int strict_root_count = 0;
    int relaxed_boundary_root_count = 0;

    int table_fallback_count = 0;
    int continuation_failure_count = 0;
    int bisection_failure_count = 0;

    double max_abs_slingshot_residual_at_root = 0.0;
    double max_abs_problem1_residual_seconds_at_root = 0.0;
};

struct Problem2GravityAssistSolverProfile {
    int initial_sample_count = 0;
    int adaptive_midpoint_sample_count = 0;
    int total_sample_build_count = 0;
    int cached_sample_hit_count = 0;

    int problem1_nearest_node_query_count = 0;
    int direct_problem1_seed_solve_count = 0;
    int direct_problem1_seed_branch_count = 0;
    int direct_problem1_seed_derivative_attach_count = 0;
    int continuation_refine_count = 0;
    int derivative_attach_count = 0;
    int bisection_iteration_count = 0;

    int nearest_query_table_loaded_branch_count = 0;
    int nearest_query_seed_branch_count = 0;
    int nearest_query_refine_attempt_count = 0;
    int nearest_query_refine_success_count = 0;
    int nearest_query_refine_fail_count = 0;
    int nearest_query_derivative_attach_attempt_count = 0;
    int nearest_query_derivative_attach_success_count = 0;
    int nearest_query_derivative_attach_fail_count = 0;
    int nearest_query_fallback_direct_solve_count = 0;
    int nearest_query_fallback_direct_branch_count = 0;
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

    int initial_nearest_query_count = 0;
    int initial_nearest_query_refine_attempt_count = 0;

    int adaptive_nearest_query_count = 0;
    int adaptive_nearest_query_refine_attempt_count = 0;

    double initial_sampling_ms = 0.0;
    double topology_adaptive_sampling_ms = 0.0;
    double candidate_collection_ms = 0.0;
    double midpoint_continuation_ms = 0.0;
    double bisection_ms = 0.0;
    double dedup_ms = 0.0;
    double total_ms = 0.0;

    double nearest_query_table_load_ms = 0.0;
    double nearest_query_table_offset_build_ms = 0.0;
    double nearest_query_node_read_ms = 0.0;
    double nearest_query_seed_generation_ms = 0.0;
    double nearest_query_refine_ms = 0.0;
    double nearest_query_derivative_attach_ms = 0.0;
    double nearest_query_dedup_ms = 0.0;
    double nearest_query_fallback_direct_solve_ms = 0.0;
    double cell_vertex_selection_ms = 0.0;
    double cell_vertex_routeA_ms = 0.0;
    double cell_vertex_validation_ms = 0.0;
    double nearest_query_total_ms = 0.0;
    double direct_problem1_seed_solve_ms = 0.0;
    double direct_problem1_seed_derivative_attach_ms = 0.0;

    double initial_nearest_query_table_load_ms = 0.0;
    double initial_nearest_query_refine_ms = 0.0;
    double initial_nearest_query_derivative_attach_ms = 0.0;
    double initial_nearest_query_total_ms = 0.0;

    double adaptive_nearest_query_table_load_ms = 0.0;
    double adaptive_nearest_query_refine_ms = 0.0;
    double adaptive_nearest_query_derivative_attach_ms = 0.0;
    double adaptive_nearest_query_total_ms = 0.0;

    int adaptive_interval_count = 0;
    int adaptive_stable_interval_count = 0;
    int adaptive_topology_change_interval_count = 0;
    int adaptive_boundary_ambiguous_interval_count = 0;
    int adaptive_subdivided_interval_count = 0;

    int adaptive_left_endpoint_query_count = 0;
    int adaptive_right_endpoint_query_count = 0;
    int adaptive_midpoint_query_count = 0;

    int adaptive_cache_hit_count = 0;
    int adaptive_cache_miss_count = 0;

    int one_sided_interval_used_count = 0;
    int one_sided_validation_success_count = 0;
    int one_sided_validation_failure_count = 0;
    int one_sided_fallback_count = 0;

    int adaptive_midpoint_routeA_attempt_count = 0;
    int adaptive_midpoint_routeA_success_count = 0;
    int adaptive_midpoint_routeA_validation_failure_count = 0;
    int adaptive_midpoint_routeA_fallback_count = 0;
    int adaptive_midpoint_nearest_query_avoided_count = 0;

    double adaptive_midpoint_routeA_ms = 0.0;
    double adaptive_midpoint_fallback_nearest_query_ms = 0.0;

    double adaptive_interval_trace_overhead_ms = 0.0;

    struct AdaptiveIntervalTrace {
        int interval_id = 0;
        int depth = 0;

        double theta_left = 0.0;
        double theta_right = 0.0;
        double theta_mid = 0.0;

        bool left_sample_built = false;
        bool right_sample_built = false;
        bool midpoint_sample_built = false;

        bool left_from_cache = false;
        bool right_from_cache = false;
        bool midpoint_from_cache = false;

        int left_branch_count = 0;
        int right_branch_count = 0;
        int midpoint_branch_count = 0;

        bool classified_stable = false;
        bool classified_topology_change = false;
        bool classified_boundary_ambiguous = false;
        bool classified_discontinuous = false;

        bool subdivided = false;
        std::string subdivision_reason;

        int nearest_query_count_for_interval = 0;
        double nearest_query_ms_for_interval = 0.0;
    };

    struct ThetaSampleTrace {
        double theta_prime = 0.0;
        double query_nu_A = 0.0;
        double query_nu_B = 0.0;
        double query_theta_A = 0.0;
        bool adaptive = false;
        bool from_cache = false;
        long long nearest_linear_index = -1;
        int nearest_nu_A_index = -1;
        int nearest_nu_B_index = -1;
        int nearest_theta_A_index = -1;
        int chunk_index = -1;

        double total_ms = 0.0;
        double table_offset_build_ms = 0.0;
        double table_node_read_ms = 0.0;
        double seed_generation_ms = 0.0;
        double refine_ms = 0.0;
        double derivative_attach_ms = 0.0;
        double dedup_ms = 0.0;
        double fallback_direct_solve_ms = 0.0;

        int loaded_branch_count = 0;
        int seed_branch_count = 0;
        int refine_attempt_count = 0;
        int refine_success_count = 0;
        int derivative_attach_attempt_count = 0;
        int fallback_direct_solve_count = 0;
    };

    struct RouteAMidpointTrace {
        double theta_mid = 0.0;
        int seed_policy = 0;
        bool seed_endpoint_left = true;
        int seed_valid_branch_count = 0;
        int routeA_attempt_count_for_midpoint = 0;
        int routeA_success_count_for_midpoint = 0;
        int routeA_validation_failure_count_for_midpoint = 0;
        bool fallback_used = false;
        double routeA_total_ms_for_midpoint = 0.0;
        double validation_ms_for_midpoint = 0.0;
        double fallback_nearest_query_ms_for_midpoint = 0.0;
    };

    std::vector<AdaptiveIntervalTrace> adaptive_interval_traces;
    std::vector<ThetaSampleTrace> theta_sample_traces;
    std::vector<RouteAMidpointTrace> routeA_midpoint_traces;
};

struct Problem2GravityAssistSolverResult {
    bool ok = false;
    std::string error_message;

    std::vector<Problem2GravityAssistSolution> solutions;
    Problem2GravityAssistSolverSummary summary;
    Problem2GravityAssistSolverProfile profile;
};

Problem2GravityAssistSolverResult solve_problem2_gravity_assist_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options
);

}  // namespace spaceship_cpp::problem2
