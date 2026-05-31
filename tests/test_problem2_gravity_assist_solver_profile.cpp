#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

int topology_change_origin_count(const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& result) {
    return static_cast<int>(std::count_if(result.solutions.begin(), result.solutions.end(), [](const auto& solution) {
        return solution.origin_was_topology_change;
    }));
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_gravity_assist_solver_profile_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto encounter_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.25 * planet_params::planet_orbital_period(encounter_planet);

    problem2::Problem2GravityAssistSolverOptions options{};
    options.theta_sample_count = 64;
    options.topology_adaptive_enabled = true;

    const auto result = problem2::solve_problem2_gravity_assist_with_table(
        loader, encounter_planet, target_planet, encounter_time, 0.4, 1.0, options);
    assert(result.ok);
    assert(result.solutions.size() == 16);
    assert(result.summary.max_abs_slingshot_residual_at_root <= 1e-8);

    std::cout << "Problem2GravityAssistSolverProfileSummary\n";
    std::cout << "solution_count=" << result.solutions.size() << '\n';
    std::cout << "strict_root_count=" << result.summary.strict_root_count << '\n';
    std::cout << "relaxed_boundary_root_count=" << result.summary.relaxed_boundary_root_count << '\n';
    std::cout << "topology_change_origin_count=" << topology_change_origin_count(result) << '\n';
    std::cout << "initial_sample_count=" << result.profile.initial_sample_count << '\n';
    std::cout << "adaptive_midpoint_sample_count=" << result.profile.adaptive_midpoint_sample_count << '\n';
    std::cout << "total_sample_build_count=" << result.profile.total_sample_build_count << '\n';
    std::cout << "cached_sample_hit_count=" << result.profile.cached_sample_hit_count << '\n';
    std::cout << "problem1_nearest_node_query_count=" << result.profile.problem1_nearest_node_query_count << '\n';
    std::cout << "continuation_refine_count=" << result.profile.continuation_refine_count << '\n';
    std::cout << "derivative_attach_count=" << result.profile.derivative_attach_count << '\n';
    std::cout << "bisection_iteration_count=" << result.profile.bisection_iteration_count << '\n';
    std::cout << "nearest_query_table_loaded_branch_count="
              << result.profile.nearest_query_table_loaded_branch_count << '\n';
    std::cout << "nearest_query_seed_branch_count=" << result.profile.nearest_query_seed_branch_count << '\n';
    std::cout << "nearest_query_refine_attempt_count="
              << result.profile.nearest_query_refine_attempt_count << '\n';
    std::cout << "nearest_query_refine_success_count="
              << result.profile.nearest_query_refine_success_count << '\n';
    std::cout << "nearest_query_refine_fail_count=" << result.profile.nearest_query_refine_fail_count << '\n';
    std::cout << "nearest_query_derivative_attach_attempt_count="
              << result.profile.nearest_query_derivative_attach_attempt_count << '\n';
    std::cout << "nearest_query_derivative_attach_success_count="
              << result.profile.nearest_query_derivative_attach_success_count << '\n';
    std::cout << "nearest_query_derivative_attach_fail_count="
              << result.profile.nearest_query_derivative_attach_fail_count << '\n';
    std::cout << "nearest_query_fallback_direct_solve_count="
              << result.profile.nearest_query_fallback_direct_solve_count << '\n';
    std::cout << "nearest_query_fallback_direct_branch_count="
              << result.profile.nearest_query_fallback_direct_branch_count << '\n';
    std::cout << "initial_sampling_ms=" << result.profile.initial_sampling_ms << '\n';
    std::cout << "topology_adaptive_sampling_ms=" << result.profile.topology_adaptive_sampling_ms << '\n';
    std::cout << "candidate_collection_ms=" << result.profile.candidate_collection_ms << '\n';
    std::cout << "midpoint_continuation_ms=" << result.profile.midpoint_continuation_ms << '\n';
    std::cout << "bisection_ms=" << result.profile.bisection_ms << '\n';
    std::cout << "dedup_ms=" << result.profile.dedup_ms << '\n';
    std::cout << "total_ms=" << result.profile.total_ms << '\n';
    std::cout << "nearest_query_table_load_ms=" << result.profile.nearest_query_table_load_ms << '\n';
    std::cout << "nearest_query_seed_generation_ms=" << result.profile.nearest_query_seed_generation_ms << '\n';
    std::cout << "nearest_query_refine_ms=" << result.profile.nearest_query_refine_ms << '\n';
    std::cout << "nearest_query_derivative_attach_ms="
              << result.profile.nearest_query_derivative_attach_ms << '\n';
    std::cout << "nearest_query_dedup_ms=" << result.profile.nearest_query_dedup_ms << '\n';
    std::cout << "nearest_query_fallback_direct_solve_ms="
              << result.profile.nearest_query_fallback_direct_solve_ms << '\n';
    std::cout << "nearest_query_total_ms=" << result.profile.nearest_query_total_ms << '\n';
    std::cout << "initial_nearest_query_count=" << result.profile.initial_nearest_query_count << '\n';
    std::cout << "initial_nearest_query_table_load_ms="
              << result.profile.initial_nearest_query_table_load_ms << '\n';
    std::cout << "initial_nearest_query_refine_ms=" << result.profile.initial_nearest_query_refine_ms << '\n';
    std::cout << "initial_nearest_query_derivative_attach_ms="
              << result.profile.initial_nearest_query_derivative_attach_ms << '\n';
    std::cout << "initial_nearest_query_total_ms=" << result.profile.initial_nearest_query_total_ms << '\n';
    std::cout << "adaptive_nearest_query_count=" << result.profile.adaptive_nearest_query_count << '\n';
    std::cout << "adaptive_nearest_query_table_load_ms="
              << result.profile.adaptive_nearest_query_table_load_ms << '\n';
    std::cout << "adaptive_nearest_query_refine_ms=" << result.profile.adaptive_nearest_query_refine_ms << '\n';
    std::cout << "adaptive_nearest_query_derivative_attach_ms="
              << result.profile.adaptive_nearest_query_derivative_attach_ms << '\n';
    std::cout << "adaptive_nearest_query_total_ms=" << result.profile.adaptive_nearest_query_total_ms << '\n';
    std::cout << "max_abs_slingshot_residual=" << result.summary.max_abs_slingshot_residual_at_root << '\n';
    std::cout << "profile_test_ok=1\n";
    return 0;
}
