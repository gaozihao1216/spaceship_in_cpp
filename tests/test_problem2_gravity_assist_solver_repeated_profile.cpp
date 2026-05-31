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
#include <vector>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

double speedup(double cold_ms, double warm_ms) {
    return warm_ms > 0.0 ? cold_ms / warm_ms : 0.0;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_gravity_assist_solver_repeated_profile_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto encounter_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.25 * planet_params::planet_orbital_period(encounter_planet);

    problem2::Problem2GravityAssistSolverOptions options{};
    options.theta_sample_count = 64;
    options.topology_adaptive_enabled = true;

    std::vector<problem2::Problem2GravityAssistSolverResult> results;
    results.reserve(3);
    for (int run_index = 0; run_index < 3; ++run_index) {
        auto result = problem2::solve_problem2_gravity_assist_with_table(
            loader, encounter_planet, target_planet, encounter_time, 0.4, 1.0, options);
        assert(result.ok);
        assert(result.solutions.size() == 16);
        assert(result.summary.max_abs_slingshot_residual_at_root <= 1e-8);

        std::cout << "Problem2RepeatedSolverProfileSummary\n";
        std::cout << "run_index=" << run_index << '\n';
        std::cout << "solution_count=" << result.solutions.size() << '\n';
        std::cout << "total_ms=" << result.profile.total_ms << '\n';
        std::cout << "initial_sampling_ms=" << result.profile.initial_sampling_ms << '\n';
        std::cout << "nearest_query_table_load_ms=" << result.profile.nearest_query_table_load_ms << '\n';
        std::cout << "nearest_query_refine_ms=" << result.profile.nearest_query_refine_ms << '\n';
        std::cout << "nearest_query_derivative_attach_ms="
                  << result.profile.nearest_query_derivative_attach_ms << '\n';
        std::cout << "nearest_query_total_ms=" << result.profile.nearest_query_total_ms << '\n';
        std::cout << "initial_nearest_query_table_load_ms="
                  << result.profile.initial_nearest_query_table_load_ms << '\n';
        std::cout << "initial_nearest_query_refine_ms="
                  << result.profile.initial_nearest_query_refine_ms << '\n';
        std::cout << "adaptive_nearest_query_table_load_ms="
                  << result.profile.adaptive_nearest_query_table_load_ms << '\n';
        std::cout << "adaptive_nearest_query_refine_ms="
                  << result.profile.adaptive_nearest_query_refine_ms << '\n';
        std::cout << "problem1_nearest_node_query_count="
                  << result.profile.problem1_nearest_node_query_count << '\n';
        std::cout << "nearest_query_refine_attempt_count="
                  << result.profile.nearest_query_refine_attempt_count << '\n';
        std::cout << "nearest_query_derivative_attach_attempt_count="
                  << result.profile.nearest_query_derivative_attach_attempt_count << '\n';
        std::cout << "nearest_query_fallback_direct_solve_count="
                  << result.profile.nearest_query_fallback_direct_solve_count << '\n';
        std::cout << "bisection_ms=" << result.profile.bisection_ms << '\n';
        std::cout << "max_abs_slingshot_residual="
                  << result.summary.max_abs_slingshot_residual_at_root << '\n';
        results.push_back(std::move(result));
    }

    const bool same_solution_count =
        results[0].solutions.size() == results[1].solutions.size() &&
        results[0].solutions.size() == results[2].solutions.size();
    assert(same_solution_count);

    std::cout << "Problem2RepeatedSolverWarmupComparison\n";
    std::cout << "cold_total_ms=" << results[0].profile.total_ms << '\n';
    std::cout << "warm1_total_ms=" << results[1].profile.total_ms << '\n';
    std::cout << "warm2_total_ms=" << results[2].profile.total_ms << '\n';
    std::cout << "cold_initial_sampling_ms=" << results[0].profile.initial_sampling_ms << '\n';
    std::cout << "warm1_initial_sampling_ms=" << results[1].profile.initial_sampling_ms << '\n';
    std::cout << "warm2_initial_sampling_ms=" << results[2].profile.initial_sampling_ms << '\n';
    std::cout << "warm_speedup_vs_cold="
              << speedup(results[0].profile.total_ms, results[1].profile.total_ms) << '\n';
    std::cout << "warm2_speedup_vs_cold="
              << speedup(results[0].profile.total_ms, results[2].profile.total_ms) << '\n';
    std::cout << "same_solution_count=" << (same_solution_count ? 1 : 0) << '\n';
    std::cout << "repeated_profile_ok=1\n";
    return 0;
}
