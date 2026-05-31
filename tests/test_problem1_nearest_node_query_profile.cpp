#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"

#include <cassert>
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

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem1_nearest_node_query_profile_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto encounter_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.25 * planet_params::planet_orbital_period(encounter_planet);
    const double theta_prime = 1.917405109402;

    const auto departure_state = planet_params::planet_state_at_time(encounter_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double theta_A = spaceship_cpp::common::normalize_angle_0_2pi(
        departure_state.theta_global - theta_prime);

    problem1::Problem1NearestNodeQueryOptions options{};
    options.residual_tolerance_seconds = 1e-2;
    options.max_newton_iterations = 80;
    options.fallback_direct_solve = true;

    const auto result = problem1::query_problem1_from_2deg_nearest_node(
        loader, encounter_planet, target_planet, phi, beta, theta_A, 1, 1, options);

    assert(result.valid);
    assert(!result.branches.empty());
    assert(!result.used_direct_solve_fallback);
    assert(result.profile.fallback_direct_solve_count == 0);

    std::cout << "Problem1NearestNodeQueryProfileSummary\n";
    std::cout << "valid=" << (result.valid ? 1 : 0) << '\n';
    std::cout << "branch_count=" << result.branches.size() << '\n';
    std::cout << "used_direct_solve_fallback=" << (result.used_direct_solve_fallback ? 1 : 0) << '\n';
    std::cout << "table_load_node_count=" << result.profile.table_load_node_count << '\n';
    std::cout << "table_loaded_branch_count=" << result.profile.table_loaded_branch_count << '\n';
    std::cout << "seed_branch_count=" << result.profile.seed_branch_count << '\n';
    std::cout << "refine_attempt_count=" << result.profile.refine_attempt_count << '\n';
    std::cout << "refine_success_count=" << result.profile.refine_success_count << '\n';
    std::cout << "refine_fail_count=" << result.profile.refine_fail_count << '\n';
    std::cout << "derivative_attach_attempt_count=" << result.profile.derivative_attach_attempt_count << '\n';
    std::cout << "derivative_attach_success_count=" << result.profile.derivative_attach_success_count << '\n';
    std::cout << "derivative_attach_fail_count=" << result.profile.derivative_attach_fail_count << '\n';
    std::cout << "dedup_input_count=" << result.profile.dedup_input_count << '\n';
    std::cout << "dedup_output_count=" << result.profile.dedup_output_count << '\n';
    std::cout << "fallback_direct_solve_count=" << result.profile.fallback_direct_solve_count << '\n';
    std::cout << "fallback_direct_branch_count=" << result.profile.fallback_direct_branch_count << '\n';
    std::cout << "table_offset_build_count=" << result.profile.table_offset_build_count << '\n';
    std::cout << "table_load_ms=" << result.profile.table_load_ms << '\n';
    std::cout << "table_offset_build_ms=" << result.profile.table_offset_build_ms << '\n';
    std::cout << "table_node_read_ms=" << result.profile.table_node_read_ms << '\n';
    std::cout << "seed_generation_ms=" << result.profile.seed_generation_ms << '\n';
    std::cout << "refine_ms=" << result.profile.refine_ms << '\n';
    std::cout << "derivative_attach_ms=" << result.profile.derivative_attach_ms << '\n';
    std::cout << "dedup_ms=" << result.profile.dedup_ms << '\n';
    std::cout << "fallback_direct_solve_ms=" << result.profile.fallback_direct_solve_ms << '\n';
    std::cout << "total_ms=" << result.profile.total_ms << '\n';
    std::cout << "profile_ok=1\n";
    return 0;
}
