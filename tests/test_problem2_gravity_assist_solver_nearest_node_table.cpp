#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <tuple>
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

double wrapped_angle_distance(double lhs, double rhs) {
    return std::abs(spaceship_cpp::common::normalize_angle_minus_pi_pi(lhs - rhs));
}

void sort_solutions(std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>* solutions) {
    std::sort(solutions->begin(), solutions->end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.transfer_revolution, lhs.theta_prime, lhs.alpha, lhs.time_of_flight_seconds) <
               std::tie(rhs.transfer_revolution, rhs.theta_prime, rhs.alpha, rhs.time_of_flight_seconds);
    });
}

bool same_root(
    const spaceship_cpp::problem2::Problem2GravityAssistSolution& lhs,
    const spaceship_cpp::problem2::Problem2GravityAssistSolution& rhs
) {
    return lhs.transfer_revolution == rhs.transfer_revolution &&
           lhs.target_revolution == rhs.target_revolution &&
           std::abs(lhs.theta_prime - rhs.theta_prime) <= 1e-6 &&
           wrapped_angle_distance(lhs.alpha, rhs.alpha) <= 1e-6 &&
           std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds) <= 1e3;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_gravity_assist_solver_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
    const double incoming_e = 0.3;
    const double incoming_theta = 0.4;

    problem2::Problem2GravityAssistSolverOptions fast_options{};
    fast_options.theta_sample_count = 64;
    fast_options.topology_adaptive_enabled = false;

    const auto fast_start = std::chrono::steady_clock::now();
    const auto fast = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta, fast_options);
    const auto fast_end = std::chrono::steady_clock::now();
    const double fast_wall_time_ms =
        std::chrono::duration<double, std::milli>(fast_end - fast_start).count();
    assert(fast.ok);
    assert(fast.solutions.size() == 16);
    assert(fast.summary.dedup_success_count == 16);
    assert(fast.summary.table_fallback_count == 0);
    assert(fast.summary.max_abs_slingshot_residual_at_root <= 1e-8);
    assert(fast.summary.strict_root_count + fast.summary.relaxed_boundary_root_count == 16);

    problem2::Problem2GravityAssistSolverOptions default_options{};
    default_options.theta_sample_count = 64;
    assert(default_options.topology_adaptive_enabled);
    const auto default_adaptive = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta, default_options);
    assert(default_adaptive.ok);
    assert(default_adaptive.solutions.size() == 16);

    auto fast_solutions = fast.solutions;
    auto default_solutions = default_adaptive.solutions;
    sort_solutions(&fast_solutions);
    sort_solutions(&default_solutions);
    for (std::size_t i = 0; i < fast_solutions.size(); ++i) {
        assert(same_root(fast_solutions[i], default_solutions[i]));
    }

    auto no_relaxed_options = fast_options;
    no_relaxed_options.allow_roundoff_boundary_relaxed_residual = false;
    const auto no_relaxed = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta, no_relaxed_options);
    assert(no_relaxed.ok);
    assert(no_relaxed.solutions.size() <= 16);
    assert(no_relaxed.summary.strict_root_count == static_cast<int>(no_relaxed.solutions.size()));
    assert(no_relaxed.summary.relaxed_boundary_root_count == 0);

    auto invalid_theta_options = default_options;
    invalid_theta_options.theta_sample_count = 3;
    const auto invalid_theta = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta, invalid_theta_options);
    assert(!invalid_theta.ok);
    assert(invalid_theta.error_message == "theta_sample_count_must_be_at_least_4");

    auto invalid_revolution_options = default_options;
    invalid_revolution_options.max_transfer_revolution = -1;
    const auto invalid_revolution = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta,
        invalid_revolution_options);
    assert(!invalid_revolution.ok);

    auto invalid_bisection_options = default_options;
    invalid_bisection_options.bisection_max_iterations = 0;
    const auto invalid_bisection = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta,
        invalid_bisection_options);
    assert(!invalid_bisection.ok);

    std::cout << "Problem2GravityAssistSolverFastModeSummary\n";
    std::cout << "solution_count=" << fast.solutions.size() << '\n';
    std::cout << "dedup_success_count=" << fast.summary.dedup_success_count << '\n';
    std::cout << "strict_root_count=" << fast.summary.strict_root_count << '\n';
    std::cout << "relaxed_boundary_root_count=" << fast.summary.relaxed_boundary_root_count << '\n';
    std::cout << "table_fallback_count=" << fast.summary.table_fallback_count << '\n';
    std::cout << "max_abs_slingshot_residual="
              << fast.summary.max_abs_slingshot_residual_at_root << '\n';

    std::cout << "Problem2GravityAssistSolverDefaultAdaptiveSummary\n";
    std::cout << "solution_count=" << default_adaptive.solutions.size() << '\n';
    std::cout << "dedup_success_count=" << default_adaptive.summary.dedup_success_count << '\n';
    std::cout << "strict_root_count=" << default_adaptive.summary.strict_root_count << '\n';
    std::cout << "relaxed_boundary_root_count=" << default_adaptive.summary.relaxed_boundary_root_count << '\n';
    std::cout << "table_fallback_count=" << default_adaptive.summary.table_fallback_count << '\n';
    std::cout << "max_abs_slingshot_residual="
              << default_adaptive.summary.max_abs_slingshot_residual_at_root << '\n';

    std::cout << "Problem2GravityAssistSolverNoRelaxedBoundarySummary\n";
    std::cout << "solution_count=" << no_relaxed.solutions.size() << '\n';
    std::cout << "strict_root_count=" << no_relaxed.summary.strict_root_count << '\n';
    std::cout << "relaxed_boundary_root_count=" << no_relaxed.summary.relaxed_boundary_root_count << '\n';
    std::cout << "max_abs_slingshot_residual="
              << no_relaxed.summary.max_abs_slingshot_residual_at_root << '\n';
    std::cout << "solver_ok=" << (no_relaxed.ok ? 1 : 0) << '\n';

    std::cout << "Problem2GravityAssistSolverOptionValidationSummary\n";
    std::cout << "invalid_theta_sample_rejected=" << (!invalid_theta.ok ? 1 : 0) << '\n';
    std::cout << "invalid_revolution_rejected=" << (!invalid_revolution.ok ? 1 : 0) << '\n';
    std::cout << "invalid_bisection_iterations_rejected=" << (!invalid_bisection.ok ? 1 : 0) << '\n';

    const double fast_solutions_per_second =
        fast_wall_time_ms > 0.0 ? 1000.0 * static_cast<double>(fast.solutions.size()) / fast_wall_time_ms
                                : 0.0;
    std::cout << "Problem2GravityAssistSolverPerformanceSummary\n";
    std::cout << "theta_sample_count=" << fast_options.theta_sample_count << '\n';
    std::cout << "topology_adaptive_enabled=" << (fast_options.topology_adaptive_enabled ? 1 : 0) << '\n';
    std::cout << "wall_time_ms=" << fast_wall_time_ms << '\n';
    std::cout << "solution_count=" << fast.solutions.size() << '\n';
    std::cout << "solutions_per_second=" << fast_solutions_per_second << '\n';

    std::cout << "test_ok=1\n";
    return 0;
}
