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

struct Scenario {
    int id = 0;
    double encounter_time_factor = 0.0;
    double incoming_e = 0.0;
    double incoming_theta = 0.0;
};

struct TimedResult {
    spaceship_cpp::problem2::Problem2GravityAssistSolverResult result;
    double wall_time_ms = 0.0;
};

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

TimedResult solve_timed(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double incoming_e,
    double incoming_theta,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverOptions& options
) {
    const auto start = std::chrono::steady_clock::now();
    auto result = spaceship_cpp::problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, incoming_e, incoming_theta, options);
    const auto end = std::chrono::steady_clock::now();

    TimedResult timed{};
    timed.result = std::move(result);
    timed.wall_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    return timed;
}

int count_matched_roots(
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& baseline,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& adaptive
) {
    int matched = 0;
    std::vector<bool> adaptive_used(adaptive.size(), false);
    for (const auto& base_root : baseline) {
        for (std::size_t i = 0; i < adaptive.size(); ++i) {
            if (!adaptive_used[i] && same_root(base_root, adaptive[i])) {
                adaptive_used[i] = true;
                matched += 1;
                break;
            }
        }
    }
    return matched;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_gravity_assist_solver_multi_scenario_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(departure_planet);

    const std::vector<Scenario> scenarios = {
        {0, 0.17, 0.3, 0.4},
        {1, 0.05, 0.2, 0.1},
        {2, 0.25, 0.4, 1.0},
        {3, 0.50, 0.3, 2.0},
        {4, 0.75, 0.6, 3.0},
    };

    problem2::Problem2GravityAssistSolverOptions fast_options{};
    fast_options.theta_sample_count = 64;
    fast_options.max_transfer_revolution = 1;
    fast_options.max_target_revolution = 1;
    fast_options.allow_roundoff_boundary_relaxed_residual = true;
    fast_options.topology_adaptive_enabled = false;

    problem2::Problem2GravityAssistSolverOptions default_options{};
    default_options.theta_sample_count = 64;
    default_options.max_transfer_revolution = 1;
    default_options.max_target_revolution = 1;
    default_options.allow_roundoff_boundary_relaxed_residual = true;
    assert(default_options.topology_adaptive_enabled);

    int ok_count = 0;
    int failed_count = 0;
    int total_solution_count = 0;
    int total_strict_root_count = 0;
    int total_relaxed_boundary_root_count = 0;
    int max_table_fallback_count = 0;
    double max_abs_slingshot_residual_over_all = 0.0;
    double max_wall_time_ms = 0.0;

    std::vector<TimedResult> default_results;
    default_results.reserve(scenarios.size());

    for (const auto& scenario : scenarios) {
        const double encounter_time = scenario.encounter_time_factor * earth_period;
        TimedResult timed = solve_timed(
            loader,
            departure_planet,
            target_planet,
            encounter_time,
            scenario.incoming_e,
            scenario.incoming_theta,
            default_options);

        const auto& result = timed.result;
        default_results.push_back(timed);

        if (result.ok) {
            ok_count += 1;
        } else {
            failed_count += 1;
        }

        const int solution_count = static_cast<int>(result.solutions.size());
        total_solution_count += solution_count;
        total_strict_root_count += result.summary.strict_root_count;
        total_relaxed_boundary_root_count += result.summary.relaxed_boundary_root_count;
        max_table_fallback_count = std::max(max_table_fallback_count, result.summary.table_fallback_count);
        max_wall_time_ms = std::max(max_wall_time_ms, timed.wall_time_ms);
        if (solution_count > 0) {
            max_abs_slingshot_residual_over_all = std::max(
                max_abs_slingshot_residual_over_all, result.summary.max_abs_slingshot_residual_at_root);
        }

        assert(result.ok);
        assert(result.summary.strict_root_count + result.summary.relaxed_boundary_root_count == solution_count);
        if (solution_count > 0) {
            assert(result.summary.max_abs_slingshot_residual_at_root <= 1e-8);
        }
        if (scenario.id == 0) {
            assert(solution_count == 16);
        }

        std::cout << "Problem2GravityAssistSolverDefaultAdaptiveScenarioSummary\n";
        std::cout << "scenario_id=" << scenario.id << '\n';
        std::cout << "encounter_time_factor=" << scenario.encounter_time_factor << '\n';
        std::cout << "incoming_e=" << scenario.incoming_e << '\n';
        std::cout << "incoming_theta=" << scenario.incoming_theta << '\n';
        std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
        std::cout << "solution_count=" << solution_count << '\n';
        std::cout << "strict_root_count=" << result.summary.strict_root_count << '\n';
        std::cout << "relaxed_boundary_root_count=" << result.summary.relaxed_boundary_root_count << '\n';
        std::cout << "sign_change_candidate_count=" << result.summary.sign_change_candidate_count << '\n';
        std::cout << "continuation_stable_candidate_count="
                  << result.summary.continuation_stable_candidate_count << '\n';
        std::cout << "bisection_attempt_count=" << result.summary.bisection_attempt_count << '\n';
        std::cout << "bisection_success_count=" << result.summary.bisection_success_count << '\n';
        std::cout << "dedup_success_count=" << result.summary.dedup_success_count << '\n';
        std::cout << "table_fallback_count=" << result.summary.table_fallback_count << '\n';
        std::cout << "max_abs_slingshot_residual="
                  << result.summary.max_abs_slingshot_residual_at_root << '\n';
        std::cout << "max_abs_problem1_residual_seconds="
                  << result.summary.max_abs_problem1_residual_seconds_at_root << '\n';
        std::cout << "wall_time_ms=" << timed.wall_time_ms << '\n';
    }

    const Scenario& scenario2 = scenarios[2];
    const double scenario2_encounter_time = scenario2.encounter_time_factor * earth_period;
    const TimedResult scenario2_fast = solve_timed(
        loader,
        departure_planet,
        target_planet,
        scenario2_encounter_time,
        scenario2.incoming_e,
        scenario2.incoming_theta,
        fast_options);
    assert(scenario2_fast.result.ok);

    const auto& scenario2_default = default_results[2].result;
    const int scenario2_matched_root_count =
        count_matched_roots(scenario2_fast.result.solutions, scenario2_default.solutions);
    const int scenario2_fast_solution_count = static_cast<int>(scenario2_fast.result.solutions.size());
    const int scenario2_default_solution_count = static_cast<int>(scenario2_default.solutions.size());
    const int scenario2_new_adaptive_root_count =
        scenario2_default_solution_count - scenario2_matched_root_count;
    const int scenario2_missing_fast_root_count =
        scenario2_fast_solution_count - scenario2_matched_root_count;

    assert(scenario2_fast_solution_count == 12);
    assert(scenario2_default_solution_count == 16);
    assert(scenario2_matched_root_count == 12);
    assert(scenario2_new_adaptive_root_count == 4);
    assert(scenario2_missing_fast_root_count == 0);

    std::cout << "Problem2GravityAssistSolverScenario2FastVsDefaultAdaptiveRegression\n";
    std::cout << "fast_solution_count=" << scenario2_fast_solution_count << '\n';
    std::cout << "default_adaptive_solution_count=" << scenario2_default_solution_count << '\n';
    std::cout << "matched_root_count=" << scenario2_matched_root_count << '\n';
    std::cout << "new_default_adaptive_root_count=" << scenario2_new_adaptive_root_count << '\n';
    std::cout << "missing_fast_root_count=" << scenario2_missing_fast_root_count << '\n';

    const bool multi_scenario_ok = failed_count == 0;
    std::cout << "Problem2GravityAssistSolverDefaultAdaptiveMultiScenarioSummary\n";
    std::cout << "scenario_count=" << scenarios.size() << '\n';
    std::cout << "ok_count=" << ok_count << '\n';
    std::cout << "failed_count=" << failed_count << '\n';
    std::cout << "total_solution_count=" << total_solution_count << '\n';
    std::cout << "total_strict_root_count=" << total_strict_root_count << '\n';
    std::cout << "total_relaxed_boundary_root_count=" << total_relaxed_boundary_root_count << '\n';
    std::cout << "max_abs_slingshot_residual_over_all=" << max_abs_slingshot_residual_over_all << '\n';
    std::cout << "max_table_fallback_count=" << max_table_fallback_count << '\n';
    std::cout << "max_wall_time_ms=" << max_wall_time_ms << '\n';
    std::cout << "multi_scenario_ok=" << (multi_scenario_ok ? 1 : 0) << '\n';

    return multi_scenario_ok ? 0 : 1;
}
