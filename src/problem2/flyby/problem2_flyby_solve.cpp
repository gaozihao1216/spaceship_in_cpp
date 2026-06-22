/*
 * 文件作用：实现 Problem 2 飞掠求解顶层封装。
 */
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"

#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <chrono>
#include <cmath>
#include <limits>

namespace spaceship_cpp::problem2 {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

using spaceship_cpp::bfs::problem2_local_periapsis_angle_to_global;
using spaceship_cpp::common::is_finite;
using spaceship_cpp::planet_params::get_solar_system_physical_params;
using spaceship_cpp::problem1::Problem1ResidualInput;
using spaceship_cpp::problem1::Problem1ResidualStatus;

double compute_relative_problem1_residual(const problem1::Problem1ResidualResult& result) {
    const double scale = std::max(
        std::abs(result.transfer_time_scale_free),
        std::abs(result.target_time_scale_free));
    if (!is_finite(scale) || !(scale > 0.0) || !is_finite(result.residual)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::abs(result.residual) / scale;
}

std::optional<Problem2FlybySolution> enrich_flyby_G_solution(
    const Problem2FlybySolveInput& input,
    const Problem2FlybyGSolution& g_solution
) {
    if (!is_finite(g_solution.theta_prime_local) || !is_finite(g_solution.G)) {
        return std::nullopt;
    }

    const double theta_prime_global = problem2_local_periapsis_angle_to_global(
        input.flyby_planet,
        g_solution.theta_prime_local);
    const auto& branch = g_solution.outgoing_branch;

    const Problem1ResidualInput residual_input{
        input.flyby_planet,
        input.target_planet,
        input.flyby_time_seconds_since_j2000,
        theta_prime_global,
        branch.encounter_global_angle,
        branch.transfer_revolution,
        branch.target_revolution,
    };
    const auto residual_result = problem1::evaluate_problem1_residual(residual_input);
    if (residual_result.status != Problem1ResidualStatus::Success ||
        !is_finite(residual_result.transfer_time_scale_free) ||
        !(residual_result.transfer_time_scale_free > 0.0)) {
        return std::nullopt;
    }

    const double mu = get_solar_system_physical_params().GM_sun;
    const double time_of_flight_seconds =
        residual_result.transfer_time_scale_free / std::sqrt(mu);
    if (!is_finite(time_of_flight_seconds) || !(time_of_flight_seconds > 0.0)) {
        return std::nullopt;
    }

    Problem2FlybySolution solution{};
    solution.source = g_solution.source;
    solution.outgoing_theta_prime_local = g_solution.theta_prime_local;
    solution.outgoing_theta_prime_global = theta_prime_global;
    solution.outgoing_eccentricity = branch.outgoing_eccentricity;
    solution.outgoing_semi_latus_rectum = branch.outgoing_semi_latus_rectum;
    solution.encounter_global_angle = branch.encounter_global_angle;
    solution.transfer_revolution = branch.transfer_revolution;
    solution.target_revolution = branch.target_revolution;
    solution.flyby_time_seconds_since_j2000 = input.flyby_time_seconds_since_j2000;
    solution.time_of_flight_seconds = time_of_flight_seconds;
    solution.arrival_time_seconds_since_j2000 =
        input.flyby_time_seconds_since_j2000 + time_of_flight_seconds;
    solution.flyby_constraint_G = g_solution.G;
    solution.problem1_relative_residual = compute_relative_problem1_residual(residual_result);
    if (!is_finite(solution.problem1_relative_residual)) {
        solution.problem1_relative_residual = branch.relative_problem1_residual;
    }
    return solution;
}

bool same_flyby_solution(
    const Problem2FlybySolution& left,
    const Problem2FlybySolution& right,
    double theta_prime_tolerance
) {
    return left.transfer_revolution == right.transfer_revolution &&
        left.target_revolution == right.target_revolution &&
        is_finite(left.outgoing_theta_prime_local) &&
        is_finite(right.outgoing_theta_prime_local) &&
        std::abs(left.outgoing_theta_prime_local - right.outgoing_theta_prime_local) <=
        theta_prime_tolerance;
}

void append_unique_flyby_solution(
    std::vector<Problem2FlybySolution>& solutions,
    Problem2FlybySolution candidate,
    double theta_prime_tolerance
) {
    for (auto& existing : solutions) {
        if (same_flyby_solution(existing, candidate, theta_prime_tolerance)) {
            if (std::abs(candidate.flyby_constraint_G) < std::abs(existing.flyby_constraint_G)) {
                existing = std::move(candidate);
            }
            return;
        }
    }
    solutions.push_back(std::move(candidate));
}

Problem2FlybyGSearchConfig bind_flyby_g_search_config(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config
) {
    Problem2FlybyGSearchConfig g_search_config = config.g_search_config;
    g_search_config.scan_config.flyby_planet = input.flyby_planet;
    g_search_config.scan_config.target_planet = input.target_planet;
    g_search_config.scan_config.flyby_time_seconds_since_j2000 =
        input.flyby_time_seconds_since_j2000;
    return g_search_config;
}

Problem2FlybySolveResult solve_problem2_flyby_with_scan_impl(
    const Problem2FlybySolveInput& input,
    const Problem2FlybyGSearchConfig& g_search_config,
    const Problem2ThetaPrimeInitialScanResult& scan,
    bool collect_profile
) {
    Problem2FlybySolveResult result{};
    if (!scan.ok || scan.nodes.size() < 2U) {
        result.error_message = "invalid_theta_prime_initial_scan";
        return result;
    }

    Problem2FlybyGSearchProfile profile{};
    Problem2FlybyGSearchConfig g_search_config_with_profile = g_search_config;
    if (collect_profile) {
        g_search_config_with_profile.profile = &profile;
    }

    const Clock::time_point cache_start = Clock::now();
    const Problem2FlybyGContext context = build_problem2_flyby_G_context(
        input.flyby_planet,
        input.flyby_time_seconds_since_j2000,
        input.incoming_eccentricity,
        input.incoming_theta,
        input.incoming_theta_is_local);
    if (collect_profile) {
        profile.incoming_cache_ms = elapsed_ms(cache_start, Clock::now());
    }
    if (!context.incoming_cache.valid) {
        result.error_message = "invalid_flyby_G_context";
        if (collect_profile) {
            result.g_search_profile = profile;
            result.has_g_search_profile = true;
        }
        return result;
    }

    const Problem2FlybyGSearchResult search = search_flyby_constraint_G_zeros_from_initial_scan(
        context,
        g_search_config_with_profile,
        scan);
    if (!search.ok) {
        result.error_message =
            search.error_message.empty() ? "flyby_G_search_failed" : search.error_message;
        if (collect_profile) {
            result.g_search_profile = profile;
            result.has_g_search_profile = true;
        }
        return result;
    }

    const Clock::time_point enrich_start = Clock::now();
    for (const auto& g_solution : search.solutions) {
        const auto enriched = enrich_flyby_G_solution(input, g_solution);
        if (!enriched.has_value()) {
            continue;
        }
        append_unique_flyby_solution(
            result.solutions,
            *enriched,
            g_search_config.solution_theta_prime_tolerance);
    }
    if (collect_profile) {
        profile.enrich_ms = elapsed_ms(enrich_start, Clock::now());
        result.g_search_profile = profile;
        result.has_g_search_profile = true;
    }

    result.ok = true;
    return result;
}

bool validate_flyby_solve_input(const Problem2FlybySolveInput& input, std::string& error_message) {
    if (!is_finite(input.flyby_time_seconds_since_j2000)) {
        error_message = "non_finite_flyby_time";
        return false;
    }
    if (!is_finite(input.incoming_eccentricity) || !is_finite(input.incoming_theta)) {
        error_message = "non_finite_incoming_orbit";
        return false;
    }
    if (input.flyby_planet == input.target_planet) {
        error_message = "flyby_and_target_planet_must_differ";
        return false;
    }
    return true;
}

}  // namespace

Problem2ThetaPrimeInitialScanResult run_problem2_flyby_theta_prime_initial_scan(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config
) {
    Problem2ThetaPrimeInitialScanResult result{};
    std::string error_message;
    if (!validate_flyby_solve_input(input, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }

    const Problem2FlybyGSearchConfig g_search_config = bind_flyby_g_search_config(input, config);
    return run_problem2_theta_prime_initial_scan(g_search_config.scan_config);
}

Problem2FlybySolveResult solve_problem2_flyby_with_scan(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config,
    const Problem2ThetaPrimeInitialScanResult& scan
) {
    Problem2FlybySolveResult result{};
    if (!validate_flyby_solve_input(input, result.error_message)) {
        return result;
    }

    const Problem2FlybyGSearchConfig g_search_config = bind_flyby_g_search_config(input, config);
    return solve_problem2_flyby_with_scan_impl(
        input,
        g_search_config,
        scan,
        config.collect_g_search_profile);
}

Problem2FlybySolveResult solve_problem2_flyby(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config
) {
    Problem2FlybySolveResult result{};
    if (!validate_flyby_solve_input(input, result.error_message)) {
        return result;
    }

    const Problem2FlybyGSearchConfig g_search_config = bind_flyby_g_search_config(input, config);
    const Problem2ThetaPrimeInitialScanResult scan =
        run_problem2_theta_prime_initial_scan(g_search_config.scan_config);
    if (!scan.ok) {
        result.error_message =
            scan.error_message.empty() ? "theta_prime_initial_scan_failed" : scan.error_message;
        return result;
    }

    return solve_problem2_flyby_with_scan_impl(
        input,
        g_search_config,
        scan,
        config.collect_g_search_profile);
}

}  // namespace spaceship_cpp::problem2
