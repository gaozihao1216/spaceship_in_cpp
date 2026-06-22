/*
 * 文件作用：验证 Problem 2 飞掠求解返回的每个解在飞掠点处 |v_∞,in| = |v_∞,out|。
 * 主要工作：多行星对、多时刻、多入射 branch 下遍历全部 G=0 解；不一致则视为管线 bug。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>
#include <vector>

namespace {

using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2FlybySolveConfig;
using spaceship_cpp::problem2::Problem2FlybySolveInput;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;
using spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult;
using spaceship_cpp::trajectory::FlybyPhysicalFeasibilityOptions;
using spaceship_cpp::trajectory::FlybyPhysicalFilterMode;

constexpr double kDayInSeconds = 86400.0;
constexpr double kRelativeSpeedTolerance = 1e-6;

struct IncomingOrbitCase {
    double incoming_eccentricity = 0.0;
    double incoming_theta_local = 0.0;
};

struct ScenarioSpec {
    PlanetId flyby_planet{};
    PlanetId target_planet{};
    double flyby_time_seconds = 0.0;
};

FlybyPhysicalFeasibilityOptions make_verification_options() {
    FlybyPhysicalFeasibilityOptions options{};
    options.enabled = true;
    options.mode = FlybyPhysicalFilterMode::Disabled;
    options.relative_speed_tolerance = kRelativeSpeedTolerance;
    return options;
}

double v_inf_mismatch_tolerance(double v_infinity_in, double v_infinity_out) {
    const double v_inf = 0.5 * (v_infinity_in + v_infinity_out);
    return kRelativeSpeedTolerance * std::max(1.0, v_inf);
}

void collect_incoming_cases_with_derivatives(
    const Problem2ThetaPrimeInitialScanResult& scan,
    std::vector<IncomingOrbitCase>& out_cases
) {
    out_cases.clear();
    if (scan.nodes.size() < 3U) {
        return;
    }
    for (std::size_t node_index = 1; node_index + 1 < scan.nodes.size(); ++node_index) {
        const auto& node = scan.nodes[node_index];
        for (const Problem2OutgoingBranchSolution& branch : node.solutions) {
            if (branch.has_dphi_dtheta_prime && branch.has_de_dtheta_prime) {
                out_cases.push_back(
                    IncomingOrbitCase{
                        .incoming_eccentricity = branch.outgoing_eccentricity,
                        .incoming_theta_local = node.theta_prime_local,
                    });
            }
        }
    }
}

std::size_t verify_all_solutions_v_infinity_match(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& solve_config,
    const Problem2ThetaPrimeInitialScanResult& scan,
    const IncomingOrbitCase& incoming_case,
    const FlybyPhysicalFeasibilityOptions& verification_options
) {
    Problem2FlybySolveInput solve_input = input;
    solve_input.incoming_eccentricity = incoming_case.incoming_eccentricity;
    solve_input.incoming_theta = incoming_case.incoming_theta_local;
    solve_input.incoming_theta_is_local = true;

    const auto context = spaceship_cpp::problem2::build_problem2_flyby_G_context(
        solve_input.flyby_planet,
        solve_input.flyby_time_seconds_since_j2000,
        solve_input.incoming_eccentricity,
        solve_input.incoming_theta,
        solve_input.incoming_theta_is_local);
    if (!context.incoming_cache.valid) {
        return 0;
    }

    const auto solved =
        spaceship_cpp::problem2::solve_problem2_flyby_with_scan(solve_input, solve_config, scan);
    if (!solved.ok || solved.solutions.empty()) {
        return 0;
    }

    for (const auto& solution : solved.solutions) {
        const auto feasibility = spaceship_cpp::trajectory::evaluate_flyby_physical_feasibility(
            solve_input.flyby_planet,
            solve_input.flyby_time_seconds_since_j2000,
            solve_input.incoming_eccentricity,
            context.incoming_cache.incoming_theta_global,
            solution.outgoing_eccentricity,
            solution.outgoing_theta_prime_global,
            verification_options);
        assert(feasibility.valid);
        assert(feasibility.v_infinity_in > 0.0);
        assert(feasibility.v_infinity_out > 0.0);
        assert(feasibility.v_infinity_mismatch <=
            v_inf_mismatch_tolerance(
                feasibility.v_infinity_in,
                feasibility.v_infinity_out));
    }
    return solved.solutions.size();
}

std::size_t verify_full_solve_v_infinity_match(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& solve_config,
    const IncomingOrbitCase& incoming_case,
    const FlybyPhysicalFeasibilityOptions& verification_options
) {
    Problem2FlybySolveInput solve_input = input;
    solve_input.incoming_eccentricity = incoming_case.incoming_eccentricity;
    solve_input.incoming_theta = incoming_case.incoming_theta_local;
    solve_input.incoming_theta_is_local = true;

    const auto context = spaceship_cpp::problem2::build_problem2_flyby_G_context(
        solve_input.flyby_planet,
        solve_input.flyby_time_seconds_since_j2000,
        solve_input.incoming_eccentricity,
        solve_input.incoming_theta,
        solve_input.incoming_theta_is_local);
    if (!context.incoming_cache.valid) {
        return 0;
    }

    const auto solved = spaceship_cpp::problem2::solve_problem2_flyby(solve_input, solve_config);
    if (!solved.ok || solved.solutions.empty()) {
        return 0;
    }

    for (const auto& solution : solved.solutions) {
        const auto feasibility = spaceship_cpp::trajectory::evaluate_flyby_physical_feasibility(
            solve_input.flyby_planet,
            solve_input.flyby_time_seconds_since_j2000,
            solve_input.incoming_eccentricity,
            context.incoming_cache.incoming_theta_global,
            solution.outgoing_eccentricity,
            solution.outgoing_theta_prime_global,
            verification_options);
        assert(feasibility.valid);
        assert(feasibility.v_infinity_in > 0.0);
        assert(feasibility.v_infinity_out > 0.0);
        assert(feasibility.v_infinity_mismatch <=
            v_inf_mismatch_tolerance(
                feasibility.v_infinity_in,
                feasibility.v_infinity_out));
    }
    return solved.solutions.size();
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;

    const auto verification_options = make_verification_options();
    const auto& defaults = config::global_config();

    const ScenarioSpec scenarios[] = {
        {.flyby_planet = PlanetId::Earth, .target_planet = PlanetId::Mars, .flyby_time_seconds = 0.0},
        {.flyby_planet = PlanetId::Mars, .target_planet = PlanetId::Earth, .flyby_time_seconds = 0.0},
        {.flyby_planet = PlanetId::Venus, .target_planet = PlanetId::Earth, .flyby_time_seconds = 0.0},
        {.flyby_planet = PlanetId::Earth, .target_planet = PlanetId::Venus, .flyby_time_seconds = 0.0},
        {
            .flyby_planet = PlanetId::Earth,
            .target_planet = PlanetId::Mars,
            .flyby_time_seconds = 100.0 * kDayInSeconds,
        },
        {
            .flyby_planet = PlanetId::Mars,
            .target_planet = PlanetId::Earth,
            .flyby_time_seconds = 100.0 * kDayInSeconds,
        },
    };

    std::size_t scenario_count_with_solutions = 0;
    std::size_t incoming_case_count = 0;
    std::size_t solution_count_with_scan = 0;
    std::size_t solution_count_full_solve = 0;

    for (const ScenarioSpec& scenario : scenarios) {
        const auto solve_config = config::make_problem2_flyby_solve_config(
            scenario.flyby_planet,
            scenario.target_planet,
            scenario.flyby_time_seconds,
            defaults.problem2_theta_prime_scan,
            defaults.problem2_route_a_newton,
            defaults.problem1_solve);

        Problem2FlybySolveInput input{};
        input.flyby_planet = scenario.flyby_planet;
        input.target_planet = scenario.target_planet;
        input.flyby_time_seconds_since_j2000 = scenario.flyby_time_seconds;

        const auto scan =
            spaceship_cpp::problem2::run_problem2_flyby_theta_prime_initial_scan(input, solve_config);
        if (!scan.ok) {
            continue;
        }

        std::vector<IncomingOrbitCase> incoming_cases;
        collect_incoming_cases_with_derivatives(scan, incoming_cases);
        if (incoming_cases.empty()) {
            continue;
        }

        bool scenario_has_solution = false;
        for (const IncomingOrbitCase& incoming_case : incoming_cases) {
            ++incoming_case_count;
            const std::size_t verified_with_scan = verify_all_solutions_v_infinity_match(
                input,
                solve_config,
                scan,
                incoming_case,
                verification_options);
            solution_count_with_scan += verified_with_scan;
            if (verified_with_scan > 0U) {
                scenario_has_solution = true;
            }
        }

        const IncomingOrbitCase& first_incoming = incoming_cases.front();
        const std::size_t verified_full = verify_full_solve_v_infinity_match(
            input,
            solve_config,
            first_incoming,
            verification_options);
        solution_count_full_solve += verified_full;
        if (verified_full > 0U) {
            scenario_has_solution = true;
        }

        if (scenario_has_solution) {
            ++scenario_count_with_solutions;
        }
    }

    assert(scenario_count_with_solutions > 0U);
    assert(incoming_case_count > 0U);
    assert(solution_count_with_scan > 0U);
    assert(solution_count_full_solve > 0U);

    std::cout << "test_problem2_flyby_solve_physical PASSED"
              << " scenarios_with_solutions=" << scenario_count_with_solutions
              << " incoming_cases=" << incoming_case_count
              << " solutions_with_scan=" << solution_count_with_scan
              << " solutions_full_solve=" << solution_count_full_solve
              << '\n';
    return 0;
}
