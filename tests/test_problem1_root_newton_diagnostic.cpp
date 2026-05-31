#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kStrictAlphaThreshold = 1e-5;
constexpr double kStrictTimeThresholdSeconds = 1e-2;
constexpr double kStrictResidualThresholdSeconds = 1e-6;

struct NewtonModeStats {
    std::string name;
    spaceship_cpp::problem1::Problem1RootDerivativeMode derivative_mode =
        spaceship_cpp::problem1::Problem1RootDerivativeMode::AnalyticOnly;
    int total_trials = 0;
    int strict_match_count = 0;
    int small_perturbation_success_count = 0;
    int wrong_root_count = 0;
    int derivative_failed_count = 0;
    int fallback_used_count = 0;
    int finite_difference_success_count = 0;
    int invalid_refinement_count = 0;
    int high_residual_count = 0;
    double max_invalid_final_residual_seconds = 0.0;
    double max_invalid_final_residual_scale_free = 0.0;
    double max_high_residual_seconds = 0.0;
    std::map<std::string, int> invalid_reason_count;
};

double alpha_error(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return std::abs(normalize_angle_minus_pi_pi(lhs.encounter_global_angle - rhs.encounter_global_angle));
}

bool is_strict_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return
        alpha_error(candidate, exact) <= kStrictAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= kStrictTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kStrictResidualThresholdSeconds;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(departure_planet);
    const std::vector<double> launch_fractions{0.0, 0.125, 0.25};
    const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};
    const std::vector<double> perturbations{1e-6, -1e-6, 1e-5, -1e-5, 1e-4, -1e-4, 1e-3, -1e-3, 1e-2, -1e-2};
    const double newton_residual_tolerance_seconds = 1e-6;
    const double newton_residual_tolerance_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(newton_residual_tolerance_seconds);

    std::vector<NewtonModeStats> modes{
        {"AnalyticOnly", problem1::Problem1RootDerivativeMode::AnalyticOnly},
        {"FiniteDifferenceOnly", problem1::Problem1RootDerivativeMode::FiniteDifferenceOnly},
        {"AnalyticWithFiniteDifferenceFallback", problem1::Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback},
    };
    int exact_root_count = 0;
    int polished_exact_count = 0;
    int unpolished_exact_count = 0;
    double max_polish_alpha_shift = 0.0;
    double max_polish_time_shift = 0.0;

    for (double launch_fraction : launch_fractions) {
        const double launch_time = earth_period * launch_fraction;
        const planet_params::PlanetState departure_state =
            planet_params::planet_state_at_time(departure_planet, launch_time);
        const planet_params::PlanetState target_state =
            planet_params::planet_state_at_time(target_planet, launch_time);
        for (double transfer_perihelion_angle : transfer_perihelion_angles) {
            const double theta_A =
                common::normalize_angle_0_2pi(departure_state.theta_global - transfer_perihelion_angle);
            const std::vector<problem1::Problem1SolutionBranch> exact_branches =
                problem1::solve_problem1_from_departure_anomalies(
                    departure_planet,
                    target_planet,
                    departure_state.varphi,
                    target_state.varphi,
                    theta_A,
                    1,
                    1);
            for (const auto& exact_branch_raw : exact_branches) {
                if (!exact_branch_raw.valid) {
                    continue;
                }

                problem1::Problem1SolutionBranch exact_branch = exact_branch_raw;

                const problem1::Problem1RootRefinementResult polished_exact =
                    problem1::refine_problem1_root_branch_newton_diagnostic(
                        departure_planet,
                        target_planet,
                        departure_state.varphi,
                        target_state.varphi,
                        theta_A,
                        exact_branch_raw.transfer_revolution,
                        exact_branch_raw.target_revolution,
                        exact_branch_raw.encounter_global_angle,
                        40,
                        newton_residual_tolerance_scale_free,
                        1e-14,
                        std::numeric_limits<double>::infinity(),
                        true,
                        8,
                        problem1::Problem1RootDerivativeMode::AnalyticOnly,
                        1e-6);

                if (polished_exact.valid) {
                    exact_branch = polished_exact.branch;
                    polished_exact_count += 1;

                    max_polish_alpha_shift = std::max(
                        max_polish_alpha_shift,
                        std::abs(normalize_angle_minus_pi_pi(
                            exact_branch.encounter_global_angle -
                            exact_branch_raw.encounter_global_angle)));

                    max_polish_time_shift = std::max(
                        max_polish_time_shift,
                        std::abs(
                            exact_branch.time_of_flight_seconds -
                            exact_branch_raw.time_of_flight_seconds));
                } else {
                    unpolished_exact_count += 1;
                }

                exact_root_count += 1;

                for (double perturbation : perturbations) {
                    const bool small_perturbation = std::abs(perturbation) <= 1e-4;
                    for (NewtonModeStats& mode : modes) {
                        mode.total_trials += 1;
                        const problem1::Problem1RootRefinementResult refined =
                            problem1::refine_problem1_root_branch_newton_diagnostic(
                                departure_planet,
                                target_planet,
                                departure_state.varphi,
                                target_state.varphi,
                                theta_A,
                                exact_branch.transfer_revolution,
                                exact_branch.target_revolution,
                                exact_branch.encounter_global_angle + perturbation,
                                40,
                                newton_residual_tolerance_scale_free,
                                1e-12,
                                std::numeric_limits<double>::infinity(),
                                false,
                                0,
                                mode.derivative_mode,
                                1e-6);
                        mode.fallback_used_count += refined.diagnostic.fallback_used_count;
                        mode.finite_difference_success_count += refined.diagnostic.finite_difference_success_count;
                        if (refined.diagnostic.derivative_failed) {
                            mode.derivative_failed_count += 1;
                        }

                        if (!refined.valid) {
                            mode.invalid_refinement_count += 1;

                            if (std::isfinite(refined.diagnostic.final_residual_seconds)) {
                                mode.max_invalid_final_residual_seconds = std::max(
                                    mode.max_invalid_final_residual_seconds,
                                    std::abs(refined.diagnostic.final_residual_seconds));
                            }

                            if (!refined.diagnostic.trace.empty()) {
                                const auto& last_step = refined.diagnostic.trace.back();
                                if (std::isfinite(last_step.residual_scale_free)) {
                                    mode.max_invalid_final_residual_scale_free = std::max(
                                        mode.max_invalid_final_residual_scale_free,
                                        std::abs(last_step.residual_scale_free));
                                }
                            }

                            const std::string reason = refined.diagnostic.invalid_reason.empty()
                                ? "invalid_without_reason"
                                : refined.diagnostic.invalid_reason;
                            mode.invalid_reason_count[reason] += 1;
                            continue;
                        }

                        if (std::abs(refined.branch.residual_seconds) > 1e-3) {
                            mode.high_residual_count += 1;
                            mode.max_high_residual_seconds = std::max(
                                mode.max_high_residual_seconds,
                                std::abs(refined.branch.residual_seconds));
                            continue;
                        }
                        if (is_strict_match(refined.branch, exact_branch)) {
                            mode.strict_match_count += 1;
                            if (small_perturbation) {
                                mode.small_perturbation_success_count += 1;
                            }
                        } else {
                            mode.wrong_root_count += 1;
                        }
                    }
                }
            }
        }
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1 Newton diagnostic by derivative mode\n";
    std::cout << "newton_residual_tolerance_seconds=" << newton_residual_tolerance_seconds << '\n';
    std::cout << "newton_residual_tolerance_scale_free=" << newton_residual_tolerance_scale_free << '\n';
    std::cout << "exact_root_count=" << exact_root_count << '\n';
    std::cout << "polished_exact_count=" << polished_exact_count << '\n';
    std::cout << "unpolished_exact_count=" << unpolished_exact_count << '\n';
    std::cout << "max_polish_alpha_shift=" << max_polish_alpha_shift << '\n';
    std::cout << "max_polish_time_shift=" << max_polish_time_shift << '\n';
    for (const NewtonModeStats& mode : modes) {
        std::cout << "mode=" << mode.name << '\n';
        std::cout << "total_trials=" << mode.total_trials << '\n';
        std::cout << "strict_match_count=" << mode.strict_match_count << '\n';
        std::cout << "small_perturbation_success_count=" << mode.small_perturbation_success_count << '\n';
        std::cout << "wrong_root_count=" << mode.wrong_root_count << '\n';
        std::cout << "derivative_failed_count=" << mode.derivative_failed_count << '\n';
        std::cout << "fallback_used_count=" << mode.fallback_used_count << '\n';
        std::cout << "finite_difference_success_count=" << mode.finite_difference_success_count << '\n';
        std::cout << "invalid_refinement_count=" << mode.invalid_refinement_count << '\n';
        std::cout << "high_residual_count=" << mode.high_residual_count << '\n';
        std::cout << "max_invalid_final_residual_seconds=" << mode.max_invalid_final_residual_seconds << '\n';
        std::cout << "max_invalid_final_residual_scale_free=" << mode.max_invalid_final_residual_scale_free << '\n';
        std::cout << "max_high_residual_seconds=" << mode.max_high_residual_seconds << '\n';
        for (const auto& [reason, count] : mode.invalid_reason_count) {
            std::cout << "invalid_reason[" << reason << "]=" << count << '\n';
        }
    }

    assert(exact_root_count > 0);
    for (const NewtonModeStats& mode : modes) {
        assert(mode.total_trials > 0);
    }
    return 0;
}
