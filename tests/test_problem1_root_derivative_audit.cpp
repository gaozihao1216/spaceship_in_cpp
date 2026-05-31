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

std::string transfer_e_bucket(double transfer_e) {
    if (!std::isfinite(transfer_e)) {
        return "non_finite_e";
    }
    if (transfer_e < 0.0) {
        return "negative_e_raw";
    }
    if (std::abs(transfer_e - 1.0) <= 1e-8) {
        return "near_parabolic";
    }
    if (transfer_e < 1.0) {
        return "elliptic";
    }
    return "hyperbolic";
}

}

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(departure_planet);
    const std::vector<double> launch_fractions{0.0, 0.125, 0.25, 0.375, 0.5};
    const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};

    int exact_root_count = 0;
    int residual_valid_count = 0;
    int derivative_valid_count = 0;
    int derivative_invalid_count = 0;
    std::map<std::string, int> invalid_reason_counts;
    std::map<int, int> invalid_by_k;
    std::map<int, int> invalid_by_q;
    std::map<std::string, int> invalid_by_transfer_e_bucket;
    int invalid_by_transfer_p_nonpositive = 0;
    int invalid_by_et_near_one = 0;
    int invalid_by_hyperbolic_domain = 0;
    int invalid_by_target_unwrap = 0;
    int invalid_by_transfer_unwrap = 0;
    int invalid_by_negative_e_normalization_branch = 0;
    int negative_e_branch_count = 0;
    int residual_valid_but_derivative_invalid = 0;

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
            for (const auto& branch : exact_branches) {
                if (!branch.valid) {
                    continue;
                }
                exact_root_count += 1;
                const problem1::Problem1RootResidualResult residual =
                    problem1::evaluate_problem1_root_residual(
                        departure_planet,
                        target_planet,
                        departure_state.varphi,
                        target_state.varphi,
                        theta_A,
                        branch.encounter_global_angle,
                        branch.transfer_revolution,
                        branch.target_revolution);
                const problem1::Problem1RootResidualDerivatives derivatives =
                    problem1::evaluate_problem1_root_residual_derivatives(
                        departure_planet,
                        target_planet,
                        departure_state.varphi,
                        target_state.varphi,
                        theta_A,
                        branch.encounter_global_angle,
                        branch.transfer_revolution,
                        branch.target_revolution);
                if (residual.valid) {
                    residual_valid_count += 1;
                }
                if (derivatives.valid) {
                    derivative_valid_count += 1;
                    if (residual.transfer_e_raw < 0.0) {
                        negative_e_branch_count += 1;
                    }
                    continue;
                }
                derivative_invalid_count += 1;
                if (residual.valid) {
                    residual_valid_but_derivative_invalid += 1;
                }
                if (residual.transfer_e_raw < 0.0) {
                    negative_e_branch_count += 1;
                    invalid_by_negative_e_normalization_branch += 1;
                }
                const std::string reason = derivatives.invalid_reason.empty() ? "none" : derivatives.invalid_reason;
                invalid_reason_counts[reason] += 1;
                invalid_by_k[branch.transfer_revolution] += 1;
                invalid_by_q[branch.target_revolution] += 1;
                invalid_by_transfer_e_bucket[transfer_e_bucket(residual.transfer_e)] += 1;
                if (!(residual.transfer_p > 0.0)) {
                    invalid_by_transfer_p_nonpositive += 1;
                }
                if (std::isfinite(residual.transfer_e) && std::abs(residual.transfer_e - 1.0) <= 1e-8) {
                    invalid_by_et_near_one += 1;
                }
                if (std::isfinite(residual.transfer_e) && residual.transfer_e >= 1.0) {
                    invalid_by_hyperbolic_domain += 1;
                }
                if (reason.find("target unwrap boundary") != std::string::npos) {
                    invalid_by_target_unwrap += 1;
                }
                if (reason.find("transfer unwrap boundary") != std::string::npos) {
                    invalid_by_transfer_unwrap += 1;
                }
            }
        }
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1 root derivative audit\n";
    std::cout << "exact_root_count=" << exact_root_count << '\n';
    std::cout << "residual_valid_count=" << residual_valid_count << '\n';
    std::cout << "derivative_valid_count=" << derivative_valid_count << '\n';
    std::cout << "derivative_invalid_count=" << derivative_invalid_count << '\n';
    std::cout << "residual_valid_but_derivative_invalid=" << residual_valid_but_derivative_invalid << '\n';
    for (const auto& [reason, count] : invalid_reason_counts) {
        std::cout << "invalid_reason[" << reason << "]=" << count << '\n';
    }
    for (const auto& [k, count] : invalid_by_k) {
        std::cout << "invalid_by_k[" << k << "]=" << count << '\n';
    }
    for (const auto& [q, count] : invalid_by_q) {
        std::cout << "invalid_by_q[" << q << "]=" << count << '\n';
    }
    for (const auto& [bucket, count] : invalid_by_transfer_e_bucket) {
        std::cout << "invalid_by_transfer_e_bucket[" << bucket << "]=" << count << '\n';
    }
    std::cout << "invalid_by_transfer_p_nonpositive=" << invalid_by_transfer_p_nonpositive << '\n';
    std::cout << "invalid_by_et_near_one=" << invalid_by_et_near_one << '\n';
    std::cout << "invalid_by_hyperbolic_domain=" << invalid_by_hyperbolic_domain << '\n';
    std::cout << "invalid_by_target_unwrap=" << invalid_by_target_unwrap << '\n';
    std::cout << "invalid_by_transfer_unwrap=" << invalid_by_transfer_unwrap << '\n';
    std::cout << "negative_e_branch_count=" << negative_e_branch_count << '\n';
    std::cout << "invalid_by_negative_e_normalization_branch=" << invalid_by_negative_e_normalization_branch << '\n';

    assert(exact_root_count > 0);
    assert(residual_valid_count > 0);
    return 0;
}
