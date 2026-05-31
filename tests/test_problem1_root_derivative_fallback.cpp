#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

double rel_error(double lhs, double rhs) {
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) / scale;
}

}

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(departure_planet);
    const std::vector<double> launch_fractions{0.0, 0.125, 0.25};
    const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};

    int analytic_valid_count = 0;
    int negative_e_analytic_valid_count = 0;
    int negative_e_skipped_count = 0;
    double max_rel_R_alpha = 0.0;
    double max_rel_R_nu_A = 0.0;
    double max_rel_R_nu_B = 0.0;
    double max_rel_R_theta_A = 0.0;
    double max_rel_dalpha_nu_A = 0.0;
    double max_rel_dalpha_nu_B = 0.0;
    double max_rel_dalpha_theta_A = 0.0;

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
                const problem1::Problem1RootResidualDerivatives analytic =
                    problem1::evaluate_problem1_root_residual_derivatives_with_mode(
                        departure_planet,
                        target_planet,
                        departure_state.varphi,
                        target_state.varphi,
                        theta_A,
                        branch.encounter_global_angle,
                        branch.transfer_revolution,
                        branch.target_revolution,
                        problem1::Problem1RootDerivativeMode::AnalyticOnly,
                        1e-6);
                if (!analytic.valid) {
                    if (problem1::evaluate_problem1_root_residual(
                            departure_planet,
                            target_planet,
                            departure_state.varphi,
                            target_state.varphi,
                            theta_A,
                            branch.encounter_global_angle,
                            branch.transfer_revolution,
                            branch.target_revolution).transfer_e_raw < 0.0) {
                        negative_e_skipped_count += 1;
                    }
                    continue;
                }
                const auto residual = problem1::evaluate_problem1_root_residual(
                    departure_planet,
                    target_planet,
                    departure_state.varphi,
                    target_state.varphi,
                    theta_A,
                    branch.encounter_global_angle,
                    branch.transfer_revolution,
                    branch.target_revolution);
                const problem1::Problem1RootResidualDerivatives finite_difference =
                    problem1::evaluate_problem1_root_residual_derivatives_with_mode(
                        departure_planet,
                        target_planet,
                        departure_state.varphi,
                        target_state.varphi,
                        theta_A,
                        branch.encounter_global_angle,
                        branch.transfer_revolution,
                        branch.target_revolution,
                        problem1::Problem1RootDerivativeMode::FiniteDifferenceOnly,
                        1e-6);
                assert(finite_difference.valid);
                analytic_valid_count += 1;
                if (residual.transfer_e_raw < 0.0) {
                    negative_e_analytic_valid_count += 1;
                }
                max_rel_R_alpha = std::max(max_rel_R_alpha, rel_error(analytic.R_alpha, finite_difference.R_alpha));
                max_rel_R_nu_A = std::max(max_rel_R_nu_A, rel_error(analytic.R_nu_A, finite_difference.R_nu_A));
                max_rel_R_nu_B = std::max(max_rel_R_nu_B, rel_error(analytic.R_nu_B, finite_difference.R_nu_B));
                max_rel_R_theta_A = std::max(
                    max_rel_R_theta_A,
                    rel_error(analytic.R_theta_A, finite_difference.R_theta_A));
                max_rel_dalpha_nu_A = std::max(
                    max_rel_dalpha_nu_A,
                    rel_error(analytic.d_alpha_d_nu_A, finite_difference.d_alpha_d_nu_A));
                max_rel_dalpha_nu_B = std::max(
                    max_rel_dalpha_nu_B,
                    rel_error(analytic.d_alpha_d_nu_B, finite_difference.d_alpha_d_nu_B));
                max_rel_dalpha_theta_A = std::max(
                    max_rel_dalpha_theta_A,
                    rel_error(analytic.d_alpha_d_theta_A, finite_difference.d_alpha_d_theta_A));
            }
        }
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1 derivative fallback comparison\n";
    std::cout << "analytic_valid_count=" << analytic_valid_count << '\n';
    std::cout << "negative_e_analytic_valid_count=" << negative_e_analytic_valid_count << '\n';
    std::cout << "negative_e_skipped_count=" << negative_e_skipped_count << '\n';
    std::cout << "max_rel_R_alpha=" << max_rel_R_alpha << '\n';
    std::cout << "max_rel_R_nu_A=" << max_rel_R_nu_A << '\n';
    std::cout << "max_rel_R_nu_B=" << max_rel_R_nu_B << '\n';
    std::cout << "max_rel_R_theta_A=" << max_rel_R_theta_A << '\n';
    std::cout << "max_rel_dalpha_nu_A=" << max_rel_dalpha_nu_A << '\n';
    std::cout << "max_rel_dalpha_nu_B=" << max_rel_dalpha_nu_B << '\n';
    std::cout << "max_rel_dalpha_theta_A=" << max_rel_dalpha_theta_A << '\n';

    assert(analytic_valid_count > 0);
    assert(negative_e_analytic_valid_count > 0);
    return 0;
}
