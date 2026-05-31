#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

namespace {

bool approx_equal(double lhs, double rhs, double abs_tol = 1e-6, double rel_tol = 1e-6) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= abs_tol + rel_tol * scale;
}

double residual_scale_free(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    const auto result = spaceship_cpp::problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        encounter_global_angle,
        transfer_revolution,
        target_revolution);
    assert(result.valid);
    return result.residual_scale_free;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace config = spaceship_cpp::config;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    {
        // 中文注释：先用有限差分验证 orbit_F 对 theta 和 e 的解析偏导。
        const std::vector<double> elliptic_e_values{0.0, 0.1, 0.5, 0.8};
        const std::vector<double> hyperbolic_e_values{1.2, 1.5, 2.0};
        const std::vector<double> theta_values{0.1, 1.0, 2.0, 4.0};
        const double theta_h = 1e-7;
        const double e_h = 1e-7;

        for (double e : elliptic_e_values) {
            for (double theta : theta_values) {
                const double theta_fd =
                    (common::orbit_F(e, theta + theta_h) - common::orbit_F(e, theta - theta_h)) / (2.0 * theta_h);
                const double theta_analytic = common::orbit_F_theta_derivative(e, theta);
                assert(approx_equal(theta_analytic, theta_fd, 1e-6, 1e-6));

                if (e - e_h < 0.0 || e + e_h >= 1.0) {
                    continue;
                }
                const double e_fd =
                    (common::orbit_F(e + e_h, theta) - common::orbit_F(e - e_h, theta)) / (2.0 * e_h);
                const double e_analytic = common::orbit_F_e_derivative(e, theta);
                assert(approx_equal(e_analytic, e_fd, 1e-5, 1e-5));
            }
        }

        for (double e : hyperbolic_e_values) {
            for (double theta : {0.1, 0.5, 1.0, 2.0}) {
                if (!common::is_orbit_F_defined(e, theta)) {
                    continue;
                }
                if (!common::is_orbit_F_defined(e, theta + theta_h) ||
                    !common::is_orbit_F_defined(e, theta - theta_h)) {
                    continue;
                }
                const double theta_fd =
                    (common::orbit_F(e, theta + theta_h) - common::orbit_F(e, theta - theta_h)) / (2.0 * theta_h);
                const double theta_analytic = common::orbit_F_theta_derivative(e, theta);
                assert(approx_equal(theta_analytic, theta_fd, 1e-6, 1e-6));

                if (e - e_h <= 1.0 ||
                    !common::is_orbit_F_defined(e + e_h, theta) ||
                    !common::is_orbit_F_defined(e - e_h, theta)) {
                    continue;
                }
                const double e_fd =
                    (common::orbit_F(e + e_h, theta) - common::orbit_F(e - e_h, theta)) / (2.0 * e_h);
                const double e_analytic = common::orbit_F_e_derivative(e, theta);
                assert(approx_equal(e_analytic, e_fd, 1e-5, 1e-5));
            }
        }
    }

    {
        // 中文注释：再验证 root residual 的解析导数与有限差分一致。
        const double earth_period = planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
        const std::vector<std::pair<planet_params::PlanetId, planet_params::PlanetId>> planet_pairs{
            {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars},
            {planet_params::PlanetId::Earth, planet_params::PlanetId::Venus},
        };
        const std::vector<double> launch_times{0.0, 0.25 * earth_period};
        const std::vector<double> transfer_perihelion_angles{0.2, 0.5};
        bool checked_any_branch = false;

        for (const auto& [departure_planet, target_planet] : planet_pairs) {
            for (double launch_time_seconds_since_j2000 : launch_times) {
                const planet_params::PlanetState departure_state =
                    planet_params::planet_state_at_time(departure_planet, launch_time_seconds_since_j2000);
                const planet_params::PlanetState target_state =
                    planet_params::planet_state_at_time(target_planet, launch_time_seconds_since_j2000);

                for (double transfer_perihelion_angle : transfer_perihelion_angles) {
                    const double theta_A = common::normalize_angle_0_2pi(
                        departure_state.theta_global - transfer_perihelion_angle);
                    const auto branches = problem1::solve_problem1_from_departure_anomalies(
                        departure_planet,
                        target_planet,
                        departure_state.varphi,
                        target_state.varphi,
                        theta_A,
                        1,
                        1);

                    for (const auto& branch : branches) {
                        const auto derivatives = problem1::evaluate_problem1_root_residual_derivatives(
                            departure_planet,
                            target_planet,
                            departure_state.varphi,
                            target_state.varphi,
                            theta_A,
                            branch.encounter_global_angle,
                            branch.transfer_revolution,
                            branch.target_revolution);
                        if (!derivatives.valid) {
                            continue;
                        }

                        const double h = 1e-7;
                        const double R_alpha_fd =
                            (residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi,
                                 target_state.varphi,
                                 theta_A,
                                 branch.encounter_global_angle + h,
                                 branch.transfer_revolution,
                                 branch.target_revolution) -
                             residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi,
                                 target_state.varphi,
                                 theta_A,
                                 branch.encounter_global_angle - h,
                                 branch.transfer_revolution,
                                 branch.target_revolution)) /
                            (2.0 * h);
                        const double R_nu_A_fd =
                            (residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi + h,
                                 target_state.varphi,
                                 theta_A,
                                 branch.encounter_global_angle,
                                 branch.transfer_revolution,
                                 branch.target_revolution) -
                             residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi - h,
                                 target_state.varphi,
                                 theta_A,
                                 branch.encounter_global_angle,
                                 branch.transfer_revolution,
                                 branch.target_revolution)) /
                            (2.0 * h);
                        const double R_nu_B_fd =
                            (residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi,
                                 target_state.varphi + h,
                                 theta_A,
                                 branch.encounter_global_angle,
                                 branch.transfer_revolution,
                                 branch.target_revolution) -
                             residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi,
                                 target_state.varphi - h,
                                 theta_A,
                                 branch.encounter_global_angle,
                                 branch.transfer_revolution,
                                 branch.target_revolution)) /
                            (2.0 * h);
                        const double R_theta_A_fd =
                            (residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi,
                                 target_state.varphi,
                                 theta_A + h,
                                 branch.encounter_global_angle,
                                 branch.transfer_revolution,
                                 branch.target_revolution) -
                             residual_scale_free(
                                 departure_planet,
                                 target_planet,
                                 departure_state.varphi,
                                 target_state.varphi,
                                 theta_A - h,
                                 branch.encounter_global_angle,
                                 branch.transfer_revolution,
                                 branch.target_revolution)) /
                            (2.0 * h);

                        assert(approx_equal(derivatives.R_alpha, R_alpha_fd, 1e-4, 1e-4));
                        assert(approx_equal(derivatives.R_nu_A, R_nu_A_fd, 1e-4, 1e-4));
                        assert(approx_equal(derivatives.R_nu_B, R_nu_B_fd, 1e-4, 1e-4));
                        assert(approx_equal(derivatives.R_theta_A, R_theta_A_fd, 1e-4, 1e-4));

                        const auto enriched = problem1::attach_problem1_root_derivatives(
                            departure_planet,
                            target_planet,
                            departure_state.varphi,
                            target_state.varphi,
                            theta_A,
                            branch);
                        assert(enriched.derivatives_available);
                        assert(std::isfinite(enriched.d_encounter_global_angle_d_nu_A));
                        assert(std::isfinite(enriched.d_encounter_global_angle_d_nu_B));
                        assert(std::isfinite(enriched.d_encounter_global_angle_d_theta_A));

                        checked_any_branch = true;
                        break;
                    }
                    if (checked_any_branch) {
                        break;
                    }
                }
                if (checked_any_branch) {
                    break;
                }
            }
            if (checked_any_branch) {
                break;
            }
        }

        assert(checked_any_branch);
    }

    return 0;
}
