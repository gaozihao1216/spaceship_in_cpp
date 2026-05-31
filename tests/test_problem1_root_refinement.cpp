#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace {

bool approx_equal(double lhs, double rhs, double abs_tol = 1e-6, double rel_tol = 1e-6) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= abs_tol + rel_tol * scale;
}

std::optional<spaceship_cpp::problem1::Problem1SolutionBranch> find_same_kq_nearest_branch(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    int transfer_revolution,
    int target_revolution,
    double encounter_global_angle
) {
    namespace common = spaceship_cpp::common;
    const spaceship_cpp::problem1::Problem1SolutionBranch* best = nullptr;
    double best_distance = std::numeric_limits<double>::infinity();
    for (const auto& branch : branches) {
        if (!branch.valid ||
            branch.transfer_revolution != transfer_revolution ||
            branch.target_revolution != target_revolution) {
            continue;
        }
        const double distance = std::abs(
            common::normalize_angle_minus_pi_pi(branch.encounter_global_angle - encounter_global_angle));
        if (distance < best_distance) {
            best_distance = distance;
            best = &branch;
        }
    }
    if (best == nullptr) {
        return std::nullopt;
    }
    return *best;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const double newton_residual_tolerance_seconds = 1e-6;
    const double newton_residual_tolerance_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(newton_residual_tolerance_seconds);

    {
        // 中文注释：从真实 root 附近的小扰动出发，Newton 应能回到同一条 root。
        const double earth_period = planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
        const std::vector<std::pair<planet_params::PlanetId, planet_params::PlanetId>> planet_pairs{
            {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars},
            {planet_params::PlanetId::Earth, planet_params::PlanetId::Venus},
        };
        const std::vector<double> launch_times{0.0, 0.25 * earth_period, 0.5 * earth_period};
        const std::vector<double> transfer_perihelion_angles{0.2, 0.5};
        const std::vector<int> max_revolution_values{0, 1};
        bool verified = false;

        for (const auto& [departure_planet, target_planet] : planet_pairs) {
            for (double launch_time_seconds_since_j2000 : launch_times) {
                const planet_params::PlanetState departure_state =
                    planet_params::planet_state_at_time(departure_planet, launch_time_seconds_since_j2000);
                const planet_params::PlanetState target_state =
                    planet_params::planet_state_at_time(target_planet, launch_time_seconds_since_j2000);

                for (double transfer_perihelion_angle : transfer_perihelion_angles) {
                    const double theta_A =
                        common::normalize_angle_0_2pi(departure_state.theta_global - transfer_perihelion_angle);

                    for (int max_revolutions : max_revolution_values) {
                        const auto branches = problem1::solve_problem1_from_departure_anomalies(
                            departure_planet,
                            target_planet,
                            departure_state.varphi,
                            target_state.varphi,
                            theta_A,
                            max_revolutions,
                            max_revolutions);

                        for (const auto& branch : branches) {
                            if (!branch.valid) {
                                continue;
                            }
                            const auto differentiated = problem1::attach_problem1_root_derivatives(
                                departure_planet,
                                target_planet,
                                departure_state.varphi,
                                target_state.varphi,
                                theta_A,
                                branch);
                            if (!differentiated.derivatives_available) {
                                continue;
                            }

                            for (double perturbation : {1e-4, -1e-4, 1e-3}) {
                                const auto refined = problem1::refine_problem1_root_branch_newton(
                                    departure_planet,
                                    target_planet,
                                    departure_state.varphi,
                                    target_state.varphi,
                                    theta_A,
                                    branch.transfer_revolution,
                                    branch.target_revolution,
                                    branch.encounter_global_angle + perturbation,
                                    40,
                                    newton_residual_tolerance_scale_free,
                                    1e-10);
                                const auto refined_seconds = problem1::refine_problem1_root_branch_newton_seconds(
                                    departure_planet,
                                    target_planet,
                                    departure_state.varphi,
                                    target_state.varphi,
                                    theta_A,
                                    branch.transfer_revolution,
                                    branch.target_revolution,
                                    branch.encounter_global_angle + perturbation,
                                    40,
                                    newton_residual_tolerance_seconds,
                                    1e-10);
                                if (!refined.valid) {
                                    continue;
                                }
                                assert(refined_seconds.valid == refined.valid);
                                assert(std::abs(
                                    common::normalize_angle_minus_pi_pi(
                                        refined_seconds.encounter_global_angle -
                                        refined.encounter_global_angle)) <= 1e-12);
                                assert(approx_equal(
                                    refined_seconds.time_of_flight_seconds,
                                    refined.time_of_flight_seconds,
                                    1e-12,
                                    1e-12));
                                assert(approx_equal(
                                    refined_seconds.residual_seconds,
                                    refined.residual_seconds,
                                    1e-12,
                                    1e-12));
                                assert(refined.transfer_revolution == branch.transfer_revolution);
                                assert(refined.target_revolution == branch.target_revolution);
                                assert(std::abs(
                                    common::normalize_angle_minus_pi_pi(
                                        refined.encounter_global_angle - branch.encounter_global_angle)) <= 1e-6);
                                assert(std::abs(refined.residual_seconds) <= 1e-5);
                                verified = true;
                                break;
                            }
                            if (verified) {
                                break;
                            }
                        }
                        if (verified) {
                            break;
                        }
                    }
                    if (verified) {
                        break;
                    }
                }
                if (verified) {
                    break;
                }
            }
            if (verified) {
                break;
            }
        }

        assert(verified);
    }

    {
        // 中文注释：最近邻一阶预测加 Newton refine 后，应与 query 点上的 exact solve 对齐。
        const problem1::Problem1RootTableConfig config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kTwoPi / 8.0,
            8,
            0.0,
            common::kTwoPi / 8.0,
            8,
            0.0,
            common::kTwoPi / 16.0,
            16,
            1,
            1,
            "problem1_root_table_draft_v0",
        };
        const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
        bool verified = false;

        for (const auto& cell : table.cells()) {
            if (cell.solutions_sorted_by_time_of_flight.empty()) {
                continue;
            }

            const auto& raw_node_branch = cell.solutions_sorted_by_time_of_flight.front();
            if (!raw_node_branch.valid) {
                continue;
            }
            const auto node_branch = problem1::attach_problem1_root_derivatives(
                config.departure_planet,
                config.target_planet,
                cell.nu_A_depart,
                cell.nu_B_depart,
                cell.theta_A,
                raw_node_branch);
            if (!node_branch.derivatives_available) {
                continue;
            }

            const double query_nu_A =
                common::normalize_angle_0_2pi(cell.nu_A_depart + 0.05 * config.nu_A_step);
            const double query_nu_B =
                common::normalize_angle_0_2pi(cell.nu_B_depart + 0.05 * config.nu_B_depart_step);
            const double query_theta_A =
                common::normalize_angle_0_2pi(cell.theta_A + 0.05 * config.theta_A_step);

            const auto prediction = problem1::predict_problem1_root_branch_linear_from_node(
                config.departure_planet,
                config.target_planet,
                cell.nu_A_depart,
                cell.nu_B_depart,
                cell.theta_A,
                node_branch,
                query_nu_A,
                query_nu_B,
                query_theta_A);
            if (!prediction.valid) {
                continue;
            }

            const auto refined = problem1::refine_problem1_root_branch_from_linear_prediction(
                config.departure_planet,
                config.target_planet,
                cell.nu_A_depart,
                cell.nu_B_depart,
                cell.theta_A,
                node_branch,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                20,
                newton_residual_tolerance_scale_free,
                1e-12);
            const auto refined_seconds =
                problem1::refine_problem1_root_branch_from_linear_prediction_seconds(
                    config.departure_planet,
                    config.target_planet,
                    cell.nu_A_depart,
                    cell.nu_B_depart,
                    cell.theta_A,
                    node_branch,
                    query_nu_A,
                    query_nu_B,
                    query_theta_A,
                    20,
                    newton_residual_tolerance_seconds,
                    1e-12);
            if (!refined.valid) {
                continue;
            }
            assert(refined_seconds.valid == refined.valid);
            assert(std::abs(
                common::normalize_angle_minus_pi_pi(
                    refined_seconds.encounter_global_angle -
                    refined.encounter_global_angle)) <= 1e-12);
            assert(approx_equal(
                refined_seconds.time_of_flight_seconds,
                refined.time_of_flight_seconds,
                1e-12,
                1e-12));
            assert(approx_equal(
                refined_seconds.residual_seconds,
                refined.residual_seconds,
                1e-12,
                1e-12));

            const auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
                config.departure_planet,
                config.target_planet,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                config.max_transfer_revolution,
                config.max_target_revolution);
            const auto exact_match = find_same_kq_nearest_branch(
                exact_branches,
                refined.transfer_revolution,
                refined.target_revolution,
                refined.encounter_global_angle);
            if (!exact_match.has_value()) {
                continue;
            }

            assert(std::abs(
                common::normalize_angle_minus_pi_pi(
                    refined.encounter_global_angle - exact_match->encounter_global_angle)) <= 1e-5);
            assert(approx_equal(
                refined.time_of_flight_seconds,
                exact_match->time_of_flight_seconds,
                1e-5,
                1e-5));
            assert(std::abs(refined.residual_seconds) <= 1e-6);
            verified = true;
            break;
        }

        assert(verified);
    }

    return 0;
}
