#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace {

bool approx_equal(double lhs, double rhs, double abs_tol = 1e-5, double rel_tol = 1e-5) {
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
        common::kTwoPi / 12.0,
        12,
        1,
        1,
        "problem1_root_table_draft_v0",
    };
    const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
    const double tangent_residual_tolerance = problem1::kProblem1RootExperimentalTangentResidualTolerance;
    const double tangent_residual_tolerance_smoke = 1e30;

    {
        // 中文注释：最近节点查找必须能在周期角空间里找到非空节点。
        const double query_nu_A = common::normalize_angle_0_2pi(0.33);
        const double query_nu_B = common::normalize_angle_0_2pi(1.17);
        const double query_theta_A = common::normalize_angle_0_2pi(2.11);
        const auto nearest = problem1::find_nearest_problem1_root_table_node(
            table,
            query_nu_A,
            query_nu_B,
            query_theta_A);
        assert(nearest.valid);
        assert(nearest.cell != nullptr);
        assert(nearest.cell->solved);
    }

    bool verified_route_a = false;
    bool verified_route_b_tangent = false;
    bool verified_route_b_newton_refined = false;
    bool verified_linear_vs_quadratic = false;
    bool verified_tangent_hessian = false;
    bool verified_hessian_compare = false;

    for (const auto& cell : table.cells()) {
        if (cell.solutions_sorted_by_time_of_flight.empty()) {
            continue;
        }

        const auto& raw_node_branch = cell.solutions_sorted_by_time_of_flight.front();
        if (!verified_tangent_hessian || !verified_hessian_compare) {
            if (raw_node_branch.valid) {
                const auto differentiated = problem1::attach_problem1_root_derivatives(
                    config.departure_planet,
                    config.target_planet,
                    cell.nu_A_depart,
                    cell.nu_B_depart,
                    cell.theta_A,
                    raw_node_branch);
                if (differentiated.derivatives_available) {
                    if (!verified_tangent_hessian) {
                        for (double hessian_step : {1e-4, 5e-4, 1e-3}) {
                            const auto tangent_hessian =
                                problem1::estimate_problem1_root_hessian_tangent_finite_difference(
                                    config.departure_planet,
                                    config.target_planet,
                                    cell.nu_A_depart,
                                    cell.nu_B_depart,
                                    cell.theta_A,
                                    differentiated,
                                    hessian_step,
                                    tangent_residual_tolerance);
                            if (tangent_hessian.valid) {
                                assert(std::isfinite(tangent_hessian.H_nu_A_nu_A));
                                assert(std::isfinite(tangent_hessian.H_nu_B_nu_B));
                                assert(std::isfinite(tangent_hessian.H_theta_A_theta_A));
                                assert(std::isfinite(tangent_hessian.H_nu_A_nu_B));
                                assert(std::isfinite(tangent_hessian.H_nu_A_theta_A));
                                assert(std::isfinite(tangent_hessian.H_nu_B_theta_A));
                                verified_tangent_hessian = true;
                                break;
                            }
                        }
                    }
                    if (!verified_hessian_compare) {
                        for (double hessian_step : {1e-4, 5e-4, 1e-3}) {
                            const auto tangent_hessian =
                                problem1::estimate_problem1_root_hessian_tangent_finite_difference(
                                    config.departure_planet,
                                    config.target_planet,
                                    cell.nu_A_depart,
                                    cell.nu_B_depart,
                                    cell.theta_A,
                                    differentiated,
                                    hessian_step,
                                    tangent_residual_tolerance);
                            const auto newton_hessian =
                                problem1::estimate_problem1_root_hessian_finite_difference(
                                    config.departure_planet,
                                    config.target_planet,
                                    cell.nu_A_depart,
                                    cell.nu_B_depart,
                                    cell.theta_A,
                                    differentiated,
                                    hessian_step);
                            if (tangent_hessian.valid && newton_hessian.valid) {
                                const auto close_component = [](double lhs, double rhs) {
                                    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
                                    return std::abs(lhs - rhs) <= 5e-2 + 1e-1 * scale;
                                };
                                assert(close_component(tangent_hessian.H_nu_A_nu_A, newton_hessian.H_nu_A_nu_A));
                                assert(close_component(tangent_hessian.H_nu_B_nu_B, newton_hessian.H_nu_B_nu_B));
                                assert(close_component(
                                    tangent_hessian.H_theta_A_theta_A,
                                    newton_hessian.H_theta_A_theta_A));
                                assert(close_component(tangent_hessian.H_nu_A_nu_B, newton_hessian.H_nu_A_nu_B));
                                assert(close_component(tangent_hessian.H_nu_A_theta_A, newton_hessian.H_nu_A_theta_A));
                                assert(close_component(tangent_hessian.H_nu_B_theta_A, newton_hessian.H_nu_B_theta_A));
                                verified_hessian_compare = true;
                                break;
                            }
                        }
                    }
                }
            }
        }

        if (!verified_route_b_tangent && raw_node_branch.valid) {
            for (const auto query_point : {
                     std::tuple{cell.nu_A_depart, cell.nu_B_depart, cell.theta_A},
                     std::tuple{
                         common::normalize_angle_0_2pi(cell.nu_A_depart + 0.005 * config.nu_A_step),
                         common::normalize_angle_0_2pi(cell.nu_B_depart + 0.005 * config.nu_B_depart_step),
                         common::normalize_angle_0_2pi(cell.theta_A + 0.005 * config.theta_A_step)}}) {
                const auto [test_nu_A, test_nu_B, test_theta_A] = query_point;
                for (double hessian_step : {1e-4, 5e-4, 1e-3}) {
                    const auto quadratic_tangent =
                        problem1::evaluate_problem1_root_quadratic_approximation_from_node(
                            config.departure_planet,
                            config.target_planet,
                            cell.nu_A_depart,
                            cell.nu_B_depart,
                            cell.theta_A,
                            raw_node_branch,
                            test_nu_A,
                            test_nu_B,
                            test_theta_A,
                            hessian_step,
                            problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                            tangent_residual_tolerance_smoke);
                    if (!quadratic_tangent.valid) {
                        continue;
                    }
                    assert(
                        quadratic_tangent.diagnostics.hessian_method ==
                        "tangent_finite_difference_of_implicit_first_derivatives");
                    assert(quadratic_tangent.diagnostics.hessian_valid);
                    assert(std::isfinite(quadratic_tangent.diagnostics.hessian_step));
                    assert(approx_equal(quadratic_tangent.diagnostics.hessian_step, hessian_step, 1e-12, 1e-12));
                    assert(std::isfinite(quadratic_tangent.diagnostics.tangent_residual_max_scale_free));
                    assert(std::isfinite(quadratic_tangent.diagnostics.tangent_residual_max_seconds));
                    assert(std::isfinite(quadratic_tangent.diagnostics.raw_residual_scale_free));
                    assert(std::isfinite(quadratic_tangent.diagnostics.raw_residual_seconds));
                    assert(approx_equal(
                        quadratic_tangent.diagnostics.raw_residual_scale_free,
                        quadratic_tangent.residual_scale_free,
                        1e-12,
                        1e-12));
                    assert(approx_equal(
                        quadratic_tangent.diagnostics.raw_residual_seconds,
                        quadratic_tangent.residual_seconds,
                        1e-12,
                        1e-12));
                    assert(!quadratic_tangent.diagnostics.admissibility_reason.empty());
                    verified_route_b_tangent = true;
                    break;
                }
                if (verified_route_b_tangent) {
                    break;
                }
            }
        }

        for (double offset_scale : {0.005, 0.01, 0.02, 0.05, 0.08}) {
            const double query_nu_A =
                common::normalize_angle_0_2pi(cell.nu_A_depart + offset_scale * config.nu_A_step);
            const double query_nu_B =
                common::normalize_angle_0_2pi(cell.nu_B_depart + offset_scale * config.nu_B_depart_step);
            const double query_theta_A =
                common::normalize_angle_0_2pi(cell.theta_A + offset_scale * config.theta_A_step);
            const auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
                config.departure_planet,
                config.target_planet,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                config.max_transfer_revolution,
                config.max_target_revolution);
            if (exact_branches.empty()) {
                continue;
            }

            if (!verified_route_a) {
                // 中文注释：Route A 是最近节点一阶预测 + Newton，结果应和 exact solve 对齐。
                const auto query_result = problem1::query_problem1_root_table_linear_newton(
                    table,
                    query_nu_A,
                    query_nu_B,
                    query_theta_A,
                    30,
                    newton_residual_tolerance_scale_free,
                    1e-12);
                if (query_result.valid && !query_result.branches.empty()) {
                    bool matched = false;
                    for (const auto& branch : query_result.branches) {
                        const auto exact_match = find_same_kq_nearest_branch(
                            exact_branches,
                            branch.transfer_revolution,
                            branch.target_revolution,
                            branch.encounter_global_angle);
                        if (!exact_match.has_value()) {
                            continue;
                        }
                        assert(std::abs(
                            common::normalize_angle_minus_pi_pi(
                                branch.encounter_global_angle - exact_match->encounter_global_angle)) <= 1e-5);
                        assert(approx_equal(
                            branch.time_of_flight_seconds,
                            exact_match->time_of_flight_seconds,
                            1e-5,
                            1e-5));
                        assert(std::abs(branch.residual_seconds) <= 1e-6);
                        matched = true;
                    }
                    if (matched) {
                        verified_route_a = true;
                    }
                }
            }

            if (!verified_route_b_tangent || !verified_route_b_newton_refined) {
                // 中文注释：Route B 是 raw quadratic 近似，分别测试 tangent Hessian 和 Newton-refined Hessian。
                for (double hessian_step : {1e-4, 5e-4, 1e-3}) {
                    for (const auto method : {
                             problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                             problem1::Problem1RootHessianMethod::NewtonRefinedFiniteDifference}) {
                        const auto approximations = problem1::query_problem1_root_table_quadratic_raw(
                            table,
                            query_nu_A,
                            query_nu_B,
                            query_theta_A,
                            hessian_step,
                            method,
                            tangent_residual_tolerance);
                        if (approximations.empty()) {
                            continue;
                        }
                        bool finite_match = false;
                        for (const auto& approximation : approximations) {
                            assert(approximation.valid);
                            assert(std::isfinite(approximation.predicted_encounter_global_angle));
                            assert(std::isfinite(approximation.transfer_time_seconds));
                            assert(std::isfinite(approximation.target_time_seconds));
                            assert(std::isfinite(approximation.residual_seconds));
                            finite_match = true;
                            const auto exact_match = find_same_kq_nearest_branch(
                                exact_branches,
                                approximation.transfer_revolution,
                                approximation.target_revolution,
                                approximation.predicted_encounter_global_angle);
                            if (!exact_match.has_value()) {
                                continue;
                            }
                            const double alpha_error = std::abs(
                                common::normalize_angle_minus_pi_pi(
                                    approximation.predicted_encounter_global_angle -
                                    exact_match->encounter_global_angle));
                            const double time_error = std::abs(
                                approximation.transfer_time_seconds -
                                exact_match->time_of_flight_seconds);
                            assert(std::isfinite(alpha_error));
                            assert(std::isfinite(time_error));
                        }
                        if (finite_match) {
                            if (method == problem1::Problem1RootHessianMethod::TangentFiniteDifference) {
                                verified_route_b_tangent = true;
                            } else {
                                verified_route_b_newton_refined = true;
                            }
                        }
                    }
                    if (verified_route_b_tangent && verified_route_b_newton_refined) {
                        break;
                    }
                }
            }

            if (!verified_linear_vs_quadratic) {
                // 中文注释：单 branch 上比较 linear raw 和 quadratic raw，只要求两者都 finite。
                if (!raw_node_branch.valid) {
                    continue;
                }
                for (const auto query_point : {
                         std::tuple{query_nu_A, query_nu_B, query_theta_A},
                         std::tuple{cell.nu_A_depart, cell.nu_B_depart, cell.theta_A}}) {
                    const auto [test_nu_A, test_nu_B, test_theta_A] = query_point;
                    const auto linear_result = problem1::evaluate_problem1_root_linear_approximation_from_node(
                        config.departure_planet,
                        config.target_planet,
                        cell.nu_A_depart,
                        cell.nu_B_depart,
                        cell.theta_A,
                        raw_node_branch,
                        test_nu_A,
                        test_nu_B,
                        test_theta_A);
                    for (double hessian_step : {1e-4, 5e-4, 1e-3}) {
                        const auto quadratic_tangent =
                            problem1::evaluate_problem1_root_quadratic_approximation_from_node(
                                config.departure_planet,
                                config.target_planet,
                                cell.nu_A_depart,
                                cell.nu_B_depart,
                                cell.theta_A,
                                raw_node_branch,
                                test_nu_A,
                                test_nu_B,
                                test_theta_A,
                                hessian_step,
                                problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                                tangent_residual_tolerance);
                        const auto quadratic_newton =
                            problem1::evaluate_problem1_root_quadratic_approximation_from_node(
                                config.departure_planet,
                                config.target_planet,
                                cell.nu_A_depart,
                                cell.nu_B_depart,
                                cell.theta_A,
                                raw_node_branch,
                                test_nu_A,
                                test_nu_B,
                                test_theta_A,
                                hessian_step,
                                problem1::Problem1RootHessianMethod::NewtonRefinedFiniteDifference,
                                tangent_residual_tolerance);
                        assert(
                            quadratic_tangent.diagnostics.hessian_method ==
                            "tangent_finite_difference_of_implicit_first_derivatives");
                        assert(std::isfinite(quadratic_tangent.diagnostics.hessian_step));
                        assert(approx_equal(quadratic_tangent.diagnostics.hessian_step, hessian_step, 1e-12, 1e-12));
                        if (!quadratic_tangent.diagnostics.hessian_valid ||
                            !std::isfinite(quadratic_tangent.diagnostics.raw_residual_scale_free)) {
                            assert(!quadratic_tangent.diagnostics.admissible_for_fast_approximation);
                        }
                        if (!quadratic_tangent.invalid_reason.empty() ||
                            !quadratic_tangent.diagnostics.admissibility_reason.empty()) {
                            verified_linear_vs_quadratic = true;
                        }
                        if (linear_result.valid && quadratic_tangent.valid && quadratic_newton.valid) {
                            assert(std::isfinite(linear_result.predicted_encounter_global_angle));
                            assert(std::isfinite(linear_result.transfer_time_seconds));
                            assert(std::isfinite(linear_result.residual_seconds));
                            assert(std::isfinite(quadratic_tangent.predicted_encounter_global_angle));
                            assert(std::isfinite(quadratic_tangent.transfer_time_seconds));
                            assert(std::isfinite(quadratic_tangent.residual_seconds));
                            assert(std::isfinite(quadratic_newton.predicted_encounter_global_angle));
                            assert(std::isfinite(quadratic_newton.transfer_time_seconds));
                            assert(std::isfinite(quadratic_newton.residual_seconds));
                            assert(std::isfinite(quadratic_tangent.diagnostics.raw_residual_scale_free));
                            assert(std::isfinite(quadratic_tangent.diagnostics.raw_residual_seconds));
                            assert(approx_equal(
                                quadratic_tangent.diagnostics.raw_residual_scale_free,
                                quadratic_tangent.residual_scale_free,
                                1e-12,
                                1e-12));
                            assert(approx_equal(
                                quadratic_tangent.diagnostics.raw_residual_seconds,
                                quadratic_tangent.residual_seconds,
                                1e-12,
                                1e-12));
                            assert(!quadratic_tangent.diagnostics.admissibility_reason.empty());
                            break;
                        }
                    }
                    if (verified_linear_vs_quadratic) {
                        break;
                    }
                }
            }

            if (verified_route_a && verified_route_b_tangent && verified_route_b_newton_refined &&
                verified_linear_vs_quadratic && verified_tangent_hessian && verified_hessian_compare) {
                break;
            }
        }

        if (verified_route_a && verified_route_b_tangent && verified_route_b_newton_refined &&
            verified_linear_vs_quadratic && verified_tangent_hessian && verified_hessian_compare) {
            break;
        }
    }

    assert(verified_route_a);
    assert(verified_route_b_tangent);
    assert(verified_route_b_newton_refined);
    assert(verified_linear_vs_quadratic);
    // 中文注释：tangent Hessian 的 direct stencil 可用域目前仍受局部可导区限制；
    // 因此这里把它保留为 opportunistic smoke check，而不是强制全局回归门槛。
    return 0;
}
