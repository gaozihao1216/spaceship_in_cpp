/*
 * 文件作用：测试 Problem 2 Route A 线性外推与 Newton 精修。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_route_a.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

namespace {

using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;
using spaceship_cpp::problem2::Problem2RouteANewtonStatus;
using spaceship_cpp::problem2::Problem2ThetaPrimeScanConfig;

std::optional<Problem2OutgoingBranchSolution> find_solution_with_derivatives(
    const spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult& scan,
    std::size_t& out_node_index
) {
    for (std::size_t node_index = 1; node_index + 1 < scan.nodes.size(); ++node_index) {
        for (const auto& solution : scan.nodes[node_index].solutions) {
            if (solution.has_dphi_dtheta_prime && solution.has_de_dtheta_prime) {
                out_node_index = node_index;
                return solution;
            }
        }
    }
    return std::nullopt;
}

bool approx_equal(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;
    namespace problem2 = spaceship_cpp::problem2;

    const auto& defaults = config::global_config();
    const auto scan_config = config::make_problem2_theta_prime_scan_config(
        PlanetId::Earth,
        PlanetId::Mars,
        0.0,
        defaults.problem2_theta_prime_scan,
        defaults.problem1_solve);
    const auto newton_options = config::make_problem2_route_a_newton_options(
        defaults.problem2_route_a_newton,
        defaults.problem1_solve);

    const auto scan = problem2::run_problem2_theta_prime_initial_scan(scan_config);
    assert(scan.ok);
    assert(scan.nodes.size() >= 3U);

    std::size_t endpoint_index = 0;
    const auto endpoint_solution = find_solution_with_derivatives(scan, endpoint_index);
    assert(endpoint_solution.has_value());

    const double endpoint_theta = scan.nodes[endpoint_index].theta_prime_local;
    const double coarse_step =
        scan.nodes[endpoint_index + 1].theta_prime_local - endpoint_theta;
    const double target_theta = endpoint_theta + 0.5 * coarse_step;

    {
        const auto linear = problem2::estimate_problem2_route_a_solution_linear(
            endpoint_theta,
            *endpoint_solution,
            target_theta);
        assert(linear.valid);
        assert(std::isfinite(linear.predicted_phi));
        assert(std::isfinite(linear.predicted_eccentricity));
        assert(approx_equal(linear.delta_theta_prime, target_theta - endpoint_theta, 1e-15));

        const auto exact_target = problem2::evaluate_problem2_theta_prime_node(scan_config, target_theta);
        const auto match_index = problem2::find_best_matching_outgoing_branch_index(
            *endpoint_solution,
            exact_target.solutions,
            scan_config.branch_phi_pairing_max_gap);
        assert(match_index.has_value());
        const auto& exact_solution = exact_target.solutions[*match_index];

        double second_phi_derivative = 0.0;
        bool has_second_phi_derivative = false;
        if (endpoint_index > 0U) {
            const auto left_match = problem2::find_best_matching_outgoing_branch_index(
                *endpoint_solution,
                scan.nodes[endpoint_index - 1].solutions,
                scan_config.branch_phi_pairing_max_gap);
            const auto right_match = problem2::find_best_matching_outgoing_branch_index(
                *endpoint_solution,
                scan.nodes[endpoint_index + 1].solutions,
                scan_config.branch_phi_pairing_max_gap);
            if (left_match.has_value() && right_match.has_value()) {
                const double left_theta = scan.nodes[endpoint_index - 1].theta_prime_local;
                const double right_theta = scan.nodes[endpoint_index + 1].theta_prime_local;
                const double left_phi =
                    scan.nodes[endpoint_index - 1].solutions[*left_match].encounter_global_angle;
                const double right_phi =
                    scan.nodes[endpoint_index + 1].solutions[*right_match].encounter_global_angle;
                second_phi_derivative = (
                    right_phi - 2.0 * endpoint_solution->encounter_global_angle + left_phi) /
                    (coarse_step * coarse_step);
                has_second_phi_derivative = std::isfinite(second_phi_derivative);
            }
        }

        const double phi_error = std::abs(
            spaceship_cpp::common::normalize_angle_minus_pi_pi(
                linear.predicted_phi - exact_solution.encounter_global_angle));
        const double e_error = std::abs(linear.predicted_eccentricity - exact_solution.outgoing_eccentricity);
        const double delta_theta = linear.delta_theta_prime;
        const double linear_phi_bound = has_second_phi_derivative
            ? (0.5 * std::abs(second_phi_derivative) * delta_theta * delta_theta + 1e-9)
            : (0.25 * coarse_step * coarse_step + 1e-6);
        const double linear_e_bound = 0.5 * std::abs(endpoint_solution->de_dtheta_prime) *
            delta_theta * delta_theta + 1e-6;
        assert(phi_error <= linear_phi_bound);
        assert(e_error <= linear_e_bound);
    }

    {
        const auto linear = problem2::estimate_problem2_route_a_solution_linear(
            endpoint_theta,
            *endpoint_solution,
            target_theta);
        assert(linear.valid);

        const auto newton = problem2::refine_problem2_route_a_phi_by_newton(
            scan_config,
            target_theta,
            linear.predicted_phi,
            endpoint_solution->transfer_revolution,
            endpoint_solution->target_revolution,
            newton_options);
        assert(newton.ok);
        assert(newton.status == Problem2RouteANewtonStatus::Converged);
        assert(newton.iterations_used > 0);
        assert(newton.solution.relative_problem1_residual <= newton_options.max_relative_residual);

        const auto exact_target = problem2::evaluate_problem2_theta_prime_node(scan_config, target_theta);
        const auto match_index = problem2::find_best_matching_outgoing_branch_index(
            *endpoint_solution,
            exact_target.solutions,
            scan_config.branch_phi_pairing_max_gap);
        assert(match_index.has_value());
        const auto& exact_solution = exact_target.solutions[*match_index];

        assert(approx_equal(
            spaceship_cpp::common::normalize_angle_minus_pi_pi(
                newton.final_phi - exact_solution.encounter_global_angle),
            0.0,
            1e-6));
        assert(approx_equal(
            newton.solution.outgoing_eccentricity,
            exact_solution.outgoing_eccentricity,
            1e-6));
    }

    {
        Problem2OutgoingBranchSolution missing_derivative = *endpoint_solution;
        missing_derivative.has_dphi_dtheta_prime = false;
        const auto linear = problem2::estimate_problem2_route_a_solution_linear(
            endpoint_theta,
            missing_derivative,
            target_theta);
        assert(!linear.valid);
        assert(linear.invalid_reason == "missing_dphi_dtheta_prime");
    }

    std::cout << "problem2_route_a_ok\n";
    return 0;
}
