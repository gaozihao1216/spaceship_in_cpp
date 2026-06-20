/*
 * 文件作用：实现 Problem 2 Route A 线性外推与 Newton 精修。
 * 主要工作：先线性预测 φ/e，再以 Newton 将 P1 残差收敛到零点。
 */
#include "spaceship_cpp/problem2/problem2_theta_prime_route_a.hpp"

#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <limits>

namespace spaceship_cpp::problem2 {
namespace {

using spaceship_cpp::bfs::problem2_local_periapsis_angle_to_global;
using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::near_zero;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::problem1::Problem1ResidualInput;
using spaceship_cpp::problem1::Problem1ResidualStatus;
using spaceship_cpp::problem1::estimate_problem1_residual_derivative_central;

double quiet_nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

Problem1ResidualInput make_route_a_residual_input(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_global,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    return Problem1ResidualInput{
        config.flyby_planet,
        config.target_planet,
        config.flyby_time_seconds_since_j2000,
        target_theta_prime_global,
        encounter_global_angle,
        transfer_revolution,
        target_revolution,
    };
}

Problem2OutgoingBranchSolution make_branch_solution_from_residual(
    const problem1::Problem1ResidualResult& residual_result,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution,
    double relative_residual
) {
    Problem2OutgoingBranchSolution solution{};
    solution.transfer_revolution = transfer_revolution;
    solution.target_revolution = target_revolution;
    solution.encounter_global_angle = encounter_global_angle;
    solution.outgoing_eccentricity = residual_result.transfer_e;
    solution.outgoing_semi_latus_rectum = residual_result.transfer_p;
    solution.relative_problem1_residual = relative_residual;
    return solution;
}

double compute_relative_residual(const problem1::Problem1ResidualResult& result) {
    const double scale = std::max(
        std::abs(result.transfer_time_scale_free),
        std::abs(result.target_time_scale_free));
    if (!is_finite(scale) || !(scale > 0.0)) {
        return std::numeric_limits<double>::infinity();
    }
    return std::abs(result.residual) / scale;
}

std::optional<problem1::Problem1ResidualResult> evaluate_route_a_residual_if_success(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_global,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    const auto result = problem1::evaluate_problem1_residual(make_route_a_residual_input(
        config,
        target_theta_prime_global,
        encounter_global_angle,
        transfer_revolution,
        target_revolution));
    if (result.status != Problem1ResidualStatus::Success || !is_finite(result.residual)) {
        return std::nullopt;
    }
    return result;
}

}  // namespace

Problem2RouteALinearEstimateResult estimate_problem2_route_a_solution_linear(
    double endpoint_theta_prime_local,
    const Problem2OutgoingBranchSolution& endpoint_solution,
    double target_theta_prime_local
) {
    Problem2RouteALinearEstimateResult result{};
    result.endpoint_theta_prime_local = endpoint_theta_prime_local;
    result.target_theta_prime_local = target_theta_prime_local;

    if (!is_finite(endpoint_theta_prime_local) || !is_finite(target_theta_prime_local)) {
        result.invalid_reason = "non_finite_theta_prime";
        return result;
    }
    if (!endpoint_solution.has_dphi_dtheta_prime) {
        result.invalid_reason = "missing_dphi_dtheta_prime";
        return result;
    }
    if (!is_finite(endpoint_solution.dphi_dtheta_prime)) {
        result.invalid_reason = "non_finite_dphi_dtheta_prime";
        return result;
    }

    result.delta_theta_prime = target_theta_prime_local - endpoint_theta_prime_local;
    result.predicted_phi = normalize_angle_0_2pi(
        endpoint_solution.encounter_global_angle +
        endpoint_solution.dphi_dtheta_prime * result.delta_theta_prime);

    if (endpoint_solution.has_de_dtheta_prime && is_finite(endpoint_solution.de_dtheta_prime)) {
        result.predicted_eccentricity =
            endpoint_solution.outgoing_eccentricity +
            endpoint_solution.de_dtheta_prime * result.delta_theta_prime;
    } else {
        result.predicted_eccentricity = endpoint_solution.outgoing_eccentricity;
    }

    result.predicted_semi_latus_rectum = endpoint_solution.outgoing_semi_latus_rectum;
    if (!is_finite(result.predicted_phi)) {
        result.invalid_reason = "non_finite_predicted_phi";
        return result;
    }

    result.valid = true;
    return result;
}

Problem2RouteANewtonResult refine_problem2_route_a_phi_by_newton(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_local,
    double initial_phi,
    int transfer_revolution,
    int target_revolution,
    const Problem2RouteANewtonOptions& options
) {
    Problem2RouteANewtonResult result{};
    result.target_theta_prime_local = target_theta_prime_local;
    result.final_phi = initial_phi;

    if (!is_finite(target_theta_prime_local) || !is_finite(initial_phi)) {
        result.status = Problem2RouteANewtonStatus::DivergedInvalidResidual;
        result.status_reason = "non_finite_newton_input";
        return result;
    }
    if (options.max_iterations <= 0 || !is_finite(options.phi_tolerance) || !(options.phi_tolerance > 0.0) ||
        !is_finite(options.phi_derivative_step) || !(options.phi_derivative_step > 0.0) ||
        !is_finite(options.max_relative_residual) || !(options.max_relative_residual > 0.0)) {
        result.status = Problem2RouteANewtonStatus::DivergedInvalidResidual;
        result.status_reason = "invalid_newton_options";
        return result;
    }

    result.target_theta_prime_global =
        problem2_local_periapsis_angle_to_global(config.flyby_planet, target_theta_prime_local);

    double phi = normalize_angle_0_2pi(initial_phi);
    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
        const auto center_result = evaluate_route_a_residual_if_success(
            config,
            result.target_theta_prime_global,
            phi,
            transfer_revolution,
            target_revolution);
        if (!center_result.has_value()) {
            result.status = Problem2RouteANewtonStatus::DivergedInvalidResidual;
            result.status_reason = "invalid_residual_evaluation";
            result.iterations_used = iteration;
            return result;
        }

        const double residual = center_result->residual;
        const double abs_residual = std::abs(residual);
        const double relative_residual = compute_relative_residual(*center_result);
        if ((options.residual_tolerance > 0.0 && abs_residual <= options.residual_tolerance) ||
            relative_residual <= options.max_relative_residual) {
            result.ok = true;
            result.status = Problem2RouteANewtonStatus::Converged;
            result.iterations_used = iteration;
            result.final_phi = phi;
            result.residual_result = *center_result;
            result.solution = make_branch_solution_from_residual(
                *center_result,
                phi,
                transfer_revolution,
                target_revolution,
                compute_relative_residual(*center_result));
            return result;
        }

        const double derivative_step = options.phi_derivative_step;
        const auto minus_result = evaluate_route_a_residual_if_success(
            config,
            result.target_theta_prime_global,
            normalize_angle_0_2pi(phi - derivative_step),
            transfer_revolution,
            target_revolution);
        const auto plus_result = evaluate_route_a_residual_if_success(
            config,
            result.target_theta_prime_global,
            normalize_angle_0_2pi(phi + derivative_step),
            transfer_revolution,
            target_revolution);
        if (!minus_result.has_value() || !plus_result.has_value()) {
            result.status = Problem2RouteANewtonStatus::DivergedInvalidResidual;
            result.status_reason = "invalid_residual_derivative_support";
            result.iterations_used = iteration;
            return result;
        }

        const double derivative = estimate_problem1_residual_derivative_central(
            phi - derivative_step,
            minus_result->residual,
            phi + derivative_step,
            plus_result->residual);
        if (!is_finite(derivative) || near_zero(derivative, kDefaultEpsilon)) {
            result.status = Problem2RouteANewtonStatus::DivergedSingularDerivative;
            result.status_reason = "singular_residual_derivative";
            result.iterations_used = iteration;
            return result;
        }

        const double phi_next = normalize_angle_0_2pi(phi - residual / derivative);
        const auto next_result = evaluate_route_a_residual_if_success(
            config,
            result.target_theta_prime_global,
            phi_next,
            transfer_revolution,
            target_revolution);
        if (!next_result.has_value()) {
            result.status = Problem2RouteANewtonStatus::DivergedInvalidResidual;
            result.status_reason = "invalid_residual_after_newton_step";
            result.iterations_used = iteration + 1;
            return result;
        }

        if (options.reject_on_residual_increase &&
            std::abs(next_result->residual) > abs_residual) {
            result.status = Problem2RouteANewtonStatus::DivergedResidualIncreased;
            result.status_reason = "residual_increased_after_newton_step";
            result.iterations_used = iteration + 1;
            result.final_phi = phi_next;
            result.residual_result = *next_result;
            return result;
        }

        if (std::abs(normalize_angle_minus_pi_pi(phi_next - phi)) <= options.phi_tolerance &&
            compute_relative_residual(*next_result) <= options.max_relative_residual) {
            result.ok = true;
            result.status = Problem2RouteANewtonStatus::Converged;
            result.iterations_used = iteration + 1;
            result.final_phi = phi_next;
            result.residual_result = *next_result;
            result.solution = make_branch_solution_from_residual(
                *next_result,
                phi_next,
                transfer_revolution,
                target_revolution,
                compute_relative_residual(*next_result));
            return result;
        }

        phi = phi_next;
        result.final_phi = phi;
        result.residual_result = *next_result;
    }

    result.status = Problem2RouteANewtonStatus::DivergedMaxIterations;
    result.status_reason = "max_newton_iterations_exceeded";
    result.iterations_used = options.max_iterations;
    return result;
}

}  // namespace spaceship_cpp::problem2
