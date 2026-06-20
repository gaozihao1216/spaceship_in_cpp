/*
 * 文件作用：实现 Problem 2 Route A 线性外推与 Newton 精修。
 * 主要工作：先线性预测 φ/e，再以 Newton 将 P1 残差收敛到零点。
 */
#include "spaceship_cpp/problem2/problem2_theta_prime_route_a.hpp"

#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <cmath>
#include <limits>

namespace spaceship_cpp::problem2 {
namespace {

using spaceship_cpp::bfs::problem2_local_periapsis_angle_to_global;
using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::near_zero;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::planet_params::planet_state_at_time;
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
    result.predicted_phi_unwrapped =
        endpoint_solution.encounter_global_angle +
        endpoint_solution.dphi_dtheta_prime * result.delta_theta_prime;
    result.predicted_phi = normalize_angle_0_2pi(result.predicted_phi_unwrapped);
    result.original_target_revolution = endpoint_solution.target_revolution;
    result.adjusted_target_revolution = endpoint_solution.target_revolution;

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

Problem2RouteATargetRevolutionAdjustment adjust_target_revolution_for_route_a_linear_prediction(
    planet_params::PlanetId flyby_planet,
    double flyby_time_seconds_since_j2000,
    const Problem2OutgoingBranchSolution& endpoint_solution,
    double delta_theta_prime,
    int max_target_revolution
) {
    Problem2RouteATargetRevolutionAdjustment adjustment{};
    adjustment.original_target_revolution = endpoint_solution.target_revolution;
    adjustment.adjusted_target_revolution = endpoint_solution.target_revolution;
    adjustment.endpoint_phi = endpoint_solution.encounter_global_angle;

    if (!is_finite(flyby_time_seconds_since_j2000) || !is_finite(delta_theta_prime)) {
        adjustment.invalid_reason = "non_finite_q_adjustment_input";
        return adjustment;
    }
    if (max_target_revolution < 0) {
        adjustment.invalid_reason = "invalid_max_target_revolution";
        return adjustment;
    }
    if (!endpoint_solution.has_dphi_dtheta_prime || !is_finite(endpoint_solution.dphi_dtheta_prime)) {
        adjustment.invalid_reason = "missing_dphi_dtheta_prime";
        return adjustment;
    }

    const auto flyby_state = planet_state_at_time(flyby_planet, flyby_time_seconds_since_j2000);
    adjustment.omega_J = flyby_state.theta_global;
    adjustment.predicted_phi_unwrapped =
        endpoint_solution.encounter_global_angle +
        endpoint_solution.dphi_dtheta_prime * delta_theta_prime;

    const int endpoint_winding = static_cast<int>(
        std::floor((adjustment.endpoint_phi - adjustment.omega_J) / kTwoPi));
    const int predicted_winding = static_cast<int>(
        std::floor((adjustment.predicted_phi_unwrapped - adjustment.omega_J) / kTwoPi));
    adjustment.target_revolution_delta = predicted_winding - endpoint_winding;

    int adjusted_q = endpoint_solution.target_revolution + adjustment.target_revolution_delta;
    adjusted_q = std::clamp(adjusted_q, 0, max_target_revolution);
    adjustment.adjusted_target_revolution = adjusted_q;
    adjustment.valid = true;
    return adjustment;
}

void apply_target_revolution_adjustment_to_route_a_linear_estimate(
    const Problem2RouteATargetRevolutionAdjustment& adjustment,
    Problem2RouteALinearEstimateResult& linear_estimate
) {
    if (!adjustment.valid) {
        return;
    }
    linear_estimate.original_target_revolution = adjustment.original_target_revolution;
    linear_estimate.adjusted_target_revolution = adjustment.adjusted_target_revolution;
    linear_estimate.target_revolution_delta = adjustment.target_revolution_delta;
    linear_estimate.target_revolution_adjusted =
        adjustment.target_revolution_delta != 0;
    linear_estimate.predicted_phi_unwrapped = adjustment.predicted_phi_unwrapped;
}

Problem2RouteANewtonResult refine_problem2_route_a_at_theta_prime(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_local,
    const Problem2OutgoingBranchSolution& linear_endpoint_branch,
    double linear_endpoint_theta_prime_local,
    int transfer_revolution,
    const Problem2RouteANewtonOptions& options
) {
    Problem2RouteANewtonResult result{};

    const Problem2RouteALinearEstimateResult linear = estimate_problem2_route_a_solution_linear(
        linear_endpoint_theta_prime_local,
        linear_endpoint_branch,
        target_theta_prime_local);
    if (!linear.valid) {
        result.status = Problem2RouteANewtonStatus::DivergedInvalidResidual;
        result.status_reason = linear.invalid_reason.empty()
            ? "route_a_linear_estimate_failed"
            : linear.invalid_reason;
        return result;
    }

    int target_revolution = linear_endpoint_branch.target_revolution;
    if (options.adjust_target_revolution_on_omega_J_crossing) {
        const Problem2RouteATargetRevolutionAdjustment adjustment =
            adjust_target_revolution_for_route_a_linear_prediction(
                config.flyby_planet,
                config.flyby_time_seconds_since_j2000,
                linear_endpoint_branch,
                linear.delta_theta_prime,
                config.problem1_solve.max_target_revolution);
        if (adjustment.valid) {
            target_revolution = adjustment.adjusted_target_revolution;
        }
    }

    return refine_problem2_route_a_phi_by_newton(
        config,
        target_theta_prime_local,
        linear.predicted_phi,
        transfer_revolution,
        target_revolution,
        options);
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
