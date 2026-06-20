/*
 * 文件作用：实现 Problem 2 统一 G(θ') 搜索管线。
 */
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"

#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace spaceship_cpp::problem2 {
namespace {

using spaceship_cpp::bfs::problem2_local_periapsis_angle_to_global;
using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::near_zero;
using spaceship_cpp::planet_params::get_planet_params;
using spaceship_cpp::planet_params::planet_state_at_time;

bool same_revolution_branch(
    const Problem2OutgoingBranchSolution& left,
    const Problem2OutgoingBranchSolution& right
) {
    return left.transfer_revolution == right.transfer_revolution &&
        left.target_revolution == right.target_revolution;
}

bool same_flyby_G_solution(
    const Problem2FlybyGSolution& left,
    const Problem2FlybyGSolution& right,
    double theta_prime_tolerance
) {
    return same_revolution_branch(left.outgoing_branch, right.outgoing_branch) &&
        is_finite(left.theta_prime_local) &&
        is_finite(right.theta_prime_local) &&
        std::abs(left.theta_prime_local - right.theta_prime_local) <= theta_prime_tolerance;
}

struct RouteAGEvaluation {
    Problem2RouteANewtonResult route_a{};
    Problem2FlybyGSample g_sample{};
};

std::optional<RouteAGEvaluation> evaluate_route_a_G_at_theta_prime(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2ThetaPrimeScanConfig& scan_config,
    const Problem2RouteANewtonOptions& route_a_options,
    double target_theta_prime_local,
    const Problem2OutgoingBranchSolution& reference_branch,
    const Problem2OutgoingBranchSolution& linear_endpoint_branch,
    double linear_endpoint_theta_prime_local
) {
    if (!linear_endpoint_branch.has_de_dtheta_prime ||
        !is_finite(linear_endpoint_branch.de_dtheta_prime) ||
        !linear_endpoint_branch.has_dphi_dtheta_prime ||
        !is_finite(linear_endpoint_branch.dphi_dtheta_prime)) {
        return std::nullopt;
    }

    const Problem2RouteANewtonResult route_a = refine_problem2_route_a_at_theta_prime(
        scan_config,
        target_theta_prime_local,
        linear_endpoint_branch,
        linear_endpoint_theta_prime_local,
        reference_branch.transfer_revolution,
        route_a_options);
    if (!route_a.ok) {
        return std::nullopt;
    }

    RouteAGEvaluation evaluation{};
    evaluation.route_a = route_a;
    evaluation.g_sample = evaluate_flyby_constraint_G_and_derivative(
        incoming_cache,
        route_a.solution.outgoing_eccentricity,
        target_theta_prime_local,
        linear_endpoint_branch.dphi_dtheta_prime,
        linear_endpoint_branch.de_dtheta_prime);
    if (!evaluation.g_sample.valid) {
        return std::nullopt;
    }
    return evaluation;
}

void append_unique_flyby_G_solution(
    std::vector<Problem2FlybyGSolution>& solutions,
    Problem2FlybyGSolution candidate,
    double theta_prime_tolerance
) {
    for (auto& existing : solutions) {
        if (same_flyby_G_solution(existing, candidate, theta_prime_tolerance)) {
            if (std::abs(candidate.G) < std::abs(existing.G)) {
                existing = std::move(candidate);
            }
            return;
        }
    }
    solutions.push_back(std::move(candidate));
}

const Problem2OutgoingBranchSolution* choose_linear_endpoint_for_theta_prime(
    double theta_prime_left,
    const Problem2OutgoingBranchSolution& left_branch,
    double theta_prime_right,
    const Problem2OutgoingBranchSolution& right_branch,
    double target_theta_prime_local
) {
    const bool left_usable =
        left_branch.has_dphi_dtheta_prime && is_finite(left_branch.dphi_dtheta_prime);
    const bool right_usable =
        right_branch.has_dphi_dtheta_prime && is_finite(right_branch.dphi_dtheta_prime);
    if (left_usable && !right_usable) {
        return &left_branch;
    }
    if (right_usable && !left_usable) {
        return &right_branch;
    }
    if (!left_usable && !right_usable) {
        return nullptr;
    }

    const double distance_to_left = std::abs(target_theta_prime_local - theta_prime_left);
    const double distance_to_right = std::abs(target_theta_prime_local - theta_prime_right);
    return distance_to_left <= distance_to_right ? &left_branch : &right_branch;
}

Problem2FlybyGSolution make_near_zero_solution(
    Problem2FlybyGSolutionSource source,
    double theta_prime_local,
    double G,
    const Problem2OutgoingBranchSolution& branch
) {
    Problem2FlybyGSolution solution{};
    solution.source = source;
    solution.theta_prime_local = theta_prime_local;
    solution.G = G;
    solution.outgoing_branch = branch;
    return solution;
}

bool interval_width_within_tolerance(
    double theta_prime_left,
    double theta_prime_right,
    double theta_prime_tolerance
) {
    return (theta_prime_right - theta_prime_left) <= theta_prime_tolerance;
}

bool should_discard_g_search_interval(
    double theta_prime_left,
    double theta_prime_right,
    double theta_prime_tolerance,
    int recursion_depth,
    int max_recursion_depth
) {
    return interval_width_within_tolerance(
               theta_prime_left,
               theta_prime_right,
               theta_prime_tolerance) ||
        recursion_depth >= max_recursion_depth;
}

bool attach_middle_branch_derivatives_from_dense_endpoint(
    double dense_theta_prime,
    const Problem2OutgoingBranchSolution& dense_branch,
    double middle_theta_prime,
    Problem2OutgoingBranchSolution& middle_branch
) {
    const double delta_theta = middle_theta_prime - dense_theta_prime;
    if (!is_finite(delta_theta) || delta_theta == 0.0) {
        return false;
    }

    const double phi_delta = spaceship_cpp::common::normalize_angle_minus_pi_pi(
        middle_branch.encounter_global_angle - dense_branch.encounter_global_angle);
    const double dphi = phi_delta / delta_theta;
    const double de =
        (middle_branch.outgoing_eccentricity - dense_branch.outgoing_eccentricity) / delta_theta;
    if (!is_finite(dphi) || !is_finite(de)) {
        return false;
    }

    middle_branch.dphi_dtheta_prime = dphi;
    middle_branch.de_dtheta_prime = de;
    middle_branch.has_dphi_dtheta_prime = true;
    middle_branch.has_de_dtheta_prime = true;
    return true;
}

struct Problem2CaseBMiddleProbeResult {
    bool ok = false;
    std::string status_reason;
    double theta_prime_middle = 0.0;
    std::vector<Problem2OutgoingBranchSolution> middle_branches;
};

Problem2CaseBMiddleProbeResult probe_case_b_middle_for_G_search(
    const Problem2ThetaPrimeScanConfig& config,
    const Problem2RouteANewtonOptions& options,
    double theta_prime_left,
    const std::vector<Problem2OutgoingBranchSolution>& left_branches,
    double theta_prime_right,
    const std::vector<Problem2OutgoingBranchSolution>& right_branches
) {
    Problem2CaseBMiddleProbeResult result{};
    if (!is_finite(theta_prime_left) || !is_finite(theta_prime_right)) {
        result.status_reason = "non_finite_theta_prime_interval";
        return result;
    }
    if (!(theta_prime_right > theta_prime_left)) {
        result.status_reason = "invalid_theta_prime_interval_order";
        return result;
    }
    if (classify_problem2_theta_prime_interval_case(left_branches.size(), right_branches.size()) !=
        Problem2ThetaPrimeIntervalCase::BranchCountDifferenceOne) {
        result.status_reason = "case_b_requires_branch_count_difference_one";
        return result;
    }

    const bool dense_is_left = left_branches.size() > right_branches.size();
    const auto& dense_branches = dense_is_left ? left_branches : right_branches;
    const double dense_theta_prime = dense_is_left ? theta_prime_left : theta_prime_right;

    result.theta_prime_middle = 0.5 * (theta_prime_left + theta_prime_right);
    result.middle_branches.reserve(dense_branches.size());

    for (const auto& dense_branch : dense_branches) {
        if (!dense_branch.has_dphi_dtheta_prime || !is_finite(dense_branch.dphi_dtheta_prime)) {
            continue;
        }

        const Problem2RouteANewtonResult newton = refine_problem2_route_a_at_theta_prime(
            config,
            result.theta_prime_middle,
            dense_branch,
            dense_theta_prime,
            dense_branch.transfer_revolution,
            options);
        if (!newton.ok) {
            continue;
        }

        Problem2OutgoingBranchSolution middle_branch = newton.solution;
        if (!attach_middle_branch_derivatives_from_dense_endpoint(
                dense_theta_prime,
                dense_branch,
                result.theta_prime_middle,
                middle_branch)) {
            continue;
        }
        result.middle_branches.push_back(std::move(middle_branch));
    }

    if (result.middle_branches.empty()) {
        result.status_reason = "case_b_middle_probe_produced_no_branches";
        return result;
    }

    result.ok = true;
    return result;
}

std::vector<Problem2FlybyGSolution> search_flyby_constraint_G_on_k_layer_interval(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2FlybyGSearchConfig& config,
    int transfer_revolution,
    double theta_prime_left,
    const std::vector<Problem2OutgoingBranchSolution>& left_branches,
    double theta_prime_right,
    const std::vector<Problem2OutgoingBranchSolution>& right_branches,
    int recursion_depth
) {
    std::vector<Problem2FlybyGSolution> solutions;
    if (!incoming_cache.valid) {
        return solutions;
    }
    if (!is_finite(theta_prime_left) || !is_finite(theta_prime_right) ||
        !(theta_prime_right > theta_prime_left)) {
        return solutions;
    }
    if (recursion_depth < 0 || recursion_depth > config.max_recursion_depth) {
        return solutions;
    }

    const Problem2ThetaPrimeIntervalCase interval_case = classify_problem2_theta_prime_interval_case(
        left_branches.size(),
        right_branches.size());

    if (interval_case == Problem2ThetaPrimeIntervalCase::EqualBranchCount) {
        const auto pairs = pair_outgoing_branch_solutions_by_phi(
            left_branches,
            right_branches,
            config.scan_config.branch_phi_pairing_max_gap);
        for (const auto& pair : pairs) {
            const Problem2FlybyGBranchIntervalInput interval{
                .theta_prime_left = theta_prime_left,
                .theta_prime_right = theta_prime_right,
                .left_branch = left_branches[pair.left_index],
                .right_branch = right_branches[pair.right_index],
            };
            const auto interval_solutions = process_flyby_constraint_G_branch_interval(
                incoming_cache,
                config,
                interval);
            for (const auto& candidate : interval_solutions) {
                append_unique_flyby_G_solution(
                    solutions,
                    candidate,
                    config.solution_theta_prime_tolerance);
            }
        }
        return solutions;
    }

    if (should_discard_g_search_interval(
            theta_prime_left,
            theta_prime_right,
            config.theta_prime_tolerance,
            recursion_depth,
            config.max_recursion_depth)) {
        return solutions;
    }

    double theta_prime_middle = 0.0;
    std::vector<Problem2OutgoingBranchSolution> middle_branches;

    if (interval_case == Problem2ThetaPrimeIntervalCase::BranchCountDifferenceOne) {
        const auto probe = probe_case_b_middle_for_G_search(
            config.scan_config,
            config.route_a_newton_options,
            theta_prime_left,
            left_branches,
            theta_prime_right,
            right_branches);
        if (!probe.ok) {
            return solutions;
        }
        theta_prime_middle = probe.theta_prime_middle;
        middle_branches = probe.middle_branches;
    } else {
        const auto case_c = solve_case_c_middle_branches_on_k_layer(
            config.scan_config,
            transfer_revolution,
            theta_prime_left,
            left_branches,
            theta_prime_right,
            right_branches);
        if (!case_c.ok) {
            return solutions;
        }
        theta_prime_middle = case_c.theta_prime_middle;
        middle_branches = case_c.middle_branches;
    }

    const auto left_solutions = search_flyby_constraint_G_on_k_layer_interval(
        incoming_cache,
        config,
        transfer_revolution,
        theta_prime_left,
        left_branches,
        theta_prime_middle,
        middle_branches,
        recursion_depth + 1);
    for (const auto& candidate : left_solutions) {
        append_unique_flyby_G_solution(solutions, candidate, config.solution_theta_prime_tolerance);
    }

    const auto right_solutions = search_flyby_constraint_G_on_k_layer_interval(
        incoming_cache,
        config,
        transfer_revolution,
        theta_prime_middle,
        middle_branches,
        theta_prime_right,
        right_branches,
        recursion_depth + 1);
    for (const auto& candidate : right_solutions) {
        append_unique_flyby_G_solution(solutions, candidate, config.solution_theta_prime_tolerance);
    }
    return solutions;
}

}  // namespace

using spaceship_cpp::common::is_finite;
using spaceship_cpp::bfs::problem2_local_periapsis_angle_to_global;

Problem2ThetaPrimeIntervalCase classify_problem2_theta_prime_interval_case(
    std::size_t left_branch_count,
    std::size_t right_branch_count
) {
    if (left_branch_count == right_branch_count) {
        return Problem2ThetaPrimeIntervalCase::EqualBranchCount;
    }
    if (left_branch_count + 1U == right_branch_count || right_branch_count + 1U == left_branch_count) {
        return Problem2ThetaPrimeIntervalCase::BranchCountDifferenceOne;
    }
    return Problem2ThetaPrimeIntervalCase::BranchCountDifferenceGreaterThanOne;
}

Problem2CaseCMiddleSolveResult solve_case_c_middle_branches_on_k_layer(
    const Problem2ThetaPrimeScanConfig& config,
    int transfer_revolution,
    double theta_prime_left,
    const std::vector<Problem2OutgoingBranchSolution>& left_branches,
    double theta_prime_right,
    const std::vector<Problem2OutgoingBranchSolution>& right_branches
) {
    Problem2CaseCMiddleSolveResult result{};
    if (!is_finite(theta_prime_left) || !is_finite(theta_prime_right)) {
        result.status_reason = "non_finite_theta_prime_interval";
        return result;
    }
    if (!(theta_prime_right > theta_prime_left)) {
        result.status_reason = "invalid_theta_prime_interval_order";
        return result;
    }
    if (transfer_revolution < 0) {
        result.status_reason = "invalid_transfer_revolution";
        return result;
    }
    if (classify_problem2_theta_prime_interval_case(left_branches.size(), right_branches.size()) !=
        Problem2ThetaPrimeIntervalCase::BranchCountDifferenceGreaterThanOne) {
        result.status_reason = "case_c_requires_branch_count_difference_greater_than_one";
        return result;
    }
    if (left_branches.empty() || right_branches.empty()) {
        result.status_reason = "case_c_requires_nonempty_endpoint_branches";
        return result;
    }

    result.theta_prime_middle = 0.5 * (theta_prime_left + theta_prime_right);
    Problem2ThetaPrimeNodeSnapshot middle_snapshot =
        evaluate_problem2_theta_prime_node(config, result.theta_prime_middle);
    const std::size_t layer_index = static_cast<std::size_t>(transfer_revolution);
    if (layer_index >= middle_snapshot.solutions_by_k.size()) {
        result.status_reason = "case_c_middle_node_missing_k_layer";
        return result;
    }

    result.middle_branches = middle_snapshot.solutions_by_k[layer_index];
    if (result.middle_branches.empty()) {
        result.status_reason = "case_c_middle_k_layer_empty";
        return result;
    }

    for (auto& middle_branch : result.middle_branches) {
        const auto left_match = find_best_matching_outgoing_branch_index(
            middle_branch,
            left_branches,
            config.branch_phi_pairing_max_gap);
        const auto right_match = find_best_matching_outgoing_branch_index(
            middle_branch,
            right_branches,
            config.branch_phi_pairing_max_gap);

        if (left_match.has_value() && right_match.has_value()) {
            const auto& left_branch = left_branches[*left_match];
            const auto& right_branch = right_branches[*right_match];
            const double dphi = estimate_theta_prime_derivative_central(
                theta_prime_left,
                left_branch.encounter_global_angle,
                theta_prime_right,
                right_branch.encounter_global_angle);
            const double de = estimate_theta_prime_derivative_central(
                theta_prime_left,
                left_branch.outgoing_eccentricity,
                theta_prime_right,
                right_branch.outgoing_eccentricity);
            if (is_finite(dphi)) {
                middle_branch.dphi_dtheta_prime = dphi;
                middle_branch.has_dphi_dtheta_prime = true;
            }
            if (is_finite(de)) {
                middle_branch.de_dtheta_prime = de;
                middle_branch.has_de_dtheta_prime = true;
            }
            continue;
        }

        if (left_match.has_value()) {
            const auto& left_branch = left_branches[*left_match];
            const double dphi = estimate_theta_prime_derivative_forward(
                theta_prime_left,
                left_branch.encounter_global_angle,
                result.theta_prime_middle,
                middle_branch.encounter_global_angle);
            const double de = estimate_theta_prime_derivative_forward(
                theta_prime_left,
                left_branch.outgoing_eccentricity,
                result.theta_prime_middle,
                middle_branch.outgoing_eccentricity);
            if (is_finite(dphi)) {
                middle_branch.dphi_dtheta_prime = dphi;
                middle_branch.has_dphi_dtheta_prime = true;
            }
            if (is_finite(de)) {
                middle_branch.de_dtheta_prime = de;
                middle_branch.has_de_dtheta_prime = true;
            }
            continue;
        }

        if (right_match.has_value()) {
            const auto& right_branch = right_branches[*right_match];
            const double dphi = estimate_theta_prime_derivative_backward(
                result.theta_prime_middle,
                middle_branch.encounter_global_angle,
                theta_prime_right,
                right_branch.encounter_global_angle);
            const double de = estimate_theta_prime_derivative_backward(
                result.theta_prime_middle,
                middle_branch.outgoing_eccentricity,
                theta_prime_right,
                right_branch.outgoing_eccentricity);
            if (is_finite(dphi)) {
                middle_branch.dphi_dtheta_prime = dphi;
                middle_branch.has_dphi_dtheta_prime = true;
            }
            if (is_finite(de)) {
                middle_branch.de_dtheta_prime = de;
                middle_branch.has_de_dtheta_prime = true;
            }
        }
    }

    result.ok = true;
    return result;
}

Problem2FlybyGContext build_problem2_flyby_G_context(
    planet_params::PlanetId flyby_planet,
    double flyby_time_seconds_since_j2000,
    double incoming_eccentricity,
    double incoming_theta,
    bool incoming_theta_is_local
) {
    Problem2FlybyGContext context{};
    const auto planet_state = planet_state_at_time(flyby_planet, flyby_time_seconds_since_j2000);
    const double flyby_planet_eccentricity = get_planet_params(flyby_planet).orbit.e;
    context.incoming_cache = build_flyby_constraint_incoming_cache(
        incoming_eccentricity,
        incoming_theta,
        planet_state.theta_global,
        flyby_planet_eccentricity,
        incoming_theta_is_local);
    return context;
}

Problem2FlybyGSample evaluate_flyby_constraint_G_and_derivative(
    const FlybyConstraintIncomingCache& incoming_cache,
    double orbit_eccentricity,
    double theta_prime_local,
    double dphi_dtheta_prime,
    double de_dtheta_prime
) {
    Problem2FlybyGSample sample{};
    if (!incoming_cache.valid) {
        sample.invalid_reason = "invalid_incoming_cache";
        return sample;
    }
    if (!is_finite(orbit_eccentricity) ||
        !is_finite(theta_prime_local) ||
        !is_finite(dphi_dtheta_prime) ||
        !is_finite(de_dtheta_prime)) {
        sample.invalid_reason = "non_finite_G_derivative_input";
        return sample;
    }

    const FlybyConstraintFResult outgoing_F = evaluate_flyby_constraint_F(
        orbit_eccentricity,
        flyby_orbit_theta_global_from_input(theta_prime_local, incoming_cache.flyby_true_anomaly_phi, true),
        incoming_cache.flyby_true_anomaly_phi,
        incoming_cache.flyby_planet_eccentricity);
    if (!outgoing_F.valid) {
        sample.invalid_reason = outgoing_F.invalid_reason;
        return sample;
    }

    const FlybyConstraintFPartialDerivatives partials = evaluate_flyby_constraint_F_partial_derivatives(
        orbit_eccentricity,
        flyby_orbit_theta_global_from_input(theta_prime_local, incoming_cache.flyby_true_anomaly_phi, true),
        incoming_cache.flyby_true_anomaly_phi,
        incoming_cache.flyby_planet_eccentricity);
    if (!partials.valid) {
        sample.invalid_reason = partials.invalid_reason;
        return sample;
    }

    const double dF_dtheta_prime =
        partials.dF_de * de_dtheta_prime +
        partials.dF_dtheta_global * (dphi_dtheta_prime - 1.0);
    sample.G = incoming_cache.incoming_F - outgoing_F.value;
    sample.dG_dtheta_prime = -dF_dtheta_prime;
    if (!is_finite(sample.G) || !is_finite(sample.dG_dtheta_prime)) {
        sample.invalid_reason = "non_finite_G_or_derivative";
        return sample;
    }

    sample.valid = true;
    return sample;
}

std::optional<double> evaluate_flyby_constraint_G_at_branch(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2OutgoingBranchSolution& branch,
    double theta_prime_local
) {
    if (!incoming_cache.valid) {
        return std::nullopt;
    }
    const FlybyConstraintResidualResult residual = evaluate_flyby_constraint_residual_from_incoming_cache(
        incoming_cache,
        branch.outgoing_eccentricity,
        theta_prime_local,
        true);
    if (!residual.valid || !is_finite(residual.residual)) {
        return std::nullopt;
    }
    return residual.residual;
}

Problem2GNewtonResult refine_flyby_constraint_G_zero_by_newton(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2ThetaPrimeScanConfig& scan_config,
    const Problem2RouteANewtonOptions& route_a_options,
    const Problem2GNewtonOptions& options,
    double initial_theta_prime_local,
    const Problem2OutgoingBranchSolution& reference_branch,
    const Problem2OutgoingBranchSolution& linear_endpoint_branch,
    double linear_endpoint_theta_prime_local
) {
    Problem2GNewtonResult result{};
    result.final_theta_prime_local = initial_theta_prime_local;

    if (!incoming_cache.valid) {
        result.status = Problem2GNewtonStatus::DivergedInvalidG;
        result.status_reason = "invalid_incoming_cache";
        return result;
    }
    if (!is_finite(initial_theta_prime_local) ||
        !is_finite(linear_endpoint_theta_prime_local)) {
        result.status = Problem2GNewtonStatus::DivergedInvalidG;
        result.status_reason = "non_finite_newton_input";
        return result;
    }
    if (options.max_iterations <= 0 ||
        !is_finite(options.G_tolerance) || !(options.G_tolerance > 0.0) ||
        !is_finite(options.theta_prime_tolerance) || !(options.theta_prime_tolerance > 0.0)) {
        result.status = Problem2GNewtonStatus::DivergedInvalidG;
        result.status_reason = "invalid_G_newton_options";
        return result;
    }

    double theta_prime = initial_theta_prime_local;
    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
        const auto evaluation = evaluate_route_a_G_at_theta_prime(
            incoming_cache,
            scan_config,
            route_a_options,
            theta_prime,
            reference_branch,
            linear_endpoint_branch,
            linear_endpoint_theta_prime_local);
        if (!evaluation.has_value()) {
            result.status = Problem2GNewtonStatus::DivergedInvalidRouteA;
            result.status_reason = "failed_to_evaluate_route_a_G_sample";
            result.iterations_used = iteration;
            return result;
        }

        result.final_G = evaluation->g_sample.G;
        if (std::abs(evaluation->g_sample.G) <= options.G_tolerance) {
            result.ok = true;
            result.status = Problem2GNewtonStatus::Converged;
            result.iterations_used = iteration;
            result.final_theta_prime_local = theta_prime;
            result.outgoing_branch = evaluation->route_a.solution;
            return result;
        }

        if (!is_finite(evaluation->g_sample.dG_dtheta_prime) ||
            near_zero(evaluation->g_sample.dG_dtheta_prime, kDefaultEpsilon)) {
            result.status = Problem2GNewtonStatus::DivergedSingularDerivative;
            result.status_reason = "singular_G_derivative";
            result.iterations_used = iteration;
            return result;
        }

        const double theta_prime_next =
            theta_prime - evaluation->g_sample.G / evaluation->g_sample.dG_dtheta_prime;
        if (!is_finite(theta_prime_next)) {
            result.status = Problem2GNewtonStatus::DivergedInvalidG;
            result.status_reason = "non_finite_theta_prime_update";
            result.iterations_used = iteration;
            return result;
        }

        const auto next_evaluation = evaluate_route_a_G_at_theta_prime(
            incoming_cache,
            scan_config,
            route_a_options,
            theta_prime_next,
            reference_branch,
            linear_endpoint_branch,
            linear_endpoint_theta_prime_local);
        if (!next_evaluation.has_value()) {
            result.status = Problem2GNewtonStatus::DivergedInvalidRouteA;
            result.status_reason = "failed_to_evaluate_route_a_after_newton_step";
            result.iterations_used = iteration + 1;
            return result;
        }

        if (options.reject_on_G_increase &&
            std::abs(next_evaluation->g_sample.G) > std::abs(evaluation->g_sample.G)) {
            result.status = Problem2GNewtonStatus::DivergedGIncreased;
            result.status_reason = "G_increased_after_newton_step";
            result.iterations_used = iteration + 1;
            result.final_theta_prime_local = theta_prime_next;
            result.final_G = next_evaluation->g_sample.G;
            return result;
        }

        if (std::abs(theta_prime_next - theta_prime) <= options.theta_prime_tolerance &&
            std::abs(next_evaluation->g_sample.G) <= options.G_tolerance) {
            result.ok = true;
            result.status = Problem2GNewtonStatus::Converged;
            result.iterations_used = iteration + 1;
            result.final_theta_prime_local = theta_prime_next;
            result.final_G = next_evaluation->g_sample.G;
            result.outgoing_branch = next_evaluation->route_a.solution;
            return result;
        }

        theta_prime = theta_prime_next;
        result.final_theta_prime_local = theta_prime;
    }

    result.status = Problem2GNewtonStatus::DivergedMaxIterations;
    result.status_reason = "max_G_newton_iterations_exceeded";
    result.iterations_used = options.max_iterations;
    return result;
}

std::vector<Problem2FlybyGSolution> process_flyby_constraint_G_branch_interval(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2FlybyGSearchConfig& config,
    const Problem2FlybyGBranchIntervalInput& interval
) {
    std::vector<Problem2FlybyGSolution> solutions;
    if (!incoming_cache.valid) {
        return solutions;
    }
    if (!is_finite(interval.theta_prime_left) ||
        !is_finite(interval.theta_prime_right) ||
        !(interval.theta_prime_right > interval.theta_prime_left)) {
        return solutions;
    }
    if (!same_revolution_branch(interval.left_branch, interval.right_branch)) {
        return solutions;
    }

    const auto maybe_G_left = evaluate_flyby_constraint_G_at_branch(
        incoming_cache,
        interval.left_branch,
        interval.theta_prime_left);
    if (maybe_G_left.has_value() &&
        std::abs(*maybe_G_left) <= config.near_zero_G_threshold) {
        append_unique_flyby_G_solution(
            solutions,
            make_near_zero_solution(
                Problem2FlybyGSolutionSource::NearZeroEndpoint,
                interval.theta_prime_left,
                *maybe_G_left,
                interval.left_branch),
            config.solution_theta_prime_tolerance);
    }

    const auto maybe_G_right = evaluate_flyby_constraint_G_at_branch(
        incoming_cache,
        interval.right_branch,
        interval.theta_prime_right);
    if (maybe_G_right.has_value() &&
        std::abs(*maybe_G_right) <= config.near_zero_G_threshold) {
        append_unique_flyby_G_solution(
            solutions,
            make_near_zero_solution(
                Problem2FlybyGSolutionSource::NearZeroEndpoint,
                interval.theta_prime_right,
                *maybe_G_right,
                interval.right_branch),
            config.solution_theta_prime_tolerance);
    }

    if (!interval.left_branch.has_dphi_dtheta_prime ||
        !interval.right_branch.has_dphi_dtheta_prime ||
        !interval.left_branch.has_de_dtheta_prime ||
        !interval.right_branch.has_de_dtheta_prime) {
        return solutions;
    }

    const Problem2FlybyGSample G_left_sample = evaluate_flyby_constraint_G_and_derivative(
        incoming_cache,
        interval.left_branch.outgoing_eccentricity,
        interval.theta_prime_left,
        interval.left_branch.dphi_dtheta_prime,
        interval.left_branch.de_dtheta_prime);
    const Problem2FlybyGSample G_right_sample = evaluate_flyby_constraint_G_and_derivative(
        incoming_cache,
        interval.right_branch.outgoing_eccentricity,
        interval.theta_prime_right,
        interval.right_branch.dphi_dtheta_prime,
        interval.right_branch.de_dtheta_prime);
    if (!G_left_sample.valid || !G_right_sample.valid) {
        return solutions;
    }

    const auto quadratic_estimate = estimate_flyby_constraint_G_quadratic_root_on_theta_prime_interval(
        interval.theta_prime_left,
        G_left_sample.G,
        G_left_sample.dG_dtheta_prime,
        interval.theta_prime_right,
        G_right_sample.G,
        G_right_sample.dG_dtheta_prime);
    if (!quadratic_estimate.has_value() || !quadratic_estimate->has_root_in_interval) {
        return solutions;
    }

    const double predicted_theta_prime = quadratic_estimate->selected_root.theta_prime;
    const Problem2OutgoingBranchSolution* linear_endpoint = choose_linear_endpoint_for_theta_prime(
        interval.theta_prime_left,
        interval.left_branch,
        interval.theta_prime_right,
        interval.right_branch,
        predicted_theta_prime);
    if (linear_endpoint == nullptr) {
        return solutions;
    }
    const double linear_endpoint_theta =
        linear_endpoint == &interval.left_branch ? interval.theta_prime_left : interval.theta_prime_right;

    const Problem2GNewtonResult newton = refine_flyby_constraint_G_zero_by_newton(
        incoming_cache,
        config.scan_config,
        config.route_a_newton_options,
        config.g_newton_options,
        predicted_theta_prime,
        interval.left_branch,
        *linear_endpoint,
        linear_endpoint_theta);
    if (!newton.ok) {
        return solutions;
    }

    append_unique_flyby_G_solution(
        solutions,
        Problem2FlybyGSolution{
            .source = Problem2FlybyGSolutionSource::QuadraticNewton,
            .theta_prime_local = newton.final_theta_prime_local,
            .G = newton.final_G,
            .outgoing_branch = newton.outgoing_branch,
        },
        config.solution_theta_prime_tolerance);
    return solutions;
}

Problem2FlybyGSearchResult search_flyby_constraint_G_zeros_on_k_layer(
    const Problem2FlybyGContext& context,
    const Problem2FlybyGSearchConfig& config,
    const Problem2ThetaPrimeInitialScanResult& scan,
    int transfer_revolution
) {
    Problem2FlybyGSearchResult result{};
    if (!context.incoming_cache.valid) {
        result.error_message = "invalid_flyby_G_context";
        return result;
    }
    if (!scan.ok || scan.nodes.size() < 2U) {
        result.error_message = "invalid_initial_scan";
        return result;
    }
    if (transfer_revolution < 0 || transfer_revolution > scan.max_transfer_revolution) {
        result.error_message = "invalid_transfer_revolution_layer";
        return result;
    }

    const std::size_t layer_index = static_cast<std::size_t>(transfer_revolution);
    for (std::size_t node_index = 0; node_index + 1 < scan.nodes.size(); ++node_index) {
        const auto& left_node = scan.nodes[node_index];
        const auto& right_node = scan.nodes[node_index + 1];
        if (layer_index >= left_node.solutions_by_k.size() ||
            layer_index >= right_node.solutions_by_k.size()) {
            continue;
        }

        const auto& left_layer = left_node.solutions_by_k[layer_index];
        const auto& right_layer = right_node.solutions_by_k[layer_index];
        const auto interval_solutions = search_flyby_constraint_G_on_k_layer_interval(
            context.incoming_cache,
            config,
            transfer_revolution,
            left_node.theta_prime_local,
            left_layer,
            right_node.theta_prime_local,
            right_layer,
            0);
        for (const auto& solution : interval_solutions) {
            append_unique_flyby_G_solution(
                result.solutions,
                solution,
                config.solution_theta_prime_tolerance);
        }
    }

    result.ok = true;
    return result;
}

Problem2FlybyGSearchResult search_flyby_constraint_G_zeros_from_initial_scan(
    const Problem2FlybyGContext& context,
    const Problem2FlybyGSearchConfig& config,
    const Problem2ThetaPrimeInitialScanResult& scan
) {
    Problem2FlybyGSearchResult merged{};
    if (!context.incoming_cache.valid) {
        merged.error_message = "invalid_flyby_G_context";
        return merged;
    }
    if (!scan.ok) {
        merged.error_message = "invalid_initial_scan";
        return merged;
    }

    for (int transfer_revolution = 0; transfer_revolution <= scan.max_transfer_revolution; ++transfer_revolution) {
        const auto layer_result = search_flyby_constraint_G_zeros_on_k_layer(
            context,
            config,
            scan,
            transfer_revolution);
        if (!layer_result.ok) {
            merged.error_message = layer_result.error_message;
            return merged;
        }
        for (const auto& solution : layer_result.solutions) {
            append_unique_flyby_G_solution(
                merged.solutions,
                solution,
                config.solution_theta_prime_tolerance);
        }
    }

    merged.ok = true;
    return merged;
}

}  // namespace spaceship_cpp::problem2
