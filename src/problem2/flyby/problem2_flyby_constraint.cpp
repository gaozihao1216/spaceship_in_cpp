/*
 * 文件作用：实现 Problem 2 弹弓约束 F 函数与可缓存残差。
 * 主要工作：提供 G = F_in - F_out 求值，以及按 branch 分组的 θ' 候选识别。
 */
#include "spaceship_cpp/problem2/problem2_flyby_constraint.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <limits>

namespace spaceship_cpp::problem2 {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kFlybyConstraintDenominatorEpsilon = 1e-12;

double quiet_nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool flyby_constraint_sign_changed(double left, double right) {
    if (!is_finite(left) || !is_finite(right)) {
        return false;
    }
    if (left == 0.0 || right == 0.0) {
        return true;
    }
    return left * right < 0.0;
}

}  // namespace

FlybyConstraintFResult evaluate_flyby_constraint_F(
    double orbit_eccentricity,
    double orbit_perihelion_angle_global,
    double flyby_true_anomaly_phi,
    double flyby_planet_eccentricity
) {
    FlybyConstraintFResult result{};
    if (!is_finite(orbit_eccentricity) ||
        !is_finite(orbit_perihelion_angle_global) ||
        !is_finite(flyby_true_anomaly_phi) ||
        !is_finite(flyby_planet_eccentricity)) {
        result.invalid_reason = "non_finite_flyby_constraint_input";
        return result;
    }

    const double cos_phi = std::cos(flyby_true_anomaly_phi);
    const double cos_delta = std::cos(flyby_true_anomaly_phi - orbit_perihelion_angle_global);
    const double cos_theta = std::cos(orbit_perihelion_angle_global);

    result.denominator = 1.0 + orbit_eccentricity * cos_delta;
    if (!is_finite(result.denominator) || result.denominator <= kFlybyConstraintDenominatorEpsilon) {
        result.invalid_reason = "non_positive_flyby_constraint_denominator";
        return result;
    }

    const double numerator =
        1.0 + flyby_planet_eccentricity * cos_phi +
        orbit_eccentricity * (cos_delta + flyby_planet_eccentricity * cos_theta);
    if (!is_finite(numerator)) {
        result.invalid_reason = "non_finite_flyby_constraint_numerator";
        return result;
    }

    result.value = numerator / std::sqrt(result.denominator);
    if (!is_finite(result.value)) {
        result.invalid_reason = "non_finite_flyby_constraint_value";
        return result;
    }

    result.valid = true;
    return result;
}

double flyby_orbit_theta_global_from_input(
    double orbit_theta,
    double flyby_true_anomaly_phi,
    bool input_theta_is_local
) {
    if (input_theta_is_local) {
        return normalize_angle_minus_pi_pi(flyby_true_anomaly_phi - orbit_theta);
    }
    return normalize_angle_minus_pi_pi(orbit_theta);
}

FlybyConstraintIncomingCache build_flyby_constraint_incoming_cache(
    double incoming_eccentricity,
    double incoming_theta,
    double flyby_true_anomaly_phi,
    double flyby_planet_eccentricity,
    bool input_theta_is_local
) {
    FlybyConstraintIncomingCache cache{};
    if (!is_finite(incoming_eccentricity) ||
        !is_finite(incoming_theta) ||
        !is_finite(flyby_true_anomaly_phi) ||
        !is_finite(flyby_planet_eccentricity)) {
        cache.invalid_reason = "non_finite_incoming_cache_input";
        return cache;
    }

    cache.flyby_true_anomaly_phi = flyby_true_anomaly_phi;
    cache.flyby_planet_eccentricity = flyby_planet_eccentricity;
    cache.incoming_eccentricity = incoming_eccentricity;
    cache.incoming_theta_global = flyby_orbit_theta_global_from_input(
        incoming_theta,
        flyby_true_anomaly_phi,
        input_theta_is_local);

    const FlybyConstraintFResult incoming_F = evaluate_flyby_constraint_F(
        cache.incoming_eccentricity,
        cache.incoming_theta_global,
        cache.flyby_true_anomaly_phi,
        cache.flyby_planet_eccentricity);
    if (!incoming_F.valid) {
        cache.invalid_reason = incoming_F.invalid_reason;
        return cache;
    }

    cache.incoming_F = incoming_F.value;
    cache.valid = true;
    return cache;
}

FlybyConstraintResidualResult evaluate_flyby_constraint_residual(
    double incoming_eccentricity,
    double incoming_theta,
    double outgoing_eccentricity,
    double outgoing_theta,
    double flyby_true_anomaly_phi,
    double flyby_planet_eccentricity,
    bool input_angles_are_local
) {
    const FlybyConstraintIncomingCache cache = build_flyby_constraint_incoming_cache(
        incoming_eccentricity,
        incoming_theta,
        flyby_true_anomaly_phi,
        flyby_planet_eccentricity,
        input_angles_are_local);
    return evaluate_flyby_constraint_residual_from_incoming_cache(
        cache,
        outgoing_eccentricity,
        outgoing_theta,
        input_angles_are_local);
}

FlybyConstraintResidualResult evaluate_flyby_constraint_residual_from_incoming_cache(
    const FlybyConstraintIncomingCache& incoming_cache,
    double outgoing_eccentricity,
    double outgoing_theta,
    bool outgoing_theta_is_local
) {
    FlybyConstraintResidualResult result{};
    if (!incoming_cache.valid) {
        result.invalid_reason = incoming_cache.invalid_reason.empty()
            ? "invalid_incoming_cache"
            : incoming_cache.invalid_reason;
        return result;
    }
    if (!is_finite(outgoing_eccentricity) || !is_finite(outgoing_theta)) {
        result.invalid_reason = "non_finite_outgoing_input";
        return result;
    }

    const double outgoing_theta_global = flyby_orbit_theta_global_from_input(
        outgoing_theta,
        incoming_cache.flyby_true_anomaly_phi,
        outgoing_theta_is_local);
    const FlybyConstraintFResult outgoing_F = evaluate_flyby_constraint_F(
        outgoing_eccentricity,
        outgoing_theta_global,
        incoming_cache.flyby_true_anomaly_phi,
        incoming_cache.flyby_planet_eccentricity);
    if (!outgoing_F.valid) {
        result.invalid_reason = outgoing_F.invalid_reason;
        return result;
    }

    result.incoming_F = incoming_cache.incoming_F;
    result.outgoing_F = outgoing_F.value;
    result.residual = result.incoming_F - result.outgoing_F;
    if (!is_finite(result.residual)) {
        result.invalid_reason = "non_finite_flyby_constraint_residual";
        return result;
    }

    result.valid = true;
    return result;
}

std::vector<FlybyThetaPrimeCandidate> detect_flyby_theta_prime_candidates_from_branch_samples(
    const std::vector<FlybyThetaPrimeBranchSample>& branch_samples,
    double near_zero_threshold
) {
    std::vector<FlybyThetaPrimeCandidate> candidates;
    if (!(near_zero_threshold > 0.0) || branch_samples.empty()) {
        return candidates;
    }

    for (const auto& sample : branch_samples) {
        if (!sample.flyby_constraint_valid || !is_finite(sample.flyby_constraint_residual)) {
            continue;
        }
        if (std::abs(sample.flyby_constraint_residual) > near_zero_threshold) {
            continue;
        }

        FlybyThetaPrimeCandidate candidate{};
        candidate.type = FlybyThetaPrimeCandidateType::NearZeroNode;
        candidate.branch_id = sample.branch_id;
        candidate.theta_prime_node = sample.theta_prime;
        candidate.G_node = sample.flyby_constraint_residual;
        candidate.outgoing_eccentricity = sample.outgoing_eccentricity;
        candidate.outgoing_semi_latus_rectum = sample.outgoing_semi_latus_rectum;
        candidates.push_back(candidate);
    }

    for (std::size_t index = 1; index < branch_samples.size(); ++index) {
        const auto& left = branch_samples[index - 1];
        const auto& right = branch_samples[index];
        if (!left.flyby_constraint_valid || !right.flyby_constraint_valid) {
            continue;
        }
        if (!is_finite(left.flyby_constraint_residual) || !is_finite(right.flyby_constraint_residual)) {
            continue;
        }
        if (!flyby_constraint_sign_changed(
                left.flyby_constraint_residual,
                right.flyby_constraint_residual)) {
            continue;
        }

        FlybyThetaPrimeCandidate candidate{};
        candidate.type = FlybyThetaPrimeCandidateType::SignChangeInterval;
        candidate.branch_id = left.branch_id;
        candidate.theta_prime_left = left.theta_prime;
        candidate.theta_prime_right = right.theta_prime;
        candidate.G_left = left.flyby_constraint_residual;
        candidate.G_right = right.flyby_constraint_residual;
        candidate.outgoing_eccentricity = left.outgoing_eccentricity;
        candidate.outgoing_semi_latus_rectum = left.outgoing_semi_latus_rectum;
        candidates.push_back(candidate);
    }

    return candidates;
}

}  // namespace spaceship_cpp::problem2
