/*
 * 文件作用：实现 Problem 2 飞掠弹弓残差方程。
 * 主要工作：根据入射轨道、出射方向和目标几何计算外向轨道与弹弓约束残差。
 */
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <limits>

namespace spaceship_cpp::problem2 {

namespace {

using spaceship_cpp::common::is_finite;

constexpr double kGeometryDenominatorEpsilon = 1e-12;

double quiet_nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

// 圆锥轨道半径 r = p / (1 + e·cos(θ-θ_p))；
// 解决给定轨道参数和相对真近点角时的位置计算。
double problem2_orbit_radius(
    double semi_latus_rectum,
    double eccentricity,
    double true_anomaly_minus_perihelion_angle
) {
    if (!is_finite(semi_latus_rectum) ||
        !is_finite(eccentricity) ||
        !is_finite(true_anomaly_minus_perihelion_angle) ||
        !(semi_latus_rectum > 0.0)) {
        return quiet_nan();
    }
    const double denominator =
        1.0 + eccentricity * std::cos(true_anomaly_minus_perihelion_angle);
    if (!is_finite(denominator) || !(denominator > 0.0)) {
        return quiet_nan();
    }
    const double radius = semi_latus_rectum / denominator;
    return is_finite(radius) ? radius : quiet_nan();
}

// 计算飞掠弹弓不变量；同一飞掠事件下入射/出射轨道的不变量应相等，
// 解决判断两条轨道能否由同一次行星飞掠连接的问题。
SlingshotInvariantResult evaluate_problem2_slingshot_invariant(
    double encounter_true_anomaly_phi,
    double planet_eccentricity_e_J,
    double orbit_eccentricity_e,
    double orbit_perihelion_angle_theta
) {
    SlingshotInvariantResult result{};
    if (!is_finite(encounter_true_anomaly_phi) ||
        !is_finite(planet_eccentricity_e_J) ||
        !is_finite(orbit_eccentricity_e) ||
        !is_finite(orbit_perihelion_angle_theta)) {
        result.invalid_reason = "non_finite_slingshot_invariant";
        return result;
    }

    const double cos_phi = std::cos(encounter_true_anomaly_phi);
    const double delta = encounter_true_anomaly_phi - orbit_perihelion_angle_theta;
    const double cos_delta = std::cos(delta);
    const double cos_theta = std::cos(orbit_perihelion_angle_theta);

    result.A = 1.0 + planet_eccentricity_e_J * cos_phi;
    result.B = 1.0 + orbit_eccentricity_e * cos_delta;

    if (!(result.A > 0.0)) {
        result.invalid_reason = "non_positive_planet_radius_factor";
        return result;
    }
    if (!(result.B > 0.0)) {
        result.invalid_reason = "non_positive_orbit_radius_factor";
        return result;
    }

    const double ratio = result.A / result.B;
    if (!is_finite(ratio) || !(ratio > 0.0)) {
        result.invalid_reason = "non_finite_slingshot_invariant";
        return result;
    }

    const double C = result.A + orbit_eccentricity_e * (cos_delta + planet_eccentricity_e_J * cos_theta);
    const double kinetic_factor =
        1.0 + 2.0 * orbit_eccentricity_e * cos_delta + orbit_eccentricity_e * orbit_eccentricity_e;
    result.invariant = ratio * kinetic_factor - 2.0 * std::sqrt(ratio) * C;
    if (!is_finite(result.invariant)) {
        result.invalid_reason = "non_finite_slingshot_invariant";
        return result;
    }

    result.valid = true;
    return result;
}

// 由飞掠点 (R_J,φ)、目标点 (R_K,α) 和出射近日点角 θ' 反推外向转移轨道 (e,p)；
// 解决 Problem 2 中"给定出射方向后出射轨道是什么"的几何问题。
SlingshotGeometryResult solve_problem2_outgoing_orbit_from_two_points(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
) {
    SlingshotGeometryResult result{};
    result.r_departure = problem2_orbit_radius(R_J, e_J, phi);
    if (!is_finite(result.r_departure) || !(result.r_departure > 0.0)) {
        result.invalid_reason = "invalid_encounter_radius";
        return result;
    }

    result.r_target = problem2_orbit_radius(R_K, e_K, alpha);
    if (!is_finite(result.r_target) || !(result.r_target > 0.0)) {
        result.invalid_reason = "invalid_target_radius";
        return result;
    }

    const double cos_departure = std::cos(phi - theta_prime);
    const double cos_target = std::cos(alpha - theta_prime);
    result.denominator = result.r_departure * cos_departure - result.r_target * cos_target;
    if (!is_finite(result.denominator) ||
        std::abs(result.denominator) <= kGeometryDenominatorEpsilon) {
        result.invalid_reason = "geometry_denominator_too_small";
        return result;
    }

    result.eccentricity = (result.r_target - result.r_departure) / result.denominator;
    if (!is_finite(result.eccentricity)) {
        result.invalid_reason = "non_finite_outgoing_eccentricity";
        return result;
    }

    result.semi_latus_rectum = result.r_departure * (1.0 + result.eccentricity * cos_departure);
    if (!is_finite(result.semi_latus_rectum) || !(result.semi_latus_rectum > 0.0)) {
        result.invalid_reason = "non_positive_outgoing_semi_latus_rectum";
        return result;
    }

    const double departure_factor = 1.0 + result.eccentricity * cos_departure;
    if (!(departure_factor > 0.0) || !is_finite(departure_factor)) {
        result.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return result;
    }
    const double target_factor = 1.0 + result.eccentricity * cos_target;
    if (!(target_factor > 0.0) || !is_finite(target_factor)) {
        result.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return result;
    }

    const double departure_back_substitute = result.semi_latus_rectum / departure_factor;
    if (!is_finite(departure_back_substitute) ||
        std::abs(departure_back_substitute - result.r_departure) > 1e-10) {
        result.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return result;
    }
    const double target_back_substitute = result.semi_latus_rectum / target_factor;
    if (!is_finite(target_back_substitute) ||
        std::abs(target_back_substitute - result.r_target) > 1e-10) {
        result.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return result;
    }

    result.valid = true;
    return result;
}

// 比较入射/出射轨道的弹弓不变量，residual = I_out - I_in；
// residual=0 表示满足弹弓约束，是 Problem 2 求根的目标函数。
SlingshotResidualResult evaluate_problem2_slingshot_residual(
    double phi,
    double e_J,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta
) {
    SlingshotResidualResult result{};
    const SlingshotInvariantResult incoming =
        evaluate_problem2_slingshot_invariant(phi, e_J, incoming_e, incoming_theta);
    if (!incoming.valid) {
        result.invalid_reason = incoming.invalid_reason;
        return result;
    }
    const SlingshotInvariantResult outgoing =
        evaluate_problem2_slingshot_invariant(phi, e_J, outgoing_e, outgoing_theta);
    if (!outgoing.valid) {
        result.invalid_reason = outgoing.invalid_reason;
        return result;
    }

    result.incoming_invariant = incoming.invariant;
    result.outgoing_invariant = outgoing.invariant;
    result.residual = outgoing.invariant - incoming.invariant;
    if (!is_finite(result.residual)) {
        result.invalid_reason = "non_finite_slingshot_invariant";
        return result;
    }
    result.valid = true;
    return result;
}

// 从候选 (θ', α) 出发：先解外向轨道，再计算弹弓残差；
// 解决搜索过程中一站式评估候选飞掠方案的问题。
SlingshotThetaAlphaResidualResult evaluate_problem2_slingshot_residual_from_theta_alpha(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    SlingshotThetaAlphaResidualResult result{};
    result.alpha = alpha;
    result.theta_prime = theta_prime;

    const SlingshotGeometryResult geometry = solve_problem2_outgoing_orbit_from_two_points(
        R_J,
        e_J,
        R_K,
        e_K,
        phi,
        alpha,
        theta_prime);
    if (!geometry.valid) {
        result.invalid_reason = geometry.invalid_reason;
        return result;
    }

    result.outgoing_eccentricity = geometry.eccentricity;
    result.outgoing_semi_latus_rectum = geometry.semi_latus_rectum;

    const SlingshotResidualResult residual = evaluate_problem2_slingshot_residual(
        phi,
        e_J,
        incoming_e,
        incoming_theta,
        geometry.eccentricity,
        theta_prime);
    if (!residual.valid) {
        result.invalid_reason = residual.invalid_reason;
        return result;
    }

    result.incoming_invariant = residual.incoming_invariant;
    result.outgoing_invariant = residual.outgoing_invariant;
    result.slingshot_residual = residual.residual;
    result.valid = true;
    return result;
}

}  // namespace spaceship_cpp::problem2
