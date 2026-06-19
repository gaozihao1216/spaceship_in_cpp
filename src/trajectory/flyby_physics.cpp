/*
 * 文件作用：实现物理飞掠可行性判断。
 * 主要工作：计算超速速度矢量、所需转角、最大可转角和近心距约束。
 */
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace spaceship_cpp::trajectory {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

double nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

// 二维速度矢量的点积。
double dot(CartesianVelocity lhs, CartesianVelocity rhs) {
    return lhs.vx * rhs.vx + lhs.vy * rhs.vy;
}

// 二维速度矢量减法。
CartesianVelocity subtract(CartesianVelocity lhs, CartesianVelocity rhs) {
    return CartesianVelocity{lhs.vx - rhs.vx, lhs.vy - rhs.vy};
}

// 二维速度矢量的模长。
double norm(CartesianVelocity velocity) {
    return std::hypot(velocity.vx, velocity.vy);
}

}  // namespace

// 规范化 (e, θ)：负偏心率等价于 e→|e|、θ→θ+π；
// 解决不同模块对轨道参数表示不一致导致的比较/搜索问题。
CanonicalOrbitETheta canonicalize_orbit_e_theta(double eccentricity, double periapsis_angle) {
    CanonicalOrbitETheta result{};
    if (!is_finite(eccentricity) || !is_finite(periapsis_angle)) {
        return result;
    }
    result.valid = true;
    if (eccentricity >= 0.0) {
        result.eccentricity = eccentricity;
        result.periapsis_angle = normalize_angle_0_2pi(periapsis_angle);
    } else {
        result.eccentricity = -eccentricity;
        result.periapsis_angle = normalize_angle_0_2pi(periapsis_angle + kPi);
    }
    return result;
}

// 给定 v_∞ 和最小近心距，计算双曲线飞掠允许的最大转角；
// 解决判断飞掠几何是否物理可行的问题。
double compute_max_flyby_turn_angle_rad(
    double planet_mu,
    double min_periapsis_radius,
    double v_inf
) {
    if (!is_finite(planet_mu) || !is_finite(min_periapsis_radius) || !is_finite(v_inf) ||
        !(planet_mu > 0.0) || !(min_periapsis_radius > 0.0) || !(v_inf > 0.0)) {
        return nan();
    }
    const double e_hyp_min = 1.0 + min_periapsis_radius * v_inf * v_inf / planet_mu;
    if (!is_finite(e_hyp_min) || !(e_hyp_min > 1.0)) {
        return nan();
    }
    return 2.0 * std::asin(std::clamp(1.0 / e_hyp_min, -1.0, 1.0));
}

// 反解达到指定转角所需的近心距 r_p；
// 与 compute_max_flyby_turn_angle_rad 互为正逆问题。
double compute_required_flyby_periapsis_radius(
    double planet_mu,
    double v_inf,
    double turn_angle_rad
) {
    if (!is_finite(planet_mu) || !is_finite(v_inf) || !is_finite(turn_angle_rad) ||
        !(planet_mu > 0.0) || !(v_inf > 0.0) || turn_angle_rad < 0.0) {
        return nan();
    }
    if (turn_angle_rad <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    const double sine = std::sin(0.5 * turn_angle_rad);
    if (!is_finite(sine) || !(sine > 0.0)) {
        return nan();
    }
    return planet_mu / (v_inf * v_inf) * (1.0 / sine - 1.0);
}

// 一站式评估飞掠物理可行性：检查 v_∞ 匹配、转角上限、近心距下界；
// 解决 BFS 搜索中过滤不可行飞掠边的问题。
FlybyPhysicalFeasibilityResult evaluate_flyby_physical_feasibility(
    planet_params::PlanetId flyby_planet,
    double flyby_time,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta,
    const FlybyPhysicalFeasibilityOptions& options
) {
    FlybyPhysicalFeasibilityResult result{};
    result.planet = flyby_planet;

    if (!options.enabled || options.mode == FlybyPhysicalFilterMode::Disabled) {
        result.valid = true;
        result.feasible = true;
        return result;
    }
    if (!is_finite(flyby_time) || !is_finite(options.min_flyby_altitude_m) ||
        options.min_flyby_altitude_m < 0.0) {
        return result;
    }
    const auto incoming = canonicalize_orbit_e_theta(incoming_e, incoming_theta);
    const auto outgoing = canonicalize_orbit_e_theta(outgoing_e, outgoing_theta);
    if (!incoming.valid || !outgoing.valid) {
        return result;
    }

    const auto& solar = planet_params::get_solar_system_physical_params();
    const auto& params = planet_params::get_planet_params(flyby_planet);
    const auto planet_state = planet_params::planet_state_at_time(flyby_planet, flyby_time);
    const auto planet_velocity = compute_planet_heliocentric_velocity(flyby_planet, flyby_time);
    const auto incoming_velocity = compute_heliocentric_velocity_on_orbit(
        solar.GM_sun,
        planet_state.radius,
        planet_state.theta_global,
        incoming.eccentricity,
        incoming.periapsis_angle);
    const auto outgoing_velocity = compute_heliocentric_velocity_on_orbit(
        solar.GM_sun,
        planet_state.radius,
        planet_state.theta_global,
        outgoing.eccentricity,
        outgoing.periapsis_angle);

    const auto u_in = subtract(incoming_velocity, planet_velocity);
    const auto u_out = subtract(outgoing_velocity, planet_velocity);
    result.v_infinity_in = norm(u_in);
    result.v_infinity_out = norm(u_out);
    result.v_infinity_mismatch = std::abs(result.v_infinity_in - result.v_infinity_out);
    if (!is_finite(result.v_infinity_in) || !is_finite(result.v_infinity_out) ||
        !(result.v_infinity_in > 0.0) || !(result.v_infinity_out > 0.0)) {
        return result;
    }

    const double v_inf = 0.5 * (result.v_infinity_in + result.v_infinity_out);
    const double mismatch_tolerance = options.relative_speed_tolerance * std::max(1.0, v_inf);
    result.rejected_by_vinf_mismatch = result.v_infinity_mismatch > mismatch_tolerance;

    const double cos_delta = dot(u_in, u_out) / (result.v_infinity_in * result.v_infinity_out);
    if (!is_finite(cos_delta)) {
        return result;
    }
    result.turn_angle_rad = std::acos(std::clamp(cos_delta, -1.0, 1.0));
    result.min_allowed_periapsis_radius_m = params.physical.radius + options.min_flyby_altitude_m;
    result.max_turn_angle_rad = compute_max_flyby_turn_angle_rad(
        params.physical.GM,
        result.min_allowed_periapsis_radius_m,
        v_inf);
    result.required_periapsis_radius_m = compute_required_flyby_periapsis_radius(
        params.physical.GM,
        v_inf,
        result.turn_angle_rad);
    if (!is_finite(result.turn_angle_rad) || !is_finite(result.max_turn_angle_rad)) {
        return result;
    }

    result.rejected_by_turn_angle =
        result.turn_angle_rad > result.max_turn_angle_rad + options.turn_angle_tolerance_rad;
    result.rejected_by_periapsis_radius =
        is_finite(result.required_periapsis_radius_m) &&
        result.required_periapsis_radius_m < result.min_allowed_periapsis_radius_m;

    result.valid = true;
    result.feasible = !result.rejected_by_vinf_mismatch &&
        !result.rejected_by_turn_angle &&
        !result.rejected_by_periapsis_radius;
    return result;
}

}  // namespace spaceship_cpp::trajectory
