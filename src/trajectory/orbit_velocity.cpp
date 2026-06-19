/*
 * 文件作用：实现轨道速度和行星相对速度计算。
 * 主要工作：计算太阳中心速度、行星速度以及发射/到达时的超速速度。
 */
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <limits>

namespace spaceship_cpp::trajectory {
namespace {

using spaceship_cpp::common::is_finite;

double nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

CartesianVelocity nan_cartesian_velocity() {
    return CartesianVelocity{nan(), nan()};
}

}  // namespace

// 由半径和轨道形状反推半通径 p = r·(1+e·cos(θ-θ_p))；
// 解决只知道位置而速度公式需要 p 的问题。
double compute_semi_latus_rectum_from_radius(
    double radius,
    double eccentricity,
    double true_anomaly,
    double periapsis_angle
) {
    if (!is_finite(radius) || !is_finite(eccentricity) ||
        !is_finite(true_anomaly) || !is_finite(periapsis_angle) ||
        !(radius > 0.0) || !(eccentricity >= 0.0)) {
        return nan();
    }
    const double factor = 1.0 + eccentricity * std::cos(true_anomaly - periapsis_angle);
    const double p = radius * factor;
    if (!is_finite(p) || !(p > 0.0)) {
        return nan();
    }
    return p;
}

// 用开普勒轨道公式计算极坐标速度分量 (v_r, v_θ)；
// 解决日心圆锥轨道上任意点的速度矢量计算。
PolarVelocity compute_orbit_polar_velocity(
    double central_mu,
    double semi_latus_rectum,
    double eccentricity,
    double true_anomaly,
    double periapsis_angle
) {
    if (!is_finite(central_mu) || !is_finite(semi_latus_rectum) ||
        !is_finite(eccentricity) || !is_finite(true_anomaly) ||
        !is_finite(periapsis_angle) || !(central_mu > 0.0) ||
        !(semi_latus_rectum > 0.0) || !(eccentricity >= 0.0)) {
        return PolarVelocity{nan(), nan()};
    }
    const double scale = std::sqrt(central_mu / semi_latus_rectum);
    if (!is_finite(scale)) {
        return PolarVelocity{nan(), nan()};
    }
    const double anomaly_from_periapsis = true_anomaly - periapsis_angle;
    return PolarVelocity{
        eccentricity * scale * std::sin(anomaly_from_periapsis),
        scale * (1.0 + eccentricity * std::cos(anomaly_from_periapsis)),
    };
}

// 将极坐标速度 (v_r, v_θ) 旋转到全局日心直角坐标 (vx, vy)。
CartesianVelocity polar_velocity_to_cartesian(
    double true_anomaly,
    const PolarVelocity& polar_velocity
) {
    if (!is_finite(true_anomaly) ||
        !is_finite(polar_velocity.v_r) ||
        !is_finite(polar_velocity.v_theta)) {
        return nan_cartesian_velocity();
    }
    const double c = std::cos(true_anomaly);
    const double s = std::sin(true_anomaly);
    return CartesianVelocity{
        polar_velocity.v_r * c - polar_velocity.v_theta * s,
        polar_velocity.v_r * s + polar_velocity.v_theta * c,
    };
}

// 一站式：由轨道参数和位置计算日心直角速度；
// 组合半通径反推、极坐标速度、坐标旋转三步。
CartesianVelocity compute_heliocentric_velocity_on_orbit(
    double central_mu,
    double radius,
    double true_anomaly,
    double eccentricity,
    double periapsis_angle
) {
    const double semi_latus_rectum = compute_semi_latus_rectum_from_radius(
        radius, eccentricity, true_anomaly, periapsis_angle);
    if (!is_finite(semi_latus_rectum) || !(semi_latus_rectum > 0.0)) {
        return nan_cartesian_velocity();
    }
    const PolarVelocity polar = compute_orbit_polar_velocity(
        central_mu, semi_latus_rectum, eccentricity, true_anomaly, periapsis_angle);
    return polar_velocity_to_cartesian(true_anomaly, polar);
}

// 计算某颗行星在指定 J2000 秒偏移时刻的日心速度。
CartesianVelocity compute_planet_heliocentric_velocity(
    planet_params::PlanetId planet,
    double time_seconds
) {
    const auto& solar = planet_params::get_solar_system_physical_params();
    const auto& params = planet_params::get_planet_params(planet);
    const auto state = planet_params::planet_state_at_time(planet, time_seconds);
    return compute_heliocentric_velocity_on_orbit(
        solar.GM_sun,
        state.radius,
        state.theta_global,
        params.orbit.e,
        params.orbit.theta_0);
}

// 计算转移轨道速度与行星速度的差值模长（发射/到达超速速度 v_∞）；
// 解决评估发射 Δv 和飞掠入射条件的问题。
double relative_speed_to_planet(
    planet_params::PlanetId planet,
    double time_seconds,
    double orbit_eccentricity,
    double orbit_periapsis_angle
) {
    const auto& solar = planet_params::get_solar_system_physical_params();
    const auto state = planet_params::planet_state_at_time(planet, time_seconds);
    const CartesianVelocity orbit_velocity = compute_heliocentric_velocity_on_orbit(
        solar.GM_sun,
        state.radius,
        state.theta_global,
        orbit_eccentricity,
        orbit_periapsis_angle);
    const CartesianVelocity planet_velocity = compute_planet_heliocentric_velocity(planet, time_seconds);
    if (!is_finite(orbit_velocity.vx) || !is_finite(orbit_velocity.vy) ||
        !is_finite(planet_velocity.vx) || !is_finite(planet_velocity.vy)) {
        return nan();
    }
    const double dx = orbit_velocity.vx - planet_velocity.vx;
    const double dy = orbit_velocity.vy - planet_velocity.vy;
    return std::hypot(dx, dy);
}

// 到达行星时的入射 v_∞；语义上等价于 relative_speed_to_planet。
double compute_arrival_v_inf_at_planet(
    planet_params::PlanetId arrival_planet,
    double arrival_time,
    double incoming_e,
    double incoming_theta
) {
    return relative_speed_to_planet(arrival_planet, arrival_time, incoming_e, incoming_theta);
}

}  // namespace spaceship_cpp::trajectory
