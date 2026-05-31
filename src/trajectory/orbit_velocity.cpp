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

double compute_arrival_v_inf_at_planet(
    planet_params::PlanetId arrival_planet,
    double arrival_time,
    double incoming_e,
    double incoming_theta
) {
    return relative_speed_to_planet(arrival_planet, arrival_time, incoming_e, incoming_theta);
}

}  // namespace spaceship_cpp::trajectory
