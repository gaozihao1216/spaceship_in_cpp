/*
 * 文件作用：声明轨道速度和行星相对速度计算接口。
 * 主要工作：计算转移轨道速度、行星速度、发射/到达超速速度等物理量。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

namespace spaceship_cpp::trajectory {

struct PolarVelocity {
    double v_r = 0.0;
    double v_theta = 0.0;
};

struct CartesianVelocity {
    double vx = 0.0;
    double vy = 0.0;
};

double compute_semi_latus_rectum_from_radius(
    double radius,
    double eccentricity,
    double true_anomaly,
    double periapsis_angle
);

PolarVelocity compute_orbit_polar_velocity(
    double central_mu,
    double semi_latus_rectum,
    double eccentricity,
    double true_anomaly,
    double periapsis_angle
);

CartesianVelocity polar_velocity_to_cartesian(
    double true_anomaly,
    const PolarVelocity& polar_velocity
);

CartesianVelocity compute_heliocentric_velocity_on_orbit(
    double central_mu,
    double radius,
    double true_anomaly,
    double eccentricity,
    double periapsis_angle
);

CartesianVelocity compute_planet_heliocentric_velocity(
    planet_params::PlanetId planet,
    double time_seconds
);

double relative_speed_to_planet(
    planet_params::PlanetId planet,
    double time_seconds,
    double orbit_eccentricity,
    double orbit_periapsis_angle
);

double compute_arrival_v_inf_at_planet(
    planet_params::PlanetId arrival_planet,
    double arrival_time,
    double incoming_e,
    double incoming_theta
);

}  // namespace spaceship_cpp::trajectory
