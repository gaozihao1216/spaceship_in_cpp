/*
 * 文件作用：声明轨道速度和行星相对速度计算接口。
 * 主要工作：计算转移轨道速度、行星速度、发射/到达超速速度等物理量。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

namespace spaceship_cpp::trajectory {

struct PolarVelocity {
    // 径向速度分量，正值表示远离太阳。
    double v_r = 0.0;
    // 横向速度分量，沿真近点角增大方向。
    double v_theta = 0.0;
};

struct CartesianVelocity {
    // 日心惯性平面坐标下的速度分量。
    double vx = 0.0;
    double vy = 0.0;
};

// 由当前位置半径和轨道形状反推半通径 p，解决只知道 e/theta 时速度公式缺 p 的问题。
double compute_semi_latus_rectum_from_radius(
    double radius,
    double eccentricity,
    double true_anomaly,
    double periapsis_angle
);

// 用开普勒轨道公式计算极坐标速度分量。
PolarVelocity compute_orbit_polar_velocity(
    double central_mu,
    double semi_latus_rectum,
    double eccentricity,
    double true_anomaly,
    double periapsis_angle
);

// 将极坐标速度转为全局平面直角坐标速度。
CartesianVelocity polar_velocity_to_cartesian(
    double true_anomaly,
    const PolarVelocity& polar_velocity
);

// 一站式计算任意圆锥轨道上某点的日心速度。
CartesianVelocity compute_heliocentric_velocity_on_orbit(
    double central_mu,
    double radius,
    double true_anomaly,
    double eccentricity,
    double periapsis_angle
);

// 计算某颗行星在指定时刻的日心速度。
CartesianVelocity compute_planet_heliocentric_velocity(
    planet_params::PlanetId planet,
    double time_seconds
);

// 计算发射时转移轨道速度相对行星速度的超速速度大小。
double relative_speed_to_planet(
    planet_params::PlanetId planet,
    double time_seconds,
    double orbit_eccentricity,
    double orbit_periapsis_angle
);

// 计算到达行星时入射轨道相对行星的 v_inf。
double compute_arrival_v_inf_at_planet(
    planet_params::PlanetId arrival_planet,
    double arrival_time,
    double incoming_e,
    double incoming_theta
);

}  // namespace spaceship_cpp::trajectory
