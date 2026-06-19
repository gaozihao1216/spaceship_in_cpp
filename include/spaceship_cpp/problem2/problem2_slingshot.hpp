/*
 * 文件作用：声明 Problem 2 飞掠弹弓残差方程。
 * 主要工作：定义轨道半径、外向转移轨道几何和弹弓约束残差计算接口。
 */
#pragma once

#include <string>

namespace spaceship_cpp::problem2 {

struct SlingshotOrbitElements {
    double semi_latus_rectum = 0.0;
    double eccentricity = 0.0;
    double perihelion_angle = 0.0;
};

struct SlingshotGeometryResult {
    bool valid = false;
    std::string invalid_reason;
    double eccentricity = 0.0;
    double semi_latus_rectum = 0.0;
    double r_departure = 0.0;
    double r_target = 0.0;
    double denominator = 0.0;
};

struct SlingshotInvariantResult {
    bool valid = false;
    std::string invalid_reason;
    double invariant = 0.0;
    double A = 0.0;
    double B = 0.0;
};

struct SlingshotResidualResult {
    bool valid = false;
    std::string invalid_reason;
    double incoming_invariant = 0.0;
    double outgoing_invariant = 0.0;
    double residual = 0.0;
};

struct SlingshotThetaAlphaResidualResult {
    bool valid = false;
    std::string invalid_reason;
    double alpha = 0.0;
    double theta_prime = 0.0;
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double slingshot_residual = 0.0;
    double incoming_invariant = 0.0;
    double outgoing_invariant = 0.0;
};

double problem2_orbit_radius(
    double semi_latus_rectum,
    double eccentricity,
    double true_anomaly_minus_perihelion_angle
);

SlingshotInvariantResult evaluate_problem2_slingshot_invariant(
    double encounter_true_anomaly_phi,
    double planet_eccentricity_e_J,
    double orbit_eccentricity_e,
    double orbit_perihelion_angle_theta
);

SlingshotGeometryResult solve_problem2_outgoing_orbit_from_two_points(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
);

SlingshotResidualResult evaluate_problem2_slingshot_residual(
    double phi,
    double e_J,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta
);

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
);

}  // namespace spaceship_cpp::problem2
