/*
 * 文件作用：声明 Problem 2 飞掠弹弓残差方程。
 * 主要工作：定义轨道半径、外向转移轨道几何和弹弓约束残差计算接口。
 */
#pragma once

#include <string>

namespace spaceship_cpp::problem2 {

struct SlingshotOrbitElements {
    // 半通径 p，用于 r = p / (1 + e cos(...)) 的轨道半径公式。
    double semi_latus_rectum = 0.0;
    // 偏心率 e，决定出射轨道形状。
    double eccentricity = 0.0;
    // 近日点角 theta，决定轨道在全局坐标中的方向。
    double perihelion_angle = 0.0;
};

struct SlingshotGeometryResult {
    // valid 表示由两个点和 theta_prime 成功反推出外向轨道。
    bool valid = false;
    std::string invalid_reason;
    // 解出的外向轨道参数。
    double eccentricity = 0.0;
    double semi_latus_rectum = 0.0;
    // 飞掠行星点和目标点的轨道半径。
    double r_departure = 0.0;
    double r_target = 0.0;
    // 反解偏心率时的分母，接近 0 时几何退化。
    double denominator = 0.0;
};

struct SlingshotInvariantResult {
    // valid 表示弹弓不变量在当前几何下可计算。
    bool valid = false;
    std::string invalid_reason;
    // 用于比较入射/出射轨道的弹弓不变量。
    double invariant = 0.0;
    // 中间量 A/B 保留给诊断，便于定位残差来源。
    double A = 0.0;
    double B = 0.0;
};

struct SlingshotResidualResult {
    // valid 表示入射和出射不变量都成功计算。
    bool valid = false;
    std::string invalid_reason;
    // 两侧不变量及其差值；residual=0 表示满足弹弓约束。
    double incoming_invariant = 0.0;
    double outgoing_invariant = 0.0;
    double residual = 0.0;
};

struct SlingshotThetaAlphaResidualResult {
    // 从 theta_prime 和 alpha 出发的一站式残差结果。
    bool valid = false;
    std::string invalid_reason;
    double alpha = 0.0;
    double theta_prime = 0.0;
    // 由飞掠点和目标点几何反推的外向轨道参数。
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double slingshot_residual = 0.0;
    double incoming_invariant = 0.0;
    double outgoing_invariant = 0.0;
};

// 解决“给定圆锥轨道参数和相对真近点角，半径是多少”的问题。
double problem2_orbit_radius(
    double semi_latus_rectum,
    double eccentricity,
    double true_anomaly_minus_perihelion_angle
);

// 计算飞掠弹弓不变量，用于比较入射轨道和出射轨道是否可由同一次飞掠连接。
SlingshotInvariantResult evaluate_problem2_slingshot_invariant(
    double encounter_true_anomaly_phi,
    double planet_eccentricity_e_J,
    double orbit_eccentricity_e,
    double orbit_perihelion_angle_theta
);

// 由飞掠点、目标点和出射近日点角反推出外向转移轨道。
SlingshotGeometryResult solve_problem2_outgoing_orbit_from_two_points(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
);

// 直接比较入射/出射轨道的不变量，得到弹弓约束残差。
SlingshotResidualResult evaluate_problem2_slingshot_residual(
    double phi,
    double e_J,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta
);

// 从候选 theta_prime/alpha 出发，先解外向轨道，再计算弹弓残差。
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
