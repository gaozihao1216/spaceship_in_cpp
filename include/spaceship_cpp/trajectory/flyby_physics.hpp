/*
 * 文件作用：声明物理飞掠可行性计算接口。
 * 主要工作：根据入射/出射超速速度、转角和近心距约束判断飞掠是否可行。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <limits>

namespace spaceship_cpp::trajectory {

enum class FlybyPhysicalFilterMode {
    // 完全不启用物理飞掠过滤。
    Disabled,
    // 只计算并记录可行性，不据此丢弃候选。
    ObserveOnly,
    // 严格执行物理可行性过滤。
    Enforce
};

struct CanonicalOrbitETheta {
    // valid 表示偏心率和近日点角已经规范化为统一表示。
    bool valid = false;
    double eccentricity = 0.0;
    double periapsis_angle = 0.0;
};

struct FlybyPhysicalFeasibilityOptions {
    // enabled 控制是否进行计算，mode 控制结果是否用于过滤。
    bool enabled = true;
    FlybyPhysicalFilterMode mode = FlybyPhysicalFilterMode::Enforce;

    double min_flyby_altitude_m = 300000.0;

    // 容差用于处理浮点误差下的速度匹配和转角边界。
    double turn_angle_tolerance_rad = 1e-10;
    double relative_speed_tolerance = 1e-6;
};

struct FlybyPhysicalFeasibilityResult {
    // valid 表示计算流程完成；feasible 表示物理约束通过。
    bool valid = false;

    bool feasible = false;

    planet_params::PlanetId planet;

    // 入射/出射超速速度大小及二者差异。
    double v_infinity_in = 0.0;
    double v_infinity_out = 0.0;
    double v_infinity_mismatch = 0.0;

    // 飞掠所需转角 δ，及在 (μ, r_p_min, v_∞) 下允许的最大转角 δ_max。
    // δ > δ_max 等价于 r_p_required < r_p_min（同一双曲线飞掠约束的两面，不独立判断）。
    double turn_angle_rad = 0.0;
    double max_turn_angle_rad = 0.0;

    // 诊断量：达到 δ 所需近心距，及 r_p_min = R_planet + h_min。
    double required_periapsis_radius_m = std::numeric_limits<double>::infinity();
    double min_allowed_periapsis_radius_m = std::numeric_limits<double>::infinity();

    bool rejected_by_vinf_mismatch = false;
    // 转角超限（含近心距不足）；max_turn_angle_rad 随本次飞掠 v_∞ 变化，非行星常数。
    bool rejected_by_turn_angle = false;
};

// 规范化偏心率和近日点角，解决负偏心率等价表示带来的比较问题。
CanonicalOrbitETheta canonicalize_orbit_e_theta(double eccentricity, double periapsis_angle);

// 给定 μ、r_p_min 与 v_∞，由双曲线飞掠公式 e_hyp = 1 + r_p v_∞²/μ 得允许的最大转角。
// δ_max 随 v_∞ 增大而减小，并非行星专属常数。
double compute_max_flyby_turn_angle_rad(
    double planet_mu,
    double min_periapsis_radius,
    double v_inf
);

// 反解达到指定转角所需的近心距。
double compute_required_flyby_periapsis_radius(
    double planet_mu,
    double v_inf,
    double turn_angle_rad
);

// 一站式评估飞掠物理可行性。
// 仅两类独立约束：(1) |v_∞,in − v_∞,out|；(2) δ ≤ δ_max(μ, r_p_min, v_∞)。
// 近心距比较与转角比较等价，只保留转角判据；required_periapsis_radius_m 仅作诊断输出。
FlybyPhysicalFeasibilityResult evaluate_flyby_physical_feasibility(
    planet_params::PlanetId flyby_planet,
    double flyby_time,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta,
    const FlybyPhysicalFeasibilityOptions& options
);

}  // namespace spaceship_cpp::trajectory
