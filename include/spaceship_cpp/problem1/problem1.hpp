/*
 * 文件作用：声明 Problem 1 的直接残差评估和求解接口。
 * 主要工作：定义输入、残差结果、候选解，并暴露扫描加二分细化的求解函数。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <vector>

namespace spaceship_cpp::problem1 {

enum class Problem1ResidualStatus {
    // 残差函数正常计算完成，可以用 residual 判断转移时间和目标时间是否一致。
    Success,
    // 输入时间、角度或圈数不满足基本有效性约束。
    InvalidInput,
    // 两点几何导致转移轨道公式分母接近 0，无法稳定求解轨道参数。
    SingularGeometry,
    // 根据两点构造出的转移轨道不物理，例如半通径非正或分母无效。
    InvalidTransferOrbit,
    // 指定的转移圈数/目标圈数分支在当前几何下不可用。
    InvalidBranch,
    // 飞行时间不是正的有限值，不能作为候选轨道。
    InvalidTimeOfFlight,
};

struct Problem1ResidualInput {
    // 出发和目标行星决定两条已知开普勒轨道。
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;
    // 发射时刻用于定位出发行星和目标行星的相位。
    double launch_time_seconds_since_j2000;
    // 转移轨道近日点全局角，是本次残差评估固定的轨道方向。
    double transfer_perihelion_angle;
    // 待检查的遇合全局角，残差函数就是关于该角度求根。
    double encounter_global_angle;
    // 多圈转移分支编号：转移轨道绕太阳的额外圈数。
    int transfer_revolution;
    // 目标行星从发射到遇合期间绕太阳的额外圈数。
    int target_revolution;
};

struct Problem1ResidualResult {
    // 说明本次残差计算是否可用，以及失败时失败在哪个阶段。
    Problem1ResidualStatus status;
    // transfer_time_scale_free - target_time_scale_free；求根目标是让它接近 0。
    double residual;

    // 出发点和目标遇合点的太阳中心半径。
    double r1;
    double r2;

    // 两点相对转移轨道近日点的角度，以及目标行星起止角。
    double xi1;
    double xi2;
    double target_theta_start;
    double target_theta_end;

    // 由两点几何解出的转移轨道偏心率、实际使用的近日点角和半通径。
    double transfer_e_raw;
    double transfer_e;
    double transfer_perihelion_angle_used;
    double transfer_p;

    // 转移轨道和目标行星轨道分别走过的偏近点角/等效异常角差。
    double deltaF_transfer;
    double deltaF_target;

    // 未除以 sqrt(mu) 的无量纲时间；两者之差构成 residual。
    double transfer_time_scale_free;
    double target_time_scale_free;
};

struct Problem1SolveInput {
    // 固定一次 Problem 1 求解的出发/目标行星组合。
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;

    // 固定发射时刻和转移轨道方向，然后对遇合角做一维求根。
    double launch_time_seconds_since_j2000;
    double transfer_perihelion_angle;

    // 控制枚举哪些多圈分支，避免无限搜索。
    int max_transfer_revolution;
    int max_target_revolution;

    // 遇合角初始扫描网格数量，用来找残差符号变化区间。
    int phi_scan_count;

    // 二分终止和候选过滤阈值。
    double phi_tolerance;
    double residual_tolerance;
    int max_bisection_iterations;
    // 相对残差过滤用于避免时间尺度很大时只看绝对残差造成误判。
    double max_candidate_relative_residual = 1e-6;
};

struct Problem1Candidate {
    // 保留完整残差结果，方便调用方追踪几何和时间来源。
    Problem1ResidualResult residual_result;

    // 求得的遇合角、发射时刻、飞行时间和到达时刻。
    double encounter_global_angle;
    double launch_time_seconds_since_j2000;
    double time_of_flight_seconds;
    double arrival_time_seconds_since_j2000;

    // 用于衡量候选解质量的绝对/相对残差指标。
    double residual_scale;
    double relative_residual;

    // 记录根所在区间和二分细化过程，便于诊断求解精度。
    double root_bracket_width;
    int bisection_iterations;
    bool refined_by_bisection;

    // 候选解所属的转移圈数和目标圈数分支。
    int transfer_revolution;
    int target_revolution;
};

// 解决“给定一个遇合角时两段时间是否匹配”的问题，是 Problem 1 求根的核心函数。
Problem1ResidualResult evaluate_problem1_residual(const Problem1ResidualInput& input);

// 根据两个轨道点和相对近日点角，解转移轨道偏心率。
double compute_transfer_e_from_two_points(
    double r1,
    double xi1,
    double r2,
    double xi2
);

// 根据出发点半径、偏心率和相对角，反推转移轨道半通径 p。
double compute_transfer_p_from_departure(
    double r1,
    double e_transfer,
    double xi1
);

// 扫描遇合角并对符号变化区间二分，返回所有通过残差过滤的候选转移轨道。
std::vector<Problem1Candidate> solve_problem1(const Problem1SolveInput& input);

}  // namespace spaceship_cpp::problem1
