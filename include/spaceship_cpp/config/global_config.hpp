/*
 * 文件作用：声明全局默认配置。
 * 主要工作：集中提供 Problem 1、行星范围和诊断程序使用的默认参数。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

namespace spaceship_cpp::config {

struct Problem1SolveDefaults {
    int max_transfer_revolution;              // 转移轨道允许的最大额外绕日圈数 k；0 表示只考虑最短前向转移分支
    int max_target_revolution;                // 目标行星允许的最大额外绕日圈数 q；0 表示目标行星只前进到最近一次相遇点
    int phi_scan_count;                       // 在 [0, 2pi) 上扫描相遇角 phi 的离散点数量；越大越不容易漏根，但越慢
    double phi_tolerance;                     // 二分求根时相遇角 phi 的收敛阈值，单位 rad
    double residual_tolerance;                // residual 的绝对容差；当前 residual 量纲较大，通常设为 0，主要依赖 phi_tolerance
    int max_bisection_iterations;             // 每个变号区间允许的最大二分迭代次数
    double max_candidate_relative_residual;   // 候选解允许的最大相对残差；用于过滤二分未真正收敛的伪候选
};

struct Problem1TableDefaults {
    double departure_true_anomaly_start;         // nu_A 起点，单位 rad
    int departure_true_anomaly_count;            // nu_A 维度离散数量
    double target_true_anomaly_start;            // nu_B 起点，单位 rad
    int target_true_anomaly_count;               // nu_B 维度离散数量
    double transfer_theta_departure_start;       // theta_A 起点，单位 rad
    int transfer_theta_departure_count;          // theta_A 维度离散数量
    int max_transfer_revolution;                 // 每个几何点预计算的椭圆多圈分支最大 k
    int max_target_revolution;                   // 每个几何点预计算的目标轨道额外绕行分支最大 q
};

struct Problem1DiagnosticsDefaults {
    int solve_phi_scan_count;                    // problem1_solve_diagnostics 中单次 solve 的 phi 扫描点数
    int table_departure_true_anomaly_count;     // problem1_table_diagnostics 中小表的 nu_A 维度数量
    int table_target_true_anomaly_count;        // problem1_table_diagnostics 中小表的 nu_B 维度数量
    int table_transfer_theta_departure_count;   // problem1_table_diagnostics 中小表的 theta_A 维度数量
    double max_candidate_relative_residual;     // diagnostics 中候选解最大相对残差阈值
};

struct GlobalConfig {
    Problem1SolveDefaults problem1_solve;              // Problem1 单次求解默认参数
    Problem1TableDefaults problem1_table_smoke;        // Problem1 小规模建表默认参数，用于快速测试
    Problem1TableDefaults problem1_table_medium;       // Problem1 中等规模建表默认参数，用于后续较认真诊断
    Problem1DiagnosticsDefaults problem1_diagnostics;  // Problem1 diagnostics app 默认参数
};

const GlobalConfig& global_config();

problem1::Problem1SolveInput make_problem1_solve_input(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    double launch_time_seconds_since_j2000,
    double transfer_perihelion_angle_global,
    const Problem1SolveDefaults& defaults
);  // 按统一默认参数构造单次 solve 输入，避免在 app 中重复手写扫描和二分配置

problem1::Problem1TableConfig make_problem1_table_config(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    const Problem1TableDefaults& defaults
);  // 按统一默认参数构造有序行星对 3D 表配置，并把三个周期角维度离散数量转换为固定角步长 2pi / count

}  // namespace spaceship_cpp::config
