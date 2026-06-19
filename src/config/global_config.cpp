/*
 * 文件作用：实现全局默认配置。
 * 主要工作：集中返回默认行星、Problem 1 扫描参数和诊断用配置值。
 */
#include "spaceship_cpp/config/global_config.hpp"

#include "spaceship_cpp/common/common.hpp"

namespace spaceship_cpp::config {

// 返回全局单例默认配置；集中管理 Problem 1 求解/建表/诊断的默认参数。
const GlobalConfig& global_config() {
    static const GlobalConfig kConfig{
        Problem1SolveDefaults{
            0,      // 默认不允许转移轨道额外绕日，先聚焦最直接的前向转移分支
            0,      // 默认不允许目标行星多绕额外圈数，先只看最近一次前向相遇
            720,    // 单次 solve 中 phi 扫描 720 个点，约等于每 0.5 度采样一次
            1e-10,  // 二分求根的角度容差
            0.0,    // residual 绝对容差当前保持为 0
            80,     // 每个变号区间最多二分 80 次
            1e-6,   // 候选解相对残差阈值
        },
        Problem1TableDefaults{
            0.0,    // nu_A 从 0 开始
            8,      // nu_A 方向取 8 个点
            0.0,    // nu_B 从 0 开始
            8,      // nu_B 方向取 8 个点
            0.0,    // theta_A 从 0 开始
            16,     // theta_A 方向取 16 个点
            0,      // 小规模表先只保留最短前向分支
            0,      // 小规模表先只保留目标轨道最近一次前向相遇分支
        },
        Problem1TableDefaults{
            0.0,    // nu_A 从 0 开始
            32,     // nu_A 方向细化
            0.0,    // nu_B 从 0 开始
            32,     // nu_B 方向细化
            0.0,    // theta_A 从 0 开始
            64,     // theta_A 方向细化
            0,      // 中等规模表仍先不扩展额外椭圆多圈分支
            0,      // 中等规模表也先不扩展目标轨道额外绕行分支
        },
        Problem1DiagnosticsDefaults{
            720,    // 单次 solve diagnostics 默认使用 720 个 phi 扫描点
            8,      // table diagnostics 默认使用 8 个 nu_A 点
            8,      // table diagnostics 默认使用 8 个 nu_B 点
            16,     // table diagnostics 默认使用 16 个 theta_A 点
            1e-6,   // diagnostics 默认沿用 1e-6 相对残差阈值
        },
    };

    return kConfig;
}

// 按统一默认参数构造 Problem1SolveInput；避免 app/测试中重复手写扫描和二分配置。
problem1::Problem1SolveInput make_problem1_solve_input(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    double launch_time_seconds_since_j2000,
    double transfer_perihelion_angle_global,
    const Problem1SolveDefaults& defaults
) {
    return problem1::Problem1SolveInput{
        departure,
        target,
        launch_time_seconds_since_j2000,
        transfer_perihelion_angle_global,
        defaults.max_transfer_revolution,
        defaults.max_target_revolution,
        defaults.phi_scan_count,
        defaults.phi_tolerance,
        defaults.residual_tolerance,
        defaults.max_bisection_iterations,
        defaults.max_candidate_relative_residual,
    };
}

// 按默认参数构造 3D 表格配置，将离散点数转换为固定角步长 2π/count。
problem1::Problem1TableConfig make_problem1_table_config(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    const Problem1TableDefaults& defaults
) {
    return problem1::Problem1TableConfig{
        departure,
        target,
        defaults.departure_true_anomaly_start,
        common::kTwoPi / static_cast<double>(defaults.departure_true_anomaly_count),
        defaults.departure_true_anomaly_count,
        defaults.target_true_anomaly_start,
        common::kTwoPi / static_cast<double>(defaults.target_true_anomaly_count),
        defaults.target_true_anomaly_count,
        defaults.transfer_theta_departure_start,
        common::kTwoPi / static_cast<double>(defaults.transfer_theta_departure_count),
        defaults.transfer_theta_departure_count,
        defaults.max_transfer_revolution,
        defaults.max_target_revolution,
    };
}

}  // namespace spaceship_cpp::config
