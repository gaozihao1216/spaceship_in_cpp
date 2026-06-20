/*
 * 文件作用：声明 Problem 2 Route A 线性外推与 Newton 精修接口。
 * 主要工作：由端点 branch 解 + 导数估计目标 θ' 处的 φ，再以 Newton 收敛 P1 零点。
 */
#pragma once

#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <algorithm>
#include <string>

namespace spaceship_cpp::problem2 {

struct Problem2RouteALinearEstimateResult {
    bool valid = false;
    std::string invalid_reason;

    double endpoint_theta_prime_local = 0.0;
    double target_theta_prime_local = 0.0;
    double delta_theta_prime = 0.0;

    double predicted_phi = 0.0;
    double predicted_phi_unwrapped = 0.0;
    double predicted_eccentricity = 0.0;
    double predicted_semi_latus_rectum = 0.0;

    int original_target_revolution = 0;
    int adjusted_target_revolution = 0;
    int target_revolution_delta = 0;
    bool target_revolution_adjusted = false;
};

struct Problem2RouteATargetRevolutionAdjustment {
    bool valid = false;
    std::string invalid_reason;

    double omega_J = 0.0;
    double endpoint_phi = 0.0;
    double predicted_phi_unwrapped = 0.0;

    int original_target_revolution = 0;
    int adjusted_target_revolution = 0;
    int target_revolution_delta = 0;
};

struct Problem2RouteANewtonOptions {
    int max_iterations = 20;
    double residual_tolerance = 0.0;
    double max_relative_residual = 1e-6;
    double phi_tolerance = 1e-10;
    double phi_derivative_step = 1e-8;
    bool reject_on_residual_increase = true;
    bool adjust_target_revolution_on_omega_J_crossing = true;
};

enum class Problem2RouteANewtonStatus {
    Converged,
    DivergedResidualIncreased,
    DivergedInvalidResidual,
    DivergedSingularDerivative,
    DivergedMaxIterations,
};

struct Problem2RouteANewtonResult {
    bool ok = false;
    Problem2RouteANewtonStatus status = Problem2RouteANewtonStatus::DivergedMaxIterations;
    std::string status_reason;

    int iterations_used = 0;
    double target_theta_prime_local = 0.0;
    double target_theta_prime_global = 0.0;
    double final_phi = 0.0;

    problem1::Problem1ResidualResult residual_result{};
    Problem2OutgoingBranchSolution solution{};
};

// 用端点解与 dφ/dθ'、de/dθ' 对目标 θ' 做线性外推。
Problem2RouteALinearEstimateResult estimate_problem2_route_a_solution_linear(
    double endpoint_theta_prime_local,
    const Problem2OutgoingBranchSolution& endpoint_solution,
    double target_theta_prime_local
);

// 线性外推 φ 后，按是否穿越飞掠行星全局角 ω_J 修正 target_revolution q。
Problem2RouteATargetRevolutionAdjustment adjust_target_revolution_for_route_a_linear_prediction(
    planet_params::PlanetId flyby_planet,
    double flyby_time_seconds_since_j2000,
    const Problem2OutgoingBranchSolution& endpoint_solution,
    double delta_theta_prime,
    int max_target_revolution
);

// 将 q 修正写入线性外推结果（不改变 predicted_phi 的 [0,2π) 归一化值）。
void apply_target_revolution_adjustment_to_route_a_linear_estimate(
    const Problem2RouteATargetRevolutionAdjustment& adjustment,
    Problem2RouteALinearEstimateResult& linear_estimate
);

// 固定目标 θ'，从 initial_phi 出发对 P1 残差做 Newton 迭代。
Problem2RouteANewtonResult refine_problem2_route_a_phi_by_newton(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_local,
    double initial_phi,
    int transfer_revolution,
    int target_revolution,
    const Problem2RouteANewtonOptions& options
);

// Route A 完整一步：线性外推 + ω_J 穿越 q 修正 + P1-Newton。
Problem2RouteANewtonResult refine_problem2_route_a_at_theta_prime(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_local,
    const Problem2OutgoingBranchSolution& linear_endpoint_branch,
    double linear_endpoint_theta_prime_local,
    int transfer_revolution,
    const Problem2RouteANewtonOptions& options
);

}  // namespace spaceship_cpp::problem2
