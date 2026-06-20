/*
 * 文件作用：声明 Problem 2 Route A 线性外推与 Newton 精修接口。
 * 主要工作：由端点 branch 解 + 导数估计目标 θ' 处的 φ，再以 Newton 收敛 P1 零点。
 */
#pragma once

#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <string>

namespace spaceship_cpp::problem2 {

struct Problem2RouteALinearEstimateResult {
    bool valid = false;
    std::string invalid_reason;

    double endpoint_theta_prime_local = 0.0;
    double target_theta_prime_local = 0.0;
    double delta_theta_prime = 0.0;

    double predicted_phi = 0.0;
    double predicted_eccentricity = 0.0;
    double predicted_semi_latus_rectum = 0.0;
};

struct Problem2RouteANewtonOptions {
    int max_iterations = 20;
    double residual_tolerance = 0.0;
    double max_relative_residual = 1e-6;
    double phi_tolerance = 1e-10;
    double phi_derivative_step = 1e-8;
    bool reject_on_residual_increase = true;
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

// 固定目标 θ'，从 initial_phi 出发对 P1 残差做 Newton 迭代。
Problem2RouteANewtonResult refine_problem2_route_a_phi_by_newton(
    const Problem2ThetaPrimeScanConfig& config,
    double target_theta_prime_local,
    double initial_phi,
    int transfer_revolution,
    int target_revolution,
    const Problem2RouteANewtonOptions& options
);

}  // namespace spaceship_cpp::problem2
