/*
 * 文件作用：声明 Problem 2 出射段 θ' 初扫与 branch 导数接口。
 * 主要工作：在固定飞掠时刻对 θ' 离散网格调用 Problem 1，并估计 dφ/dθ'、de/dθ'。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace spaceship_cpp::problem2 {

struct Problem2OutgoingBranchSolution {
    int transfer_revolution = 0;
    int target_revolution = 0;

    double encounter_global_angle = 0.0;
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double relative_problem1_residual = 0.0;

    double dphi_dtheta_prime = 0.0;
    double de_dtheta_prime = 0.0;
    bool has_dphi_dtheta_prime = false;
    bool has_de_dtheta_prime = false;
};

struct Problem2ThetaPrimeNodeSnapshot {
    double theta_prime_local = 0.0;
    double theta_prime_global = 0.0;
    std::vector<Problem2OutgoingBranchSolution> solutions;
};

struct Problem2ThetaPrimeScanConfig {
    planet_params::PlanetId flyby_planet = planet_params::PlanetId::Mars;
    planet_params::PlanetId target_planet = planet_params::PlanetId::Earth;

    double flyby_time_seconds_since_j2000 = 0.0;
    int theta_prime_count = 64;

    problem1::Problem1SolveInput problem1_solve{};

    // 相邻 θ' 节点间按 (k,q) 与 φ 配对时的最大允许角差，超过则不估计导数。
    double branch_phi_pairing_max_gap = 0.75;
};

struct Problem2ThetaPrimeInitialScanResult {
    bool ok = false;
    std::string error_message;
    std::vector<Problem2ThetaPrimeNodeSnapshot> nodes;
};

struct Problem2OutgoingBranchPair {
    std::size_t left_index = 0;
    std::size_t right_index = 0;
    double phi_gap = 0.0;
};

// 在 [-pi, pi) 上均匀生成 theta_prime_count 个局部角网格点。
std::vector<double> build_uniform_theta_prime_local_grid(int theta_prime_count);

// 标量对 theta' 的中心/单侧差分。
double estimate_theta_prime_derivative_central(
    double left_theta_prime,
    double left_value,
    double right_theta_prime,
    double right_value
);

double estimate_theta_prime_derivative_forward(
    double left_theta_prime,
    double left_value,
    double right_theta_prime,
    double right_value
);

double estimate_theta_prime_derivative_backward(
    double left_theta_prime,
    double left_value,
    double right_theta_prime,
    double right_value
);

// 在两组解之间按 (k,q) 与 φ 近邻做一一配对。
std::vector<Problem2OutgoingBranchPair> pair_outgoing_branch_solutions_by_phi(
    const std::vector<Problem2OutgoingBranchSolution>& left_solutions,
    const std::vector<Problem2OutgoingBranchSolution>& right_solutions,
    double max_phi_gap
);

// 在 right_solutions 中找与 reference 同属 (k,q) 且 φ 最近的解。
std::optional<std::size_t> find_best_matching_outgoing_branch_index(
    const Problem2OutgoingBranchSolution& reference,
    const std::vector<Problem2OutgoingBranchSolution>& candidates,
    double max_phi_gap
);

// 单个 θ' 节点：调用 Problem 1 收集全部出射根。
Problem2ThetaPrimeNodeSnapshot evaluate_problem2_theta_prime_node(
    const Problem2ThetaPrimeScanConfig& config,
    double theta_prime_local
);

// 完整初扫：默认 64 个 θ' 节点。
Problem2ThetaPrimeInitialScanResult run_problem2_theta_prime_initial_scan(
    const Problem2ThetaPrimeScanConfig& config
);

// 在初扫结果上为每条 branch 估计 dφ/dθ' 与 de/dθ'。
void attach_problem2_theta_prime_solution_derivatives(
    std::vector<Problem2ThetaPrimeNodeSnapshot>& nodes,
    double branch_phi_pairing_max_gap
);

}  // namespace spaceship_cpp::problem2
