/*
 * 文件作用：声明 Problem 2 统一 G(θ') 搜索管线。
 * 主要工作：解析 dG/dθ'、二次求根初值、G-Newton 精化，并在 k 层配对 branch 上搜 G=0。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_constraint.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_route_a.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::problem2 {

struct Problem2FlybyGContext {
    FlybyConstraintIncomingCache incoming_cache{};
};

struct Problem2FlybyGSample {
    bool valid = false;
    std::string invalid_reason;
    double G = 0.0;
    double dG_dtheta_prime = 0.0;
};

struct Problem2GNewtonOptions {
    int max_iterations = 20;
    double G_tolerance = 1e-8;
    double theta_prime_tolerance = 1e-10;
    bool reject_on_G_increase = true;
};

enum class Problem2GNewtonStatus {
    Converged,
    DivergedInvalidRouteA,
    DivergedInvalidG,
    DivergedSingularDerivative,
    DivergedGIncreased,
    DivergedMaxIterations,
};

struct Problem2GNewtonResult {
    bool ok = false;
    Problem2GNewtonStatus status = Problem2GNewtonStatus::DivergedMaxIterations;
    std::string status_reason;

    int iterations_used = 0;
    double final_theta_prime_local = 0.0;
    double final_G = 0.0;
    Problem2OutgoingBranchSolution outgoing_branch{};
};

enum class Problem2FlybyGSolutionSource {
    NearZeroEndpoint,
    QuadraticNewton,
};

struct Problem2FlybyGSolution {
    Problem2FlybyGSolutionSource source = Problem2FlybyGSolutionSource::QuadraticNewton;
    double theta_prime_local = 0.0;
    double G = 0.0;
    Problem2OutgoingBranchSolution outgoing_branch{};
};

struct Problem2FlybyGBranchIntervalInput {
    double theta_prime_left = 0.0;
    double theta_prime_right = 0.0;
    Problem2OutgoingBranchSolution left_branch{};
    Problem2OutgoingBranchSolution right_branch{};
};

// G 搜索 profiling：计数 + 叶节点耗时（Case C 中点已走 Route A，case_c_middle_ms 恒为 0）。
struct Problem2FlybyGSearchProfile {
    double incoming_cache_ms = 0.0;
    double route_a_ms = 0.0;
    double case_c_middle_ms = 0.0;
    double enrich_ms = 0.0;

    std::size_t scan_node_count = 0;
    int max_transfer_revolution = 0;

    std::size_t interval_visits = 0;
    std::size_t interval_equal = 0;
    std::size_t interval_case_b = 0;
    std::size_t interval_case_c = 0;
    std::size_t interval_discarded = 0;

    std::size_t branch_interval_process_calls = 0;
    std::size_t case_b_probe_calls = 0;
    std::size_t case_c_middle_calls = 0;

    std::size_t route_a_calls = 0;
    std::size_t route_a_iterations = 0;
    std::size_t g_newton_calls = 0;
    std::size_t g_newton_iterations = 0;

    // 顶层 63 个 θ′ 区间（recursion_depth==0）分类计数。
    std::size_t top_interval_equal = 0;
    std::size_t top_interval_case_b = 0;
    std::size_t top_interval_case_c = 0;

    // Case C 成功取到中点 branch 后：n_middle 相对 n_left/n_right 与左右子区间类型。
    std::size_t case_c_split_samples = 0;
    std::size_t case_c_endpoint_gap_sum = 0;
    std::size_t case_c_middle_gt_max_endpoints = 0;
    std::size_t case_c_middle_lt_min_endpoints = 0;
    std::size_t case_c_middle_in_endpoint_range = 0;
    std::size_t case_c_child_left_equal = 0;
    std::size_t case_c_child_left_case_b = 0;
    std::size_t case_c_child_left_case_c = 0;
    std::size_t case_c_child_right_equal = 0;
    std::size_t case_c_child_right_case_b = 0;
    std::size_t case_c_child_right_case_c = 0;

    // Case B 同理（Route A 中点）。
    std::size_t case_b_split_samples = 0;
    std::size_t case_b_middle_gt_max_endpoints = 0;
    std::size_t case_b_middle_lt_min_endpoints = 0;
    std::size_t case_b_middle_in_endpoint_range = 0;
    std::size_t case_b_child_left_case_c = 0;
    std::size_t case_b_child_right_case_c = 0;
};

void merge_problem2_flyby_G_search_profile(
    Problem2FlybyGSearchProfile& into,
    const Problem2FlybyGSearchProfile& from
);

struct Problem2FlybyGSearchConfig {
    Problem2ThetaPrimeScanConfig scan_config{};
    Problem2RouteANewtonOptions route_a_newton_options{};
    Problem2GNewtonOptions g_newton_options{};
    double near_zero_G_threshold = 1e-4;
    double solution_theta_prime_tolerance = 1e-8;
    // 情形 B/C 递归细分：区间宽度低于此值时舍去（端点大概率不是 G=0）。
    double theta_prime_tolerance = 1e-8;
    int max_recursion_depth = 32;

    // 非空时记录 G 搜索各路径计数与叶节点耗时（测试/profiling 用）。
    Problem2FlybyGSearchProfile* profile = nullptr;
};

enum class Problem2ThetaPrimeIntervalCase {
    EqualBranchCount,
    BranchCountDifferenceOne,
    BranchCountDifferenceGreaterThanOne,
};

struct Problem2CaseCMiddleSolveResult {
    bool ok = false;
    std::string status_reason;

    double theta_prime_middle = 0.0;
    std::vector<Problem2OutgoingBranchSolution> middle_branches;
};

// 按同一 k 层端点 branch 个数关系分类区间。
Problem2ThetaPrimeIntervalCase classify_problem2_theta_prime_interval_case(
    std::size_t left_branch_count,
    std::size_t right_branch_count
);

// 情形 C：|n_L - n_R| > 1 时，在 θ'_mid 处用 Route A（同情形 B）估计中点 branch。
Problem2CaseCMiddleSolveResult solve_case_c_middle_branches_on_k_layer(
    const Problem2ThetaPrimeScanConfig& config,
    int transfer_revolution,
    double theta_prime_left,
    const std::vector<Problem2OutgoingBranchSolution>& left_branches,
    double theta_prime_right,
    const std::vector<Problem2OutgoingBranchSolution>& right_branches
);

struct Problem2FlybyGSearchResult {
    bool ok = false;
    std::string error_message;
    std::vector<Problem2FlybyGSolution> solutions;
};

// 由入射轨道与飞掠时刻构造 G 上下文（φ 取行星全局方向）。
Problem2FlybyGContext build_problem2_flyby_G_context(
    planet_params::PlanetId flyby_planet,
    double flyby_time_seconds_since_j2000,
    double incoming_eccentricity,
    double incoming_theta,
    bool incoming_theta_is_local
);

// 在固定 branch 轨迹上计算 G 与 dG/dθ'（链式法则 + ∂F/∂e、∂F/∂θ）。
Problem2FlybyGSample evaluate_flyby_constraint_G_and_derivative(
    const FlybyConstraintIncomingCache& incoming_cache,
    double orbit_eccentricity,
    double theta_prime_local,
    double dphi_dtheta_prime,
    double de_dtheta_prime
);

// 仅求 G，用于端点 near-zero 判断。
std::optional<double> evaluate_flyby_constraint_G_at_branch(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2OutgoingBranchSolution& branch,
    double theta_prime_local
);

// 固定 (k,q) branch，从二次预测初值出发对 G 做 Newton 精化。
Problem2GNewtonResult refine_flyby_constraint_G_zero_by_newton(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2ThetaPrimeScanConfig& scan_config,
    const Problem2RouteANewtonOptions& route_a_options,
    const Problem2GNewtonOptions& options,
    double initial_theta_prime_local,
    const Problem2OutgoingBranchSolution& reference_branch,
    const Problem2OutgoingBranchSolution& linear_endpoint_branch,
    double linear_endpoint_theta_prime_local,
    Problem2FlybyGSearchProfile* profile = nullptr
);

// 单条配对 branch 区间：端点 near-zero 收录 + 一律二次求根 + G-Newton。
std::vector<Problem2FlybyGSolution> process_flyby_constraint_G_branch_interval(
    const FlybyConstraintIncomingCache& incoming_cache,
    const Problem2FlybyGSearchConfig& config,
    const Problem2FlybyGBranchIntervalInput& interval
);

// 在初扫结果的单个 k 层上搜索全部 G=0 解。
Problem2FlybyGSearchResult search_flyby_constraint_G_zeros_on_k_layer(
    const Problem2FlybyGContext& context,
    const Problem2FlybyGSearchConfig& config,
    const Problem2ThetaPrimeInitialScanResult& scan,
    int transfer_revolution
);

// 在初扫结果的全部 k 层上搜索 G=0 解。
Problem2FlybyGSearchResult search_flyby_constraint_G_zeros_from_initial_scan(
    const Problem2FlybyGContext& context,
    const Problem2FlybyGSearchConfig& config,
    const Problem2ThetaPrimeInitialScanResult& scan
);

}  // namespace spaceship_cpp::problem2
