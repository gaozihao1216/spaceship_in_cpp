/*
 * 文件作用：声明固定发射时刻下 leg0 可行 θ 区间扫描。
 * 主要工作：在 [θ_start, θ_end) 上离散采样，找出至少有一个 Problem 1 候选通过 leg0 剪枝的 θ 范围。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct Leg0PruningLimits {
    // 与 FixedSequenceBfsConfig 中 leg0 剪枝语义一致。
    double max_total_time_seconds = std::numeric_limits<double>::infinity();
    double max_launch_v_inf_mps = 7500.0;
};

struct Leg0Problem1Options {
    int max_transfer_revolution = 0;
    int max_target_revolution = 0;

    // 小于 0 表示沿用 config::global_config().problem1_solve 中的默认值。
    int phi_scan_count_override = -1;
};

struct Leg0CandidateFilterStats {
    std::size_t raw_candidates = 0;
    std::size_t dead_by_non_positive_transfer_e_raw = 0;
    std::size_t dead_by_time_limit = 0;
    std::size_t dead_by_invalid_v_inf = 0;
    std::size_t dead_by_launch_v_inf_limit = 0;
    std::size_t passed = 0;
};

struct FeasibleThetaInterval {
    // 半开区间 [start_rad, end_rad)，由扫描网格合并得到。
    double start_rad = 0.0;
    double end_rad = 0.0;
};

struct Leg0ThetaScanConfig {
    double launch_time_seconds_since_j2000 = 0.0;
    planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;

    Leg0PruningLimits pruning{};
    Leg0Problem1Options problem1{};

    // 在 [theta_start_rad, theta_end_rad) 上均匀采样 theta_scan_count 个点。
    double theta_start_rad = 0.0;
    double theta_end_rad = 0.0;
    int theta_scan_count = 60;
};

struct Leg0ThetaFeasibilityResult {
    bool ok = false;
    std::string error_message;

    std::vector<FeasibleThetaInterval> feasible_intervals;
    std::size_t theta_samples_tested = 0;
    std::size_t theta_samples_feasible = 0;
};

// 判断单个 Problem 1 候选是否通过 leg0 剪枝（transfer_e_raw>0、TOF、发射 v_∞ 上限等）。
// transfer_e_raw>0 用于排除 (-e, θ+π) 与 (e, θ) 的等价重复分支。
bool leg0_candidate_passes_pruning(
    const problem1::Problem1Candidate& candidate,
    double launch_time_seconds_since_j2000,
    const Leg0PruningLimits& pruning,
    Leg0CandidateFilterStats* stats = nullptr
);

// 固定发射时刻与 θ 下，是否存在至少一个通过 leg0 剪枝的 Problem 1 候选。
bool leg0_has_feasible_candidate_at_theta(
    double launch_time_seconds_since_j2000,
    planet_params::PlanetId target_planet,
    double transfer_theta_global,
    const Leg0PruningLimits& pruning,
    const Leg0Problem1Options& problem1 = {},
    Leg0CandidateFilterStats* stats = nullptr
);

// 扫描 θ 并合并相邻可行采样点为区间。
Leg0ThetaFeasibilityResult find_leg0_feasible_theta_intervals(const Leg0ThetaScanConfig& config);

}  // namespace spaceship_cpp::bfs
