/*
 * 文件作用：声明轨迹搜索管线的全局参量。
 * 主要工作：统一发射时刻、首段目标、物理约束、Problem 1 圈数限制与离散化策略，并生成各阶段子配置。
 */
#pragma once

#include "spaceship_cpp/bfs/fixed_sequence_bfs.hpp"
#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <vector>

namespace spaceship_cpp::bfs {

constexpr double kDefaultMaxTotalTimeSeconds = 40.0 * 365.25 * 86400.0;

// 任务级：固定发射时刻与首段 Problem 1 候选目标行星。
struct TrajectorySearchMission {
    double launch_time_seconds_since_j2000 = 0.0;

    // 第一步按首段目标分次扫描时使用的候选行星（默认内行星三目标）。
    std::vector<planet_params::PlanetId> first_leg_target_planets{
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };

    // Step 2 自由路径每段可选的飞掠/目标行星集合（不含 Earth；可重复访问同一行星）。
    std::vector<planet_params::PlanetId> visit_planets{
        planet_params::PlanetId::Mercury,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mars,
    };
};

// 物理与任务硬约束：全程 BFS / leg0 剪枝共用。
struct TrajectorySearchConstraints {
    double max_total_time_seconds = kDefaultMaxTotalTimeSeconds;
    double max_launch_v_inf_mps = 7500.0;
    trajectory::FlybyPhysicalFeasibilityOptions flyby_physical_filter{};
};

// Problem 1 预言机参数：圈数上限与可选扫描精度覆盖（全程搜索共用；默认 k=q=1）。
using TrajectorySearchProblem1Options = Leg0Problem1Options;

// 离散化与 Top-K 策略（由 global_config().trajectory_search 填充，不在此写死默认值）。
struct TrajectorySearchDiscretization {
    int leg0_theta_coarse_scan_count = 0;
    int leg0_theta_fine_scan_count = 0;
    int leg0_theta_seeds_per_first_leg_target = 0;
    int top_k_sequences = 0;
};

struct TrajectorySearchGlobalConfig {
    TrajectorySearchMission mission{};
    TrajectorySearchConstraints constraints{};
    TrajectorySearchProblem1Options problem1{};
    TrajectorySearchDiscretization discretization{};
    int max_search_legs = 0;

    // 为 true 时 Step 2 记录每次 P2 扩展的 scan/solve 耗时与入射状态（测试/profiling 用）。
    bool collect_p2_expansion_timings = false;

    // 为 true 时 Step 2 聚合每次 P2 solve 内的 G-search profile（route_a、区间分类等）。
    bool collect_g_search_profile = false;
};

struct Leg0TargetThetaFeasibility {
    planet_params::PlanetId target_planet{};
    Leg0ThetaFeasibilityResult result{};
};

struct Leg0MultiTargetThetaResult {
    bool ok = false;
    std::string error_message;
    std::vector<Leg0TargetThetaFeasibility> by_target;
};

// Step 2 种子：固定 (t_launch, 首段目标 P₁, θ)。
struct Leg0ThetaSeed {
    planet_params::PlanetId first_leg_target_planet{};
    double transfer_theta_global = 0.0;
};

// Step 1→2：在单段 [L,R) 内以间距 d 从 mid 向两侧扩展取点：mid±d/2, mid±3d/2, …
std::vector<double> discretize_feasible_theta_interval(
    const FeasibleThetaInterval& interval,
    double spacing_d_rad
);

std::vector<Leg0ThetaSeed> discretize_leg0_theta_seeds(
    const Leg0MultiTargetThetaResult& leg0_result,
    const TrajectorySearchDiscretization& discretization
);

std::vector<Leg0ThetaSeed> discretize_leg0_theta_seeds(
    const TrajectorySearchGlobalConfig& config,
    const Leg0MultiTargetThetaResult& leg0_result
);

// Step 4：在单个首段目标的可行 θ 区间内按 sample_budget 细扫取点。
std::vector<double> discretize_leg0_theta_samples_for_target(
    const Leg0TargetThetaFeasibility& entry,
    int sample_budget
);

std::vector<double> discretize_leg0_theta_samples_for_target(
    const Leg0MultiTargetThetaResult& leg0_result,
    planet_params::PlanetId first_leg_target_planet,
    int sample_budget
);

TrajectorySearchGlobalConfig default_trajectory_search_global_config();

TrajectorySearchDiscretization make_trajectory_search_discretization(
    const config::TrajectorySearchDefaults& defaults
);

Leg0PruningLimits make_leg0_pruning_limits(const TrajectorySearchConstraints& constraints);

config::Problem1SolveDefaults make_problem1_solve_defaults(
    const Leg0Problem1Options& problem1
);

config::Problem1SolveDefaults make_problem1_solve_defaults(
    const TrajectorySearchGlobalConfig& config
);

Leg0ThetaScanConfig make_leg0_theta_scan_config(
    const TrajectorySearchGlobalConfig& config,
    planet_params::PlanetId first_leg_target_planet
);

FixedSequenceBfsConfig make_fixed_sequence_bfs_config(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<planet_params::PlanetId>& planet_sequence,
    double launch_transfer_theta_global
);

// Step 1：对 mission.first_leg_target_planets 中每个首段目标分别扫描 leg0 可行 θ 区间。
Leg0MultiTargetThetaResult find_leg0_feasible_theta_for_first_leg_targets(
    const TrajectorySearchGlobalConfig& config
);

}  // namespace spaceship_cpp::bfs
