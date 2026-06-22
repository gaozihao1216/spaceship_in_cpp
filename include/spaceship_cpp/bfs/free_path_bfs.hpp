/*
 * 文件作用：声明从 leg0 种子出发的自由路径 BFS（Step 2）。
 * 主要工作：固定 (t_launch, P₁, θ) 后枚举遭遇序列（允许重复访问同一行星），收集可行轨迹供 Step 3 筛选。
 */
#pragma once

#include "spaceship_cpp/bfs/trajectory_search_config.hpp"
#include "spaceship_cpp/bfs/trajectory_search_state.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"

#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct P2ExpansionTimingRecord {
    planet_params::PlanetId flyby_planet{};
    planet_params::PlanetId target_planet{};
    double flyby_time_seconds = 0.0;
    double incoming_e = 0.0;
    double incoming_theta = 0.0;
    double scan_ms = 0.0;
    double solve_ms = 0.0;
    bool scan_ok = false;
    bool solve_ok = false;
    std::size_t solution_count = 0;
};

struct FreePathBfsStats {
    std::size_t leg0_candidates = 0;
    std::size_t expanded_nodes = 0;
    std::size_t enqueued_nodes = 0;
    std::size_t leaf_nodes = 0;
    std::size_t recorded_solutions = 0;

    std::size_t dead_by_time_limit = 0;
    std::size_t dead_by_no_p2_solution = 0;
    std::size_t dead_by_turn_angle = 0;
    std::size_t dead_by_invalid_v_inf = 0;
    std::size_t dead_by_launch_v_inf_limit = 0;
    std::size_t dead_by_leg0_no_candidate = 0;

    std::size_t p2_solve_calls = 0;

    // 每次 expand_flyby_to_target（一次 P2 求解尝试）产生的子状态数分布。
    std::size_t p2_expansion_attempts = 0;
    std::map<std::size_t, std::size_t> p2_expansion_child_count_hist;

    // Step 2 分段时间（毫秒），用于 profiling。
    double leg0_solve_ms = 0.0;
    double p2_scan_ms = 0.0;
    double p2_solve_ms = 0.0;

    // 聚合每次 P2 solve 的 G-search profile（需 collect_g_search_profile）。
    std::size_t g_search_profile_samples = 0;
    problem2::Problem2FlybyGSearchProfile g_search_profile{};

    std::vector<P2ExpansionTimingRecord> p2_expansion_timings;
};

void record_p2_expansion_branching(FreePathBfsStats& stats, std::size_t child_count);

void merge_free_path_bfs_stats(FreePathBfsStats& into, const FreePathBfsStats& from);


void print_p2_expansion_branching_histogram(
    std::ostream& output,
    const FreePathBfsStats& stats
);

// 打印 Step 2 计数与分段时间占比。
void print_step2_timing_breakdown(
    std::ostream& output,
    const FreePathBfsStats& stats,
    double total_wall_ms
);

// 打印 Step 2 聚合 G-search profile（route_a / interval 计数等）。
void print_step2_g_search_profile_breakdown(
    std::ostream& output,
    const FreePathBfsStats& stats
);

struct FreePathBfsSolution {
    std::vector<planet_params::PlanetId> planet_sequence;
    double score = 0.0;
    double launch_v_inf = 0.0;
    double arrival_v_inf = 0.0;
    double total_time_seconds = 0.0;
    std::vector<TrajectorySearchEdge> edges;
};

struct FreePathBfsResult {
    bool ok = false;
    std::string error_message;

    Leg0ThetaSeed seed{};
    std::vector<FreePathBfsSolution> solutions;
    FreePathBfsStats stats{};
};

FreePathBfsStats aggregate_free_path_bfs_stats(
    const std::vector<FreePathBfsResult>& by_seed
);

std::vector<P2ExpansionTimingRecord> aggregate_p2_expansion_timings(
    const std::vector<FreePathBfsResult>& by_seed
);

struct Step2FreePathSearchResult {
    bool ok = false;
    std::string error_message;
    std::vector<Leg0ThetaSeed> seeds_processed;
    std::vector<FreePathBfsResult> by_seed;
};

// 单个 leg0 种子上的自由路径 BFS。
FreePathBfsResult search_free_path_from_leg0_seed(
    const TrajectorySearchGlobalConfig& config,
    const Leg0ThetaSeed& seed
);

// Step 2：对全部 leg0 种子运行自由路径 BFS。
Step2FreePathSearchResult run_step2_free_path_search(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<Leg0ThetaSeed>& seeds
);

// Step 1 + 离散化 + Step 2 一条龙。
Step2FreePathSearchResult run_step2_free_path_search(
    const TrajectorySearchGlobalConfig& config
);

}  // namespace spaceship_cpp::bfs
