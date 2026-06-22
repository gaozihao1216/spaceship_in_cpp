/*
 * 文件作用：声明固定行星序列上的 BFS 轨迹搜索。
 * 主要工作：在固定发射时刻、初始 θ 与访问序列下，最小化 v_∞,launch + v_∞,arrival，并满足总用时上限。
 */
#pragma once

#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"
#include "spaceship_cpp/bfs/trajectory_search_state.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct FixedSequenceBfsStats {
    std::size_t leg0_candidates = 0;
    std::size_t expanded_nodes = 0;
    std::size_t enqueued_nodes = 0;
    std::size_t leaf_nodes = 0;

    std::size_t dead_by_time_limit = 0;
    std::size_t dead_by_no_p2_solution = 0;
    std::size_t dead_by_turn_angle = 0;
    std::size_t dead_by_invalid_v_inf = 0;
    std::size_t dead_by_launch_v_inf_limit = 0;

    std::size_t p2_solve_calls = 0;
};

struct FixedSequenceBfsConfig {
    double launch_time_seconds_since_j2000 = 0.0;
    double launch_transfer_theta_global = 0.0;

    std::vector<planet_params::PlanetId> planet_sequence;

    // 从发射到最终到达 sequence.back() 的累计 TOF 上限（秒）。
    double max_total_time_seconds = std::numeric_limits<double>::infinity();

    // 首段 Earth 发射相对行星速度 v_∞ 上限（m/s）；默认 7500。设为 +inf 表示不剪枝。
    double max_launch_v_inf_mps = 7500.0;

    Leg0Problem1Options problem1{};

    trajectory::FlybyPhysicalFeasibilityOptions flyby_physical_filter{};
};

struct FixedSequenceBfsResult {
    bool ok = false;
    std::string error_message;

    bool found_solution = false;
    double best_score = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = 0.0;
    double best_arrival_v_inf = 0.0;
    double best_total_time_seconds = 0.0;

    std::vector<TrajectorySearchEdge> best_edges;
    FixedSequenceBfsStats stats{};
};

// 固定序列 BFS：逐层展开部分路径，在用时约束下最小化 launch_v_inf + arrival_v_inf。
FixedSequenceBfsResult search_fixed_sequence_bfs(const FixedSequenceBfsConfig& config);

}  // namespace spaceship_cpp::bfs
