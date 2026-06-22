/*
 * 文件作用：BFS 轨迹搜索模块的公共入口。
 * 主要工作：
 *   - search_best_trajectory：给定发射时刻、出发/目标行星，返回最佳转移路径描述
 *   - 内部 Steps 1–4 管线头文件供测试与高级调用
 */
#pragma once

#include "spaceship_cpp/bfs/trajectory_solution.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"

namespace spaceship_cpp::bfs {

// 运行 Steps 1–4 完整管线，筛选 visit_sequence 终点为 destination 的最优解。
// 当前要求 departure_planet == Earth（与固定序列 BFS 一致）。
TrajectorySearchOutput search_best_trajectory(
    const TrajectorySearchInput& input,
    const TrajectorySearchGlobalConfig& config = default_trajectory_search_global_config()
);

}  // namespace spaceship_cpp::bfs

// --- 内部管线（Steps 1–4 分步 API，测试与调试使用）---
#include "spaceship_cpp/bfs/fixed_sequence_bfs.hpp"
#include "spaceship_cpp/bfs/free_path_bfs.hpp"
#include "spaceship_cpp/bfs/step3_top_k_sequences.hpp"
#include "spaceship_cpp/bfs/step4_fixed_sequence_fine_search.hpp"
#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"
#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/bfs/trajectory_search_state.hpp"
