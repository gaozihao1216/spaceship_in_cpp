/*
 * 文件作用：定义轨迹搜索使用的基础状态和边数据结构。
 * 主要工作：保存当前行星、时间、转移轨道参数和搜索扩展结果。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct TrajectorySearchState {
    // valid=false 表示该状态是占位或被过滤掉的无效节点。
    bool valid = false;

    // 当前轨迹已经到达并准备继续飞掠/转移的行星。
    planet_params::PlanetId current_planet = planet_params::PlanetId::Mercury;

    // 当前节点对应的绝对时间。
    double current_time = 0.0;

    // 到达当前行星前的入射轨道参数。
    double incoming_e = 0.0;
    double incoming_theta = 0.0;

    // 搜索深度和累计代价。
    int depth = 0;

    double accumulated_time_seconds = 0.0;
    double launch_v_inf = 0.0;
    double accumulated_score = 0.0;
};

struct TrajectorySearchEdge {
    // 表示一次从 from_planet 到 to_planet 的转移边。
    bool valid = false;

    planet_params::PlanetId from_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId to_planet = planet_params::PlanetId::Mercury;

    // 边的时间信息。
    double departure_time = 0.0;
    double arrival_time = 0.0;
    double transfer_time_seconds = 0.0;

    // 出射转移轨道参数。
    double outgoing_e = 0.0;
    double outgoing_p = 0.0;
    double outgoing_theta = 0.0;

    double theta_prime = 0.0;
    double alpha = 0.0;

    // 多圈分支和残差信息用于诊断这条边来自哪个数学分支。
    int transfer_revolution = 0;
    int target_revolution = 0;

    double slingshot_residual = 0.0;
    double problem1_residual_seconds = 0.0;

    bool boundary_ambiguous = false;
    bool origin_was_topology_change = false;
};

struct TrajectorySearchExpansionResult {
    // 一次扩展是否成功，以及失败时的简短原因。
    bool ok = false;
    std::string error_message;

    std::vector<TrajectorySearchEdge> edges;
    std::vector<TrajectorySearchState> next_states;
};

}  // namespace spaceship_cpp::bfs
