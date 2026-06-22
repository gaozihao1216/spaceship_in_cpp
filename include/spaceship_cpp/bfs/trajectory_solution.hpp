/*
 * 文件作用：声明 BFS 轨迹搜索的输入与输出描述类型。
 * 主要工作：给定发射时刻与起终点行星，返回最佳转移路径及各段轨道参数 (e, p, theta)。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

// 单段转移轨道的可可视化描述（Heliocentric 出射椭圆参数）。
struct TransferLegDescriptor {
    planet_params::PlanetId from_planet{};
    planet_params::PlanetId to_planet{};

    double departure_time_seconds_since_j2000 = 0.0;
    double arrival_time_seconds_since_j2000 = 0.0;
    double time_of_flight_seconds = 0.0;

    // 出射转移轨道：偏心率 e、半通径 p (AU)、全局近日点幅角 theta (rad)。
    double eccentricity = 0.0;
    double semi_latus_rectum_au = 0.0;
    double perihelion_angle_global_rad = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    // 飞掠段（from != Earth）附加参数；直达段可忽略。
    double flyby_theta_prime_local_rad = 0.0;
    double encounter_angle_rad = 0.0;
    double flyby_constraint_G = 0.0;
};

// 任务输入：发射时刻 + 出发行星 + 目标行星。
struct TrajectorySearchInput {
    double launch_time_seconds_since_j2000 = 0.0;
    planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    planet_params::PlanetId destination_planet = planet_params::PlanetId::Mercury;
};

// 搜索输出：完整访问序列、各段转移描述与代价。
struct TrajectorySearchOutput {
    bool ok = false;
    bool found_solution = false;
    std::string error_message;

    planet_params::PlanetId departure_planet{};
    planet_params::PlanetId destination_planet{};

    std::vector<planet_params::PlanetId> visit_sequence;
    std::vector<TransferLegDescriptor> legs;

    double score = std::numeric_limits<double>::infinity();
    double launch_v_inf_mps = 0.0;
    double arrival_v_inf_mps = 0.0;
    double total_time_seconds = 0.0;
    double leg0_theta_global_rad = 0.0;
};

}  // namespace spaceship_cpp::bfs
