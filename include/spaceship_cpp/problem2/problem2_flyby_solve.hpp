/*
 * 文件作用：声明 Problem 2 飞掠求解顶层封装。
 * 主要工作：由入射轨道与行星对出发，自动初扫并搜索 G=0，返回出射轨道与转移时间。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::problem2 {

struct Problem2FlybySolveInput {
    planet_params::PlanetId flyby_planet = planet_params::PlanetId::Earth;
    planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;

    // 飞掠出发时刻（从飞掠行星出发段 Problem 1 的 launch_time）。
    double flyby_time_seconds_since_j2000 = 0.0;

    // 入射轨道偏心率与近日点角；默认局部角（以飞掠行星近日点为零点）。
    double incoming_eccentricity = 0.0;
    double incoming_theta = 0.0;
    bool incoming_theta_is_local = true;
};

struct Problem2FlybySolveConfig {
    Problem2FlybyGSearchConfig g_search_config{};
    bool collect_g_search_profile = false;
};

struct Problem2FlybySolution {
    Problem2FlybyGSolutionSource source = Problem2FlybyGSolutionSource::QuadraticNewton;

    // 出射段近日点角 θ'（局部 / 全局）。
    double outgoing_theta_prime_local = 0.0;
    double outgoing_theta_prime_global = 0.0;

    // 出射轨道几何与 Problem 1 遇合角 φ（全局）。
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double encounter_global_angle = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    // 飞掠段出发时刻、飞行时间与到达目标行星时刻。
    double flyby_time_seconds_since_j2000 = 0.0;
    double time_of_flight_seconds = 0.0;
    double arrival_time_seconds_since_j2000 = 0.0;

    // 弹弓约束残差 G = F_in - F_out；Problem 1 相对残差。
    double flyby_constraint_G = 0.0;
    double problem1_relative_residual = 0.0;
};

struct Problem2FlybySolveResult {
    bool ok = false;
    std::string error_message;
    std::vector<Problem2FlybySolution> solutions;

    bool has_g_search_profile = false;
    Problem2FlybyGSearchProfile g_search_profile{};
};

// 给定入射轨道与行星对，执行 θ' 初扫 + G=0 搜索，返回全部飞掠出射候选。
Problem2FlybySolveResult solve_problem2_flyby(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config
);

// 在已有 θ' 初扫结果上搜索 G=0（同一 flyby/target/t_J 下可复用 scan，换入射轨道不必重扫）。
Problem2FlybySolveResult solve_problem2_flyby_with_scan(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config,
    const Problem2ThetaPrimeInitialScanResult& scan
);

// 为飞掠求解执行 θ' 初扫；与 solve 内嵌扫描参数一致，供单次 expansion 内 scan→solve 流水线使用。
Problem2ThetaPrimeInitialScanResult run_problem2_flyby_theta_prime_initial_scan(
    const Problem2FlybySolveInput& input,
    const Problem2FlybySolveConfig& config
);

}  // namespace spaceship_cpp::problem2
