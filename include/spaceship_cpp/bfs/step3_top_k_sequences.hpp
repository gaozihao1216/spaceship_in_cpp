/*
 * 文件作用：声明 Step 3 行星序列 Top-K 筛选。
 * 主要工作：汇总 Step 2 自由路径结果，按唯一行星序列取最优 score，保留 Top-K 供 Step 4 细搜。
 */
#pragma once

#include "spaceship_cpp/bfs/free_path_bfs.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct Step3TopKStats {
    std::size_t solutions_pooled = 0;
    std::size_t unique_sequences = 0;
    int top_k_requested = 0;
    std::size_t sequences_returned = 0;
};

struct RankedPlanetSequence {
    std::vector<planet_params::PlanetId> planet_sequence;
    double best_score = 0.0;
    Leg0ThetaSeed best_seed{};
    FreePathBfsSolution best_solution{};
};

struct Step3TopKSequencesResult {
    bool ok = false;
    std::string error_message;
    std::vector<RankedPlanetSequence> sequences;
    Step3TopKStats stats{};
};

// 从 Step 2 结果中按 score（launch_v_inf + arrival_v_inf）升序保留 Top-K 条唯一行星序列。
Step3TopKSequencesResult select_top_k_planet_sequences(
    const Step2FreePathSearchResult& step2,
    int top_k
);

// 同上，但只考虑 visit_sequence 终点为 destination_planet 的序列（任务导向 Top-K）。
Step3TopKSequencesResult select_top_k_planet_sequences_reaching(
    const Step2FreePathSearchResult& step2,
    int top_k,
    planet_params::PlanetId destination_planet
);

// 使用 global.discretization.top_k_sequences 的 Step 3 封装。
Step3TopKSequencesResult run_step3_select_top_k_sequences(
    const TrajectorySearchGlobalConfig& config,
    const Step2FreePathSearchResult& step2
);

}  // namespace spaceship_cpp::bfs
