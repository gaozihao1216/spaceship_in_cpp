/*
 * 文件作用：声明 Step 4 固定序列 θ 细扫 + BFS 精搜。
 * 主要工作：对 Step 3 Top-K 行星序列，在对应首段目标的可行 θ 区间内细扫并调用 search_fixed_sequence_bfs。
 */
#pragma once

#include "spaceship_cpp/bfs/fixed_sequence_bfs.hpp"
#include "spaceship_cpp/bfs/step3_top_k_sequences.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct Step4SequenceFineSearchResult {
    std::vector<planet_params::PlanetId> planet_sequence;
    bool found_solution = false;
    std::string error_message;

    double best_theta_global = 0.0;
    double best_score = std::numeric_limits<double>::infinity();
    FixedSequenceBfsResult best_bfs{};

    std::size_t theta_samples_tried = 0;
    std::size_t theta_samples_with_solution = 0;
};

struct Step4FineSearchStats {
    std::size_t sequences_processed = 0;
    std::size_t sequences_with_solution = 0;
    std::size_t total_theta_samples = 0;
    std::size_t total_bfs_calls = 0;
};

struct Step4FineSearchResult {
    bool ok = false;
    std::string error_message;

    std::vector<Step4SequenceFineSearchResult> by_sequence;
    int best_sequence_index = -1;
    double global_best_score = std::numeric_limits<double>::infinity();

    Step4FineSearchStats stats{};
};

// 在给定 θ 样本上对单条固定序列做 BFS，返回最优结果。
Step4SequenceFineSearchResult search_fixed_sequence_over_theta_samples(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<planet_params::PlanetId>& planet_sequence,
    const std::vector<double>& theta_samples
);

// Step 4：对 Step 3 每条序列，在其首段目标 P₁ 的 leg0 可行区间内细扫。
Step4FineSearchResult run_step4_fixed_sequence_fine_search(
    const TrajectorySearchGlobalConfig& config,
    const Step3TopKSequencesResult& step3,
    const Leg0MultiTargetThetaResult& leg0
);

// Step 4 封装：内部先运行 Step 1 获取 leg0 可行区间。
Step4FineSearchResult run_step4_fixed_sequence_fine_search(
    const TrajectorySearchGlobalConfig& config,
    const Step3TopKSequencesResult& step3
);

}  // namespace spaceship_cpp::bfs
