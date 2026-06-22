/*
 * 文件作用：统计 Step 2 自由路径 BFS 中 Problem 2 扩展分支分布。
 */
#include "spaceship_cpp/bfs/free_path_bfs.hpp"

#include <chrono>
#include <iostream>
#include <limits>

namespace {

constexpr double kDayInSeconds = 86400.0;

using spaceship_cpp::bfs::aggregate_free_path_bfs_stats;
using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::discretize_leg0_theta_seeds;
using spaceship_cpp::bfs::find_leg0_feasible_theta_for_first_leg_targets;
using spaceship_cpp::bfs::print_p2_expansion_branching_histogram;
using spaceship_cpp::bfs::run_step2_free_path_search;

double elapsed_ms(
    const std::chrono::steady_clock::time_point& start,
    const std::chrono::steady_clock::time_point& end
) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

}  // namespace

int main() {
    auto global = default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 40.0 * 365.25 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    for (const int max_search_legs : {3}) {
        auto case_global = global;
        case_global.max_search_legs = max_search_legs;

        const auto leg0 = find_leg0_feasible_theta_for_first_leg_targets(case_global);
        const auto seeds = discretize_leg0_theta_seeds(case_global, leg0);

        const auto start = std::chrono::steady_clock::now();
        const auto step2 = run_step2_free_path_search(case_global, seeds);
        const auto end = std::chrono::steady_clock::now();

        const auto stats = aggregate_free_path_bfs_stats(step2.by_seed);
        std::cout << "=== default k=q=1 max_search_legs=" << max_search_legs
                  << " seeds=" << seeds.size()
                  << " ms=" << elapsed_ms(start, end)
                  << " expanded=" << stats.expanded_nodes
                  << " solutions=" << stats.recorded_solutions << " ===\n";
        print_p2_expansion_branching_histogram(std::cout, stats);
    }

    return 0;
}
