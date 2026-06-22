/*
 * 文件作用：分步统计轨迹搜索管线 Step 1–4 耗时。
 */
#include "spaceship_cpp/bfs/free_path_bfs.hpp"
#include "spaceship_cpp/bfs/step3_top_k_sequences.hpp"
#include "spaceship_cpp/bfs/step4_fixed_sequence_fine_search.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"

#include <chrono>
#include <iomanip>
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
using spaceship_cpp::bfs::run_step3_select_top_k_sequences;
using spaceship_cpp::bfs::run_step4_fixed_sequence_fine_search;

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

void print_ms(const char* label, double ms) {
    std::cout << std::fixed << std::setprecision(1)
              << label << ": " << ms << " ms (" << (ms / 1000.0) << " s)\n";
}

}  // namespace

int main() {
    auto global = default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 40.0 * 365.25 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    std::cout << "Pipeline config: k=q=1"
              << " coarse_scan=" << global.discretization.leg0_theta_coarse_scan_count
              << " seeds_per_target=" << global.discretization.leg0_theta_seeds_per_first_leg_target
              << " max_search_legs=" << global.max_search_legs
              << " top_k=" << global.discretization.top_k_sequences
              << " fine_scan=" << global.discretization.leg0_theta_fine_scan_count
              << '\n';

    const Clock::time_point total_start = Clock::now();

    const Clock::time_point t0 = Clock::now();
    const auto leg0 = find_leg0_feasible_theta_for_first_leg_targets(global);
    const Clock::time_point t1 = Clock::now();
    if (!leg0.ok) {
        std::cerr << "Step 1 failed: " << leg0.error_message << '\n';
        return 1;
    }

    const auto seeds = discretize_leg0_theta_seeds(global, leg0);
    const Clock::time_point t2 = Clock::now();

    const auto step2 = run_step2_free_path_search(global, seeds);
    const Clock::time_point t3 = Clock::now();
    if (!step2.ok) {
        std::cerr << "Step 2 failed: " << step2.error_message << '\n';
        return 1;
    }

    const auto step3 = run_step3_select_top_k_sequences(global, step2);
    const Clock::time_point t4 = Clock::now();

    const Clock::time_point t5 = Clock::now();
    spaceship_cpp::bfs::Step4FineSearchResult step4{};
    if (step3.sequences.empty()) {
        std::cout << "Step 4 skipped: no sequences from Step 3\n";
    } else {
        step4 = run_step4_fixed_sequence_fine_search(global, step3, leg0);
    }
    const Clock::time_point t6 = Clock::now();

    const double ms_step1 = elapsed_ms(t0, t1);
    const double ms_discretize = elapsed_ms(t1, t2);
    const double ms_step2 = elapsed_ms(t2, t3);
    const double ms_step3 = elapsed_ms(t3, t4);
    const double ms_step4 = step3.sequences.empty() ? 0.0 : elapsed_ms(t5, t6);
    const double ms_total = elapsed_ms(total_start, t6);

    std::cout << "\n--- Timing ---\n";
    print_ms("Step 1 (leg0 theta scan)", ms_step1);
    print_ms("Discretize (theta seeds)", ms_discretize);
    print_ms("Step 2 (free path BFS)", ms_step2);
    print_ms("Step 3 (Top-K sequences)", ms_step3);
    if (step3.sequences.empty()) {
        std::cout << "Step 4 (fixed seq + fine theta): skipped (0 ms)\n";
    } else {
        print_ms("Step 4 (fixed seq + fine theta)", ms_step4);
    }
    print_ms("Total", ms_total);

    std::cout << "\n--- Step 1 summary ---\n";
    for (const auto& entry : leg0.by_target) {
        std::cout << "  target=" << static_cast<int>(entry.target_planet)
                  << " feasible_samples=" << entry.result.theta_samples_feasible
                  << " intervals=" << entry.result.feasible_intervals.size() << '\n';
    }

    std::cout << "\n--- Step 2 seeds ---\n";
    std::cout << "  count=" << seeds.size() << '\n';
    for (const auto& seed : seeds) {
        std::cout << "  target=" << static_cast<int>(seed.first_leg_target_planet)
                  << " theta=" << seed.transfer_theta_global << '\n';
    }

    const auto step2_stats = aggregate_free_path_bfs_stats(step2.by_seed);
    std::cout << "\n--- Step 2 stats ---\n";
    std::cout << "  expanded=" << step2_stats.expanded_nodes
              << " solutions=" << step2_stats.recorded_solutions
              << " p2_solve_calls=" << step2_stats.p2_solve_calls << '\n';
    print_p2_expansion_branching_histogram(std::cout, step2_stats);

    std::cout << "\n--- Step 3 ---\n";
    std::cout << "  ok=" << step3.ok
              << " sequences_returned=" << step3.sequences.size();
    if (!step3.error_message.empty()) {
        std::cout << " note=" << step3.error_message;
    }
    std::cout << '\n';
    for (std::size_t i = 0; i < step3.sequences.size(); ++i) {
        const auto& seq = step3.sequences[i];
        std::cout << "  [" << i << "] score=" << seq.best_score
                  << " len=" << seq.planet_sequence.size() << '\n';
    }

    if (!step3.sequences.empty()) {
        std::cout << "\n--- Step 4 ---\n";
        std::cout << "  ok=" << step4.ok;
        if (!step4.error_message.empty()) {
            std::cout << " note=" << step4.error_message;
        }
        std::cout << "\n  sequences_with_solution=" << step4.stats.sequences_with_solution
                  << " total_bfs_calls=" << step4.stats.total_bfs_calls
                  << " best_index=" << step4.best_sequence_index << '\n';
        if (step4.best_sequence_index >= 0) {
            std::cout << "  global_best_score=" << step4.global_best_score << '\n';
        }
    }

    return 0;
}
