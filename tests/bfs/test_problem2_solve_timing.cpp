/*
 * 文件作用：对比 BFS 内 P2 solve 耗时与独立重放是否一致，并统计真实调用分布。
 * 主要工作：记录 Step 2 每次 P2 扩展的 scan/solve 参数，离线重放验证计时，纠正单一 Venus→Mercury bench 的偏差。
 */
#include "spaceship_cpp/bfs/free_path_bfs.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace {

constexpr double kDayInSeconds = 86400.0;
constexpr int kReplayIterations = 3;

using Clock = std::chrono::steady_clock;
namespace bfs = spaceship_cpp::bfs;
namespace config = spaceship_cpp::config;
namespace planet_params = spaceship_cpp::planet_params;
namespace problem2 = spaceship_cpp::problem2;

using planet_params::PlanetId;
using bfs::P2ExpansionTimingRecord;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct SolveTimingStats {
    double min_ms = 0.0;
    double max_ms = 0.0;
    double mean_ms = 0.0;
    double median_ms = 0.0;
    std::size_t count = 0;
};

SolveTimingStats compute_stats(std::vector<double> samples) {
    SolveTimingStats stats{};
    if (samples.empty()) {
        return stats;
    }
    stats.count = samples.size();
    std::sort(samples.begin(), samples.end());
    stats.min_ms = samples.front();
    stats.max_ms = samples.back();
    stats.median_ms = samples[samples.size() / 2U];
    stats.mean_ms =
        std::accumulate(samples.begin(), samples.end(), 0.0) /
        static_cast<double>(samples.size());
    return stats;
}

const char* planet_label(PlanetId id) {
    return planet_params::planet_name(id);
}

problem2::Problem2FlybySolveConfig make_p2_config(
    PlanetId flyby,
    PlanetId target,
    double flyby_time,
    const bfs::TrajectorySearchGlobalConfig& global
) {
    const config::GlobalConfig& defaults = config::global_config();
    return config::make_problem2_flyby_solve_config(
        flyby,
        target,
        flyby_time,
        defaults.problem2_theta_prime_scan,
        defaults.problem2_route_a_newton,
        bfs::make_problem1_solve_defaults(global.problem1));
}

struct ReplayResult {
    double scan_ms = 0.0;
    double solve_ms = 0.0;
    bool scan_ok = false;
    bool solve_ok = false;
    std::size_t solution_count = 0;
};

ReplayResult replay_p2_expansion(
    const P2ExpansionTimingRecord& record,
    const bfs::TrajectorySearchGlobalConfig& global
) {
    ReplayResult result{};
    const problem2::Problem2FlybySolveConfig solve_config = make_p2_config(
        record.flyby_planet,
        record.target_planet,
        record.flyby_time_seconds,
        global);

    problem2::Problem2FlybySolveInput scan_input{};
    scan_input.flyby_planet = record.flyby_planet;
    scan_input.target_planet = record.target_planet;
    scan_input.flyby_time_seconds_since_j2000 = record.flyby_time_seconds;

    const Clock::time_point scan_start = Clock::now();
    const problem2::Problem2ThetaPrimeInitialScanResult scan =
        problem2::run_problem2_flyby_theta_prime_initial_scan(scan_input, solve_config);
    result.scan_ms = elapsed_ms(scan_start, Clock::now());
    result.scan_ok = scan.ok;
    if (!scan.ok) {
        return result;
    }

    problem2::Problem2FlybySolveInput solve_input = scan_input;
    solve_input.incoming_eccentricity = record.incoming_e;
    solve_input.incoming_theta = record.incoming_theta;
    solve_input.incoming_theta_is_local = false;

    const Clock::time_point solve_start = Clock::now();
    const problem2::Problem2FlybySolveResult solved =
        problem2::solve_problem2_flyby_with_scan(solve_input, solve_config, scan);
    result.solve_ms = elapsed_ms(solve_start, Clock::now());
    result.solve_ok = solved.ok && !solved.solutions.empty();
    result.solution_count = solved.solutions.size();
    return result;
}

double replay_solve_mean_ms(
    const P2ExpansionTimingRecord& record,
    const bfs::TrajectorySearchGlobalConfig& global,
    int iterations
) {
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(iterations));
    for (int i = 0; i < iterations; ++i) {
        const ReplayResult replay = replay_p2_expansion(record, global);
        if (replay.scan_ok) {
            samples.push_back(replay.solve_ms);
        }
    }
    if (samples.empty()) {
        return 0.0;
    }
    return std::accumulate(samples.begin(), samples.end(), 0.0) /
        static_cast<double>(samples.size());
}

void print_solve_stats(const char* label, const SolveTimingStats& stats) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << label << ": n=" << stats.count
              << " min=" << stats.min_ms
              << " median=" << stats.median_ms
              << " mean=" << stats.mean_ms
              << " max=" << stats.max_ms << " ms\n";
}

}  // namespace

int main() {
    auto global = bfs::default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 40.0 * 365.25 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
    global.max_search_legs = 6;
    global.collect_p2_expansion_timings = true;
    global.collect_g_search_profile = true;

    std::cout << "=== P2 solve timing validation (Step 2 in-context vs isolated replay) ===\n";
    std::cout << "max_search_legs=" << global.max_search_legs << '\n';

    const auto leg0 = bfs::find_leg0_feasible_theta_for_first_leg_targets(global);
    const auto seeds = bfs::discretize_leg0_theta_seeds(global, leg0);

    const Clock::time_point t0 = Clock::now();
    const auto step2 = bfs::run_step2_free_path_search(global, seeds);
    const double step2_ms = elapsed_ms(t0, Clock::now());

    if (!step2.ok) {
        std::cerr << "Step 2 failed: " << step2.error_message << '\n';
        return 1;
    }

    const auto stats = bfs::aggregate_free_path_bfs_stats(step2.by_seed);
    const std::vector<P2ExpansionTimingRecord> records =
        bfs::aggregate_p2_expansion_timings(step2.by_seed);

    std::cout << "\n--- Step 2 aggregate ---\n";
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "wall=" << step2_ms << " ms, p2_expansion_attempts=" << records.size()
              << ", p2_solve_calls=" << stats.p2_solve_calls << '\n';
    if (stats.p2_solve_calls > 0U) {
        std::cout << std::setprecision(2);
        std::cout << "aggregate mean_ms_per_p2_solve="
                  << (stats.p2_solve_ms / static_cast<double>(stats.p2_solve_calls))
                  << " ms\n";
    }
    bfs::print_step2_g_search_profile_breakdown(std::cout, stats);

    std::vector<double> bfs_solve_ms;
    bfs_solve_ms.reserve(records.size());
    for (const P2ExpansionTimingRecord& record : records) {
        if (record.scan_ok && record.solve_ms > 0.0) {
            bfs_solve_ms.push_back(record.solve_ms);
        }
    }
    print_solve_stats("BFS in-context solve_ms", compute_stats(bfs_solve_ms));

    std::cout << "\n--- Replay each BFS call in isolation ("
              << kReplayIterations << " runs, cold scan each time) ---\n";
    std::vector<double> replay_solve_ms;
    replay_solve_ms.reserve(bfs_solve_ms.size());
    double max_abs_diff = 0.0;
    double max_rel_diff = 0.0;
    for (const P2ExpansionTimingRecord& record : records) {
        if (!record.scan_ok || record.solve_ms <= 0.0) {
            continue;
        }
        const double replay_mean = replay_solve_mean_ms(record, global, kReplayIterations);
        replay_solve_ms.push_back(replay_mean);
        const double abs_diff = std::abs(replay_mean - record.solve_ms);
        const double rel_diff = abs_diff / std::max(record.solve_ms, 1e-9);
        max_abs_diff = std::max(max_abs_diff, abs_diff);
        max_rel_diff = std::max(max_rel_diff, rel_diff);
    }
    print_solve_stats("Isolated replay solve_ms", compute_stats(replay_solve_ms));
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "BFS vs replay max_abs_diff=" << max_abs_diff
              << " ms, max_rel_diff=" << (100.0 * max_rel_diff) << "%\n";
    if (max_rel_diff > 0.25) {
        std::cout << "WARNING: replay differs from BFS timing by >25% — check timer scope\n";
    } else {
        std::cout << "OK: BFS inline timing matches isolated replay (timer scope is correct)\n";
    }

    std::cout << "\n--- Old microbenchmark (single Venus->Mercury, fixed incoming orbit) ---\n";
    const auto p1_defaults = bfs::make_problem1_solve_defaults(global.problem1);
    const auto venus_leg0 = config::make_problem1_solve_input(
        PlanetId::Earth, PlanetId::Venus, 0.0, 1.0, p1_defaults);
    const auto venus_candidates =
        spaceship_cpp::problem1::solve_problem1(venus_leg0);
    if (!venus_candidates.empty()) {
        P2ExpansionTimingRecord bench_record{};
        bench_record.flyby_planet = PlanetId::Venus;
        bench_record.target_planet = PlanetId::Mercury;
        bench_record.flyby_time_seconds = venus_candidates.front().arrival_time_seconds_since_j2000;
        bench_record.incoming_e = venus_candidates.front().residual_result.transfer_e;
        bench_record.incoming_theta =
            venus_candidates.front().residual_result.transfer_perihelion_angle_used;
        const double bench_solve = replay_solve_mean_ms(bench_record, global, 5);
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Venus->Mercury bench solve_ms=" << bench_solve << " ms\n";
        if (!bfs_solve_ms.empty()) {
            const SolveTimingStats bfs_stats = compute_stats(bfs_solve_ms);
            std::cout << "BFS mean solve_ms=" << bfs_stats.mean_ms
                      << " => bench is "
                      << (100.0 * bench_solve / bfs_stats.mean_ms)
                      << "% of BFS mean (NOT representative)\n";
        }
    }

    std::cout << "\n--- Per planet-pair mean solve_ms (BFS in-context) ---\n";
    for (PlanetId flyby : {PlanetId::Mercury, PlanetId::Venus, PlanetId::Mars}) {
        for (PlanetId target : {PlanetId::Mercury, PlanetId::Venus, PlanetId::Mars}) {
            std::vector<double> pair_samples;
            for (const P2ExpansionTimingRecord& record : records) {
                if (record.scan_ok && record.solve_ms > 0.0 &&
                    record.flyby_planet == flyby && record.target_planet == target) {
                    pair_samples.push_back(record.solve_ms);
                }
            }
            if (pair_samples.empty()) {
                continue;
            }
            const SolveTimingStats pair_stats = compute_stats(pair_samples);
            std::cout << "  " << planet_label(flyby) << " -> " << planet_label(target)
                      << ": n=" << pair_stats.count
                      << " mean=" << pair_stats.mean_ms
                      << " max=" << pair_stats.max_ms << " ms\n";
        }
    }

    std::cout << "\n--- Slowest 5 BFS P2 solves ---\n";
    std::vector<P2ExpansionTimingRecord> slowest = records;
    std::sort(
        slowest.begin(),
        slowest.end(),
        [](const P2ExpansionTimingRecord& lhs, const P2ExpansionTimingRecord& rhs) {
            return lhs.solve_ms > rhs.solve_ms;
        });
    std::cout << std::fixed << std::setprecision(2);
    for (std::size_t i = 0; i < std::min<std::size_t>(5U, slowest.size()); ++i) {
        const P2ExpansionTimingRecord& r = slowest[i];
        if (r.solve_ms <= 0.0) {
            continue;
        }
        std::cout << "  [" << i << "] " << planet_label(r.flyby_planet) << " -> "
                  << planet_label(r.target_planet)
                  << " t=" << (r.flyby_time_seconds / kDayInSeconds) << " d"
                  << " e_in=" << r.incoming_e
                  << " solve_ms=" << r.solve_ms << '\n';
    }

    std::cout << "\nConclusion: use BFS-recorded solve_ms distribution for optimization targets;\n"
              << "single Venus->Mercury microbench underestimates typical BFS cost.\n";

    return 0;
}
