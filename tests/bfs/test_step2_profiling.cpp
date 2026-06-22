/*
 * 文件作用：剖析 Step 2 耗时：BFS 分支规模 vs 单次 P1/P2 求解开销。
 */
#include "spaceship_cpp/bfs/free_path_bfs.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace {

constexpr double kDayInSeconds = 86400.0;
constexpr int kMicrobenchIterations = 5;

using Clock = std::chrono::steady_clock;
namespace bfs = spaceship_cpp::bfs;
namespace config = spaceship_cpp::config;
namespace planet_params = spaceship_cpp::planet_params;
namespace problem1 = spaceship_cpp::problem1;
namespace problem2 = spaceship_cpp::problem2;

using planet_params::PlanetId;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

template <typename Fn>
double bench_mean_ms(Fn&& fn, int warmup, int iterations) {
    for (int i = 0; i < warmup; ++i) {
        fn();
    }
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(iterations));
    for (int i = 0; i < iterations; ++i) {
        const Clock::time_point t0 = Clock::now();
        fn();
        samples.push_back(elapsed_ms(t0, Clock::now()));
    }
    return std::accumulate(samples.begin(), samples.end(), 0.0) /
        static_cast<double>(samples.size());
}

void run_step2_case(bfs::TrajectorySearchGlobalConfig global, int max_search_legs) {
    global.max_search_legs = max_search_legs;
    global.collect_g_search_profile = true;

    const auto leg0 = bfs::find_leg0_feasible_theta_for_first_leg_targets(global);
    const auto seeds = bfs::discretize_leg0_theta_seeds(global, leg0);

    const Clock::time_point t0 = Clock::now();
    const auto step2 = bfs::run_step2_free_path_search(global, seeds);
    const Clock::time_point t1 = Clock::now();
    const double wall_ms = elapsed_ms(t0, t1);

    const auto stats = bfs::aggregate_free_path_bfs_stats(step2.by_seed);

    std::cout << "\n======== max_search_legs=" << max_search_legs
              << ", seeds=" << seeds.size() << " ========\n";
    bfs::print_step2_timing_breakdown(std::cout, stats, wall_ms);
    bfs::print_step2_g_search_profile_breakdown(std::cout, stats);
    bfs::print_p2_expansion_branching_histogram(std::cout, stats);

    const double oracle_ms = stats.leg0_solve_ms + stats.p2_scan_ms + stats.p2_solve_ms;
    std::cout << std::fixed << std::setprecision(1);
    if (stats.p2_solve_ms >= stats.leg0_solve_ms &&
        stats.p2_solve_ms >= stats.p2_scan_ms) {
        std::cout << "Diagnosis: dominated by P2 solve ("
                  << (100.0 * stats.p2_solve_ms / std::max(wall_ms, 1.0))
                  << "%), calls=" << stats.p2_solve_calls << '\n';
    } else if (stats.expanded_nodes <= seeds.size() + 2U) {
        std::cout << "Diagnosis: tree barely expands; cost is mostly leg0 P1 per seed\n";
    } else if (stats.p2_expansion_attempts > stats.expanded_nodes * 2U) {
        std::cout << "Diagnosis: many P2 attempts per expanded node (branching), "
                  << "oracle_ms=" << oracle_ms << " vs wall=" << wall_ms << '\n';
    }
}

void run_microbenchmarks(const bfs::TrajectorySearchGlobalConfig& global) {
    const config::GlobalConfig& cfg = config::global_config();
    const auto p1_defaults = bfs::make_problem1_solve_defaults(global.problem1);

    std::cout << "\n======== isolated oracle microbenchmarks (mean of "
              << kMicrobenchIterations << ") ========\n";

    const double p1_mercury_ms = bench_mean_ms(
        [&]() {
            const auto input = config::make_problem1_solve_input(
                PlanetId::Earth,
                PlanetId::Mercury,
                global.mission.launch_time_seconds_since_j2000,
                4.1,
                p1_defaults);
            (void)problem1::solve_problem1(input);
        },
        1,
        kMicrobenchIterations);
    std::cout << std::fixed << std::setprecision(2)
              << "  P1 solve Earth->Mercury: " << p1_mercury_ms << " ms\n";

    const auto leg0_input = config::make_problem1_solve_input(
        PlanetId::Earth,
        PlanetId::Venus,
        global.mission.launch_time_seconds_since_j2000,
        1.0,
        p1_defaults);
    const auto leg0_candidates = problem1::solve_problem1(leg0_input);
    if (leg0_candidates.empty()) {
        std::cout << "  P2 microbench skipped: no Earth->Venus leg0 candidate\n";
        return;
    }
    const problem1::Problem1Candidate& sample = leg0_candidates.front();

    const auto p2_config = config::make_problem2_flyby_solve_config(
        PlanetId::Venus,
        PlanetId::Mercury,
        sample.arrival_time_seconds_since_j2000,
        cfg.problem2_theta_prime_scan,
        cfg.problem2_route_a_newton,
        p1_defaults);

    problem2::Problem2FlybySolveInput scan_input{};
    scan_input.flyby_planet = PlanetId::Venus;
    scan_input.target_planet = PlanetId::Mercury;
    scan_input.flyby_time_seconds_since_j2000 = sample.arrival_time_seconds_since_j2000;

    const double p2_scan_ms = bench_mean_ms(
        [&]() {
            (void)problem2::run_problem2_flyby_theta_prime_initial_scan(scan_input, p2_config);
        },
        1,
        kMicrobenchIterations);

    const auto scan =
        problem2::run_problem2_flyby_theta_prime_initial_scan(scan_input, p2_config);

    problem2::Problem2FlybySolveInput solve_input = scan_input;
    solve_input.incoming_eccentricity = sample.residual_result.transfer_e;
    solve_input.incoming_theta = sample.residual_result.transfer_perihelion_angle_used;
    solve_input.incoming_theta_is_local = false;

    const double p2_solve_ms = bench_mean_ms(
        [&]() {
            (void)problem2::solve_problem2_flyby_with_scan(solve_input, p2_config, scan);
        },
        1,
        kMicrobenchIterations);

    std::cout << "  P2 theta' initial scan Venus->Mercury: " << p2_scan_ms << " ms\n";
    std::cout << "  P2 flyby solve Venus->Mercury (reuse scan from bench): "
              << p2_solve_ms << " ms\n";
    std::cout << "  P2 scan+solve (cold): " << (p2_scan_ms + p2_solve_ms) << " ms\n";
}

}  // namespace

int main() {
    auto global = bfs::default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 40.0 * 365.25 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    std::cout << "=== Step 2 profiling (Earth departure, default seeds) ===\n";
    std::cout << "visit_planets=" << global.mission.visit_planets.size()
              << " (P2 attempts per expanded node)\n";

    run_microbenchmarks(global);

    for (const int legs : {1, 3, 6}) {
        run_step2_case(global, legs);
    }

    std::cout << "\n======== interpretation ========\n";
    std::cout << "Compare max_search_legs sweep: if wall_ms grows ~linearly with\n"
              << "p2_solve_calls while mean_ms_per_p2_solve stays flat,\n"
              << "bottleneck is BFS branching volume (many P2 calls).\n"
              << "If mean_ms_per_p2_solve is large vs microbench, overhead is inside P2.\n"
              << "If expanded_nodes ~ seeds and p2_expansion child_count=0 dominates,\n"
              << "tree dies early (not branch explosion).\n";
    std::cout << "\nNOTE: Venus->Mercury microbench uses ONE easy incoming orbit from leg0.\n"
              << "Run test_problem2_solve_timing (BFS replay) and test_problem2_random_timing\n"
              << "(random e,theta + inner-planet sequences) for representative distributions.\n";

    return 0;
}
