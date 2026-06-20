/*
 * 文件作用：以 ms 为单位测量 Problem 2 飞掠求解耗时，对比是否复用初扫。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <optional>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using spaceship_cpp::planet_params::PlanetId;
using Clock = std::chrono::steady_clock;

struct TimingMs {
    double mean = 0.0;
    double min = 0.0;
    double max = 0.0;
    int samples = 0;
};

TimingMs bench_ms(int warmup, int runs, const auto& fn) {
    for (int i = 0; i < warmup; ++i) {
        fn();
    }

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(runs));
    for (int i = 0; i < runs; ++i) {
        const auto t0 = Clock::now();
        fn();
        const auto t1 = Clock::now();
        samples.push_back(
            std::chrono::duration<double, std::milli>(t1 - t0).count());
    }

    TimingMs stats{};
    stats.samples = runs;
    stats.min = samples.front();
    stats.max = samples.front();
    double sum = 0.0;
    for (const double sample : samples) {
        sum += sample;
        stats.min = std::min(stats.min, sample);
        stats.max = std::max(stats.max, sample);
    }
    stats.mean = sum / static_cast<double>(runs);
    return stats;
}

std::optional<spaceship_cpp::problem2::Problem2OutgoingBranchSolution>
find_solution_with_derivatives(
    const spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult& scan,
    std::size_t& out_node_index
) {
    for (std::size_t node_index = 1; node_index + 1 < scan.nodes.size(); ++node_index) {
        for (const auto& solution : scan.nodes[node_index].solutions) {
            if (solution.has_dphi_dtheta_prime && solution.has_de_dtheta_prime) {
                out_node_index = node_index;
                return solution;
            }
        }
    }
    return std::nullopt;
}

void print_timing(const char* label, const TimingMs& stats) {
    std::printf(
        "  %-32s mean=%8.2f ms  min=%8.2f ms  max=%8.2f ms  (n=%d)\n",
        label,
        stats.mean,
        stats.min,
        stats.max,
        stats.samples);
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;
    namespace problem2 = spaceship_cpp::problem2;

    constexpr int kWarmup = 2;
    constexpr int kRuns = 8;

    const auto& defaults = config::global_config();
    const auto solve_config = config::make_problem2_flyby_solve_config(
        PlanetId::Earth,
        PlanetId::Mars,
        0.0,
        defaults.problem2_theta_prime_scan,
        defaults.problem2_route_a_newton,
        defaults.problem1_solve);

    problem2::Problem2FlybySolveInput input{};
    input.flyby_planet = PlanetId::Earth;
    input.target_planet = PlanetId::Mars;
    input.flyby_time_seconds_since_j2000 = 0.0;

    const auto scan = problem2::run_problem2_flyby_theta_prime_initial_scan(input, solve_config);
    if (!scan.ok) {
        std::cerr << "initial scan failed\n";
        return 1;
    }

    std::size_t node_index = 0;
    const auto branch = find_solution_with_derivatives(scan, node_index);
    if (!branch.has_value()) {
        std::cerr << "no branch with derivatives\n";
        return 1;
    }

    input.incoming_eccentricity = branch->outgoing_eccentricity;
    input.incoming_theta = scan.nodes[node_index].theta_prime_local;
    input.incoming_theta_is_local = true;

    const TimingMs scan_only = bench_ms(kWarmup, kRuns, [&]() {
        (void)problem2::run_problem2_flyby_theta_prime_initial_scan(input, solve_config);
    });

    const TimingMs solve_with_scan_only = bench_ms(kWarmup, kRuns, [&]() {
        (void)problem2::solve_problem2_flyby_with_scan(input, solve_config, scan);
    });

    const TimingMs full_solve = bench_ms(kWarmup, kRuns, [&]() {
        (void)problem2::solve_problem2_flyby(input, solve_config);
    });

    const TimingMs reused_pipeline = bench_ms(kWarmup, kRuns, [&]() {
        const auto local_scan =
            problem2::run_problem2_flyby_theta_prime_initial_scan(input, solve_config);
        (void)problem2::solve_problem2_flyby_with_scan(input, solve_config, local_scan);
    });

    std::printf("problem2_flyby_solve timing (Earth->Mars)\n");
    std::printf(
        "  config: theta_prime_count=%d scan_phi_count=%d solve_phi_count=%d max_k=%d max_q=%d",
        defaults.problem2_theta_prime_scan.theta_prime_count,
        defaults.problem2_theta_prime_scan.phi_scan_count,
        defaults.problem1_solve.phi_scan_count,
        defaults.problem1_solve.max_transfer_revolution,
        defaults.problem1_solve.max_target_revolution);
    std::printf("\n");
#ifdef _OPENMP
    std::printf("  openmp_threads=%d\n", omp_get_max_threads());
#endif
    print_timing("theta_prime_initial_scan", scan_only);
    print_timing("solve_with_scan (reuse)", solve_with_scan_only);
    print_timing("solve_problem2_flyby (full)", full_solve);
    print_timing("scan + solve_with_scan", reused_pipeline);

    const double expected_full = scan_only.mean + solve_with_scan_only.mean;
    std::printf(
        "  expected_full_from_parts=%8.2f ms  full_solve_ratio=%.2f\n",
        expected_full,
        full_solve.mean / expected_full);

    if (!(solve_with_scan_only.mean < full_solve.mean)) {
        std::cerr << "timing regression: reused solve should be faster than full solve\n";
        return 1;
    }

    std::cout << "test_problem2_flyby_solve_timing PASSED\n";
    return 0;
}
