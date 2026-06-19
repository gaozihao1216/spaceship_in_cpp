/*
 * 文件作用：对比不�?phi 粗扫步长�?Problem 1 求解的正确率与耗时�? * 主要工作：以 phi_scan_count=2880 为参考解，评�?240 及更粗步长是否漏�?多根，并统计求解时间�? */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
namespace common = spaceship_cpp::common;
namespace config = spaceship_cpp::config;
namespace planet_params = spaceship_cpp::planet_params;
namespace problem1 = spaceship_cpp::problem1;

constexpr int kReferencePhiScanCount = 2880;
constexpr int kWarmupIterations = 2;
constexpr int kTimedIterations = 12;
constexpr double kLaunchTimeSeconds = 0.0;
constexpr double kTransferTheta = 0.5;

struct Scenario {
    planet_params::PlanetId departure{};
    planet_params::PlanetId target{};
    int max_k = 0;
    int max_q = 0;
};

struct SolveMetrics {
    int phi_scan_count = 0;
    int candidate_count = 0;
    double mean_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
};

struct CorrectnessMetrics {
    int matched = 0;
    int missed = 0;
    int extra = 0;
    double max_matched_angle_error_rad = 0.0;
    double max_missed_angle_error_rad = 0.0;
};

const char* planet_label(planet_params::PlanetId id) {
    return planet_params::planet_name(id);
}

problem1::Problem1SolveInput make_input(
    const Scenario& scenario,
    int phi_scan_count,
    const config::Problem1SolveDefaults& defaults
) {
    auto input = config::make_problem1_solve_input(
        scenario.departure,
        scenario.target,
        kLaunchTimeSeconds,
        kTransferTheta,
        defaults);
    input.max_transfer_revolution = scenario.max_k;
    input.max_target_revolution = scenario.max_q;
    input.phi_scan_count = phi_scan_count;
    return input;
}

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

SolveMetrics benchmark_solve(const problem1::Problem1SolveInput& input) {
    for (int i = 0; i < kWarmupIterations; ++i) {
        (void)problem1::solve_problem1(input);
    }

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(kTimedIterations));
    for (int i = 0; i < kTimedIterations; ++i) {
        const auto start = Clock::now();
        (void)problem1::solve_problem1(input);
        samples.push_back(elapsed_ms(start, Clock::now()));
    }

    SolveMetrics metrics{};
    metrics.phi_scan_count = input.phi_scan_count;
    metrics.candidate_count = static_cast<int>(problem1::solve_problem1(input).size());
    metrics.mean_ms =
        std::accumulate(samples.begin(), samples.end(), 0.0) /
        static_cast<double>(samples.size());
    metrics.min_ms = *std::min_element(samples.begin(), samples.end());
    metrics.max_ms = *std::max_element(samples.begin(), samples.end());
    return metrics;
}

bool same_branch(const problem1::Problem1Candidate& a, const problem1::Problem1Candidate& b) {
    return a.transfer_revolution == b.transfer_revolution &&
           a.target_revolution == b.target_revolution;
}

double angle_match_tolerance_rad(int reference_scan_count, int test_scan_count) {
    const int finer_scan = std::max(reference_scan_count, test_scan_count);
    return std::max(1e-6, common::kTwoPi / static_cast<double>(finer_scan));
}

const problem1::Problem1Candidate* find_best_match(
    const problem1::Problem1Candidate& reference,
    const std::vector<problem1::Problem1Candidate>& test_candidates,
    double tolerance_rad,
    std::vector<bool>& used
) {
    const problem1::Problem1Candidate* best = nullptr;
    double best_angle_error = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < test_candidates.size(); ++index) {
        if (used[index] || !same_branch(reference, test_candidates[index])) {
            continue;
        }
        const double angle_error = std::abs(common::normalize_angle_minus_pi_pi(
            reference.encounter_global_angle - test_candidates[index].encounter_global_angle));
        if (angle_error <= tolerance_rad && angle_error < best_angle_error) {
            best = &test_candidates[index];
            best_angle_error = angle_error;
        }
    }
    return best;
}

CorrectnessMetrics compare_against_reference(
    const std::vector<problem1::Problem1Candidate>& reference_candidates,
    const std::vector<problem1::Problem1Candidate>& test_candidates,
    int reference_scan_count,
    int test_scan_count
) {
    CorrectnessMetrics metrics{};
    const double tolerance_rad = angle_match_tolerance_rad(reference_scan_count, test_scan_count);
    std::vector<bool> used(test_candidates.size(), false);

    for (const problem1::Problem1Candidate& reference : reference_candidates) {
        const problem1::Problem1Candidate* match = find_best_match(
            reference,
            test_candidates,
            tolerance_rad,
            used);
        if (match != nullptr) {
            ++metrics.matched;
            const std::size_t index =
                static_cast<std::size_t>(match - test_candidates.data());
            used[index] = true;
            metrics.max_matched_angle_error_rad = std::max(
                metrics.max_matched_angle_error_rad,
                std::abs(common::normalize_angle_minus_pi_pi(
                    reference.encounter_global_angle - match->encounter_global_angle)));
        } else {
            ++metrics.missed;
            metrics.max_missed_angle_error_rad = std::max(
                metrics.max_missed_angle_error_rad,
                tolerance_rad);
        }
    }

    for (bool is_used : used) {
        if (!is_used) {
            ++metrics.extra;
        }
    }
    return metrics;
}

void print_correctness(const CorrectnessMetrics& metrics, int reference_count, int test_count) {
    std::cout
        << "ref=" << reference_count
        << ", test=" << test_count
        << ", matched=" << metrics.matched
        << ", missed=" << metrics.missed
        << ", extra=" << metrics.extra
        << ", max_match_err=" << metrics.max_matched_angle_error_rad << " rad";
    if (metrics.missed == 0 && metrics.extra == 0) {
        std::cout << " [OK]";
    }
    std::cout << '\n';
}

void run_scenario(
    const Scenario& scenario,
    const config::Problem1SolveDefaults& defaults,
    const std::vector<int>& scan_counts
) {
    const auto reference_input = make_input(scenario, kReferencePhiScanCount, defaults);
    const auto reference_candidates = problem1::solve_problem1(reference_input);
    const auto reference_timing = benchmark_solve(reference_input);

    std::cout << "\n=== "
              << planet_label(scenario.departure) << " -> "
              << planet_label(scenario.target)
              << " (k<=" << scenario.max_k << ", q<=" << scenario.max_q << ") ===\n";
    std::cout << "reference phi_scan_count=" << kReferencePhiScanCount
              << ", candidates=" << reference_candidates.size()
              << ", mean=" << reference_timing.mean_ms << " ms\n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "scan  step(deg)  mean_ms  speedup  candidates  matched  missed  extra\n";

    for (int scan_count : scan_counts) {
        const auto input = make_input(scenario, scan_count, defaults);
        const auto candidates = problem1::solve_problem1(input);
        const auto timing = benchmark_solve(input);
        const auto correctness =
            compare_against_reference(reference_candidates, candidates, kReferencePhiScanCount, scan_count);
        const double step_deg = 360.0 / static_cast<double>(scan_count);
        const double speedup = reference_timing.mean_ms / std::max(timing.mean_ms, 1e-12);

        std::cout
            << std::setw(4) << scan_count
            << std::setw(9) << step_deg
            << std::setw(9) << timing.mean_ms
            << std::setw(8) << speedup << "x"
            << std::setw(12) << candidates.size()
            << std::setw(9) << correctness.matched
            << std::setw(8) << correctness.missed
            << std::setw(8) << correctness.extra;

        if (correctness.missed == 0 && correctness.extra == 0) {
            std::cout << "  OK";
        } else if (correctness.missed > 0) {
            std::cout << "  MISS";
        } else {
            std::cout << "  EXTRA";
        }
        std::cout << '\n';
    }
}

}  // namespace

int main() {
    const config::GlobalConfig& cfg = config::global_config();
    const auto& defaults = cfg.problem1_solve;

    std::cout << "Problem 1 scan compare (reference phi_scan_count=" << kReferencePhiScanCount
              << ", current default=" << defaults.phi_scan_count << ")\n";
    std::cout << "Correctness: each scan count vs reference root set (same k,q, angle within one fine step).\n";
    std::cout << "Timed iterations per point: " << kTimedIterations << '\n';

    const std::vector<int> scan_counts{
        96,
        120,
        180,
        240,
        360,
        720,
        2880,
    };

    const std::vector<Scenario> scenarios{
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars, 0, 0},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars, 1, 1},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Venus, 0, 0},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Venus, 1, 1},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury, 0, 0},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury, 1, 1},
    };

    for (const Scenario& scenario : scenarios) {
        run_scenario(scenario, defaults, scan_counts);
    }

    return 0;
}

