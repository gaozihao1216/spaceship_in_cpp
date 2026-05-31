#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

int sample_count_from_environment() {
    const char* perf_env = std::getenv("PROBLEM1_PERF_SAMPLE_COUNT");
    if (perf_env != nullptr) {
        const int parsed = std::atoi(perf_env);
        if (parsed > 0) {
            return parsed;
        }
    }
    const char* env = std::getenv("PROBLEM1_SOLVE_PERF_SAMPLE_COUNT");
    if (env != nullptr) {
        const int parsed = std::atoi(env);
        if (parsed > 0) {
            return parsed;
        }
    }
    return 20;
}

struct ModeStats {
    std::string name;
    spaceship_cpp::problem1::Problem1SolveMode mode;
    std::vector<double> solve_times;
    double solve_time_sum = 0.0;
    int branch_count_sum = 0;
    int branch_count_max = 0;
    long long residual_evaluations_sum = 0;
    long long adaptive_residual_evaluations_sum = 0;
    long long adaptive_interval_total = 0;
    long long adaptive_sign_change_interval_total = 0;
    long long adaptive_local_min_interval_total = 0;
    long long adaptive_local_max_interval_total = 0;
    long long adaptive_valid_boundary_interval_total = 0;
    long long adaptive_near_zero_interval_total = 0;
    long long adaptive_wrap_interval_total = 0;
    long long adaptive_rapid_change_interval_total = 0;
    long long adaptive_bisection_interval_total = 0;
    long long adaptive_ternary_interval_total = 0;
    long long adaptive_local_fine_scan_interval_total = 0;
    long long adaptive_candidate_count_before_dedup = 0;
    long long adaptive_candidate_count_after_dedup = 0;
    double adaptive_coarse_scan_seconds = 0.0;
    double adaptive_interval_collection_seconds = 0.0;
    double adaptive_interval_refine_seconds = 0.0;
    double adaptive_local_fine_scan_seconds = 0.0;
    double adaptive_ternary_seconds = 0.0;
    double adaptive_bisection_seconds = 0.0;
    double adaptive_candidate_dedup_seconds = 0.0;
    double adaptive_sorting_seconds = 0.0;
    int fallback_count = 0;
};

double percentile(std::vector<double> values, double q) {
    assert(!values.empty());
    std::sort(values.begin(), values.end());
    const double position = q * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    if (lower == upper) {
        return values[lower];
    }
    const double weight = position - static_cast<double>(lower);
    return values[lower] * (1.0 - weight) + values[upper] * weight;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    const int sample_count = sample_count_from_environment();
    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;

    std::vector<ModeStats> mode_stats{
        {"FullScan2880", problem1::Problem1SolveMode::FullScan2880},
        {"AdaptiveScanWithFallback", problem1::Problem1SolveMode::AdaptiveScanWithFallback},
    };
    for (auto& stats : mode_stats) {
        stats.solve_times.reserve(static_cast<std::size_t>(sample_count));
    }

    for (int sample_index = 0; sample_index < sample_count; ++sample_index) {
        const double nu_A = common::normalize_angle_0_2pi(0.37 + static_cast<double>(sample_index) * 2.3999632297);
        const double nu_B = common::normalize_angle_0_2pi(1.11 + static_cast<double>(sample_index) * 1.7548776662);
        const double theta_A = common::normalize_angle_0_2pi(0.23 + static_cast<double>(sample_index) * 0.9182736455);

        for (auto& stats : mode_stats) {
            problem1::Problem1SolveOptions options{};
            options.mode = stats.mode;
            options.max_transfer_revolution = 1;
            options.max_target_revolution = 1;

            const auto start = clock::now();
            const problem1::Problem1SolveWithModeResult solved =
                problem1::solve_problem1_from_departure_anomalies_with_mode(
                departure_planet,
                target_planet,
                nu_A,
                nu_B,
                theta_A,
                options);
            const auto end = clock::now();
            const double solve_time_seconds = std::chrono::duration<double>(end - start).count();
            stats.solve_times.push_back(solve_time_seconds);
            stats.solve_time_sum += solve_time_seconds;
            stats.branch_count_sum += static_cast<int>(solved.branches.size());
            stats.branch_count_max = std::max(stats.branch_count_max, static_cast<int>(solved.branches.size()));
            stats.residual_evaluations_sum += solved.mode_profile.residual_evaluation_count;
            stats.adaptive_residual_evaluations_sum += solved.mode_profile.adaptive_residual_eval_count;
            stats.adaptive_interval_total += solved.mode_profile.adaptive_interval_total;
            stats.adaptive_sign_change_interval_total += solved.mode_profile.adaptive_sign_change_interval_count;
            stats.adaptive_local_min_interval_total += solved.mode_profile.adaptive_local_min_interval_count;
            stats.adaptive_local_max_interval_total += solved.mode_profile.adaptive_local_max_interval_count;
            stats.adaptive_valid_boundary_interval_total += solved.mode_profile.adaptive_valid_boundary_interval_count;
            stats.adaptive_near_zero_interval_total += solved.mode_profile.adaptive_near_zero_interval_count;
            stats.adaptive_wrap_interval_total += solved.mode_profile.adaptive_wrap_interval_count;
            stats.adaptive_rapid_change_interval_total += solved.mode_profile.adaptive_rapid_change_interval_count;
            stats.adaptive_bisection_interval_total += solved.mode_profile.adaptive_bisection_interval_count;
            stats.adaptive_ternary_interval_total += solved.mode_profile.adaptive_ternary_interval_count;
            stats.adaptive_local_fine_scan_interval_total += solved.mode_profile.adaptive_local_fine_scan_interval_count;
            stats.adaptive_candidate_count_before_dedup += solved.mode_profile.adaptive_candidate_count_before_dedup;
            stats.adaptive_candidate_count_after_dedup += solved.mode_profile.adaptive_candidate_count_after_dedup;
            stats.adaptive_coarse_scan_seconds += solved.mode_profile.adaptive_coarse_scan_seconds;
            stats.adaptive_interval_collection_seconds += solved.mode_profile.adaptive_interval_collection_seconds;
            stats.adaptive_interval_refine_seconds += solved.mode_profile.adaptive_interval_refine_seconds;
            stats.adaptive_local_fine_scan_seconds += solved.mode_profile.adaptive_local_fine_scan_seconds;
            stats.adaptive_ternary_seconds += solved.mode_profile.adaptive_ternary_seconds;
            stats.adaptive_bisection_seconds += solved.mode_profile.adaptive_bisection_seconds;
            stats.adaptive_candidate_dedup_seconds += solved.mode_profile.adaptive_candidate_dedup_seconds;
            stats.adaptive_sorting_seconds += solved.mode_profile.adaptive_sorting_seconds;
            if (solved.mode_profile.fallback_used || solved.mode_profile.adaptive_fallback_to_fullscan) {
                stats.fallback_count += 1;
            }
        }
    }

    const unsigned int hardware_threads = std::thread::hardware_concurrency();

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1 solve mode performance\n";
    std::cout << "sample_count=" << sample_count << '\n';
    std::cout << "hardware_threads_reported=" << hardware_threads << '\n';
    for (const auto& stats : mode_stats) {
        const double avg_time = stats.solve_time_sum / static_cast<double>(sample_count);
        std::cout << "Problem1SolveModePerformanceSummary\n";
        std::cout << "mode=" << stats.name << '\n';
        std::cout << "solve_count=" << sample_count << '\n';
        std::cout << "total_time_ms=" << stats.solve_time_sum * 1000.0 << '\n';
        std::cout << "average_time_ms_per_solve=" << avg_time * 1000.0 << '\n';
        std::cout << "median_time_ms_per_solve=" << percentile(stats.solve_times, 0.50) * 1000.0 << '\n';
        std::cout << "min_time_ms=" << *std::min_element(stats.solve_times.begin(), stats.solve_times.end()) * 1000.0 << '\n';
        std::cout << "max_time_ms=" << *std::max_element(stats.solve_times.begin(), stats.solve_times.end()) * 1000.0 << '\n';
        std::cout << "total_residual_eval_count=" << stats.residual_evaluations_sum << '\n';
        std::cout << "average_residual_eval_count=" <<
            static_cast<double>(stats.residual_evaluations_sum) / static_cast<double>(sample_count) << '\n';
        std::cout << "total_adaptive_coarse_residual_eval_count=" << stats.adaptive_residual_evaluations_sum << '\n';
        std::cout << "average_adaptive_coarse_residual_eval_count=" <<
            static_cast<double>(stats.adaptive_residual_evaluations_sum) / static_cast<double>(sample_count) << '\n';
        const long long effective_residual_eval_count =
            stats.residual_evaluations_sum + stats.adaptive_residual_evaluations_sum;
        std::cout << "total_effective_residual_eval_count=" << effective_residual_eval_count << '\n';
        std::cout << "average_effective_residual_eval_count=" <<
            static_cast<double>(effective_residual_eval_count) / static_cast<double>(sample_count) << '\n';
        std::cout << "total_candidate_count=" << stats.branch_count_sum << '\n';
        std::cout << "average_candidate_count=" <<
            static_cast<double>(stats.branch_count_sum) / static_cast<double>(sample_count) << '\n';
        std::cout << "max_candidate_count=" << stats.branch_count_max << '\n';
        std::cout << "fallback_count=" << stats.fallback_count << '\n';
        std::cout << "adaptive_interval_total=" << stats.adaptive_interval_total << '\n';
        std::cout << "adaptive_sign_change_interval_total=" << stats.adaptive_sign_change_interval_total << '\n';
        std::cout << "adaptive_local_min_interval_total=" << stats.adaptive_local_min_interval_total << '\n';
        std::cout << "adaptive_local_max_interval_total=" << stats.adaptive_local_max_interval_total << '\n';
        std::cout << "adaptive_valid_boundary_interval_total=" << stats.adaptive_valid_boundary_interval_total << '\n';
        std::cout << "adaptive_near_zero_interval_total=" << stats.adaptive_near_zero_interval_total << '\n';
        std::cout << "adaptive_wrap_interval_total=" << stats.adaptive_wrap_interval_total << '\n';
        std::cout << "adaptive_rapid_change_interval_total=" << stats.adaptive_rapid_change_interval_total << '\n';
        std::cout << "adaptive_bisection_interval_total=" << stats.adaptive_bisection_interval_total << '\n';
        std::cout << "adaptive_ternary_interval_total=" << stats.adaptive_ternary_interval_total << '\n';
        std::cout << "adaptive_local_fine_scan_interval_total=" << stats.adaptive_local_fine_scan_interval_total << '\n';
        std::cout << "adaptive_candidate_count_before_dedup=" << stats.adaptive_candidate_count_before_dedup << '\n';
        std::cout << "adaptive_candidate_count_after_dedup=" << stats.adaptive_candidate_count_after_dedup << '\n';
        std::cout << "adaptive_coarse_scan_ms=" << stats.adaptive_coarse_scan_seconds * 1000.0 << '\n';
        std::cout << "adaptive_interval_collection_ms=" << stats.adaptive_interval_collection_seconds * 1000.0 << '\n';
        std::cout << "adaptive_interval_refine_ms=" << stats.adaptive_interval_refine_seconds * 1000.0 << '\n';
        std::cout << "adaptive_local_fine_scan_ms=" << stats.adaptive_local_fine_scan_seconds * 1000.0 << '\n';
        std::cout << "adaptive_ternary_ms=" << stats.adaptive_ternary_seconds * 1000.0 << '\n';
        std::cout << "adaptive_bisection_ms=" << stats.adaptive_bisection_seconds * 1000.0 << '\n';
        std::cout << "adaptive_candidate_dedup_ms=" << stats.adaptive_candidate_dedup_seconds * 1000.0 << '\n';
        std::cout << "adaptive_sorting_ms=" << stats.adaptive_sorting_seconds * 1000.0 << '\n';
    }
    const double full_avg = mode_stats[0].solve_time_sum / static_cast<double>(sample_count);
    const double adaptive_avg = mode_stats[1].solve_time_sum / static_cast<double>(sample_count);
    std::cout << "speedup_fullscan_over_adaptive=" << full_avg / adaptive_avg << '\n';

    assert(sample_count > 0);
    assert(mode_stats[0].branch_count_sum > 0);
    assert(mode_stats[0].residual_evaluations_sum > 0);
    assert(mode_stats[1].residual_evaluations_sum > 0);

    return 0;
}
