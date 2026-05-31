#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr double kGridStepDegrees = 2.0;
constexpr double kGridStepRadians = kPi / 90.0;
constexpr int kDimensionCount = 3;
constexpr int kSamplesPerDimension = 180;
constexpr long long kTotalNodes =
    1LL * kSamplesPerDimension * kSamplesPerDimension * kSamplesPerDimension;

struct SampleNode {
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
};

struct NodeTimingRecord {
    double exact_solve_ms = 0.0;
    double derivative_attach_ms = 0.0;
    double total_ms = 0.0;
    int branch_count = 0;
    int valid_branch_count = 0;
    int derivative_success_count = 0;
    int derivative_fail_count = 0;
    bool invalid_node = false;
};

struct SummaryStats {
    int sample_count = 0;
    int invalid_node_count = 0;
    long long branch_count_sum = 0;
    long long valid_branch_count_sum = 0;
    long long derivative_success_count = 0;
    long long derivative_fail_count = 0;
    std::vector<double> exact_ms;
    std::vector<double> derivative_ms;
    std::vector<double> total_ms;
};

struct ParallelRunStats {
    int threads = 1;
    int sample_count = 0;
    double wall_time_seconds = 0.0;
};

int get_env_int(const char* name, int default_value) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') {
        return default_value;
    }
    const int value = std::atoi(raw);
    return value > 0 ? value : default_value;
}

double percentile_ms(std::vector<double> values, double p) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const double index = p * static_cast<double>(values.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(index));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(index));
    if (lo == hi) {
        return values[lo];
    }
    const double t = index - static_cast<double>(lo);
    return values[lo] * (1.0 - t) + values[hi] * t;
}

double mean_ms(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / static_cast<double>(values.size());
}

SampleNode sample_node_from_linear_index(long long linear_index) {
    SampleNode node{};
    node.theta_A_index = static_cast<int>(linear_index % kSamplesPerDimension);
    linear_index /= kSamplesPerDimension;
    node.nu_B_index = static_cast<int>(linear_index % kSamplesPerDimension);
    linear_index /= kSamplesPerDimension;
    node.nu_A_index = static_cast<int>(linear_index % kSamplesPerDimension);
    node.nu_A = normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(node.nu_A_index));
    node.nu_B = normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(node.nu_B_index));
    node.theta_A = normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(node.theta_A_index));
    return node;
}

std::vector<SampleNode> build_sample_nodes(int sample_count) {
    std::vector<SampleNode> nodes;
    nodes.reserve(static_cast<std::size_t>(sample_count));
    for (int i = 0; i < sample_count; ++i) {
        const long long linear_index =
            std::min<long long>(kTotalNodes - 1, (static_cast<long long>(i) * kTotalNodes) / sample_count);
        nodes.push_back(sample_node_from_linear_index(linear_index));
    }
    return nodes;
}

NodeTimingRecord evaluate_node_cost(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const SampleNode& node,
    int max_transfer_revolution,
    int max_target_revolution
) {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    NodeTimingRecord record{};
    const auto solve_begin = clock::now();
    const std::vector<problem1::Problem1SolutionBranch> branches =
        problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            node.nu_A,
            node.nu_B,
            node.theta_A,
            max_transfer_revolution,
            max_target_revolution);
    const auto solve_end = clock::now();
    record.exact_solve_ms = std::chrono::duration<double, std::milli>(solve_end - solve_begin).count();
    record.branch_count = static_cast<int>(branches.size());

    const auto derivative_begin = clock::now();
    for (const auto& branch : branches) {
        if (!branch.valid) {
            continue;
        }
        record.valid_branch_count += 1;
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            node.nu_A,
            node.nu_B,
            node.theta_A,
            branch,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (attached.valid && attached.derivatives_available) {
            record.derivative_success_count += 1;
        } else {
            record.derivative_fail_count += 1;
        }
    }
    const auto derivative_end = clock::now();
    record.derivative_attach_ms =
        std::chrono::duration<double, std::milli>(derivative_end - derivative_begin).count();
    record.total_ms = record.exact_solve_ms + record.derivative_attach_ms;
    record.invalid_node = record.valid_branch_count == 0;
    return record;
}

SummaryStats run_serial_summary(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const std::vector<SampleNode>& nodes,
    int max_transfer_revolution,
    int max_target_revolution
) {
    SummaryStats stats{};
    stats.sample_count = static_cast<int>(nodes.size());
    stats.exact_ms.reserve(nodes.size());
    stats.derivative_ms.reserve(nodes.size());
    stats.total_ms.reserve(nodes.size());
    for (const auto& node : nodes) {
        const auto record = evaluate_node_cost(
            departure_planet, target_planet, node, max_transfer_revolution, max_target_revolution);
        stats.exact_ms.push_back(record.exact_solve_ms);
        stats.derivative_ms.push_back(record.derivative_attach_ms);
        stats.total_ms.push_back(record.total_ms);
        stats.branch_count_sum += record.branch_count;
        stats.valid_branch_count_sum += record.valid_branch_count;
        stats.derivative_success_count += record.derivative_success_count;
        stats.derivative_fail_count += record.derivative_fail_count;
        if (record.invalid_node) {
            stats.invalid_node_count += 1;
        }
    }
    return stats;
}

ParallelRunStats run_parallel_scaling(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const std::vector<SampleNode>& nodes,
    int max_transfer_revolution,
    int max_target_revolution,
    int thread_count
) {
    using clock = std::chrono::steady_clock;
    ParallelRunStats stats{};
    stats.threads = thread_count;
    stats.sample_count = static_cast<int>(nodes.size());

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));
    const auto begin = clock::now();
    for (int t = 0; t < thread_count; ++t) {
        workers.emplace_back([&, t]() {
            for (std::size_t i = static_cast<std::size_t>(t); i < nodes.size();
                 i += static_cast<std::size_t>(thread_count)) {
                (void)evaluate_node_cost(
                    departure_planet, target_planet, nodes[i], max_transfer_revolution, max_target_revolution);
            }
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }
    const auto end = clock::now();
    stats.wall_time_seconds = std::chrono::duration<double>(end - begin).count();
    return stats;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;

    std::cout << std::setprecision(6) << std::scientific;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const int sample_count = get_env_int("ROOT_TABLE_2DEG_COST_SAMPLE_COUNT", 1000);
    const int scaling_sample_count = get_env_int("ROOT_TABLE_2DEG_SCALING_SAMPLE_COUNT", 1000);
    const unsigned hardware_threads = std::max(1u, std::thread::hardware_concurrency());

    const auto sample_nodes = build_sample_nodes(sample_count);
    const auto scaling_nodes = build_sample_nodes(scaling_sample_count);
    const auto summary = run_serial_summary(
        departure_planet, target_planet, sample_nodes, max_transfer_revolution, max_target_revolution);

    std::cout << "RootTable2DegGridSize\n";
    std::cout << "dimension_count=" << kDimensionCount << '\n';
    std::cout << "samples_per_dimension=" << kSamplesPerDimension << '\n';
    std::cout << "total_nodes=" << kTotalNodes << '\n';
    std::cout << "grid_step_rad=" << kGridStepRadians << '\n';
    std::cout << "grid_step_deg=" << kGridStepDegrees << '\n';

    const double mean_exact_solve_ms = mean_ms(summary.exact_ms);
    const double mean_derivative_attach_ms = mean_ms(summary.derivative_ms);
    const double mean_total_ms = mean_ms(summary.total_ms);
    const double derivative_success_ratio =
        summary.valid_branch_count_sum > 0
            ? static_cast<double>(summary.derivative_success_count) /
                  static_cast<double>(summary.valid_branch_count_sum)
            : 0.0;

    std::cout << "RootTable2DegCostEstimateSummary\n";
    std::cout << "sample_count=" << sample_count << '\n';
    std::cout << "total_nodes=" << kTotalNodes << '\n';
    std::cout << "mean_exact_solve_ms=" << mean_exact_solve_ms << '\n';
    std::cout << "median_exact_solve_ms=" << percentile_ms(summary.exact_ms, 0.5) << '\n';
    std::cout << "p90_exact_solve_ms=" << percentile_ms(summary.exact_ms, 0.9) << '\n';
    std::cout << "p99_exact_solve_ms=" << percentile_ms(summary.exact_ms, 0.99) << '\n';
    std::cout << "mean_derivative_attach_ms=" << mean_derivative_attach_ms << '\n';
    std::cout << "median_derivative_attach_ms=" << percentile_ms(summary.derivative_ms, 0.5) << '\n';
    std::cout << "p90_derivative_attach_ms=" << percentile_ms(summary.derivative_ms, 0.9) << '\n';
    std::cout << "p99_derivative_attach_ms=" << percentile_ms(summary.derivative_ms, 0.99) << '\n';
    std::cout << "mean_total_ms=" << mean_total_ms << '\n';
    std::cout << "median_total_ms=" << percentile_ms(summary.total_ms, 0.5) << '\n';
    std::cout << "p90_total_ms=" << percentile_ms(summary.total_ms, 0.9) << '\n';
    std::cout << "p99_total_ms=" << percentile_ms(summary.total_ms, 0.99) << '\n';
    std::cout << "mean_branch_count="
              << (static_cast<double>(summary.branch_count_sum) / static_cast<double>(summary.sample_count)) << '\n';
    std::cout << "mean_valid_branch_count="
              << (static_cast<double>(summary.valid_branch_count_sum) / static_cast<double>(summary.sample_count))
              << '\n';
    std::cout << "derivative_success_ratio=" << derivative_success_ratio << '\n';
    std::cout << "invalid_node_count=" << summary.invalid_node_count << '\n';

    const double single_thread_seconds = mean_total_ms * 1e-3 * static_cast<double>(kTotalNodes);
    const double single_thread_hours = single_thread_seconds / 3600.0;
    auto estimate_hours = [&](double threads) {
        return single_thread_hours / threads;
    };
    std::cout << "Estimated2DegBuildTime\n";
    std::cout << "single_thread_seconds=" << single_thread_seconds << '\n';
    std::cout << "single_thread_hours=" << single_thread_hours << '\n';
    std::cout << "estimated_8_threads_hours=" << estimate_hours(8.0) << '\n';
    std::cout << "estimated_16_threads_hours=" << estimate_hours(16.0) << '\n';
    std::cout << "estimated_32_threads_hours=" << estimate_hours(32.0) << '\n';
    std::cout << "hardware_concurrency=" << hardware_threads << '\n';
    std::cout << "estimated_hardware_threads_hours=" << estimate_hours(static_cast<double>(hardware_threads)) << '\n';

    std::vector<int> thread_counts{1, 2, 4, 8, static_cast<int>(std::min(16u, hardware_threads))};
    std::sort(thread_counts.begin(), thread_counts.end());
    thread_counts.erase(std::unique(thread_counts.begin(), thread_counts.end()), thread_counts.end());
    double baseline_wall_time = 0.0;
    for (int threads : thread_counts) {
        const auto scaling = run_parallel_scaling(
            departure_planet, target_planet, scaling_nodes, max_transfer_revolution, max_target_revolution, threads);
        if (threads == 1) {
            baseline_wall_time = scaling.wall_time_seconds;
        }
        const double nodes_per_second =
            scaling.wall_time_seconds > 0.0
                ? static_cast<double>(scaling.sample_count) / scaling.wall_time_seconds
                : 0.0;
        const double speedup = (baseline_wall_time > 0.0) ? baseline_wall_time / scaling.wall_time_seconds : 0.0;
        const double efficiency = threads > 0 ? speedup / static_cast<double>(threads) : 0.0;
        std::cout << "RootTable2DegParallelScalingSummary\n";
        std::cout << "threads=" << threads << '\n';
        std::cout << "sample_count=" << scaling.sample_count << '\n';
        std::cout << "wall_time_seconds=" << scaling.wall_time_seconds << '\n';
        std::cout << "nodes_per_second=" << nodes_per_second << '\n';
        std::cout << "speedup_vs_1_thread=" << speedup << '\n';
        std::cout << "parallel_efficiency=" << efficiency << '\n';
    }

    std::cout << "RootTable2DegOptimizationNotes\n";
    std::cout << "repeated_trig_candidates=solve_problem1_from_departure_anomalies,attach_problem1_root_derivatives_with_mode,"
                 "residual_and_derivative_inner_loops\n";
    std::cout << "precomputable_grid_values=sin_cos_nu_A,sin_cos_nu_B,sin_cos_theta_A,planet_radius_factors,"
                 "planet_velocity_prefactors\n";
    std::cout << "inner_loop_trig_values_not_precomputable=alpha_dependent_transfer_geometry_and_newton_iteration_terms\n";
    std::cout << "possible_sincos_pairing=yes_for_angle_difference_terms\n";
    std::cout << "normalize_angle_hotspot=likely_in_residual_root_scans_and_derivative_attach_but_not_directly_profiled_here\n";
    std::cout << "recommended_next_optimization=precompute_grid_trig_and_planet_terms_before_any_full_2deg_build,"
                 "then_parallelize_chunked_node_build_with_checkpointing\n";

    std::cout << "problem1_root_table_2deg_build_cost_estimate_ok\n";
    return 0;
}
