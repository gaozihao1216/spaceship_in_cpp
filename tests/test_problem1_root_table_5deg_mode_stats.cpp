#include "problem1_5deg_test_utils.hpp"

#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct SampleNode {
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
};

struct SampleStats {
    int sampled_node_count = 0;
    int sampled_valid_node_count = 0;
    long long branch_count_total = 0;
    long long derivative_attach_success_total = 0;
    double solve_ms_total = 0.0;
    double derivative_attach_ms_total = 0.0;
    double wall_ms = 0.0;
    std::vector<double> branch_counts;
    std::vector<double> build_ms;
};

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

std::vector<SampleNode> make_sample_nodes(int samples_per_dimension, double step_rad, int sample_count) {
    const long long total_nodes =
        1LL * samples_per_dimension * samples_per_dimension * samples_per_dimension;
    std::vector<SampleNode> nodes;
    nodes.reserve(static_cast<std::size_t>(sample_count));
    for (int i = 0; i < sample_count; ++i) {
        long long linear = std::min<long long>(
            total_nodes - 1, (static_cast<long long>(i) * total_nodes) / std::max(1, sample_count));
        SampleNode node{};
        node.theta_A_index = static_cast<int>(linear % samples_per_dimension);
        linear /= samples_per_dimension;
        node.nu_B_index = static_cast<int>(linear % samples_per_dimension);
        linear /= samples_per_dimension;
        node.nu_A_index = static_cast<int>(linear % samples_per_dimension);
        node.nu_A = spaceship_cpp::common::normalize_angle_0_2pi(step_rad * static_cast<double>(node.nu_A_index));
        node.nu_B = spaceship_cpp::common::normalize_angle_0_2pi(step_rad * static_cast<double>(node.nu_B_index));
        node.theta_A = spaceship_cpp::common::normalize_angle_0_2pi(step_rad * static_cast<double>(node.theta_A_index));
        nodes.push_back(node);
    }
    return nodes;
}

SampleStats run_sample_build(
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    const std::vector<SampleNode>& nodes
) {
    namespace problem1 = spaceship_cpp::problem1;
    SampleStats stats{};
    stats.sampled_node_count = static_cast<int>(nodes.size());
    stats.branch_counts.reserve(nodes.size());
    stats.build_ms.reserve(nodes.size());
    const auto wall_start = Clock::now();
    for (const auto& node : nodes) {
        const auto node_start = Clock::now();
        const auto solve_start = Clock::now();
        auto branches = problem1::solve_problem1_from_departure_anomalies(
            from, to, node.nu_A, node.nu_B, node.theta_A, 1, 1);
        stats.solve_ms_total += elapsed_ms(solve_start, Clock::now());
        stats.branch_count_total += static_cast<long long>(branches.size());
        if (!branches.empty()) {
            stats.sampled_valid_node_count += 1;
        }

        const auto attach_start = Clock::now();
        for (const auto& branch : branches) {
            if (!branch.valid) continue;
            const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
                from,
                to,
                node.nu_A,
                node.nu_B,
                node.theta_A,
                branch,
                problem1::Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
            if (attached.valid && attached.derivatives_available) {
                stats.derivative_attach_success_total += 1;
            }
        }
        stats.derivative_attach_ms_total += elapsed_ms(attach_start, Clock::now());
        stats.branch_counts.push_back(static_cast<double>(branches.size()));
        stats.build_ms.push_back(elapsed_ms(node_start, Clock::now()));
    }
    stats.wall_ms = elapsed_ms(wall_start, Clock::now());
    return stats;
}

void print_stats(
    const std::string& mode_name,
    double grid_step_degrees,
    const SampleStats& stats,
    int samples_per_dimension,
    long long total_nodes
) {
    const double mean_branches_per_node =
        static_cast<double>(stats.branch_count_total) / static_cast<double>(std::max(1, stats.sampled_node_count));
    const double mean_ms_per_node = stats.wall_ms / static_cast<double>(std::max(1, stats.sampled_node_count));
    const double estimated_wall_ms = mean_ms_per_node * static_cast<double>(total_nodes);
    const double branch_storage_bytes =
        static_cast<double>(total_nodes) * mean_branches_per_node *
        static_cast<double>(sizeof(spaceship_cpp::problem1::Problem1SolutionBranch));
    const double vector_object_bytes =
        static_cast<double>(total_nodes) * static_cast<double>(sizeof(std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>));

    std::cout << "Problem1RootTableGridModeStats\n";
    std::cout << "mode_name=" << mode_name << '\n';
    std::cout << "grid_step_degrees=" << grid_step_degrees << '\n';
    std::cout << "nu_A_sample_count=" << samples_per_dimension << '\n';
    std::cout << "nu_B_sample_count=" << samples_per_dimension << '\n';
    std::cout << "theta_A_sample_count=" << samples_per_dimension << '\n';
    std::cout << "total_sample_point_count=" << total_nodes << '\n';
    std::cout << "sizeof_problem1_solution_branch=" << sizeof(spaceship_cpp::problem1::Problem1SolutionBranch) << '\n';
    std::cout << "sizeof_branch_vector_object=" << sizeof(std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>) << '\n';
    std::cout << "sampled_node_count=" << stats.sampled_node_count << '\n';
    std::cout << "sampled_valid_node_count=" << stats.sampled_valid_node_count << '\n';
    std::cout << "mean_branch_count_per_point=" << mean_branches_per_node << '\n';
    std::cout << "p50_branch_count_per_point=" << problem1_5deg_test::percentile(stats.branch_counts, 0.50) << '\n';
    std::cout << "p90_branch_count_per_point=" << problem1_5deg_test::percentile(stats.branch_counts, 0.90) << '\n';
    std::cout << "p99_branch_count_per_point=" << problem1_5deg_test::percentile(stats.branch_counts, 0.99) << '\n';
    std::cout << "max_branch_count_per_point=" << problem1_5deg_test::percentile(stats.branch_counts, 1.0) << '\n';
    std::cout << "mean_build_ms_per_point=" << mean_ms_per_node << '\n';
    std::cout << "p50_build_ms_per_point=" << problem1_5deg_test::percentile(stats.build_ms, 0.50) << '\n';
    std::cout << "p90_build_ms_per_point=" << problem1_5deg_test::percentile(stats.build_ms, 0.90) << '\n';
    std::cout << "p99_build_ms_per_point=" << problem1_5deg_test::percentile(stats.build_ms, 0.99) << '\n';
    std::cout << "max_build_ms_per_point=" << problem1_5deg_test::percentile(stats.build_ms, 1.0) << '\n';
    std::cout << "sample_wall_ms=" << stats.wall_ms << '\n';
    std::cout << "sample_solve_ms=" << stats.solve_ms_total << '\n';
    std::cout << "sample_derivative_attach_ms=" << stats.derivative_attach_ms_total << '\n';
    std::cout << "estimated_total_build_ms=" << estimated_wall_ms << '\n';
    std::cout << "estimated_total_build_seconds=" << estimated_wall_ms / 1000.0 << '\n';
    std::cout << "estimated_branch_storage_bytes=" << branch_storage_bytes << '\n';
    std::cout << "estimated_vector_object_bytes=" << vector_object_bytes << '\n';
    std::cout << "estimated_total_memory_bytes=" << branch_storage_bytes + vector_object_bytes << '\n';
    std::cout << "estimated_total_memory_mib=" << (branch_storage_bytes + vector_object_bytes) / (1024.0 * 1024.0) << '\n';
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;

    std::cout << std::setprecision(12) << std::scientific;
    constexpr double grid_step_degrees = 5.0;
    constexpr int samples_per_dimension = 72;
    constexpr long long total_nodes =
        1LL * samples_per_dimension * samples_per_dimension * samples_per_dimension;
    const double step_rad = grid_step_degrees * common::kPi / 180.0;
    const int sample_count = problem1_5deg_test::env_int(
        "ROOT_TABLE_5DEG_TEST_SAMPLE_COUNT",
        problem1_5deg_test::env_int("ROOT_TABLE_5DEG_STATS_SAMPLE_COUNT", 100));
    const auto nodes = make_sample_nodes(samples_per_dimension, step_rad, sample_count);

    const auto pairs = problem1_5deg_test::selected_pairs_from_env();
    for (const auto& pair : pairs) {
        const auto stats = run_sample_build(pair.from, pair.to, nodes);
        std::cout << "Problem1RootTable5DegBuildCostPair\n";
        std::cout << "from_planet=" << planet_params::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(pair.to) << '\n';
        print_stats("root_table_5deg_estimate", grid_step_degrees, stats, samples_per_dimension, total_nodes);
    }

    std::cout << "Problem1RootTable5DegModeConfig\n";
    std::cout << "nu_A_start=0\n";
    std::cout << "nu_A_step=" << step_rad << '\n';
    std::cout << "nu_A_count=" << samples_per_dimension << '\n';
    std::cout << "nu_B_depart_start=0\n";
    std::cout << "nu_B_depart_step=" << step_rad << '\n';
    std::cout << "nu_B_depart_count=" << samples_per_dimension << '\n';
    std::cout << "theta_A_start=0\n";
    std::cout << "theta_A_step=" << step_rad << '\n';
    std::cout << "theta_A_count=" << samples_per_dimension << '\n';
    std::cout << "problem1_root_table_5deg_mode_stats_ok=1\n";
    return 0;
}
