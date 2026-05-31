#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = kPi / 90.0;
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

struct TestSummary {
    std::size_t chunk_count = 0;
    long long total_nodes = 0;
    int fixed_node_test_count = 0;
    int direct_compare_sample_count = 0;
    int direct_compare_success_count = 0;
    long long loader_derivative_missing_count = 0;
    double max_angle_diff = 0.0;
    double max_time_diff_seconds = 0.0;
    double max_residual_diff_seconds = 0.0;
    double mean_loader_query_ms = 0.0;
    double max_loader_query_ms = 0.0;
    bool loader_ok = true;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

SampleNode sample_node_from_linear_index(long long linear_index) {
    SampleNode node{};
    spaceship_cpp::problem1::Problem1RootTable2DegLoader::indices_from_linear_index(
        linear_index, &node.nu_A_index, &node.nu_B_index, &node.theta_A_index);
    node.nu_A = normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(node.nu_A_index));
    node.nu_B = normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(node.nu_B_index));
    node.theta_A = normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(node.theta_A_index));
    return node;
}

double angle_distance(double lhs, double rhs) {
    double delta = std::abs(normalize_angle_0_2pi(lhs) - normalize_angle_0_2pi(rhs));
    return std::min(delta, 2.0 * kPi - delta);
}

bool branch_less(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.time_of_flight_seconds) <
           std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.time_of_flight_seconds);
}

std::vector<long long> direct_compare_indices() {
    std::vector<long long> indices;
    indices.reserve(20);
    std::uint64_t state = 0x5eed1234ULL;
    for (int i = 0; i < 20; ++i) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        indices.push_back(static_cast<long long>(state % static_cast<std::uint64_t>(kTotalNodes)));
    }
    return indices;
}

bool all_valid_branches_have_derivatives(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoadedNode& node,
    long long* missing_count
) {
    bool ok = true;
    for (const auto& branch : node.branches) {
        if (!branch.valid) {
            continue;
        }
        if (!branch.derivatives_available ||
            !std::isfinite(branch.d_encounter_global_angle_d_nu_A) ||
            !std::isfinite(branch.d_encounter_global_angle_d_nu_B) ||
            !std::isfinite(branch.d_encounter_global_angle_d_theta_A)) {
            ok = false;
            *missing_count += 1;
        }
    }
    return ok;
}

bool compare_loaded_with_direct_solve(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoadedNode& loaded,
    long long linear_index,
    TestSummary* summary
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    if (!loaded.valid || loaded.branches.empty()) {
        std::cout << "direct_compare_error=loader_failed linear_index=" << linear_index
                  << " reason=" << loaded.invalid_reason << '\n';
        return false;
    }
    all_valid_branches_have_derivatives(loaded, &summary->loader_derivative_missing_count);

    const SampleNode sample = sample_node_from_linear_index(linear_index);
    auto direct = problem1::solve_problem1_from_departure_anomalies(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        sample.nu_A,
        sample.nu_B,
        sample.theta_A,
        1,
        1);
    for (auto& branch : direct) {
        if (!branch.valid) {
            continue;
        }
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            sample.nu_A,
            sample.nu_B,
            sample.theta_A,
            branch,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (attached.valid && attached.derivatives_available) {
            branch = attached;
        }
    }

    auto loaded_branches = loaded.branches;
    std::sort(loaded_branches.begin(), loaded_branches.end(), branch_less);
    std::sort(direct.begin(), direct.end(), branch_less);
    if (loaded_branches.size() != direct.size()) {
        std::cout << "direct_compare_error=branch_count_mismatch linear_index=" << linear_index
                  << " loader=" << loaded_branches.size() << " direct=" << direct.size() << '\n';
        return false;
    }

    bool ok = true;
    for (std::size_t i = 0; i < loaded_branches.size(); ++i) {
        const auto& lhs = loaded_branches[i];
        const auto& rhs = direct[i];
        if (lhs.transfer_revolution != rhs.transfer_revolution ||
            lhs.target_revolution != rhs.target_revolution) {
            std::cout << "direct_compare_error=kq_mismatch linear_index=" << linear_index << '\n';
            ok = false;
            continue;
        }
        const double alpha_diff = angle_distance(lhs.encounter_global_angle, rhs.encounter_global_angle);
        const double target_diff = angle_distance(lhs.target_arrival_true_anomaly, rhs.target_arrival_true_anomaly);
        const double tof_diff = std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds);
        const double residual_diff = std::abs(lhs.residual_seconds - rhs.residual_seconds);
        summary->max_angle_diff = std::max(summary->max_angle_diff, std::max(alpha_diff, target_diff));
        summary->max_time_diff_seconds = std::max(summary->max_time_diff_seconds, tof_diff);
        summary->max_residual_diff_seconds = std::max(summary->max_residual_diff_seconds, residual_diff);
        if (alpha_diff >= 1e-9 || target_diff >= 1e-9 || tof_diff >= 1e-6 || residual_diff >= 1e-6) {
            std::cout << "direct_compare_error=value_mismatch linear_index=" << linear_index
                      << " alpha_diff=" << alpha_diff
                      << " target_diff=" << target_diff
                      << " tof_diff=" << tof_diff
                      << " residual_diff=" << residual_diff << '\n';
            ok = false;
        }
        if (lhs.valid && !lhs.derivatives_available) {
            ok = false;
        }
    }
    return ok;
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    const std::filesystem::path table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "root_table_2deg_loader_skipped_missing_table\n";
        return 0;
    }

    TestSummary summary{};
    problem1::Problem1RootTable2DegLoader loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    summary.chunk_count = loader.chunk_count();
    summary.total_nodes = loader.total_nodes();
    if (summary.chunk_count != 59 || summary.total_nodes != kTotalNodes) {
        summary.loader_ok = false;
    }

    const std::vector<long long> fixed_indices = {0, 500000, 800000, 2900000, 5831999};
    double loader_query_ms_sum = 0.0;
    int loader_query_count = 0;
    for (const long long linear_index : fixed_indices) {
        const auto begin = clock::now();
        const auto node = loader.load_node_by_linear_index(linear_index);
        const auto end = clock::now();
        const double elapsed_ms = std::chrono::duration<double, std::milli>(end - begin).count();
        loader_query_ms_sum += elapsed_ms;
        summary.max_loader_query_ms = std::max(summary.max_loader_query_ms, elapsed_ms);
        loader_query_count += 1;

        summary.fixed_node_test_count += 1;
        if (!node.valid || node.branches.empty()) {
            std::cout << "fixed_node_error=loader_failed linear_index=" << linear_index
                      << " reason=" << node.invalid_reason << '\n';
            summary.loader_ok = false;
            continue;
        }
        if (!all_valid_branches_have_derivatives(node, &summary.loader_derivative_missing_count)) {
            summary.loader_ok = false;
        }
    }

    const std::vector<long long> compare_indices = direct_compare_indices();
    summary.direct_compare_sample_count = static_cast<int>(compare_indices.size());
    for (const long long linear_index : compare_indices) {
        const auto begin = clock::now();
        const auto node = loader.load_node_by_linear_index(linear_index);
        const auto end = clock::now();
        const double elapsed_ms = std::chrono::duration<double, std::milli>(end - begin).count();
        loader_query_ms_sum += elapsed_ms;
        summary.max_loader_query_ms = std::max(summary.max_loader_query_ms, elapsed_ms);
        loader_query_count += 1;

        const bool ok = compare_loaded_with_direct_solve(node, linear_index, &summary);
        if (ok) {
            summary.direct_compare_success_count += 1;
        } else {
            summary.loader_ok = false;
        }
    }
    summary.mean_loader_query_ms =
        loader_query_count > 0 ? loader_query_ms_sum / static_cast<double>(loader_query_count) : 0.0;
    if (summary.loader_derivative_missing_count != 0) {
        summary.loader_ok = false;
    }

    std::cout << std::setprecision(12) << std::scientific;
    std::cout << "Problem1RootTable2DegLoaderSummary\n";
    std::cout << "table_dir=" << table_dir.string() << '\n';
    std::cout << "chunk_count=" << summary.chunk_count << '\n';
    std::cout << "total_nodes=" << summary.total_nodes << '\n';
    std::cout << "fixed_node_test_count=" << summary.fixed_node_test_count << '\n';
    std::cout << "direct_compare_sample_count=" << summary.direct_compare_sample_count << '\n';
    std::cout << "direct_compare_success_count=" << summary.direct_compare_success_count << '\n';
    std::cout << "loader_derivative_missing_count=" << summary.loader_derivative_missing_count << '\n';
    std::cout << "max_angle_diff=" << summary.max_angle_diff << '\n';
    std::cout << "max_time_diff_seconds=" << summary.max_time_diff_seconds << '\n';
    std::cout << "max_residual_diff_seconds=" << summary.max_residual_diff_seconds << '\n';
    std::cout << "mean_loader_query_ms=" << summary.mean_loader_query_ms << '\n';
    std::cout << "max_loader_query_ms=" << summary.max_loader_query_ms << '\n';
    std::cout << "loader_ok=" << (summary.loader_ok ? 1 : 0) << '\n';

    return summary.loader_ok ? 0 : 1;
}
