#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = kPi / 90.0;
constexpr long long kTotalNodes =
    1LL * kSamplesPerDimension * kSamplesPerDimension * kSamplesPerDimension;

struct NodeAngles {
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
};

struct Summary {
    int exact_grid_sample_count = 0;
    int near_node_sample_count = 0;
    int query_success_count = 0;
    int fallback_count = 0;
    double fallback_ratio = 0.0;
    double mean_branch_count = 0.0;
    double coverage_vs_direct = 0.0;
    double same_root_success_ratio = 0.0;
    double max_abs_residual_seconds = 0.0;
    double max_time_diff_seconds = 0.0;
    long long derivative_missing_count = 0;
    bool nearest_node_query_ok = true;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

NodeAngles node_angles(long long linear_index) {
    NodeAngles node{};
    spaceship_cpp::problem1::Problem1RootTable2DegLoader::indices_from_linear_index(
        linear_index, &node.nu_A_index, &node.nu_B_index, &node.theta_A_index);
    node.nu_A = normalize_angle_0_2pi(kGridStepRadians * node.nu_A_index);
    node.nu_B = normalize_angle_0_2pi(kGridStepRadians * node.nu_B_index);
    node.theta_A = normalize_angle_0_2pi(kGridStepRadians * node.theta_A_index);
    return node;
}

bool branch_less(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.time_of_flight_seconds) <
           std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.time_of_flight_seconds);
}

std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> direct_solve(
    double nu_A,
    double nu_B,
    double theta_A
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    auto branches = problem1::solve_problem1_from_departure_anomalies(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        nu_A,
        nu_B,
        theta_A,
        1,
        1);
    std::sort(branches.begin(), branches.end(), branch_less);
    return branches;
}

void accumulate_query_result(
    const spaceship_cpp::problem1::Problem1NearestNodeQueryResult& query,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& direct,
    Summary* summary
) {
    if (query.valid && !query.branches.empty()) {
        summary->query_success_count += 1;
    } else {
        summary->nearest_node_query_ok = false;
    }
    if (query.used_direct_solve_fallback) {
        summary->fallback_count += 1;
    }

    summary->mean_branch_count += static_cast<double>(query.branches.size());
    int valid_query_count = 0;
    for (const auto& branch : query.branches) {
        if (!branch.valid) {
            continue;
        }
        valid_query_count += 1;
        summary->max_abs_residual_seconds =
            std::max(summary->max_abs_residual_seconds, std::abs(branch.residual_seconds));
        if (!branch.derivatives_available ||
            !std::isfinite(branch.d_encounter_global_angle_d_nu_A) ||
            !std::isfinite(branch.d_encounter_global_angle_d_nu_B) ||
            !std::isfinite(branch.d_encounter_global_angle_d_theta_A)) {
            summary->derivative_missing_count += 1;
            summary->nearest_node_query_ok = false;
        }
    }
    if (summary->max_abs_residual_seconds > 1e-2) {
        summary->nearest_node_query_ok = false;
    }

    auto query_branches = query.branches;
    std::sort(query_branches.begin(), query_branches.end(), branch_less);
    const int pair_count = std::min<int>(query_branches.size(), direct.size());
    int same_root_count = 0;
    for (int i = 0; i < pair_count; ++i) {
        const auto& lhs = query_branches[static_cast<std::size_t>(i)];
        const auto& rhs = direct[static_cast<std::size_t>(i)];
        const double time_diff = std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds);
        summary->max_time_diff_seconds = std::max(summary->max_time_diff_seconds, time_diff);
        if (lhs.transfer_revolution == rhs.transfer_revolution &&
            lhs.target_revolution == rhs.target_revolution &&
            time_diff < 1e-3 &&
            std::abs(normalize_angle_minus_pi_pi(
                lhs.encounter_global_angle - rhs.encounter_global_angle)) < 1e-9) {
            same_root_count += 1;
        }
    }
    if (!direct.empty()) {
        summary->coverage_vs_direct +=
            static_cast<double>(valid_query_count) / static_cast<double>(direct.size());
        summary->same_root_success_ratio +=
            static_cast<double>(same_root_count) / static_cast<double>(direct.size());
    }
}

std::vector<long long> near_node_indices() {
    std::vector<long long> indices;
    std::uint64_t state = 0x2d2d2d2dULL;
    for (int i = 0; i < 50; ++i) {
        state = state * 2862933555777941757ULL + 3037000493ULL;
        indices.push_back(static_cast<long long>(state % static_cast<std::uint64_t>(kTotalNodes)));
    }
    return indices;
}

double unit_random(std::uint64_t* state) {
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>((*state >> 11) & ((1ULL << 53) - 1)) /
           static_cast<double>(1ULL << 53);
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "root_table_2deg_nearest_node_query_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    problem1::Problem1NearestNodeQueryOptions options{};
    options.residual_tolerance_seconds = 1e-2;
    options.max_newton_iterations = 80;
    options.fallback_direct_solve = true;

    Summary summary{};
    for (const long long linear_index : std::vector<long long>{0, 500000, 800000, 2900000, 5831999}) {
        const auto node = node_angles(linear_index);
        const auto direct = direct_solve(node.nu_A, node.nu_B, node.theta_A);
        const auto query = problem1::query_problem1_from_2deg_nearest_node(
            loader,
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            node.nu_A,
            node.nu_B,
            node.theta_A,
            1,
            1,
            options);
        summary.exact_grid_sample_count += 1;
        accumulate_query_result(query, direct, &summary);
    }

    std::uint64_t rng = 0xabcdefULL;
    for (const long long linear_index : near_node_indices()) {
        const auto node = node_angles(linear_index);
        const double delta_A = (2.0 * unit_random(&rng) - 1.0) * (kPi / 180.0);
        const double delta_B = (2.0 * unit_random(&rng) - 1.0) * (kPi / 180.0);
        const double delta_theta = (2.0 * unit_random(&rng) - 1.0) * (kPi / 180.0);
        const double nu_A = normalize_angle_0_2pi(node.nu_A + delta_A);
        const double nu_B = normalize_angle_0_2pi(node.nu_B + delta_B);
        const double theta_A = normalize_angle_0_2pi(node.theta_A + delta_theta);
        const auto direct = direct_solve(nu_A, nu_B, theta_A);
        const auto query = problem1::query_problem1_from_2deg_nearest_node(
            loader,
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            nu_A,
            nu_B,
            theta_A,
            1,
            1,
            options);
        summary.near_node_sample_count += 1;
        accumulate_query_result(query, direct, &summary);
    }

    const int total_samples = summary.exact_grid_sample_count + summary.near_node_sample_count;
    if (total_samples > 0) {
        summary.fallback_ratio = static_cast<double>(summary.fallback_count) / static_cast<double>(total_samples);
        summary.mean_branch_count /= static_cast<double>(total_samples);
        summary.coverage_vs_direct /= static_cast<double>(total_samples);
        summary.same_root_success_ratio /= static_cast<double>(total_samples);
    }
    if (summary.fallback_ratio >= 0.2) {
        summary.nearest_node_query_ok = false;
    }
    if (summary.derivative_missing_count != 0) {
        summary.nearest_node_query_ok = false;
    }

    std::cout << std::setprecision(12) << std::scientific;
    std::cout << "Problem1RootTable2DegNearestNodeQuerySummary\n";
    std::cout << "exact_grid_sample_count=" << summary.exact_grid_sample_count << '\n';
    std::cout << "near_node_sample_count=" << summary.near_node_sample_count << '\n';
    std::cout << "query_success_count=" << summary.query_success_count << '\n';
    std::cout << "fallback_count=" << summary.fallback_count << '\n';
    std::cout << "fallback_ratio=" << summary.fallback_ratio << '\n';
    std::cout << "mean_branch_count=" << summary.mean_branch_count << '\n';
    std::cout << "coverage_vs_direct=" << summary.coverage_vs_direct << '\n';
    std::cout << "same_root_success_ratio=" << summary.same_root_success_ratio << '\n';
    std::cout << "max_abs_residual_seconds=" << summary.max_abs_residual_seconds << '\n';
    std::cout << "max_time_diff_seconds=" << summary.max_time_diff_seconds << '\n';
    std::cout << "derivative_missing_count=" << summary.derivative_missing_count << '\n';
    std::cout << "nearest_node_query_ok=" << (summary.nearest_node_query_ok ? 1 : 0) << '\n';

    return summary.nearest_node_query_ok ? 0 : 1;
}
