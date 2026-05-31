#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kMidCellOffsetDegrees = 1.0;

struct QuerySample {
    int sample_index = -1;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

struct Problem1RootVirtualGridNode {
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double nu_A_node = 0.0;
    double nu_B_node = 0.0;
    double theta_A_node = 0.0;
    double dnu_A = 0.0;
    double dnu_B = 0.0;
    double dtheta_A = 0.0;
};

Problem1RootVirtualGridNode find_nearest_virtual_root_table_node(
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    const int axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));
    const auto nearest_index = [&](double angle) {
        long long index = std::llround(normalize_angle_0_2pi(angle) / kVirtualGridStepRadians);
        index %= axis_count;
        if (index < 0) {
            index += axis_count;
        }
        return static_cast<int>(index);
    };
    Problem1RootVirtualGridNode node{};
    node.nu_A_index = nearest_index(query_nu_A);
    node.nu_B_index = nearest_index(query_nu_B);
    node.theta_A_index = nearest_index(query_theta_A);
    node.nu_A_node = normalize_angle_0_2pi(static_cast<double>(node.nu_A_index) * kVirtualGridStepRadians);
    node.nu_B_node = normalize_angle_0_2pi(static_cast<double>(node.nu_B_index) * kVirtualGridStepRadians);
    node.theta_A_node = normalize_angle_0_2pi(static_cast<double>(node.theta_A_index) * kVirtualGridStepRadians);
    node.dnu_A = spaceship_cpp::common::normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = spaceship_cpp::common::normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A =
        spaceship_cpp::common::normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

std::vector<QuerySample> build_midcell_samples(int samples_per_group) {
    std::vector<QuerySample> samples;
    const double mid_offset = kMidCellOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < samples_per_group; ++i) {
        const double base_nu_A = normalize_angle_0_2pi(0.91 + static_cast<double>(i) * 1.8123456789);
        const double base_nu_B = normalize_angle_0_2pi(1.41 + static_cast<double>(i) * 2.2718281828);
        const double base_theta_A = normalize_angle_0_2pi(0.63 + static_cast<double>(i) * 1.1447298860);
        const auto nearest = find_nearest_virtual_root_table_node(base_nu_A, base_nu_B, base_theta_A);
        samples.push_back({
            i,
            normalize_angle_0_2pi(nearest.nu_A_node + mid_offset),
            normalize_angle_0_2pi(nearest.nu_B_node + mid_offset),
            normalize_angle_0_2pi(nearest.theta_A_node + mid_offset),
        });
    }
    return samples;
}

int lower_corner_index(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
}

spaceship_cpp::problem1::Problem1RootTable build_local_cell_table_for_query(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    int max_transfer_revolution,
    int max_target_revolution,
    int* lower_i,
    int* lower_j,
    int* lower_k
) {
    const auto nearest = find_nearest_virtual_root_table_node(
        sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
    *lower_i = lower_corner_index(nearest.nu_A_index, nearest.dnu_A);
    *lower_j = lower_corner_index(nearest.nu_B_index, nearest.dnu_B);
    *lower_k = lower_corner_index(nearest.theta_A_index, nearest.dtheta_A);

    spaceship_cpp::problem1::Problem1RootTableConfig config{};
    config.departure_planet = departure_planet;
    config.target_planet = target_planet;
    config.nu_A_start = static_cast<double>(*lower_i) * kVirtualGridStepRadians;
    config.nu_A_step = kVirtualGridStepRadians;
    config.nu_A_count = 2;
    config.nu_B_depart_start = static_cast<double>(*lower_j) * kVirtualGridStepRadians;
    config.nu_B_depart_step = kVirtualGridStepRadians;
    config.nu_B_depart_count = 2;
    config.theta_A_start = static_cast<double>(*lower_k) * kVirtualGridStepRadians;
    config.theta_A_step = kVirtualGridStepRadians;
    config.theta_A_count = 2;
    config.max_transfer_revolution = max_transfer_revolution;
    config.max_target_revolution = max_target_revolution;
    return spaceship_cpp::problem1::build_problem1_root_table(config);
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    std::cout << std::setprecision(6) << std::scientific;
    bool ok = true;
    const auto expect = [&](bool condition, const char* message) {
        if (!condition) {
            std::cerr << "EXPECTATION_FAILED: " << message << '\n';
            ok = false;
        }
    };
    const auto map_value = [](const std::map<int, int>& values, int key) {
        const auto it = values.find(key);
        return it == values.end() ? 0 : it->second;
    };

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_target_revolution = 1;
    const auto midcell_samples = build_midcell_samples(12);

    int lower_i = 0;
    int lower_j = 0;
    int lower_k = 0;

    const auto sample3_table = build_local_cell_table_for_query(
        departure_planet,
        target_planet,
        midcell_samples.at(3),
        0,
        max_target_revolution,
        &lower_i,
        &lower_j,
        &lower_k);
    const auto sample3_safe = problem1::query_problem1_root_table_route_b_linear_safe(
        sample3_table,
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A);
    std::cout << "sample3_linear_safe_fallback_required=" << sample3_safe.fallback_required << '\n';
    std::cout << "sample3_linear_safe_reason=" << sample3_safe.reason << '\n';
    std::cout << "sample3_linear_safe_approximation_count=" << sample3_safe.approximations.size() << '\n';
    std::cout << "sample3_linear_safe_branch_count_complete=" << sample3_safe.branch_count_complete << '\n';
    std::cout << "sample3_expected_count_k0=" << map_value(sample3_safe.expected_count_by_k, 0) << '\n';
    std::cout << "sample3_candidate_count_k0=" << map_value(sample3_safe.candidate_count_by_k, 0) << '\n';
    std::cout << "sample3_missing_count_k0=" << map_value(sample3_safe.missing_count_by_k, 0) << '\n';
    std::cout << "sample3_extra_count_k0=" << map_value(sample3_safe.extra_count_by_k, 0) << '\n';
    expect(!sample3_safe.fallback_required, "sample3 admissible cell should not fallback in linear Route B");
    expect(sample3_safe.branch_count_complete, "sample3 admissible cell should have complete branch count");
    expect(!sample3_safe.approximations.empty(), "sample3 admissible cell should produce linear approximations");
    expect(map_value(sample3_safe.expected_count_by_k, 0) == 3, "sample3 expected_count_by_k[0] should be 3");
    expect(map_value(sample3_safe.candidate_count_by_k, 0) == 3, "sample3 candidate_count_by_k[0] should be 3");
    expect(map_value(sample3_safe.missing_count_by_k, 0) == 0, "sample3 missing_count_by_k[0] should be 0");
    expect(map_value(sample3_safe.extra_count_by_k, 0) == 0, "sample3 extra_count_by_k[0] should be 0");
    for (const auto& approximation : sample3_safe.approximations) {
        expect(approximation.method == "route_b_linear_no_newton", "linear Route B method label mismatch");
    }

    const auto sample5_table = build_local_cell_table_for_query(
        departure_planet,
        target_planet,
        midcell_samples.at(5),
        0,
        max_target_revolution,
        &lower_i,
        &lower_j,
        &lower_k);
    const auto sample5_safe = problem1::query_problem1_root_table_route_b_linear_safe(
        sample5_table,
        midcell_samples.at(5).query_nu_A,
        midcell_samples.at(5).query_nu_B,
        midcell_samples.at(5).query_theta_A);
    std::cout << "sample5_linear_safe_fallback_required=" << sample5_safe.fallback_required << '\n';
    std::cout << "sample5_linear_safe_reason=" << sample5_safe.reason << '\n';
    expect(sample5_safe.fallback_required, "sample5 non-admissible cell should fallback in linear Route B");
    expect(
        sample5_safe.reason == "cell_non_admissible_branch_count_inconsistent",
        "sample5 non-admissible cell should use admissibility fallback reason");

    const auto sample3_nearest_global = find_nearest_virtual_root_table_node(
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A);
    const auto sample3_node_branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        sample3_nearest_global.nu_A_node,
        sample3_nearest_global.nu_B_node,
        sample3_nearest_global.theta_A_node,
        0,
        max_target_revolution);
    std::vector<const problem1::Problem1SolutionBranch*> k0_branches;
    for (const auto& branch : sample3_node_branches) {
        if (branch.valid && branch.transfer_revolution == 0) {
            k0_branches.push_back(&branch);
        }
    }
    expect(k0_branches.size() == 3, "sample3 should have 3 k=0 node branches");
    const auto* middle_k0_source = k0_branches.size() >= 2 ? k0_branches.at(1) : nullptr;
    expect(middle_k0_source != nullptr, "sample3 should have middle k=0 source branch");
    if (middle_k0_source != nullptr) {
        const auto approximation = problem1::evaluate_problem1_root_linear_route_b_approximation_from_node(
            departure_planet,
            target_planet,
            sample3_nearest_global.nu_A_node,
            sample3_nearest_global.nu_B_node,
            sample3_nearest_global.theta_A_node,
            *middle_k0_source,
            midcell_samples.at(3).query_nu_A,
            midcell_samples.at(3).query_nu_B,
            midcell_samples.at(3).query_theta_A,
            max_target_revolution);
        std::cout << "sample3_linear_source_q=" << middle_k0_source->target_revolution << '\n';
        std::cout << "sample3_linear_selected_q=" << approximation.diagnostics.selected_target_revolution << '\n';
        std::cout << "sample3_linear_q_changed=" << approximation.diagnostics.q_sheet_selection_changed << '\n';
        std::cout << "sample3_linear_q_failed=" << approximation.diagnostics.q_sheet_selection_failed << '\n';
        expect(middle_k0_source->target_revolution == 0, "sample3 middle source should start with q=0");
        expect(
            approximation.diagnostics.selected_target_revolution == 1,
            "sample3 middle source should select q=1 in linear Route B");
        expect(approximation.diagnostics.q_sheet_selection_changed, "sample3 linear Route B should mark q changed");
        expect(!approximation.diagnostics.q_sheet_selection_failed, "sample3 linear Route B q selection should succeed");
        expect(approximation.method == "route_b_linear_no_newton", "linear helper method label mismatch");
    }

    std::cout << "route_b_linear_safety_ok\n";
    return ok ? 0 : 1;
}
