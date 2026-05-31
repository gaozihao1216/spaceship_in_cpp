#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
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

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const double hessian_step = 5e-4;
    const double tangent_residual_tolerance = problem1::kProblem1RootExperimentalTangentResidualTolerance;
    bool ok = true;
    const auto expect = [&](bool condition, const char* message) {
        if (!condition) {
            std::cerr << "EXPECTATION_FAILED: " << message << '\n';
            ok = false;
        }
    };

    const auto midcell_samples = build_midcell_samples(12);

    int lower_i = 0;
    int lower_j = 0;
    int lower_k = 0;

    const auto sample3_admissible_table = build_local_cell_table_for_query(
        departure_planet,
        target_planet,
        midcell_samples.at(3),
        0,
        max_target_revolution,
        &lower_i,
        &lower_j,
        &lower_k);
    const auto sample3_safe = problem1::query_problem1_root_table_route_b_safe(
        sample3_admissible_table,
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A,
        hessian_step,
        problem1::Problem1RootHessianMethod::TangentFiniteDifference,
        tangent_residual_tolerance);
    std::cout << "sample3_route_b_safe_fallback_required=" << sample3_safe.fallback_required << '\n';
    std::cout << "sample3_route_b_safe_reason=" << sample3_safe.reason << '\n';
    std::cout << "sample3_route_b_safe_approximation_count=" << sample3_safe.approximations.size() << '\n';
    expect(!sample3_safe.fallback_required, "sample3 admissible cell should not require fallback");
    expect(sample3_safe.cell_admissibility.admissible, "sample3 admissible cell should be admissible");
    const auto sample5_table = build_local_cell_table_for_query(
        departure_planet,
        target_planet,
        midcell_samples.at(5),
        0,
        max_target_revolution,
        &lower_i,
        &lower_j,
        &lower_k);
    const auto sample5_safe = problem1::query_problem1_root_table_route_b_safe(
        sample5_table,
        midcell_samples.at(5).query_nu_A,
        midcell_samples.at(5).query_nu_B,
        midcell_samples.at(5).query_theta_A,
        hessian_step,
        problem1::Problem1RootHessianMethod::TangentFiniteDifference,
        tangent_residual_tolerance);
    std::cout << "sample5_route_b_safe_fallback_required=" << sample5_safe.fallback_required << '\n';
    std::cout << "sample5_route_b_safe_reason=" << sample5_safe.reason << '\n';
    expect(sample5_safe.fallback_required, "sample5 non-admissible cell should require fallback");
    expect(
        sample5_safe.reason == "cell_non_admissible_branch_count_inconsistent",
        "sample5 non-admissible cell should use branch_count_inconsistent reason");

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
    expect(k0_branches.size() == 3, "sample3 fixed q-sheet case should have 3 k=0 node branches");
    const problem1::Problem1SolutionBranch* middle_k0_source = k0_branches.size() >= 2 ? k0_branches.at(1) : nullptr;
    expect(middle_k0_source != nullptr, "sample3 fixed q-sheet case should expose middle source branch");
    if (middle_k0_source == nullptr) {
        return 1;
    }
    const auto linear_prediction = problem1::predict_problem1_root_branch_linear_from_node(
        departure_planet,
        target_planet,
        sample3_nearest_global.nu_A_node,
        sample3_nearest_global.nu_B_node,
        sample3_nearest_global.theta_A_node,
        *middle_k0_source,
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A);
    expect(linear_prediction.valid, "sample3 fixed source branch should have valid linear prediction");
    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        departure_planet,
        target_planet,
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A,
        middle_k0_source->transfer_revolution,
        linear_prediction.predicted_encounter_global_angle,
        *middle_k0_source,
        max_target_revolution);
    std::cout << "sample3_helper_selected_q=" << q_selection.selected_q << '\n';
    std::cout << "sample3_helper_q_changed=" << q_selection.q_changed << '\n';
    std::cout << "sample3_helper_selection_failed=" << q_selection.selection_failed << '\n';
    expect(q_selection.selected_q == 1, "sample3 helper q-sheet selection should choose q=1");
    expect(q_selection.q_changed, "sample3 helper q-sheet selection should mark q_changed");
    expect(!q_selection.selection_failed, "sample3 helper q-sheet selection should not fail");

    const auto quadratic_qsheet = problem1::evaluate_problem1_root_quadratic_approximation_from_node(
        departure_planet,
        target_planet,
        sample3_nearest_global.nu_A_node,
        sample3_nearest_global.nu_B_node,
        sample3_nearest_global.theta_A_node,
        *middle_k0_source,
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A,
        hessian_step,
        problem1::Problem1RootHessianMethod::TangentFiniteDifference,
        tangent_residual_tolerance,
        max_target_revolution);
    std::cout << "sample3_qsheet_source_q=" << middle_k0_source->target_revolution << '\n';
    std::cout << "sample3_qsheet_selected_q=" << quadratic_qsheet.diagnostics.selected_target_revolution << '\n';
    std::cout << "sample3_qsheet_q_changed=" << quadratic_qsheet.diagnostics.q_sheet_selection_changed << '\n';
    std::cout << "sample3_qsheet_selection_failed=" << quadratic_qsheet.diagnostics.q_sheet_selection_failed << '\n';
    expect(middle_k0_source->target_revolution == 0, "sample3 fixed source branch should start with source_q=0");
    expect(!quadratic_qsheet.diagnostics.q_sheet_selection_failed, "sample3 q-sheet should not fail");

    const auto disabled_newton_hessian = problem1::query_problem1_root_table_route_b_safe(
        sample3_admissible_table,
        midcell_samples.at(3).query_nu_A,
        midcell_samples.at(3).query_nu_B,
        midcell_samples.at(3).query_theta_A,
        hessian_step,
        problem1::Problem1RootHessianMethod::NewtonRefinedFiniteDifference,
        tangent_residual_tolerance);
    std::cout << "newton_hessian_guard_fallback_required=" << disabled_newton_hessian.fallback_required << '\n';
    std::cout << "newton_hessian_guard_reason=" << disabled_newton_hessian.reason << '\n';
    expect(disabled_newton_hessian.fallback_required, "Route B safe wrapper should reject Newton Hessian mode");
    expect(
        disabled_newton_hessian.reason == "route_b_newton_hessian_method_disabled",
        "Route B safe wrapper should use Newton Hessian disabled reason");

    std::cout << "route_b_safety_ok\n";
    return ok ? 0 : 1;
}
