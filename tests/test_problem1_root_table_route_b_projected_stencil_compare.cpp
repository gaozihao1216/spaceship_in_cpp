#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kNearNodeOffsetDegrees = 0.1;
constexpr double kMidCellOffsetDegrees = 1.0;
constexpr double kEngineeringResidualToleranceSeconds = 1e-2;
constexpr double kSameRootTimeThresholdSeconds = 1.0;
constexpr double kSameRootAlphaThreshold = 1e-3;
constexpr int kSamplesPerGroup = 12;
constexpr std::array<double, 4> kAdaptiveHessianSteps{{1e-5, 5e-6, 2e-6, 1e-6}};

struct QuerySample {
    std::string group_name;
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

struct RouteAEvalResult {
    bool valid = false;
    std::string invalid_reason;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct OldStencilEval {
    std::string stencil_id;
    bool residual_valid = false;
    double residual_seconds = std::numeric_limits<double>::quiet_NaN();
    bool residual_gate_passed = false;
    bool attach_success = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct ProjectedStencilEval {
    std::string stencil_id;
    bool valid = false;
    std::string invalid_reason;
    double residual_before_seconds = std::numeric_limits<double>::quiet_NaN();
    double residual_after_seconds = std::numeric_limits<double>::quiet_NaN();
    double delta_alpha = std::numeric_limits<double>::quiet_NaN();
    double alpha_tangent = std::numeric_limits<double>::quiet_NaN();
    double alpha_projected = std::numeric_limits<double>::quiet_NaN();
    bool derivative_before_valid = false;
    bool attach_success = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct StepEvalResult {
    double hessian_step = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
    std::string invalid_reason;
    int stencil_pass_count = 0;
    int attach_success_count = 0;
    double max_abs_stencil_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    std::string worst_stencil_id;
    double max_projected_stencil_delta_alpha = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_residual_before = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_residual_after = std::numeric_limits<double>::quiet_NaN();
    double residual_after_correction = std::numeric_limits<double>::quiet_NaN();
    double time_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
    double alpha_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
};

struct AdaptiveEvalResult {
    bool valid = false;
    std::string invalid_reason;
    double selected_step = std::numeric_limits<double>::quiet_NaN();
    double residual_after_correction = std::numeric_limits<double>::quiet_NaN();
    double time_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
    double alpha_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_delta_alpha = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_residual_before = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_residual_after = std::numeric_limits<double>::quiet_NaN();
};

struct RecoveredCase {
    std::string group;
    int sample_index = -1;
    int source_index = -1;
    int source_rank_in_k = -1;
    int k = 0;
    int source_q = 0;
    double selected_step = std::numeric_limits<double>::quiet_NaN();
    std::string old_invalid_reason;
    double projected_residual_after = std::numeric_limits<double>::quiet_NaN();
    double projected_time_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
    double projected_alpha_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_delta_alpha = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_residual_before = std::numeric_limits<double>::quiet_NaN();
    double max_projected_stencil_residual_after = std::numeric_limits<double>::quiet_NaN();
};

struct SuspiciousCase {
    std::string group;
    int sample_index = -1;
    int source_index = -1;
    int source_rank_in_k = -1;
    int k = 0;
    double selected_step = std::numeric_limits<double>::quiet_NaN();
    double residual_after = std::numeric_limits<double>::quiet_NaN();
    double time_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
    double alpha_diff_vs_route_a = std::numeric_limits<double>::quiet_NaN();
};

struct GroupSummary {
    int route_a_valid_count = 0;
    int old_route_b_valid_count = 0;
    int projected_route_b_valid_count = 0;
    int newly_recovered_count = 0;
    int projected_lost_count = 0;
    int projected_same_root_success_count = 0;
    int projected_same_root_valid_count = 0;
    double sum_time_diff_vs_route_a = 0.0;
    double max_time_diff_vs_route_a = 0.0;
    double sum_alpha_diff_vs_route_a = 0.0;
    double max_alpha_diff_vs_route_a = 0.0;
    double sum_abs_residual_seconds = 0.0;
    double max_abs_residual_seconds = 0.0;
    std::map<std::string, int> old_invalid_reasons;
    std::map<std::string, int> projected_invalid_reasons;
    std::map<std::string, int> projected_selected_step_distribution;
    std::vector<RecoveredCase> recovered_cases;
    std::vector<SuspiciousCase> suspicious_cases;
};

double wrapped_alpha_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

int lower_corner_index(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
}

std::string step_label(double step) {
    if (std::abs(step - 1e-5) < 1e-12) {
        return "1e-5";
    }
    if (std::abs(step - 5e-6) < 1e-12) {
        return "5e-6";
    }
    if (std::abs(step - 2e-6) < 1e-12) {
        return "2e-6";
    }
    if (std::abs(step - 1e-6) < 1e-12) {
        return "1e-6";
    }
    return "other";
}

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
    node.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

std::vector<QuerySample> build_samples(int samples_per_group) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::vector<QuerySample> samples;
    const double near_offset = kNearNodeOffsetDegrees * kPi / 180.0;
    const double mid_offset = kMidCellOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < samples_per_group; ++i) {
        const double base_nu_A = normalize_angle_0_2pi(0.37 + static_cast<double>(i) * 2.3999632297);
        const double base_nu_B = normalize_angle_0_2pi(1.11 + static_cast<double>(i) * 1.7548776662);
        const double base_theta_A = normalize_angle_0_2pi(0.23 + static_cast<double>(i) * 0.9182736455);
        const auto near_node = find_nearest_virtual_root_table_node(base_nu_A, base_nu_B, base_theta_A);
        samples.push_back({"near_node",
                           normalize_angle_0_2pi(near_node.nu_A_node + near_offset),
                           normalize_angle_0_2pi(near_node.nu_B_node + near_offset),
                           normalize_angle_0_2pi(near_node.theta_A_node + near_offset)});

        const double base2_nu_A = normalize_angle_0_2pi(0.91 + static_cast<double>(i) * 1.8123456789);
        const double base2_nu_B = normalize_angle_0_2pi(1.41 + static_cast<double>(i) * 2.2718281828);
        const double base2_theta_A = normalize_angle_0_2pi(0.63 + static_cast<double>(i) * 1.1447298860);
        const auto mid_node = find_nearest_virtual_root_table_node(base2_nu_A, base2_nu_B, base2_theta_A);
        samples.push_back({"mid_cell",
                           normalize_angle_0_2pi(mid_node.nu_A_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.nu_B_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.theta_A_node + mid_offset)});
    }
    const double earth_period = planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
    const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};
    for (int i = 0; i < samples_per_group; ++i) {
        const double launch_fraction = std::fmod(0.17 + 0.31 * static_cast<double>(i), 1.0);
        const double launch_time = launch_fraction * earth_period;
        const auto departure_state = planet_params::planet_state_at_time(planet_params::PlanetId::Earth, launch_time);
        const auto target_state = planet_params::planet_state_at_time(planet_params::PlanetId::Mars, launch_time);
        const double transfer_perihelion_angle =
            transfer_perihelion_angles[static_cast<std::size_t>(i) % transfer_perihelion_angles.size()];
        samples.push_back({"physical_launch",
                           departure_state.varphi,
                           target_state.varphi,
                           normalize_angle_0_2pi(departure_state.theta_global - transfer_perihelion_angle)});
    }
    return samples;
}

spaceship_cpp::problem1::Problem1RootTable build_local_cell_table_for_query(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    int max_transfer_revolution,
    int max_target_revolution
) {
    const auto nearest = find_nearest_virtual_root_table_node(
        sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
    const int lower_i = lower_corner_index(nearest.nu_A_index, nearest.dnu_A);
    const int lower_j = lower_corner_index(nearest.nu_B_index, nearest.dnu_B);
    const int lower_k = lower_corner_index(nearest.theta_A_index, nearest.dtheta_A);

    spaceship_cpp::problem1::Problem1RootTableConfig config{};
    config.departure_planet = departure_planet;
    config.target_planet = target_planet;
    config.nu_A_start = static_cast<double>(lower_i) * kVirtualGridStepRadians;
    config.nu_A_step = kVirtualGridStepRadians;
    config.nu_A_count = 2;
    config.nu_B_depart_start = static_cast<double>(lower_j) * kVirtualGridStepRadians;
    config.nu_B_depart_step = kVirtualGridStepRadians;
    config.nu_B_depart_count = 2;
    config.theta_A_start = static_cast<double>(lower_k) * kVirtualGridStepRadians;
    config.theta_A_step = kVirtualGridStepRadians;
    config.theta_A_count = 2;
    config.max_transfer_revolution = max_transfer_revolution;
    config.max_target_revolution = max_target_revolution;
    return spaceship_cpp::problem1::build_problem1_root_table(config);
}

RouteAEvalResult evaluate_route_a_from_source_branch_test_only(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    RouteAEvalResult result{};
    if (!source_branch.valid) {
        result.invalid_reason = "source_branch_invalid";
        return result;
    }
    const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
        table.config().departure_planet,
        table.config().target_planet,
        nearest.node_nu_A,
        nearest.node_nu_B,
        nearest.node_theta_A,
        source_branch,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    if (!attached.valid || !attached.derivatives_available) {
        result.invalid_reason = "attach_problem1_root_derivatives_failed";
        return result;
    }
    const double dnu_A = normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A);
    const double dnu_B = normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B);
    const double dtheta_A = normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A);
    const double alpha_linear = normalize_angle_0_2pi(
        attached.encounter_global_angle +
        attached.d_encounter_global_angle_d_nu_A * dnu_A +
        attached.d_encounter_global_angle_d_nu_B * dnu_B +
        attached.d_encounter_global_angle_d_theta_A * dtheta_A);
    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        source_branch.transfer_revolution,
        alpha_linear,
        source_branch,
        table.config().max_target_revolution);
    if (q_selection.selection_failed) {
        result.invalid_reason = "q_sheet_selection_failed";
        return result;
    }
    const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        source_branch.transfer_revolution,
        q_selection.selected_q,
        alpha_linear,
        80,
        1e-2,
        1e-12,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    if (!refined.valid) {
        result.invalid_reason =
            refined.diagnostic.invalid_reason.empty() ? "route_a_refinement_failed" : refined.diagnostic.invalid_reason;
        return result;
    }
    result.valid = true;
    result.branch = refined.branch;
    return result;
}

OldStencilEval evaluate_old_stencil(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& differentiated_source_branch,
    double stencil_nu_A,
    double stencil_nu_B,
    double stencil_theta_A,
    std::string stencil_id
) {
    namespace problem1 = spaceship_cpp::problem1;
    OldStencilEval eval{};
    eval.stencil_id = std::move(stencil_id);
    const double dx1 = normalize_angle_minus_pi_pi(stencil_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(stencil_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(stencil_theta_A - node_theta_A);
    const double alpha_tangent = normalize_angle_0_2pi(
        differentiated_source_branch.encounter_global_angle +
        differentiated_source_branch.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated_source_branch.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated_source_branch.d_encounter_global_angle_d_theta_A * dx3);
    const auto residual = problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        alpha_tangent,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    eval.residual_valid = residual.valid && std::isfinite(residual.residual_seconds);
    eval.residual_seconds = residual.residual_seconds;
    eval.residual_gate_passed = eval.residual_valid &&
        std::abs(residual.residual_seconds) <= kEngineeringResidualToleranceSeconds;
    if (!eval.residual_gate_passed) {
        return eval;
    }
    spaceship_cpp::problem1::Problem1SolutionBranch tangent_branch{};
    tangent_branch.valid = true;
    tangent_branch.encounter_global_angle = residual.encounter_global_angle;
    tangent_branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    tangent_branch.transfer_revolution = differentiated_source_branch.transfer_revolution;
    tangent_branch.target_revolution = differentiated_source_branch.target_revolution;
    tangent_branch.time_of_flight_seconds = residual.transfer_time_seconds;
    tangent_branch.target_time_seconds = residual.target_time_seconds;
    tangent_branch.residual_seconds = residual.residual_seconds;
    tangent_branch.transfer_e = residual.transfer_e;
    tangent_branch.transfer_p = residual.transfer_p;
    tangent_branch.transfer_a = residual.transfer_a;
    tangent_branch.theta_B = residual.theta_B;
    tangent_branch = problem1::attach_problem1_root_derivatives(
        departure_planet, target_planet, stencil_nu_A, stencil_nu_B, stencil_theta_A, tangent_branch);
    eval.attach_success = tangent_branch.valid && tangent_branch.derivatives_available;
    eval.branch = tangent_branch;
    return eval;
}

ProjectedStencilEval evaluate_projected_stencil_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& differentiated_source_branch,
    double stencil_nu_A,
    double stencil_nu_B,
    double stencil_theta_A,
    std::string stencil_id
) {
    namespace problem1 = spaceship_cpp::problem1;
    ProjectedStencilEval eval{};
    eval.stencil_id = std::move(stencil_id);
    const double dx1 = normalize_angle_minus_pi_pi(stencil_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(stencil_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(stencil_theta_A - node_theta_A);
    eval.alpha_tangent = normalize_angle_0_2pi(
        differentiated_source_branch.encounter_global_angle +
        differentiated_source_branch.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated_source_branch.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated_source_branch.d_encounter_global_angle_d_theta_A * dx3);
    const auto residual_before = problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        eval.alpha_tangent,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    const auto derivative_before = problem1::evaluate_problem1_root_residual_derivatives(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        eval.alpha_tangent,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    eval.derivative_before_valid = derivative_before.valid && std::isfinite(derivative_before.R_alpha);
    eval.residual_before_seconds = residual_before.residual_seconds;
    if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds) ||
        !eval.derivative_before_valid) {
        eval.invalid_reason = "projected_stencil_residual_or_derivative_invalid";
        return eval;
    }
    if (std::abs(derivative_before.R_alpha) <= 1e-12) {
        eval.invalid_reason = "projected_stencil_derivative_too_small";
        return eval;
    }
    const double residual_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    eval.delta_alpha = -residual_scale_free / derivative_before.R_alpha;
    eval.alpha_projected = normalize_angle_0_2pi(eval.alpha_tangent + eval.delta_alpha);
    const auto residual_after = problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        eval.alpha_projected,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
        eval.invalid_reason = "projected_stencil_residual_after_invalid";
        return eval;
    }
    eval.residual_after_seconds = residual_after.residual_seconds;
    if (std::abs(residual_after.residual_seconds) > kEngineeringResidualToleranceSeconds) {
        eval.invalid_reason = "projected_stencil_residual_after_too_large";
        return eval;
    }
    spaceship_cpp::problem1::Problem1SolutionBranch projected_branch{};
    projected_branch.valid = true;
    projected_branch.encounter_global_angle = residual_after.encounter_global_angle;
    projected_branch.target_arrival_true_anomaly = residual_after.target_arrival_true_anomaly;
    projected_branch.transfer_revolution = differentiated_source_branch.transfer_revolution;
    projected_branch.target_revolution = differentiated_source_branch.target_revolution;
    projected_branch.time_of_flight_seconds = residual_after.transfer_time_seconds;
    projected_branch.target_time_seconds = residual_after.target_time_seconds;
    projected_branch.residual_seconds = residual_after.residual_seconds;
    projected_branch.transfer_e = residual_after.transfer_e;
    projected_branch.transfer_p = residual_after.transfer_p;
    projected_branch.transfer_a = residual_after.transfer_a;
    projected_branch.theta_B = residual_after.theta_B;
    projected_branch = problem1::attach_problem1_root_derivatives(
        departure_planet, target_planet, stencil_nu_A, stencil_nu_B, stencil_theta_A, projected_branch);
    eval.attach_success = projected_branch.valid && projected_branch.derivatives_available;
    if (!eval.attach_success) {
        eval.invalid_reason = "hessian_attach_failed";
        return eval;
    }
    eval.valid = true;
    eval.branch = projected_branch;
    return eval;
}

StepEvalResult evaluate_old_step(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    const spaceship_cpp::problem1::Problem1SolutionBranch& route_a_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step
) {
    namespace problem1 = spaceship_cpp::problem1;
    StepEvalResult result{};
    result.hessian_step = hessian_step;
    auto differentiated = source_branch;
    if (!differentiated.derivatives_available) {
        differentiated = problem1::attach_problem1_root_derivatives(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            source_branch);
    }
    if (!differentiated.derivatives_available) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }
    const std::array<OldStencilEval, 6> stencils{{
        evaluate_old_stencil(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, normalize_angle_0_2pi(nearest.node_nu_A + hessian_step), nearest.node_nu_B, nearest.node_theta_A, "nu_A_plus"),
        evaluate_old_stencil(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, normalize_angle_0_2pi(nearest.node_nu_A - hessian_step), nearest.node_nu_B, nearest.node_theta_A, "nu_A_minus"),
        evaluate_old_stencil(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, normalize_angle_0_2pi(nearest.node_nu_B + hessian_step), nearest.node_theta_A, "nu_B_plus"),
        evaluate_old_stencil(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, normalize_angle_0_2pi(nearest.node_nu_B - hessian_step), nearest.node_theta_A, "nu_B_minus"),
        evaluate_old_stencil(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, nearest.node_nu_B, normalize_angle_0_2pi(nearest.node_theta_A + hessian_step), "theta_A_plus"),
        evaluate_old_stencil(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, nearest.node_nu_B, normalize_angle_0_2pi(nearest.node_theta_A - hessian_step), "theta_A_minus"),
    }};
    double max_abs = -1.0;
    for (const auto& stencil : stencils) {
        if (stencil.residual_gate_passed) {
            result.stencil_pass_count += 1;
        }
        if (stencil.attach_success) {
            result.attach_success_count += 1;
        }
        const double metric = stencil.residual_valid ? std::abs(stencil.residual_seconds) : std::numeric_limits<double>::infinity();
        if (metric > max_abs) {
            max_abs = metric;
            result.max_abs_stencil_residual_seconds = stencil.residual_seconds;
            result.worst_stencil_id = stencil.stencil_id;
        }
    }
    if (result.stencil_pass_count != 6) {
        result.invalid_reason = "hessian_stencil_seconds_gate_failed";
        return result;
    }
    if (result.attach_success_count != 6) {
        result.invalid_reason = "hessian_attach_failed";
        return result;
    }
    const double inv_2h = 1.0 / (2.0 * hessian_step);
    const double H_aa = (stencils[0].branch.d_encounter_global_angle_d_nu_A - stencils[1].branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H_bb = (stencils[2].branch.d_encounter_global_angle_d_nu_B - stencils[3].branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H_cc = (stencils[4].branch.d_encounter_global_angle_d_theta_A - stencils[5].branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H12_from_g1 = (stencils[2].branch.d_encounter_global_angle_d_nu_A - stencils[3].branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 = (stencils[0].branch.d_encounter_global_angle_d_nu_B - stencils[1].branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 = (stencils[4].branch.d_encounter_global_angle_d_nu_A - stencils[5].branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 = (stencils[0].branch.d_encounter_global_angle_d_theta_A - stencils[1].branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H23_from_g2 = (stencils[4].branch.d_encounter_global_angle_d_nu_B - stencils[5].branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H32_from_g3 = (stencils[2].branch.d_encounter_global_angle_d_theta_A - stencils[3].branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H_ab = 0.5 * (H12_from_g1 + H21_from_g2);
    const double H_ac = 0.5 * (H13_from_g1 + H31_from_g3);
    const double H_bc = 0.5 * (H23_from_g2 + H32_from_g3);
    if (!(std::isfinite(H_aa) && std::isfinite(H_bb) && std::isfinite(H_cc) && std::isfinite(H_ab) && std::isfinite(H_ac) && std::isfinite(H_bc))) {
        result.invalid_reason = "hessian_non_finite";
        return result;
    }
    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A);
    const double alpha_linear = normalize_angle_0_2pi(
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated.d_encounter_global_angle_d_theta_A * dx3);
    const double alpha_quadratic = normalize_angle_0_2pi(
        alpha_linear + 0.5 * (H_aa * dx1 * dx1 + H_bb * dx2 * dx2 + H_cc * dx3 * dx3 +
                              2.0 * H_ab * dx1 * dx2 + 2.0 * H_ac * dx1 * dx3 + 2.0 * H_bc * dx2 * dx3));
    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        source_branch.transfer_revolution, alpha_quadratic, source_branch, table.config().max_target_revolution);
    if (q_selection.selection_failed) {
        result.invalid_reason = "q_sheet_selection_failed";
        return result;
    }
    const auto residual_before = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        alpha_quadratic, source_branch.transfer_revolution, q_selection.selected_q);
    if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
        result.invalid_reason = "residual_before_correction_invalid";
        return result;
    }
    const auto derivatives_before = problem1::evaluate_problem1_root_residual_derivatives(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        alpha_quadratic, source_branch.transfer_revolution, q_selection.selected_q);
    if (!derivatives_before.valid || !std::isfinite(derivatives_before.R_alpha)) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }
    if (std::abs(derivatives_before.R_alpha) <= 1e-12) {
        result.invalid_reason = "correction_derivative_too_small";
        return result;
    }
    const double residual_scale_free = problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    const double delta_alpha = -residual_scale_free / derivatives_before.R_alpha;
    const double alpha_corrected = normalize_angle_0_2pi(alpha_quadratic + delta_alpha);
    const auto residual_after = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        alpha_corrected, source_branch.transfer_revolution, q_selection.selected_q);
    if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
        result.invalid_reason = "residual_after_correction_invalid";
        return result;
    }
    result.residual_after_correction = residual_after.residual_seconds;
    if (std::abs(residual_after.residual_seconds) > kEngineeringResidualToleranceSeconds) {
        result.invalid_reason = "residual_after_correction_too_large";
        return result;
    }
    result.valid = true;
    result.time_diff_vs_route_a = std::abs(residual_after.transfer_time_seconds - route_a_branch.time_of_flight_seconds);
    result.alpha_diff_vs_route_a = wrapped_alpha_distance(residual_after.encounter_global_angle, route_a_branch.encounter_global_angle);
    return result;
}

StepEvalResult evaluate_projected_step(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    const spaceship_cpp::problem1::Problem1SolutionBranch& route_a_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step
) {
    namespace problem1 = spaceship_cpp::problem1;
    StepEvalResult result{};
    result.hessian_step = hessian_step;
    auto differentiated = source_branch;
    if (!differentiated.derivatives_available) {
        differentiated = problem1::attach_problem1_root_derivatives(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            source_branch);
    }
    if (!differentiated.derivatives_available) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }
    const std::array<ProjectedStencilEval, 6> stencils{{
        evaluate_projected_stencil_test_only(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, normalize_angle_0_2pi(nearest.node_nu_A + hessian_step), nearest.node_nu_B, nearest.node_theta_A, "nu_A_plus"),
        evaluate_projected_stencil_test_only(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, normalize_angle_0_2pi(nearest.node_nu_A - hessian_step), nearest.node_nu_B, nearest.node_theta_A, "nu_A_minus"),
        evaluate_projected_stencil_test_only(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, normalize_angle_0_2pi(nearest.node_nu_B + hessian_step), nearest.node_theta_A, "nu_B_plus"),
        evaluate_projected_stencil_test_only(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, normalize_angle_0_2pi(nearest.node_nu_B - hessian_step), nearest.node_theta_A, "nu_B_minus"),
        evaluate_projected_stencil_test_only(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, nearest.node_nu_B, normalize_angle_0_2pi(nearest.node_theta_A + hessian_step), "theta_A_plus"),
        evaluate_projected_stencil_test_only(table.config().departure_planet, table.config().target_planet, nearest.node_nu_A, nearest.node_nu_B, nearest.node_theta_A, differentiated, nearest.node_nu_A, nearest.node_nu_B, normalize_angle_0_2pi(nearest.node_theta_A - hessian_step), "theta_A_minus"),
    }};
    double max_res_before = 0.0;
    double max_res_after = 0.0;
    double max_delta_alpha = 0.0;
    for (const auto& stencil : stencils) {
        if (stencil.valid) {
            result.stencil_pass_count += 1;
            result.attach_success_count += 1;
        }
        if (!std::isnan(stencil.residual_before_seconds) && std::abs(stencil.residual_before_seconds) >= max_res_before) {
            max_res_before = std::abs(stencil.residual_before_seconds);
        }
        if (!std::isnan(stencil.residual_after_seconds) && std::abs(stencil.residual_after_seconds) >= max_res_after) {
            max_res_after = std::abs(stencil.residual_after_seconds);
        }
        if (!std::isnan(stencil.delta_alpha) && std::abs(stencil.delta_alpha) >= max_delta_alpha) {
            max_delta_alpha = std::abs(stencil.delta_alpha);
        }
        const double metric = std::isnan(stencil.residual_after_seconds)
            ? (std::isnan(stencil.residual_before_seconds) ? std::numeric_limits<double>::infinity() : std::abs(stencil.residual_before_seconds))
            : std::abs(stencil.residual_after_seconds);
        if (std::isnan(result.max_abs_stencil_residual_seconds) || metric >= std::abs(result.max_abs_stencil_residual_seconds)) {
            result.max_abs_stencil_residual_seconds = std::isnan(stencil.residual_after_seconds) ? stencil.residual_before_seconds : stencil.residual_after_seconds;
            result.worst_stencil_id = stencil.stencil_id;
        }
    }
    result.max_projected_stencil_delta_alpha = max_delta_alpha;
    result.max_projected_stencil_residual_before = max_res_before;
    result.max_projected_stencil_residual_after = max_res_after;
    if (result.stencil_pass_count != 6) {
        bool any_derivative = false;
        bool any_after_invalid = false;
        bool any_after_large = false;
        bool any_attach_failed = false;
        for (const auto& stencil : stencils) {
            any_derivative = any_derivative || stencil.invalid_reason == "projected_stencil_derivative_too_small" ||
                stencil.invalid_reason == "projected_stencil_residual_or_derivative_invalid";
            any_after_invalid = any_after_invalid || stencil.invalid_reason == "projected_stencil_residual_after_invalid";
            any_after_large = any_after_large || stencil.invalid_reason == "projected_stencil_residual_after_too_large";
            any_attach_failed = any_attach_failed || stencil.invalid_reason == "hessian_attach_failed";
        }
        if (any_attach_failed) {
            result.invalid_reason = "hessian_attach_failed";
        } else if (any_after_large) {
            result.invalid_reason = "projected_stencil_residual_after_too_large";
        } else if (any_after_invalid) {
            result.invalid_reason = "projected_stencil_residual_after_invalid";
        } else if (any_derivative) {
            result.invalid_reason = "projected_stencil_residual_or_derivative_invalid";
        } else {
            result.invalid_reason = "projected_stencil_failed";
        }
        return result;
    }
    const double inv_2h = 1.0 / (2.0 * hessian_step);
    const double H_aa = (stencils[0].branch.d_encounter_global_angle_d_nu_A - stencils[1].branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H_bb = (stencils[2].branch.d_encounter_global_angle_d_nu_B - stencils[3].branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H_cc = (stencils[4].branch.d_encounter_global_angle_d_theta_A - stencils[5].branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H12_from_g1 = (stencils[2].branch.d_encounter_global_angle_d_nu_A - stencils[3].branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 = (stencils[0].branch.d_encounter_global_angle_d_nu_B - stencils[1].branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 = (stencils[4].branch.d_encounter_global_angle_d_nu_A - stencils[5].branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 = (stencils[0].branch.d_encounter_global_angle_d_theta_A - stencils[1].branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H23_from_g2 = (stencils[4].branch.d_encounter_global_angle_d_nu_B - stencils[5].branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H32_from_g3 = (stencils[2].branch.d_encounter_global_angle_d_theta_A - stencils[3].branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H_ab = 0.5 * (H12_from_g1 + H21_from_g2);
    const double H_ac = 0.5 * (H13_from_g1 + H31_from_g3);
    const double H_bc = 0.5 * (H23_from_g2 + H32_from_g3);
    if (!(std::isfinite(H_aa) && std::isfinite(H_bb) && std::isfinite(H_cc) && std::isfinite(H_ab) && std::isfinite(H_ac) && std::isfinite(H_bc))) {
        result.invalid_reason = "hessian_non_finite";
        return result;
    }
    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A);
    const double alpha_linear = normalize_angle_0_2pi(
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated.d_encounter_global_angle_d_theta_A * dx3);
    const double alpha_quadratic = normalize_angle_0_2pi(
        alpha_linear + 0.5 * (H_aa * dx1 * dx1 + H_bb * dx2 * dx2 + H_cc * dx3 * dx3 +
                              2.0 * H_ab * dx1 * dx2 + 2.0 * H_ac * dx1 * dx3 + 2.0 * H_bc * dx2 * dx3));
    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        source_branch.transfer_revolution, alpha_quadratic, source_branch, table.config().max_target_revolution);
    if (q_selection.selection_failed) {
        result.invalid_reason = "q_sheet_selection_failed";
        return result;
    }
    const auto residual_before = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        alpha_quadratic, source_branch.transfer_revolution, q_selection.selected_q);
    if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
        result.invalid_reason = "residual_before_correction_invalid";
        return result;
    }
    const auto derivatives_before = problem1::evaluate_problem1_root_residual_derivatives(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        alpha_quadratic, source_branch.transfer_revolution, q_selection.selected_q);
    if (!derivatives_before.valid || !std::isfinite(derivatives_before.R_alpha)) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }
    if (std::abs(derivatives_before.R_alpha) <= 1e-12) {
        result.invalid_reason = "correction_derivative_too_small";
        return result;
    }
    const double residual_scale_free = problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    const double delta_alpha = -residual_scale_free / derivatives_before.R_alpha;
    const double alpha_corrected = normalize_angle_0_2pi(alpha_quadratic + delta_alpha);
    const auto residual_after = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet, table.config().target_planet, query_nu_A, query_nu_B, query_theta_A,
        alpha_corrected, source_branch.transfer_revolution, q_selection.selected_q);
    if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
        result.invalid_reason = "residual_after_correction_invalid";
        return result;
    }
    result.residual_after_correction = residual_after.residual_seconds;
    if (std::abs(residual_after.residual_seconds) > kEngineeringResidualToleranceSeconds) {
        result.invalid_reason = "residual_after_correction_too_large";
        return result;
    }
    result.valid = true;
    result.time_diff_vs_route_a = std::abs(residual_after.transfer_time_seconds - route_a_branch.time_of_flight_seconds);
    result.alpha_diff_vs_route_a = wrapped_alpha_distance(residual_after.encounter_global_angle, route_a_branch.encounter_global_angle);
    return result;
}

AdaptiveEvalResult evaluate_adaptive_old(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    const spaceship_cpp::problem1::Problem1SolutionBranch& route_a_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    AdaptiveEvalResult result{};
    for (const double step : kAdaptiveHessianSteps) {
        const auto step_result = evaluate_old_step(table, nearest, source_branch, route_a_branch, query_nu_A, query_nu_B, query_theta_A, step);
        if (step_result.valid) {
            result.valid = true;
            result.selected_step = step;
            result.residual_after_correction = step_result.residual_after_correction;
            result.time_diff_vs_route_a = step_result.time_diff_vs_route_a;
            result.alpha_diff_vs_route_a = step_result.alpha_diff_vs_route_a;
            return result;
        }
        result.invalid_reason = step_result.invalid_reason;
    }
    return result;
}

AdaptiveEvalResult evaluate_adaptive_projected(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    const spaceship_cpp::problem1::Problem1SolutionBranch& route_a_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    AdaptiveEvalResult result{};
    for (const double step : kAdaptiveHessianSteps) {
        const auto step_result = evaluate_projected_step(table, nearest, source_branch, route_a_branch, query_nu_A, query_nu_B, query_theta_A, step);
        if (step_result.valid) {
            result.valid = true;
            result.selected_step = step;
            result.residual_after_correction = step_result.residual_after_correction;
            result.time_diff_vs_route_a = step_result.time_diff_vs_route_a;
            result.alpha_diff_vs_route_a = step_result.alpha_diff_vs_route_a;
            result.max_projected_stencil_delta_alpha = step_result.max_projected_stencil_delta_alpha;
            result.max_projected_stencil_residual_before = step_result.max_projected_stencil_residual_before;
            result.max_projected_stencil_residual_after = step_result.max_projected_stencil_residual_after;
            return result;
        }
        result.invalid_reason = step_result.invalid_reason;
        result.max_projected_stencil_delta_alpha = step_result.max_projected_stencil_delta_alpha;
        result.max_projected_stencil_residual_before = step_result.max_projected_stencil_residual_before;
        result.max_projected_stencil_residual_after = step_result.max_projected_stencil_residual_after;
    }
    return result;
}

void print_summary(const std::string& group, const GroupSummary& s) {
    std::cout << "ProjectedStencilCompareSummary\n";
    std::cout << "group=" << group << '\n';
    std::cout << "route_a_valid_count=" << s.route_a_valid_count << '\n';
    std::cout << "old_route_b_valid_count=" << s.old_route_b_valid_count << '\n';
    std::cout << "old_route_b_valid_ratio="
              << (s.route_a_valid_count > 0 ? static_cast<double>(s.old_route_b_valid_count) / s.route_a_valid_count : 0.0)
              << '\n';
    std::cout << "old_missing_count=" << (s.route_a_valid_count - s.old_route_b_valid_count) << '\n';
    for (const auto& [reason, count] : s.old_invalid_reasons) {
        std::cout << "old_invalid_reason[" << reason << "]=" << count << '\n';
    }
    std::cout << "projected_route_b_valid_count=" << s.projected_route_b_valid_count << '\n';
    std::cout << "projected_route_b_valid_ratio="
              << (s.route_a_valid_count > 0 ? static_cast<double>(s.projected_route_b_valid_count) / s.route_a_valid_count : 0.0)
              << '\n';
    std::cout << "projected_missing_count=" << (s.route_a_valid_count - s.projected_route_b_valid_count) << '\n';
    for (const auto& [reason, count] : s.projected_invalid_reasons) {
        std::cout << "projected_invalid_reason[" << reason << "]=" << count << '\n';
    }
    std::cout << "newly_recovered_count=" << s.newly_recovered_count << '\n';
    std::cout << "newly_recovered_ratio_over_old_missing="
              << ((s.route_a_valid_count - s.old_route_b_valid_count) > 0
                      ? static_cast<double>(s.newly_recovered_count) / (s.route_a_valid_count - s.old_route_b_valid_count)
                      : 0.0)
              << '\n';
    std::cout << "projected_lost_count=" << s.projected_lost_count << '\n';
    std::cout << "projected_same_root_success_count=" << s.projected_same_root_success_count << '\n';
    std::cout << "projected_same_root_success_ratio="
              << (s.projected_same_root_valid_count > 0
                      ? static_cast<double>(s.projected_same_root_success_count) / s.projected_same_root_valid_count
                      : 0.0)
              << '\n';
    std::cout << "mean_time_diff_vs_route_a="
              << (s.projected_same_root_valid_count > 0 ? s.sum_time_diff_vs_route_a / s.projected_same_root_valid_count : 0.0)
              << '\n';
    std::cout << "max_time_diff_vs_route_a=" << s.max_time_diff_vs_route_a << '\n';
    std::cout << "mean_alpha_diff_vs_route_a="
              << (s.projected_same_root_valid_count > 0 ? s.sum_alpha_diff_vs_route_a / s.projected_same_root_valid_count : 0.0)
              << '\n';
    std::cout << "max_alpha_diff_vs_route_a=" << s.max_alpha_diff_vs_route_a << '\n';
    std::cout << "mean_abs_residual_seconds="
              << (s.projected_same_root_valid_count > 0 ? s.sum_abs_residual_seconds / s.projected_same_root_valid_count : 0.0)
              << '\n';
    std::cout << "max_abs_residual_seconds=" << s.max_abs_residual_seconds << '\n';
    for (const auto& key : std::array<std::string, 4>{"1e-5", "5e-6", "2e-6", "1e-6"}) {
        const auto it = s.projected_selected_step_distribution.find(key);
        std::cout << "projected_selected_step[" << key << "]=" << (it != s.projected_selected_step_distribution.end() ? it->second : 0) << '\n';
    }
    const int recovered_to_print = std::min<int>(5, s.recovered_cases.size());
    for (int i = 0; i < recovered_to_print; ++i) {
        const auto& c = s.recovered_cases[static_cast<std::size_t>(i)];
        std::cout << "NewlyRecoveredProjectedStencilCase\n";
        std::cout << "group=" << c.group << '\n';
        std::cout << "sample_index=" << c.sample_index << '\n';
        std::cout << "source_index=" << c.source_index << '\n';
        std::cout << "source_rank_in_k=" << c.source_rank_in_k << '\n';
        std::cout << "k=" << c.k << '\n';
        std::cout << "source_q=" << c.source_q << '\n';
        std::cout << "selected_step=" << c.selected_step << '\n';
        std::cout << "old_invalid_reason=" << c.old_invalid_reason << '\n';
        std::cout << "projected_residual_after=" << c.projected_residual_after << '\n';
        std::cout << "projected_time_diff_vs_route_a=" << c.projected_time_diff_vs_route_a << '\n';
        std::cout << "projected_alpha_diff_vs_route_a=" << c.projected_alpha_diff_vs_route_a << '\n';
        std::cout << "max_projected_stencil_delta_alpha=" << c.max_projected_stencil_delta_alpha << '\n';
        std::cout << "max_projected_stencil_residual_before=" << c.max_projected_stencil_residual_before << '\n';
        std::cout << "max_projected_stencil_residual_after=" << c.max_projected_stencil_residual_after << '\n';
    }
    const int suspicious_to_print = std::min<int>(5, s.suspicious_cases.size());
    for (int i = 0; i < suspicious_to_print; ++i) {
        const auto& c = s.suspicious_cases[static_cast<std::size_t>(i)];
        std::cout << "SuspiciousProjectedStencilCase\n";
        std::cout << "group=" << c.group << '\n';
        std::cout << "sample_index=" << c.sample_index << '\n';
        std::cout << "source_index=" << c.source_index << '\n';
        std::cout << "source_rank_in_k=" << c.source_rank_in_k << '\n';
        std::cout << "k=" << c.k << '\n';
        std::cout << "selected_step=" << c.selected_step << '\n';
        std::cout << "residual_after=" << c.residual_after << '\n';
        std::cout << "time_diff_vs_route_a=" << c.time_diff_vs_route_a << '\n';
        std::cout << "alpha_diff_vs_route_a=" << c.alpha_diff_vs_route_a << '\n';
    }
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
    const auto samples = build_samples(kSamplesPerGroup);
    std::map<std::string, GroupSummary> summaries;

    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const auto& sample = samples[sample_index];
        auto table = build_local_cell_table_for_query(
            departure_planet, target_planet, sample, max_transfer_revolution, max_target_revolution);
        const auto nearest = problem1::find_nearest_problem1_root_table_node(
            table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
        if (!nearest.valid || nearest.cell == nullptr) {
            continue;
        }
        std::map<int, int> source_rank_by_k;
        for (std::size_t source_index = 0; source_index < nearest.cell->solutions_sorted_by_time_of_flight.size(); ++source_index) {
            const auto& source_branch = nearest.cell->solutions_sorted_by_time_of_flight[source_index];
            if (!source_branch.valid) {
                continue;
            }
            const int source_rank_in_k = source_rank_by_k[source_branch.transfer_revolution]++;
            auto& summary = summaries[sample.group_name];
            const auto route_a = evaluate_route_a_from_source_branch_test_only(
                table, nearest, source_branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            if (!route_a.valid) {
                continue;
            }
            summary.route_a_valid_count += 1;
            const auto old_result = evaluate_adaptive_old(
                table, nearest, source_branch, route_a.branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            const auto projected_result = evaluate_adaptive_projected(
                table, nearest, source_branch, route_a.branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);

            if (old_result.valid) {
                summary.old_route_b_valid_count += 1;
            } else {
                summary.old_invalid_reasons[old_result.invalid_reason.empty() ? "other" : old_result.invalid_reason] += 1;
            }
            if (projected_result.valid) {
                summary.projected_route_b_valid_count += 1;
                summary.projected_same_root_valid_count += 1;
                summary.projected_selected_step_distribution[step_label(projected_result.selected_step)] += 1;
                const double time_diff = projected_result.time_diff_vs_route_a;
                const double alpha_diff = projected_result.alpha_diff_vs_route_a;
                const double abs_residual = std::abs(projected_result.residual_after_correction);
                summary.sum_time_diff_vs_route_a += time_diff;
                summary.max_time_diff_vs_route_a = std::max(summary.max_time_diff_vs_route_a, time_diff);
                summary.sum_alpha_diff_vs_route_a += alpha_diff;
                summary.max_alpha_diff_vs_route_a = std::max(summary.max_alpha_diff_vs_route_a, alpha_diff);
                summary.sum_abs_residual_seconds += abs_residual;
                summary.max_abs_residual_seconds = std::max(summary.max_abs_residual_seconds, abs_residual);
                const bool same_root = time_diff <= kSameRootTimeThresholdSeconds && alpha_diff <= kSameRootAlphaThreshold;
                if (same_root) {
                    summary.projected_same_root_success_count += 1;
                }
                if (abs_residual <= kEngineeringResidualToleranceSeconds &&
                    (time_diff > kSameRootTimeThresholdSeconds || alpha_diff > kSameRootAlphaThreshold)) {
                    summary.suspicious_cases.push_back(SuspiciousCase{
                        sample.group_name,
                        static_cast<int>(sample_index),
                        static_cast<int>(source_index),
                        source_rank_in_k,
                        source_branch.transfer_revolution,
                        projected_result.selected_step,
                        projected_result.residual_after_correction,
                        time_diff,
                        alpha_diff,
                    });
                }
            } else {
                summary.projected_invalid_reasons[projected_result.invalid_reason.empty() ? "other" : projected_result.invalid_reason] += 1;
            }
            if (!old_result.valid && projected_result.valid) {
                summary.newly_recovered_count += 1;
                summary.recovered_cases.push_back(RecoveredCase{
                    sample.group_name,
                    static_cast<int>(sample_index),
                    static_cast<int>(source_index),
                    source_rank_in_k,
                    source_branch.transfer_revolution,
                    source_branch.target_revolution,
                    projected_result.selected_step,
                    old_result.invalid_reason,
                    projected_result.residual_after_correction,
                    projected_result.time_diff_vs_route_a,
                    projected_result.alpha_diff_vs_route_a,
                    projected_result.max_projected_stencil_delta_alpha,
                    projected_result.max_projected_stencil_residual_before,
                    projected_result.max_projected_stencil_residual_after,
                });
            }
            if (old_result.valid && !projected_result.valid) {
                summary.projected_lost_count += 1;
            }
        }
    }

    for (const auto& [group, summary] : summaries) {
        print_summary(group, summary);
    }
    return 0;
}
