#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
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
constexpr double kResidualToleranceSeconds = 1e-2;
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
    double linear_seed_ms = 0.0;
    double q_sheet_selection_ms = 0.0;
    double newton_refinement_ms = 0.0;
    double total_ms = 0.0;
};

struct PrecomputedProjectedData {
    bool valid = false;
    std::string invalid_reason;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double source_alpha = 0.0;
    double d_alpha_d_nu_A = 0.0;
    double d_alpha_d_nu_B = 0.0;
    double d_alpha_d_theta_A = 0.0;
    spaceship_cpp::problem1::Problem1RootHessian hessian;
    double precompute_hessian_ms = 0.0;
};

struct RouteBProjectedEvalResult {
    bool valid = false;
    std::string invalid_reason;
    double selected_step = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
    double projected_hessian_total_ms = 0.0;
    double quadratic_prediction_eval_ms = 0.0;
    double q_sheet_selection_ms = 0.0;
    double residual_before_correction_ms = 0.0;
    double derivative_before_correction_ms = 0.0;
    double correction_and_residual_after_ms = 0.0;
    double total_route_b_projected_ms = 0.0;
};

struct FallbackCase {
    std::string group;
    int sample_index = -1;
    int source_index = -1;
    int source_rank_in_k = -1;
    int k = 0;
    int source_q = 0;
    std::string fallback_reason;
};

struct GroupSummary {
    int route_a_valid_count = 0;
    int route_b_valid_count = 0;
    int precomputed_route_b_valid_count = 0;
    int precomputed_same_root_success_count = 0;
    double route_a_total_ms = 0.0;
    double route_b_total_ms = 0.0;
    double projected_hessian_total_ms = 0.0;
    double quadratic_prediction_eval_ms = 0.0;
    double q_sheet_selection_ms = 0.0;
    double residual_before_correction_ms = 0.0;
    double derivative_before_correction_ms = 0.0;
    double correction_and_residual_after_ms = 0.0;
    double precompute_hessian_ms = 0.0;
    double precomputed_query_ms = 0.0;
    std::map<std::string, int> fallback_reasons;
    std::map<std::string, int> selected_step_distribution;
};

double wrapped_alpha_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

int lower_corner_index(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
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

RouteAEvalResult evaluate_route_a_from_source_branch_timed(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    RouteAEvalResult result{};
    const auto total_start = clock::now();
    if (!source_branch.valid) {
        result.invalid_reason = "source_branch_invalid";
        return result;
    }

    const auto seed_start = clock::now();
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
    const auto seed_end = clock::now();
    result.linear_seed_ms = std::chrono::duration<double, std::milli>(seed_end - seed_start).count();

    const auto q_start = clock::now();
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
    const auto q_end = clock::now();
    result.q_sheet_selection_ms = std::chrono::duration<double, std::milli>(q_end - q_start).count();
    if (q_selection.selection_failed) {
        result.invalid_reason = "q_sheet_selection_failed";
        return result;
    }

    const auto newton_start = clock::now();
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
    const auto newton_end = clock::now();
    result.newton_refinement_ms = std::chrono::duration<double, std::milli>(newton_end - newton_start).count();
    if (!refined.valid) {
        result.invalid_reason =
            refined.diagnostic.invalid_reason.empty() ? "route_a_refinement_failed" : refined.diagnostic.invalid_reason;
        return result;
    }
    result.valid = true;
    result.branch = refined.branch;
    const auto total_end = clock::now();
    result.total_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    return result;
}

PrecomputedProjectedData precompute_projected_hessian_test_only(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    double hessian_step
) {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    PrecomputedProjectedData result{};
    const auto start = clock::now();
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
        result.invalid_reason = "node_branch_derivatives_unavailable";
        return result;
    }
    const auto hessian = problem1::estimate_problem1_root_hessian_projected_tangent_finite_difference(
        table.config().departure_planet,
        table.config().target_planet,
        nearest.node_nu_A,
        nearest.node_nu_B,
        nearest.node_theta_A,
        differentiated,
        hessian_step,
        1e-2,
        1e-12);
    const auto end = clock::now();
    result.precompute_hessian_ms = std::chrono::duration<double, std::milli>(end - start).count();
    if (!hessian.valid) {
        result.invalid_reason = hessian.invalid_reason.empty() ? "projected_hessian_invalid" : hessian.invalid_reason;
        return result;
    }
    result.valid = true;
    result.transfer_revolution = differentiated.transfer_revolution;
    result.target_revolution = differentiated.target_revolution;
    result.source_alpha = differentiated.encounter_global_angle;
    result.d_alpha_d_nu_A = differentiated.d_encounter_global_angle_d_nu_A;
    result.d_alpha_d_nu_B = differentiated.d_encounter_global_angle_d_nu_B;
    result.d_alpha_d_theta_A = differentiated.d_encounter_global_angle_d_theta_A;
    result.hessian = hessian;
    return result;
}

RouteBProjectedEvalResult evaluate_route_b_projected_with_timing(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    RouteBProjectedEvalResult result{};
    const auto total_start = clock::now();
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
        result.invalid_reason = "node_branch_derivatives_unavailable";
        return result;
    }

    for (const double step : kAdaptiveHessianSteps) {
        const auto hessian_start = clock::now();
        const auto hessian = problem1::estimate_problem1_root_hessian_projected_tangent_finite_difference(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            differentiated,
            step,
            1e-2,
            1e-12);
        const auto hessian_end = clock::now();
        result.projected_hessian_total_ms +=
            std::chrono::duration<double, std::milli>(hessian_end - hessian_start).count();
        if (!hessian.valid) {
            result.invalid_reason = hessian.invalid_reason.empty() ? "projected_hessian_invalid" : hessian.invalid_reason;
            continue;
        }

        const auto quad_start = clock::now();
        const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A);
        const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B);
        const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A);
        const double alpha_linear = normalize_angle_0_2pi(
            differentiated.encounter_global_angle +
            differentiated.d_encounter_global_angle_d_nu_A * dx1 +
            differentiated.d_encounter_global_angle_d_nu_B * dx2 +
            differentiated.d_encounter_global_angle_d_theta_A * dx3);
        const double alpha_quadratic = normalize_angle_0_2pi(
            alpha_linear + 0.5 * (
                hessian.H_nu_A_nu_A * dx1 * dx1 +
                hessian.H_nu_B_nu_B * dx2 * dx2 +
                hessian.H_theta_A_theta_A * dx3 * dx3 +
                2.0 * hessian.H_nu_A_nu_B * dx1 * dx2 +
                2.0 * hessian.H_nu_A_theta_A * dx1 * dx3 +
                2.0 * hessian.H_nu_B_theta_A * dx2 * dx3));
        const auto quad_end = clock::now();
        result.quadratic_prediction_eval_ms +=
            std::chrono::duration<double, std::milli>(quad_end - quad_start).count();

        const auto q_start = clock::now();
        const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            source_branch.transfer_revolution,
            alpha_quadratic,
            source_branch,
            table.config().max_target_revolution);
        const auto q_end = clock::now();
        result.q_sheet_selection_ms += std::chrono::duration<double, std::milli>(q_end - q_start).count();
        if (q_selection.selection_failed) {
            result.invalid_reason = "q_sheet_selection_failed";
            continue;
        }

        const auto residual_before_start = clock::now();
        const auto residual_before = problem1::evaluate_problem1_root_residual(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_quadratic,
            source_branch.transfer_revolution,
            q_selection.selected_q);
        const auto residual_before_end = clock::now();
        result.residual_before_correction_ms +=
            std::chrono::duration<double, std::milli>(residual_before_end - residual_before_start).count();
        if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
            result.invalid_reason = "residual_before_correction_invalid";
            continue;
        }

        const auto derivative_start = clock::now();
        const auto derivative_before = problem1::evaluate_problem1_root_residual_derivatives(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_quadratic,
            source_branch.transfer_revolution,
            q_selection.selected_q);
        const auto derivative_end = clock::now();
        result.derivative_before_correction_ms +=
            std::chrono::duration<double, std::milli>(derivative_end - derivative_start).count();
        if (!derivative_before.valid || !std::isfinite(derivative_before.R_alpha)) {
            result.invalid_reason = "derivative_before_correction_invalid";
            continue;
        }
        if (std::abs(derivative_before.R_alpha) <= 1e-12) {
            result.invalid_reason = "correction_derivative_too_small";
            continue;
        }

        const auto correction_start = clock::now();
        const double residual_scale_free =
            problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
        const double delta_alpha = -residual_scale_free / derivative_before.R_alpha;
        const double alpha_corrected = normalize_angle_0_2pi(alpha_quadratic + delta_alpha);
        const auto residual_after = problem1::evaluate_problem1_root_residual(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_corrected,
            source_branch.transfer_revolution,
            q_selection.selected_q);
        const auto correction_end = clock::now();
        result.correction_and_residual_after_ms +=
            std::chrono::duration<double, std::milli>(correction_end - correction_start).count();
        if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
            result.invalid_reason = "residual_after_correction_invalid";
            continue;
        }
        if (std::abs(residual_after.residual_seconds) > kResidualToleranceSeconds) {
            result.invalid_reason = "residual_after_correction_too_large";
            continue;
        }

        result.valid = true;
        result.selected_step = step;
        result.branch.valid = true;
        result.branch.encounter_global_angle = residual_after.encounter_global_angle;
        result.branch.target_arrival_true_anomaly = residual_after.target_arrival_true_anomaly;
        result.branch.transfer_revolution = source_branch.transfer_revolution;
        result.branch.target_revolution = q_selection.selected_q;
        result.branch.time_of_flight_seconds = residual_after.transfer_time_seconds;
        result.branch.target_time_seconds = residual_after.target_time_seconds;
        result.branch.residual_seconds = residual_after.residual_seconds;
        result.branch.transfer_e = residual_after.transfer_e;
        result.branch.transfer_p = residual_after.transfer_p;
        result.branch.transfer_a = residual_after.transfer_a;
        result.branch.theta_B = residual_after.theta_B;
        break;
    }

    const auto total_end = clock::now();
    result.total_route_b_projected_ms =
        std::chrono::duration<double, std::milli>(total_end - total_start).count();
    return result;
}

RouteBProjectedEvalResult evaluate_route_b_from_precomputed_hessian(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    const PrecomputedProjectedData& precomputed,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;

    RouteBProjectedEvalResult result{};
    const auto total_start = clock::now();
    if (!precomputed.valid) {
        result.invalid_reason = precomputed.invalid_reason;
        return result;
    }
    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A);

    const auto quad_start = clock::now();
    const double alpha_linear = normalize_angle_0_2pi(
        precomputed.source_alpha +
        precomputed.d_alpha_d_nu_A * dx1 +
        precomputed.d_alpha_d_nu_B * dx2 +
        precomputed.d_alpha_d_theta_A * dx3);
    const double alpha_quadratic = normalize_angle_0_2pi(
        alpha_linear + 0.5 * (
            precomputed.hessian.H_nu_A_nu_A * dx1 * dx1 +
            precomputed.hessian.H_nu_B_nu_B * dx2 * dx2 +
            precomputed.hessian.H_theta_A_theta_A * dx3 * dx3 +
            2.0 * precomputed.hessian.H_nu_A_nu_B * dx1 * dx2 +
            2.0 * precomputed.hessian.H_nu_A_theta_A * dx1 * dx3 +
            2.0 * precomputed.hessian.H_nu_B_theta_A * dx2 * dx3));
    const auto quad_end = clock::now();
    result.quadratic_prediction_eval_ms =
        std::chrono::duration<double, std::milli>(quad_end - quad_start).count();

    const auto q_start = clock::now();
    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        source_branch.transfer_revolution,
        alpha_quadratic,
        source_branch,
        table.config().max_target_revolution);
    const auto q_end = clock::now();
    result.q_sheet_selection_ms = std::chrono::duration<double, std::milli>(q_end - q_start).count();
    if (q_selection.selection_failed) {
        result.invalid_reason = "q_sheet_selection_failed";
        return result;
    }

    const auto residual_before_start = clock::now();
    const auto residual_before = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        source_branch.transfer_revolution,
        q_selection.selected_q);
    const auto residual_before_end = clock::now();
    result.residual_before_correction_ms =
        std::chrono::duration<double, std::milli>(residual_before_end - residual_before_start).count();
    if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
        result.invalid_reason = "residual_before_correction_invalid";
        return result;
    }

    const auto derivative_start = clock::now();
    const auto derivative_before = problem1::evaluate_problem1_root_residual_derivatives(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        source_branch.transfer_revolution,
        q_selection.selected_q);
    const auto derivative_end = clock::now();
    result.derivative_before_correction_ms =
        std::chrono::duration<double, std::milli>(derivative_end - derivative_start).count();
    if (!derivative_before.valid || !std::isfinite(derivative_before.R_alpha)) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }
    if (std::abs(derivative_before.R_alpha) <= 1e-12) {
        result.invalid_reason = "correction_derivative_too_small";
        return result;
    }

    const auto correction_start = clock::now();
    const double residual_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    const double delta_alpha = -residual_scale_free / derivative_before.R_alpha;
    const double alpha_corrected = normalize_angle_0_2pi(alpha_quadratic + delta_alpha);
    const auto residual_after = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_corrected,
        source_branch.transfer_revolution,
        q_selection.selected_q);
    const auto correction_end = clock::now();
    result.correction_and_residual_after_ms =
        std::chrono::duration<double, std::milli>(correction_end - correction_start).count();
    if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
        result.invalid_reason = "residual_after_correction_invalid";
        return result;
    }
    if (std::abs(residual_after.residual_seconds) > kResidualToleranceSeconds) {
        result.invalid_reason = "residual_after_correction_too_large";
        return result;
    }

    result.valid = true;
    result.selected_step = 1e-5;
    result.branch.valid = true;
    result.branch.encounter_global_angle = residual_after.encounter_global_angle;
    result.branch.target_arrival_true_anomaly = residual_after.target_arrival_true_anomaly;
    result.branch.transfer_revolution = source_branch.transfer_revolution;
    result.branch.target_revolution = q_selection.selected_q;
    result.branch.time_of_flight_seconds = residual_after.transfer_time_seconds;
    result.branch.target_time_seconds = residual_after.target_time_seconds;
    result.branch.residual_seconds = residual_after.residual_seconds;
    result.branch.transfer_e = residual_after.transfer_e;
    result.branch.transfer_p = residual_after.transfer_p;
    result.branch.transfer_a = residual_after.transfer_a;
    result.branch.theta_B = residual_after.theta_B;
    const auto total_end = clock::now();
    result.total_route_b_projected_ms =
        std::chrono::duration<double, std::milli>(total_end - total_start).count();
    return result;
}

void print_summary(const std::string& group, const GroupSummary& s) {
    const double avg_route_a_total_ms =
        s.route_a_valid_count > 0 ? s.route_a_total_ms / s.route_a_valid_count : 0.0;
    const double avg_route_b_total_ms =
        s.route_a_valid_count > 0 ? s.route_b_total_ms / s.route_a_valid_count : 0.0;
    const double avg_projected_hessian_ms =
        s.route_a_valid_count > 0 ? s.projected_hessian_total_ms / s.route_a_valid_count : 0.0;
    const double avg_quadratic_prediction_eval_ms =
        s.route_a_valid_count > 0 ? s.quadratic_prediction_eval_ms / s.route_a_valid_count : 0.0;
    const double avg_q_sheet_selection_ms =
        s.route_a_valid_count > 0 ? s.q_sheet_selection_ms / s.route_a_valid_count : 0.0;
    const double avg_residual_before_correction_ms =
        s.route_a_valid_count > 0 ? s.residual_before_correction_ms / s.route_a_valid_count : 0.0;
    const double avg_derivative_before_correction_ms =
        s.route_a_valid_count > 0 ? s.derivative_before_correction_ms / s.route_a_valid_count : 0.0;
    const double avg_correction_and_residual_after_ms =
        s.route_a_valid_count > 0 ? s.correction_and_residual_after_ms / s.route_a_valid_count : 0.0;
    const double fraction_projected_hessian =
        s.route_b_total_ms > 0.0 ? s.projected_hessian_total_ms / s.route_b_total_ms : 0.0;
    const double fraction_derivative_before_correction =
        s.route_b_total_ms > 0.0 ? s.derivative_before_correction_ms / s.route_b_total_ms : 0.0;
    const double fraction_residual_eval =
        s.route_b_total_ms > 0.0
            ? (s.residual_before_correction_ms + s.correction_and_residual_after_ms) / s.route_b_total_ms
            : 0.0;
    const double avg_precompute_hessian_ms_per_branch =
        s.route_a_valid_count > 0 ? s.precompute_hessian_ms / s.route_a_valid_count : 0.0;
    const double avg_query_time_with_precomputed_hessian_ms =
        s.route_a_valid_count > 0 ? s.precomputed_query_ms / s.route_a_valid_count : 0.0;
    const double query_speedup_vs_route_a =
        s.precomputed_query_ms > 0.0 ? s.route_a_total_ms / s.precomputed_query_ms : 0.0;
    const double precompute_plus_query_speedup_vs_route_a =
        (s.precompute_hessian_ms + s.precomputed_query_ms) > 0.0
            ? s.route_a_total_ms / (s.precompute_hessian_ms + s.precomputed_query_ms)
            : 0.0;
    const double same_root_success_ratio =
        s.precomputed_route_b_valid_count > 0
            ? static_cast<double>(s.precomputed_same_root_success_count) / s.precomputed_route_b_valid_count
            : 0.0;

    std::cout << "TimingBreakdownSummary\n";
    std::cout << "group=" << group << '\n';
    std::cout << "route_a_valid_count=" << s.route_a_valid_count << '\n';
    std::cout << "route_b_valid_count=" << s.route_b_valid_count << '\n';
    std::cout << "avg_route_a_total_ms=" << avg_route_a_total_ms << '\n';
    std::cout << "avg_route_b_total_ms=" << avg_route_b_total_ms << '\n';
    std::cout << "avg_projected_hessian_ms=" << avg_projected_hessian_ms << '\n';
    std::cout << "avg_quadratic_prediction_eval_ms=" << avg_quadratic_prediction_eval_ms << '\n';
    std::cout << "avg_q_sheet_selection_ms=" << avg_q_sheet_selection_ms << '\n';
    std::cout << "avg_residual_before_correction_ms=" << avg_residual_before_correction_ms << '\n';
    std::cout << "avg_derivative_before_correction_ms=" << avg_derivative_before_correction_ms << '\n';
    std::cout << "avg_correction_and_residual_after_ms=" << avg_correction_and_residual_after_ms << '\n';
    std::cout << "fraction_projected_hessian=" << fraction_projected_hessian << '\n';
    std::cout << "fraction_derivative_before_correction=" << fraction_derivative_before_correction << '\n';
    std::cout << "fraction_residual_eval=" << fraction_residual_eval << '\n';

    std::cout << "PrecomputedHessianQuerySummary\n";
    std::cout << "group=" << group << '\n';
    std::cout << "route_a_valid_count=" << s.route_a_valid_count << '\n';
    std::cout << "precomputed_route_b_valid_count=" << s.precomputed_route_b_valid_count << '\n';
    std::cout << "same_root_success_ratio=" << same_root_success_ratio << '\n';
    std::cout << "avg_precompute_hessian_ms_per_branch=" << avg_precompute_hessian_ms_per_branch << '\n';
    std::cout << "avg_query_time_with_precomputed_hessian_ms=" << avg_query_time_with_precomputed_hessian_ms << '\n';
    std::cout << "avg_pure_route_a_ms=" << avg_route_a_total_ms << '\n';
    std::cout << "query_speedup_vs_route_a=" << query_speedup_vs_route_a << '\n';
    std::cout << "precompute_plus_query_speedup_vs_route_a=" << precompute_plus_query_speedup_vs_route_a << '\n';
}

}  // namespace

int main() {
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
        const auto nearest = spaceship_cpp::problem1::find_nearest_problem1_root_table_node(
            table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
        if (!nearest.valid || nearest.cell == nullptr) {
            continue;
        }

        for (const auto& source_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
            if (!source_branch.valid) {
                continue;
            }
            auto& summary = summaries[sample.group_name];

            const auto route_a = evaluate_route_a_from_source_branch_timed(
                table, nearest, source_branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            if (!route_a.valid) {
                continue;
            }
            summary.route_a_valid_count += 1;
            summary.route_a_total_ms += route_a.total_ms;

            const auto route_b = evaluate_route_b_projected_with_timing(
                table, nearest, source_branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            summary.route_b_total_ms += route_b.total_route_b_projected_ms;
            summary.projected_hessian_total_ms += route_b.projected_hessian_total_ms;
            summary.quadratic_prediction_eval_ms += route_b.quadratic_prediction_eval_ms;
            summary.q_sheet_selection_ms += route_b.q_sheet_selection_ms;
            summary.residual_before_correction_ms += route_b.residual_before_correction_ms;
            summary.derivative_before_correction_ms += route_b.derivative_before_correction_ms;
            summary.correction_and_residual_after_ms += route_b.correction_and_residual_after_ms;
            if (route_b.valid) {
                summary.route_b_valid_count += 1;
                summary.selected_step_distribution["1e-5"] += (std::abs(route_b.selected_step - 1e-5) < 1e-12) ? 1 : 0;
            } else {
                summary.fallback_reasons[route_b.invalid_reason] += 1;
            }

            const auto precomputed = precompute_projected_hessian_test_only(
                table, nearest, source_branch, 1e-5);
            summary.precompute_hessian_ms += precomputed.precompute_hessian_ms;
            const auto precomputed_query = evaluate_route_b_from_precomputed_hessian(
                table, nearest, source_branch, precomputed, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            summary.precomputed_query_ms += precomputed_query.total_route_b_projected_ms;
            if (precomputed_query.valid) {
                summary.precomputed_route_b_valid_count += 1;
                const double time_diff =
                    std::abs(precomputed_query.branch.time_of_flight_seconds - route_a.branch.time_of_flight_seconds);
                const double alpha_diff =
                    wrapped_alpha_distance(precomputed_query.branch.encounter_global_angle, route_a.branch.encounter_global_angle);
                if (time_diff <= kSameRootTimeThresholdSeconds &&
                    alpha_diff <= kSameRootAlphaThreshold &&
                    std::abs(precomputed_query.branch.residual_seconds) <= kResidualToleranceSeconds) {
                    summary.precomputed_same_root_success_count += 1;
                }
            }
        }
    }

    for (const auto& [group, summary] : summaries) {
        print_summary(group, summary);
    }

    std::cout << "route_b_projected_timing_breakdown_ok\n";
    return EXIT_SUCCESS;
}
