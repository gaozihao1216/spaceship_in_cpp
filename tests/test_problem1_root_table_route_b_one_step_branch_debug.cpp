#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
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
constexpr double kQuadraticOneStepHessianStep = 1e-5;
constexpr double kHugeTimeErrorThresholdSeconds = 1e6;

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

struct StencilAttachedBranch {
    bool residual_gate_passed = false;
    bool attach_success = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct TestOnlyHessian {
    bool valid = false;
    std::string invalid_reason;
    int stencil_pass_count = 0;
    double H_aa = 0.0;
    double H_bb = 0.0;
    double H_cc = 0.0;
    double H_ab = 0.0;
    double H_ac = 0.0;
    double H_bc = 0.0;
};

struct CandidateDebug {
    int source_index = -1;
    int source_rank_in_k = -1;
    int k = 0;
    int source_q = 0;
    int selected_q = 0;
    double hessian_step = kQuadraticOneStepHessianStep;
    bool hessian_valid = false;
    std::string hessian_invalid_reason;
    int hessian_stencil_seconds_gate_pass_count = 0;
    double alpha_linear = std::numeric_limits<double>::quiet_NaN();
    double alpha_quadratic = std::numeric_limits<double>::quiet_NaN();
    double alpha_corrected = std::numeric_limits<double>::quiet_NaN();
    double residual_before_correction = std::numeric_limits<double>::quiet_NaN();
    double residual_after_correction = std::numeric_limits<double>::quiet_NaN();
    double transfer_time_seconds_after = std::numeric_limits<double>::quiet_NaN();
    double target_time_seconds_after = std::numeric_limits<double>::quiet_NaN();
    bool valid_candidate = false;
    std::string invalid_reason;
    spaceship_cpp::problem1::Problem1RootApproximationResult approximation;
};

struct SampleDebug {
    int sample_index = -1;
    QuerySample sample;
    spaceship_cpp::problem1::Problem1RootTable table;
    spaceship_cpp::problem1::Problem1RootNearestNode nearest;
    spaceship_cpp::problem1::Problem1RouteBSafeQueryResult safe_result;
    std::map<std::string, int> invalid_reason_distribution;
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> exact_roots;
    std::vector<CandidateDebug> candidate_debugs;
    double matching_a_max_time_error = 0.0;
};

struct MatchPair {
    int pair_rank = -1;
    int exact_rank = -1;
    int candidate_rank = -1;
    int k = 0;
    int exact_q = 0;
    int candidate_q = 0;
    double exact_time = std::numeric_limits<double>::quiet_NaN();
    double candidate_time = std::numeric_limits<double>::quiet_NaN();
    double time_error = std::numeric_limits<double>::quiet_NaN();
    double exact_alpha = std::numeric_limits<double>::quiet_NaN();
    double candidate_alpha = std::numeric_limits<double>::quiet_NaN();
    double alpha_error = std::numeric_limits<double>::quiet_NaN();
    double candidate_residual = std::numeric_limits<double>::quiet_NaN();
};

double wrapped_alpha_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

double branch_time_duplicate_tolerance_seconds(double t1, double t2) {
    return std::max(1e-4, 1e-10 * std::max(std::abs(t1), std::abs(t2)));
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

StencilAttachedBranch evaluate_stencil_branch_under_seconds_gate(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& differentiated_source_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - node_theta_A);
    const double alpha_tangent = normalize_angle_0_2pi(
        differentiated_source_branch.encounter_global_angle +
        differentiated_source_branch.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated_source_branch.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated_source_branch.d_encounter_global_angle_d_theta_A * dx3);
    const auto residual = problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_tangent,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    StencilAttachedBranch result{};
    result.residual_gate_passed =
        residual.valid &&
        std::isfinite(residual.residual_seconds) &&
        std::abs(residual.residual_seconds) <= kEngineeringResidualToleranceSeconds;
    if (!result.residual_gate_passed) {
        return result;
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
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        tangent_branch);
    result.attach_success = tangent_branch.valid && tangent_branch.derivatives_available;
    result.branch = tangent_branch;
    return result;
}

TestOnlyHessian build_test_only_hessian(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& differentiated_source_branch,
    double hessian_step
) {
    TestOnlyHessian result{};
    const auto nu_A_plus = evaluate_stencil_branch_under_seconds_gate(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated_source_branch,
        normalize_angle_0_2pi(node_nu_A + hessian_step),
        node_nu_B,
        node_theta_A);
    const auto nu_A_minus = evaluate_stencil_branch_under_seconds_gate(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated_source_branch,
        normalize_angle_0_2pi(node_nu_A - hessian_step),
        node_nu_B,
        node_theta_A);
    const auto nu_B_plus = evaluate_stencil_branch_under_seconds_gate(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated_source_branch,
        node_nu_A,
        normalize_angle_0_2pi(node_nu_B + hessian_step),
        node_theta_A);
    const auto nu_B_minus = evaluate_stencil_branch_under_seconds_gate(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated_source_branch,
        node_nu_A,
        normalize_angle_0_2pi(node_nu_B - hessian_step),
        node_theta_A);
    const auto theta_A_plus = evaluate_stencil_branch_under_seconds_gate(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated_source_branch,
        node_nu_A,
        node_nu_B,
        normalize_angle_0_2pi(node_theta_A + hessian_step));
    const auto theta_A_minus = evaluate_stencil_branch_under_seconds_gate(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated_source_branch,
        node_nu_A,
        node_nu_B,
        normalize_angle_0_2pi(node_theta_A - hessian_step));

    result.stencil_pass_count =
        static_cast<int>(nu_A_plus.residual_gate_passed) +
        static_cast<int>(nu_A_minus.residual_gate_passed) +
        static_cast<int>(nu_B_plus.residual_gate_passed) +
        static_cast<int>(nu_B_minus.residual_gate_passed) +
        static_cast<int>(theta_A_plus.residual_gate_passed) +
        static_cast<int>(theta_A_minus.residual_gate_passed);

    const bool all_gate_pass =
        nu_A_plus.residual_gate_passed && nu_A_minus.residual_gate_passed &&
        nu_B_plus.residual_gate_passed && nu_B_minus.residual_gate_passed &&
        theta_A_plus.residual_gate_passed && theta_A_minus.residual_gate_passed;
    if (!all_gate_pass) {
        result.invalid_reason = "hessian_stencil_seconds_gate_failed";
        return result;
    }
    const bool all_attach_success =
        nu_A_plus.attach_success && nu_A_minus.attach_success &&
        nu_B_plus.attach_success && nu_B_minus.attach_success &&
        theta_A_plus.attach_success && theta_A_minus.attach_success;
    if (!all_attach_success) {
        result.invalid_reason = "hessian_attach_failed";
        return result;
    }

    const double inv_2h = 1.0 / (2.0 * hessian_step);
    result.H_aa =
        (nu_A_plus.branch.d_encounter_global_angle_d_nu_A - nu_A_minus.branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    result.H_bb =
        (nu_B_plus.branch.d_encounter_global_angle_d_nu_B - nu_B_minus.branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    result.H_cc =
        (theta_A_plus.branch.d_encounter_global_angle_d_theta_A -
         theta_A_minus.branch.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H12_from_g1 =
        (nu_B_plus.branch.d_encounter_global_angle_d_nu_A - nu_B_minus.branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 =
        (nu_A_plus.branch.d_encounter_global_angle_d_nu_B - nu_A_minus.branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 =
        (theta_A_plus.branch.d_encounter_global_angle_d_nu_A - theta_A_minus.branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 =
        (nu_A_plus.branch.d_encounter_global_angle_d_theta_A - nu_A_minus.branch.d_encounter_global_angle_d_theta_A) *
        inv_2h;
    const double H23_from_g2 =
        (theta_A_plus.branch.d_encounter_global_angle_d_nu_B - theta_A_minus.branch.d_encounter_global_angle_d_nu_B) *
        inv_2h;
    const double H32_from_g3 =
        (nu_B_plus.branch.d_encounter_global_angle_d_theta_A - nu_B_minus.branch.d_encounter_global_angle_d_theta_A) *
        inv_2h;
    result.H_ab = 0.5 * (H12_from_g1 + H21_from_g2);
    result.H_ac = 0.5 * (H13_from_g1 + H31_from_g3);
    result.H_bc = 0.5 * (H23_from_g2 + H32_from_g3);
    result.valid =
        std::isfinite(result.H_aa) &&
        std::isfinite(result.H_bb) &&
        std::isfinite(result.H_cc) &&
        std::isfinite(result.H_ab) &&
        std::isfinite(result.H_ac) &&
        std::isfinite(result.H_bc);
    if (!result.valid) {
        result.invalid_reason = "hessian_non_finite";
    }
    return result;
}

CandidateDebug evaluate_candidate_debug(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& node_branch,
    int source_index,
    int source_rank_in_k
) {
    namespace problem1 = spaceship_cpp::problem1;
    CandidateDebug debug{};
    debug.source_index = source_index;
    debug.source_rank_in_k = source_rank_in_k;
    debug.k = node_branch.transfer_revolution;
    debug.source_q = node_branch.target_revolution;

    auto differentiated = node_branch;
    if (!differentiated.derivatives_available) {
        differentiated = problem1::attach_problem1_root_derivatives(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            node_branch);
    }
    if (!differentiated.derivatives_available) {
        debug.invalid_reason = "derivative_before_correction_invalid";
        return debug;
    }

    const auto hessian = build_test_only_hessian(
        table.config().departure_planet,
        table.config().target_planet,
        nearest.node_nu_A,
        nearest.node_nu_B,
        nearest.node_theta_A,
        differentiated,
        kQuadraticOneStepHessianStep);
    debug.hessian_valid = hessian.valid;
    debug.hessian_invalid_reason = hessian.invalid_reason;
    debug.hessian_stencil_seconds_gate_pass_count = hessian.stencil_pass_count;
    if (!hessian.valid) {
        debug.invalid_reason = hessian.invalid_reason;
        return debug;
    }

    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A);
    const double alpha_linear = normalize_angle_0_2pi(
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated.d_encounter_global_angle_d_theta_A * dx3);
    debug.alpha_linear = alpha_linear;
    const double alpha_quadratic = normalize_angle_0_2pi(
        alpha_linear +
        0.5 * (
            hessian.H_aa * dx1 * dx1 +
            hessian.H_bb * dx2 * dx2 +
            hessian.H_cc * dx3 * dx3 +
            2.0 * hessian.H_ab * dx1 * dx2 +
            2.0 * hessian.H_ac * dx1 * dx3 +
            2.0 * hessian.H_bc * dx2 * dx3));
    debug.alpha_quadratic = alpha_quadratic;

    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        node_branch.transfer_revolution,
        alpha_quadratic,
        node_branch,
        table.config().max_target_revolution);
    debug.selected_q = q_selection.selected_q;
    if (q_selection.selection_failed) {
        debug.invalid_reason = "q_sheet_selection_failed";
        return debug;
    }

    const auto residual_before = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        node_branch.transfer_revolution,
        q_selection.selected_q);
    if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
        debug.invalid_reason = "residual_before_correction_invalid";
        return debug;
    }
    debug.residual_before_correction = residual_before.residual_seconds;

    const auto derivatives_before = problem1::evaluate_problem1_root_residual_derivatives(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        node_branch.transfer_revolution,
        q_selection.selected_q);
    if (!derivatives_before.valid || !std::isfinite(derivatives_before.R_alpha)) {
        debug.invalid_reason = "derivative_before_correction_invalid";
        return debug;
    }
    if (std::abs(derivatives_before.R_alpha) <= 1e-12) {
        debug.invalid_reason = "correction_derivative_too_small";
        return debug;
    }

    const double residual_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    const double delta_alpha = -residual_scale_free / derivatives_before.R_alpha;
    const double alpha_corrected = normalize_angle_0_2pi(alpha_quadratic + delta_alpha);
    debug.alpha_corrected = alpha_corrected;

    const auto residual_after = problem1::evaluate_problem1_root_residual(
        table.config().departure_planet,
        table.config().target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_corrected,
        node_branch.transfer_revolution,
        q_selection.selected_q);
    if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
        debug.invalid_reason = "residual_after_correction_invalid";
        return debug;
    }
    debug.residual_after_correction = residual_after.residual_seconds;
    debug.transfer_time_seconds_after = residual_after.transfer_time_seconds;
    debug.target_time_seconds_after = residual_after.target_time_seconds;
    if (std::abs(residual_after.residual_seconds) > kEngineeringResidualToleranceSeconds) {
        debug.invalid_reason = "residual_after_correction_too_large";
        return debug;
    }

    debug.valid_candidate = true;
    debug.approximation.valid = true;
    debug.approximation.method = "route_b_quadratic_one_step_corrected_test_only";
    debug.approximation.transfer_revolution = node_branch.transfer_revolution;
    debug.approximation.target_revolution = q_selection.selected_q;
    debug.approximation.predicted_encounter_global_angle = residual_after.encounter_global_angle;
    debug.approximation.target_arrival_true_anomaly = residual_after.target_arrival_true_anomaly;
    debug.approximation.residual_scale_free = residual_after.residual_scale_free;
    debug.approximation.residual_seconds = residual_after.residual_seconds;
    debug.approximation.transfer_time_seconds = residual_after.transfer_time_seconds;
    debug.approximation.target_time_seconds = residual_after.target_time_seconds;
    debug.approximation.transfer_e = residual_after.transfer_e;
    debug.approximation.transfer_p = residual_after.transfer_p;
    debug.approximation.transfer_a = residual_after.transfer_a;
    debug.approximation.theta_B = residual_after.theta_B;
    debug.approximation.diagnostics.hessian_valid = true;
    return debug;
}

std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> exact_sorted_by_k_time(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& exact_roots
) {
    auto sorted = exact_roots;
    std::sort(sorted.begin(), sorted.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.transfer_revolution != rhs.transfer_revolution) {
            return lhs.transfer_revolution < rhs.transfer_revolution;
        }
        if (lhs.time_of_flight_seconds != rhs.time_of_flight_seconds) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        }
        return lhs.encounter_global_angle < rhs.encounter_global_angle;
    });
    return sorted;
}

std::vector<MatchPair> build_matching_a(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& exact_roots,
    const std::vector<CandidateDebug>& candidate_debugs
) {
    std::map<int, std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>> exact_by_k;
    std::map<int, std::vector<CandidateDebug>> candidate_by_k;
    for (const auto& root : exact_roots) {
        if (root.valid) {
            exact_by_k[root.transfer_revolution].push_back(root);
        }
    }
    for (const auto& candidate : candidate_debugs) {
        if (candidate.valid_candidate) {
            candidate_by_k[candidate.k].push_back(candidate);
        }
    }
    for (auto& [k, branches] : exact_by_k) {
        std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    }
    for (auto& [k, branches] : candidate_by_k) {
        std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.approximation.transfer_time_seconds < rhs.approximation.transfer_time_seconds;
        });
    }

    std::vector<MatchPair> pairs;
    for (const auto& [k, exact_group] : exact_by_k) {
        const auto it = candidate_by_k.find(k);
        if (it == candidate_by_k.end()) {
            continue;
        }
        const auto& candidate_group = it->second;
        const std::size_t pair_count = std::min(exact_group.size(), candidate_group.size());
        for (std::size_t i = 0; i < pair_count; ++i) {
            MatchPair pair{};
            pair.pair_rank = static_cast<int>(i);
            pair.exact_rank = static_cast<int>(i);
            pair.candidate_rank = static_cast<int>(i);
            pair.k = k;
            pair.exact_q = exact_group[i].target_revolution;
            pair.candidate_q = candidate_group[i].approximation.target_revolution;
            pair.exact_time = exact_group[i].time_of_flight_seconds;
            pair.candidate_time = candidate_group[i].approximation.transfer_time_seconds;
            pair.time_error = std::abs(pair.candidate_time - pair.exact_time);
            pair.exact_alpha = exact_group[i].encounter_global_angle;
            pair.candidate_alpha = candidate_group[i].approximation.predicted_encounter_global_angle;
            pair.alpha_error = wrapped_alpha_distance(pair.candidate_alpha, pair.exact_alpha);
            pair.candidate_residual = candidate_group[i].approximation.residual_seconds;
            pairs.push_back(pair);
        }
    }
    return pairs;
}

std::optional<SampleDebug> find_typical_failure_sample() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const auto samples = build_samples(12);

    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const auto& sample = samples[sample_index];
        if (sample.group_name != "near_node") {
            continue;
        }

        auto table = build_local_cell_table_for_query(departure_planet, target_planet, sample, 1, 1);
        const auto nearest = problem1::find_nearest_problem1_root_table_node(
            table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
        if (!nearest.valid || nearest.cell == nullptr) {
            continue;
        }

        std::vector<CandidateDebug> candidate_debugs;
        std::map<int, int> seen_rank_by_k;
        for (std::size_t source_index = 0; source_index < nearest.cell->solutions_sorted_by_time_of_flight.size(); ++source_index) {
            const auto& node_branch = nearest.cell->solutions_sorted_by_time_of_flight[source_index];
            if (!node_branch.valid) {
                continue;
            }
            const int source_rank_in_k = seen_rank_by_k[node_branch.transfer_revolution]++;
            candidate_debugs.push_back(evaluate_candidate_debug(
                table,
                nearest,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                node_branch,
                static_cast<int>(source_index),
                source_rank_in_k));
        }

        auto exact_roots = problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            1,
            1);
        exact_roots.erase(
            std::remove_if(exact_roots.begin(), exact_roots.end(), [](const auto& branch) { return !branch.valid; }),
            exact_roots.end());

        spaceship_cpp::problem1::Problem1RouteBSafeQueryResult safe_result{};
        const int cell_i = lower_corner_index(nearest.i, normalize_angle_minus_pi_pi(sample.query_nu_A - nearest.node_nu_A));
        const int cell_j = lower_corner_index(nearest.j, normalize_angle_minus_pi_pi(sample.query_nu_B - nearest.node_nu_B));
        const int cell_k = lower_corner_index(nearest.k, normalize_angle_minus_pi_pi(sample.query_theta_A - nearest.node_theta_A));
        safe_result.cell_admissibility = problem1::evaluate_problem1_root_cell_admissibility(
            table,
            cell_i,
            cell_j,
            cell_k,
            table.config().max_transfer_revolution,
            table.config().max_target_revolution);
        safe_result.expected_count_by_k = safe_result.cell_admissibility.reference_root_count_by_k;
        for (const auto& candidate : candidate_debugs) {
            if (candidate.valid_candidate) {
                safe_result.candidate_count_by_k[candidate.k] += 1;
            }
        }
        safe_result.branch_count_complete = true;
        for (const auto& [k, expected] : safe_result.expected_count_by_k) {
            const int actual = safe_result.candidate_count_by_k[k];
            if (actual != expected) {
                safe_result.branch_count_complete = false;
                safe_result.missing_count_by_k[k] = std::max(0, expected - actual);
                safe_result.extra_count_by_k[k] = std::max(0, actual - expected);
                safe_result.incomplete_by_k[k] = true;
            } else {
                safe_result.incomplete_by_k[k] = false;
            }
        }

        std::map<std::string, int> invalid_reason_distribution;
        for (const auto& candidate : candidate_debugs) {
            if (!candidate.valid_candidate) {
                invalid_reason_distribution[candidate.invalid_reason] += 1;
            }
        }

        const auto matching_a = build_matching_a(exact_roots, candidate_debugs);
        double max_time_error = 0.0;
        for (const auto& pair : matching_a) {
            max_time_error = std::max(max_time_error, pair.time_error);
        }
        bool has_small_residual = false;
        for (const auto& candidate : candidate_debugs) {
            if (candidate.valid_candidate &&
                std::abs(candidate.approximation.residual_seconds) <= kEngineeringResidualToleranceSeconds) {
                has_small_residual = true;
                break;
            }
        }

        if (has_small_residual && max_time_error > kHugeTimeErrorThresholdSeconds) {
            return SampleDebug{
                static_cast<int>(sample_index),
                sample,
                std::move(table),
                nearest,
                safe_result,
                std::move(invalid_reason_distribution),
                exact_sorted_by_k_time(exact_roots),
                std::move(candidate_debugs),
                max_time_error,
            };
        }
    }

    return std::nullopt;
}

void print_count_map(const std::string& label, const std::map<int, int>& values) {
    std::cout << label;
    bool first = true;
    for (const auto& [k, count] : values) {
        std::cout << (first ? "" : ",") << "k" << k << ":" << count;
        first = false;
    }
    if (first) {
        std::cout << "(empty)";
    }
    std::cout << '\n';
}

std::map<int, int> build_source_rank_map(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    std::vector<int>* rank_by_index
) {
    std::map<int, int> seen_rank_by_k;
    rank_by_index->assign(branches.size(), -1);
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (!branches[i].valid) {
            continue;
        }
        (*rank_by_index)[i] = seen_rank_by_k[branches[i].transfer_revolution]++;
    }
    return seen_rank_by_k;
}

std::vector<MatchPair> build_matching_b_or_c(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& exact_roots,
    const std::vector<CandidateDebug>& candidate_debugs,
    bool by_time
) {
    std::map<int, std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>> exact_by_k;
    for (const auto& root : exact_roots) {
        if (root.valid) {
            exact_by_k[root.transfer_revolution].push_back(root);
        }
    }
    for (auto& [k, group] : exact_by_k) {
        std::sort(group.begin(), group.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    }

    std::vector<MatchPair> matches;
    std::map<int, int> candidate_rank_by_k;
    for (const auto& candidate : candidate_debugs) {
        if (!candidate.valid_candidate) {
            continue;
        }
        const auto it = exact_by_k.find(candidate.k);
        if (it == exact_by_k.end() || it->second.empty()) {
            continue;
        }
        const auto& exact_group = it->second;
        double best_metric = std::numeric_limits<double>::infinity();
        int best_index = -1;
        for (std::size_t i = 0; i < exact_group.size(); ++i) {
            const double metric = by_time
                ? std::abs(candidate.approximation.transfer_time_seconds - exact_group[i].time_of_flight_seconds)
                : wrapped_alpha_distance(candidate.approximation.predicted_encounter_global_angle,
                                         exact_group[i].encounter_global_angle);
            if (metric < best_metric) {
                best_metric = metric;
                best_index = static_cast<int>(i);
            }
        }
        MatchPair pair{};
        pair.candidate_rank = candidate_rank_by_k[candidate.k]++;
        pair.exact_rank = best_index;
        pair.k = candidate.k;
        pair.exact_q = exact_group[static_cast<std::size_t>(best_index)].target_revolution;
        pair.candidate_q = candidate.approximation.target_revolution;
        pair.exact_time = exact_group[static_cast<std::size_t>(best_index)].time_of_flight_seconds;
        pair.candidate_time = candidate.approximation.transfer_time_seconds;
        pair.time_error = std::abs(pair.candidate_time - pair.exact_time);
        pair.exact_alpha = exact_group[static_cast<std::size_t>(best_index)].encounter_global_angle;
        pair.candidate_alpha = candidate.approximation.predicted_encounter_global_angle;
        pair.alpha_error = wrapped_alpha_distance(pair.candidate_alpha, pair.exact_alpha);
        pair.candidate_residual = candidate.approximation.residual_seconds;
        matches.push_back(pair);
    }
    return matches;
}

std::string classify_failure(
    const SampleDebug& debug,
    const std::vector<MatchPair>& matching_a,
    const std::vector<MatchPair>& matching_b,
    const std::vector<MatchPair>& matching_c
) {
    bool count_mismatch = false;
    for (const auto& [k, expected] : debug.safe_result.expected_count_by_k) {
        const int actual = debug.safe_result.candidate_count_by_k.contains(k)
            ? debug.safe_result.candidate_count_by_k.at(k)
            : 0;
        if (actual != expected) {
            count_mismatch = true;
            break;
        }
    }
    double max_a = 0.0;
    double max_b = 0.0;
    double max_c = 0.0;
    for (const auto& pair : matching_a) {
        max_a = std::max(max_a, pair.time_error);
    }
    for (const auto& pair : matching_b) {
        max_b = std::max(max_b, pair.time_error);
    }
    for (const auto& pair : matching_c) {
        max_c = std::max(max_c, pair.time_error);
    }
    bool q_mismatch_found = false;
    for (const auto& pair : matching_b) {
        if (pair.exact_q != pair.candidate_q) {
            q_mismatch_found = true;
            break;
        }
    }
    int hessian_fail_count = 0;
    for (const auto& candidate : debug.candidate_debugs) {
        if (candidate.invalid_reason == "hessian_stencil_seconds_gate_failed") {
            hessian_fail_count += 1;
        }
    }

    if (count_mismatch && max_b < 1.0 && max_a > kHugeTimeErrorThresholdSeconds) {
        return "count_mismatch_causes_order_shift";
    }
    if (q_mismatch_found && max_b < 1.0) {
        return "selected_q_wrong";
    }
    if (max_b < 1.0 || max_c < 1.0) {
        return "matching_metric_wrong";
    }
    if (hessian_fail_count > 0 && count_mismatch) {
        return "hessian_failed_for_missing_branch";
    }
    return "other";
}

}  // namespace

int main() {
    std::cout << std::setprecision(15) << std::scientific;

    const auto debug_opt = find_typical_failure_sample();
    if (!debug_opt.has_value()) {
        std::cout << "failed_to_find_typical_failure_sample\n";
        return 1;
    }
    const auto& debug = *debug_opt;
    const auto& nearest = debug.nearest;

    std::cout << "OneStepBranchDebug\n";
    std::cout << "group=" << debug.sample.group_name << '\n';
    std::cout << "sample_index=" << debug.sample_index << '\n';
    std::cout << "query_nu_A=" << debug.sample.query_nu_A << '\n';
    std::cout << "query_nu_B=" << debug.sample.query_nu_B << '\n';
    std::cout << "query_theta_A=" << debug.sample.query_theta_A << '\n';
    std::cout << "cell_admissible=" << debug.safe_result.cell_admissibility.admissible << '\n';
    std::cout << "branch_count_complete=" << debug.safe_result.branch_count_complete << '\n';
    print_count_map("expected_count_by_k=", debug.safe_result.expected_count_by_k);
    print_count_map("candidate_count_by_k=", debug.safe_result.candidate_count_by_k);
    print_count_map("missing_count_by_k=", debug.safe_result.missing_count_by_k);
    print_count_map("extra_count_by_k=", debug.safe_result.extra_count_by_k);
    for (const auto& [reason, count] : debug.invalid_reason_distribution) {
        std::cout << "invalid_reason[" << reason << "]=" << count << '\n';
    }

    std::cout << "ExactRootsBegin\n";
    std::map<int, int> exact_rank_by_k;
    for (const auto& root : debug.exact_roots) {
        const int rank = exact_rank_by_k[root.transfer_revolution]++;
        std::cout << "exact_rank_in_k=" << rank
                  << " k=" << root.transfer_revolution
                  << " q=" << root.target_revolution
                  << " alpha=" << root.encounter_global_angle
                  << " time_of_flight_seconds=" << root.time_of_flight_seconds
                  << " target_time_seconds=" << root.target_time_seconds
                  << " residual_seconds=" << root.residual_seconds << '\n';
    }
    std::cout << "ExactRootsEnd\n";

    std::cout << "SourceBranchesBegin\n";
    std::vector<int> source_rank_by_index;
    build_source_rank_map(nearest.cell->solutions_sorted_by_time_of_flight, &source_rank_by_index);
    for (std::size_t i = 0; i < nearest.cell->solutions_sorted_by_time_of_flight.size(); ++i) {
        const auto& root = nearest.cell->solutions_sorted_by_time_of_flight[i];
        if (!root.valid) {
            continue;
        }
        std::cout << "source_rank_in_k=" << source_rank_by_index[i]
                  << " source_index=" << i
                  << " k=" << root.transfer_revolution
                  << " q=" << root.target_revolution
                  << " alpha=" << root.encounter_global_angle
                  << " time_of_flight_seconds=" << root.time_of_flight_seconds
                  << " target_time_seconds=" << root.target_time_seconds
                  << " derivatives_available=" << root.derivatives_available << '\n';
    }
    std::cout << "SourceBranchesEnd\n";

    std::cout << "OneStepCandidatesBegin\n";
    for (const auto& candidate : debug.candidate_debugs) {
        std::cout << "source_index=" << candidate.source_index
                  << " source_rank_in_k=" << candidate.source_rank_in_k
                  << " k=" << candidate.k
                  << " source_q=" << candidate.source_q
                  << " selected_q=" << candidate.selected_q
                  << " hessian_step=" << candidate.hessian_step
                  << " hessian_valid=" << candidate.hessian_valid
                  << " hessian_invalid_reason=" << candidate.hessian_invalid_reason
                  << " hessian_stencil_seconds_gate_pass_count=" << candidate.hessian_stencil_seconds_gate_pass_count
                  << " alpha_linear=" << candidate.alpha_linear
                  << " alpha_quadratic=" << candidate.alpha_quadratic
                  << " alpha_corrected=" << candidate.alpha_corrected
                  << " residual_before_correction=" << candidate.residual_before_correction
                  << " residual_after_correction=" << candidate.residual_after_correction
                  << " transfer_time_seconds_after=" << candidate.transfer_time_seconds_after
                  << " target_time_seconds_after=" << candidate.target_time_seconds_after
                  << " valid_candidate=" << candidate.valid_candidate
                  << " invalid_reason=" << candidate.invalid_reason << '\n';
    }
    std::cout << "OneStepCandidatesEnd\n";

    const auto matching_a = build_matching_a(debug.exact_roots, debug.candidate_debugs);
    const auto matching_b = build_matching_b_or_c(debug.exact_roots, debug.candidate_debugs, true);
    const auto matching_c = build_matching_b_or_c(debug.exact_roots, debug.candidate_debugs, false);

    std::cout << "MatchingABegin\n";
    for (const auto& pair : matching_a) {
        std::cout << "pair_rank=" << pair.pair_rank
                  << " exact_rank=" << pair.exact_rank
                  << " candidate_rank=" << pair.candidate_rank
                  << " k=" << pair.k
                  << " exact_time=" << pair.exact_time
                  << " candidate_time=" << pair.candidate_time
                  << " time_error=" << pair.time_error
                  << " exact_alpha=" << pair.exact_alpha
                  << " candidate_alpha=" << pair.candidate_alpha
                  << " alpha_error=" << pair.alpha_error
                  << " exact_q=" << pair.exact_q
                  << " candidate_q=" << pair.candidate_q
                  << " candidate_residual=" << pair.candidate_residual << '\n';
    }
    std::cout << "MatchingAEnd\n";

    std::cout << "MatchingBBegin\n";
    for (const auto& pair : matching_b) {
        std::cout << "candidate_rank=" << pair.candidate_rank
                  << " nearest_exact_rank_by_time=" << pair.exact_rank
                  << " k=" << pair.k
                  << " time_error=" << pair.time_error
                  << " alpha_error=" << pair.alpha_error
                  << " exact_q=" << pair.exact_q
                  << " candidate_q=" << pair.candidate_q << '\n';
    }
    std::cout << "MatchingBEnd\n";

    std::cout << "MatchingCBegin\n";
    for (const auto& pair : matching_c) {
        std::cout << "candidate_rank=" << pair.candidate_rank
                  << " nearest_exact_rank_by_alpha=" << pair.exact_rank
                  << " k=" << pair.k
                  << " time_error=" << pair.time_error
                  << " alpha_error=" << pair.alpha_error
                  << " exact_q=" << pair.exact_q
                  << " candidate_q=" << pair.candidate_q << '\n';
    }
    std::cout << "MatchingCEnd\n";

    std::cout << "failure_type=" << classify_failure(debug, matching_a, matching_b, matching_c) << '\n';
    return 0;
}
