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
};

struct RouteBEvalResult {
    bool valid = false;
    std::string invalid_reason;
    double selected_step = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct FallbackCase {
    std::string group;
    int sample_index = -1;
    int source_index = -1;
    int source_rank_in_k = -1;
    int k = 0;
    int source_q = 0;
    std::string fallback_reason;
    double route_a_time_of_flight = std::numeric_limits<double>::quiet_NaN();
    double route_a_alpha = std::numeric_limits<double>::quiet_NaN();
    double route_a_residual_seconds = std::numeric_limits<double>::quiet_NaN();
};

struct GroupSummary {
    int route_a_valid_count = 0;
    int route_b_used_count = 0;
    int route_a_fallback_count = 0;
    int final_result_count = 0;
    int route_b_same_root_success_count = 0;
    double sum_time_diff_route_b_vs_route_a = 0.0;
    double max_time_diff_route_b_vs_route_a = 0.0;
    double sum_alpha_diff_route_b_vs_route_a = 0.0;
    double max_alpha_diff_route_b_vs_route_a = 0.0;
    double sum_abs_route_b_residual_seconds = 0.0;
    double max_abs_route_b_residual_seconds = 0.0;
    double pure_route_a_total_ms = 0.0;
    double hybrid_route_b_attempt_total_ms = 0.0;
    double hybrid_fallback_route_a_total_ms = 0.0;
    std::map<std::string, int> fallback_reasons;
    std::map<std::string, int> route_b_selected_step_distribution;
    std::vector<FallbackCase> fallback_cases;
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

RouteBEvalResult evaluate_route_b_projected_one_step_from_source_branch_test_only(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const spaceship_cpp::problem1::Problem1RootNearestNode& nearest,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    RouteBEvalResult result{};
    for (const double step : kAdaptiveHessianSteps) {
        const auto prediction = problem1::predict_problem1_root_branch_quadratic_from_node(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            source_branch,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            step);
        if (!prediction.valid) {
            result.invalid_reason =
                prediction.invalid_reason.empty() ? "quadratic_prediction_invalid" : prediction.invalid_reason;
            continue;
        }
        if (prediction.hessian_method != "projected_tangent_finite_difference_of_implicit_first_derivatives") {
            result.invalid_reason = "default_hessian_method_is_not_projected";
            return result;
        }
        const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            source_branch.transfer_revolution,
            prediction.predicted_encounter_global_angle,
            source_branch,
            table.config().max_target_revolution);
        if (q_selection.selection_failed) {
            result.invalid_reason = "q_sheet_selection_failed";
            continue;
        }
        const auto residual_before = problem1::evaluate_problem1_root_residual(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            prediction.predicted_encounter_global_angle,
            source_branch.transfer_revolution,
            q_selection.selected_q);
        if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
            result.invalid_reason = "residual_before_correction_invalid";
            continue;
        }
        const auto derivative_before = problem1::evaluate_problem1_root_residual_derivatives(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            prediction.predicted_encounter_global_angle,
            source_branch.transfer_revolution,
            q_selection.selected_q);
        if (!derivative_before.valid || !std::isfinite(derivative_before.R_alpha)) {
            result.invalid_reason = "derivative_before_correction_invalid";
            continue;
        }
        if (std::abs(derivative_before.R_alpha) <= 1e-12) {
            result.invalid_reason = "correction_derivative_too_small";
            continue;
        }
        const double residual_scale_free =
            problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
        const double delta_alpha = -residual_scale_free / derivative_before.R_alpha;
        const double alpha_corrected =
            normalize_angle_0_2pi(prediction.predicted_encounter_global_angle + delta_alpha);
        const auto residual_after = problem1::evaluate_problem1_root_residual(
            table.config().departure_planet,
            table.config().target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_corrected,
            source_branch.transfer_revolution,
            q_selection.selected_q);
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
        return result;
    }
    return result;
}

void print_summary(const std::string& group, const GroupSummary& s) {
    const double route_b_usage_ratio =
        s.route_a_valid_count > 0 ? static_cast<double>(s.route_b_used_count) / s.route_a_valid_count : 0.0;
    const double route_a_fallback_ratio =
        s.route_a_valid_count > 0 ? static_cast<double>(s.route_a_fallback_count) / s.route_a_valid_count : 0.0;
    const double final_coverage =
        s.route_a_valid_count > 0 ? static_cast<double>(s.final_result_count) / s.route_a_valid_count : 0.0;
    const double route_b_same_root_ratio =
        s.route_b_used_count > 0 ? static_cast<double>(s.route_b_same_root_success_count) / s.route_b_used_count : 0.0;
    const double pure_avg = s.route_a_valid_count > 0 ? s.pure_route_a_total_ms / s.route_a_valid_count : 0.0;
    const double hybrid_total = s.hybrid_route_b_attempt_total_ms + s.hybrid_fallback_route_a_total_ms;
    const double hybrid_avg = s.route_a_valid_count > 0 ? hybrid_total / s.route_a_valid_count : 0.0;
    const double hybrid_speedup = hybrid_total > 0.0 ? s.pure_route_a_total_ms / hybrid_total : 0.0;

    std::cout << "HybridProjectedRouteBRouteAFallbackSummary\n";
    std::cout << "group=" << group << '\n';
    std::cout << "route_a_valid_count=" << s.route_a_valid_count << '\n';
    std::cout << "route_b_used_count=" << s.route_b_used_count << '\n';
    std::cout << "route_a_fallback_count=" << s.route_a_fallback_count << '\n';
    std::cout << "route_b_usage_ratio=" << route_b_usage_ratio << '\n';
    std::cout << "route_a_fallback_ratio=" << route_a_fallback_ratio << '\n';
    std::cout << "final_result_count=" << s.final_result_count << '\n';
    std::cout << "final_coverage_vs_route_a_valid=" << final_coverage << '\n';
    std::cout << "route_b_same_root_success_count=" << s.route_b_same_root_success_count << '\n';
    std::cout << "route_b_same_root_success_ratio=" << route_b_same_root_ratio << '\n';
    std::cout << "mean_time_diff_route_b_vs_route_a="
              << (s.route_b_used_count > 0 ? s.sum_time_diff_route_b_vs_route_a / s.route_b_used_count : 0.0) << '\n';
    std::cout << "max_time_diff_route_b_vs_route_a=" << s.max_time_diff_route_b_vs_route_a << '\n';
    std::cout << "mean_alpha_diff_route_b_vs_route_a="
              << (s.route_b_used_count > 0 ? s.sum_alpha_diff_route_b_vs_route_a / s.route_b_used_count : 0.0) << '\n';
    std::cout << "max_alpha_diff_route_b_vs_route_a=" << s.max_alpha_diff_route_b_vs_route_a << '\n';
    std::cout << "mean_abs_route_b_residual_seconds="
              << (s.route_b_used_count > 0 ? s.sum_abs_route_b_residual_seconds / s.route_b_used_count : 0.0) << '\n';
    std::cout << "max_abs_route_b_residual_seconds=" << s.max_abs_route_b_residual_seconds << '\n';
    for (const auto& [reason, count] : s.fallback_reasons) {
        std::cout << "fallback_reason[" << reason << "]=" << count << '\n';
    }
    for (const auto& key : std::array<std::string, 4>{"1e-5", "5e-6", "2e-6", "1e-6"}) {
        const auto it = s.route_b_selected_step_distribution.find(key);
        std::cout << "route_b_selected_step[" << key << "]="
                  << (it != s.route_b_selected_step_distribution.end() ? it->second : 0) << '\n';
    }
    std::cout << "pure_route_a_total_ms=" << s.pure_route_a_total_ms << '\n';
    std::cout << "hybrid_route_b_attempt_total_ms=" << s.hybrid_route_b_attempt_total_ms << '\n';
    std::cout << "hybrid_fallback_route_a_total_ms=" << s.hybrid_fallback_route_a_total_ms << '\n';
    std::cout << "hybrid_total_ms=" << hybrid_total << '\n';
    std::cout << "pure_route_a_avg_per_branch_ms=" << pure_avg << '\n';
    std::cout << "hybrid_avg_per_branch_ms=" << hybrid_avg << '\n';
    std::cout << "hybrid_speedup_vs_pure_route_a=" << hybrid_speedup << '\n';
    std::cout << "hybrid_covers_all_route_a_valid=" << (final_coverage == 1.0 ? 1 : 0) << '\n';
    std::cout << "route_b_accuracy_ok=" << (route_b_same_root_ratio == 1.0 ? 1 : 0) << '\n';
    std::cout << "recommended_next_step="
              << ((final_coverage == 1.0 && route_b_same_root_ratio == 1.0 && hybrid_speedup > 1.1)
                      ? "consider_production_high_level_wrapper"
                      : "keep_projected_route_b_as_low_level_capability_and_optimize_fallback_or_speed")
              << '\n';

    const int fallback_to_print = std::min<int>(5, s.fallback_cases.size());
    for (int i = 0; i < fallback_to_print; ++i) {
        const auto& c = s.fallback_cases[static_cast<std::size_t>(i)];
        std::cout << "HybridFallbackCase\n";
        std::cout << "group=" << c.group << '\n';
        std::cout << "sample_index=" << c.sample_index << '\n';
        std::cout << "source_index=" << c.source_index << '\n';
        std::cout << "source_rank_in_k=" << c.source_rank_in_k << '\n';
        std::cout << "k=" << c.k << '\n';
        std::cout << "source_q=" << c.source_q << '\n';
        std::cout << "fallback_reason=" << c.fallback_reason << '\n';
        std::cout << "route_a_time_of_flight=" << c.route_a_time_of_flight << '\n';
        std::cout << "route_a_alpha=" << c.route_a_alpha << '\n';
        std::cout << "route_a_residual_seconds=" << c.route_a_residual_seconds << '\n';
    }
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    using clock = std::chrono::steady_clock;

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
        std::map<int, int> source_rank_by_k;
        for (std::size_t source_index = 0; source_index < nearest.cell->solutions_sorted_by_time_of_flight.size(); ++source_index) {
            const auto& source_branch = nearest.cell->solutions_sorted_by_time_of_flight[source_index];
            if (!source_branch.valid) {
                continue;
            }
            auto& summary = summaries[sample.group_name];
            const int source_rank_in_k = source_rank_by_k[source_branch.transfer_revolution]++;

            const auto pure_start = clock::now();
            const auto route_a_reference = evaluate_route_a_from_source_branch_test_only(
                table, nearest, source_branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            const auto pure_end = clock::now();
            summary.pure_route_a_total_ms +=
                std::chrono::duration<double, std::milli>(pure_end - pure_start).count();
            if (!route_a_reference.valid) {
                continue;
            }
            summary.route_a_valid_count += 1;

            const auto route_b_start = clock::now();
            const auto route_b = evaluate_route_b_projected_one_step_from_source_branch_test_only(
                table, nearest, source_branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            const auto route_b_end = clock::now();
            summary.hybrid_route_b_attempt_total_ms +=
                std::chrono::duration<double, std::milli>(route_b_end - route_b_start).count();

            if (route_b.valid) {
                summary.route_b_used_count += 1;
                summary.final_result_count += 1;
                summary.route_b_selected_step_distribution[step_label(route_b.selected_step)] += 1;
                const double time_diff =
                    std::abs(route_b.branch.time_of_flight_seconds - route_a_reference.branch.time_of_flight_seconds);
                const double alpha_diff =
                    wrapped_alpha_distance(route_b.branch.encounter_global_angle, route_a_reference.branch.encounter_global_angle);
                const double abs_residual = std::abs(route_b.branch.residual_seconds);
                summary.sum_time_diff_route_b_vs_route_a += time_diff;
                summary.max_time_diff_route_b_vs_route_a = std::max(summary.max_time_diff_route_b_vs_route_a, time_diff);
                summary.sum_alpha_diff_route_b_vs_route_a += alpha_diff;
                summary.max_alpha_diff_route_b_vs_route_a = std::max(summary.max_alpha_diff_route_b_vs_route_a, alpha_diff);
                summary.sum_abs_route_b_residual_seconds += abs_residual;
                summary.max_abs_route_b_residual_seconds = std::max(summary.max_abs_route_b_residual_seconds, abs_residual);
                if (time_diff <= kSameRootTimeThresholdSeconds &&
                    alpha_diff <= kSameRootAlphaThreshold &&
                    abs_residual <= kResidualToleranceSeconds) {
                    summary.route_b_same_root_success_count += 1;
                }
                continue;
            }

            summary.route_a_fallback_count += 1;
            summary.final_result_count += 1;
            summary.fallback_reasons[route_b.invalid_reason] += 1;
            if (summary.fallback_cases.size() < 5) {
                summary.fallback_cases.push_back(FallbackCase{
                    sample.group_name,
                    static_cast<int>(sample_index),
                    static_cast<int>(source_index),
                    source_rank_in_k,
                    source_branch.transfer_revolution,
                    source_branch.target_revolution,
                    route_b.invalid_reason,
                    route_a_reference.branch.time_of_flight_seconds,
                    route_a_reference.branch.encounter_global_angle,
                    route_a_reference.branch.residual_seconds,
                });
            }

            const auto fallback_start = clock::now();
            const auto fallback_route_a = evaluate_route_a_from_source_branch_test_only(
                table, nearest, source_branch, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            const auto fallback_end = clock::now();
            summary.hybrid_fallback_route_a_total_ms +=
                std::chrono::duration<double, std::milli>(fallback_end - fallback_start).count();
            if (!fallback_route_a.valid) {
                std::cerr << "hybrid_fallback_route_a_failed\n";
                return EXIT_FAILURE;
            }
        }
    }

    for (const auto& [group, summary] : summaries) {
        print_summary(group, summary);
    }

    std::cout << "route_b_projected_with_route_a_fallback_compare_ok\n";
    return EXIT_SUCCESS;
}
