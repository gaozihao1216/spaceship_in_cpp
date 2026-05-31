#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <array>
#include <chrono>
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
constexpr double kMatchTimeThresholdSeconds = 1.0;
constexpr double kMatchAlphaThreshold = 1e-3;
constexpr double kEngineeringResidualToleranceSeconds = 1e-2;
constexpr double kQuadraticOneStepHessianStep = 1e-5;

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

struct MethodStats {
    int total_queries = 0;
    int success_count = 0;
    int fallback_required_count = 0;
    int branch_count_complete_count = 0;
    int no_valid_approximation_count = 0;
    int non_admissible_count = 0;

    int matched_branch_count = 0;
    int total_exact_branch_count = 0;
    double sum_abs_time_error = 0.0;
    double max_abs_time_error = 0.0;
    double sum_alpha_wrapped_error = 0.0;
    double max_alpha_wrapped_error = 0.0;
    double sum_abs_residual_seconds = 0.0;
    double max_abs_residual_seconds = 0.0;
    int paired_branch_count = 0;

    double total_time_seconds = 0.0;
    std::map<std::string, int> invalid_reason_distribution;
};

struct StencilAttachedBranch {
    bool residual_gate_passed = false;
    bool attach_success = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct TestOnlyHessian {
    bool valid = false;
    std::string invalid_reason;
    double H_aa = 0.0;
    double H_bb = 0.0;
    double H_cc = 0.0;
    double H_ab = 0.0;
    double H_ac = 0.0;
    double H_bc = 0.0;
};

struct TestOnlyRouteBQueryResult {
    spaceship_cpp::problem1::Problem1RouteBSafeQueryResult safe_result;
    std::map<std::string, int> invalid_reason_distribution;
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
    node.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

int lower_corner_index(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
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

double wrapped_alpha_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

double branch_time_duplicate_tolerance_seconds(double t1, double t2) {
    return std::max(1e-4, 1e-10 * std::max(std::abs(t1), std::abs(t2)));
}

bool add_branch_if_not_duplicate(
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* branches,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    for (auto& existing : *branches) {
        if (existing.transfer_revolution != branch.transfer_revolution) {
            continue;
        }
        const double time_tol =
            branch_time_duplicate_tolerance_seconds(existing.time_of_flight_seconds, branch.time_of_flight_seconds);
        const double time_diff = std::abs(existing.time_of_flight_seconds - branch.time_of_flight_seconds);
        const double angle_diff = wrapped_alpha_distance(existing.encounter_global_angle, branch.encounter_global_angle);
        if (time_diff <= time_tol && angle_diff <= 1e-4) {
            if (std::abs(branch.residual_seconds) < std::abs(existing.residual_seconds)) {
                existing = branch;
            }
            return true;
        }
    }
    branches->push_back(branch);
    return false;
}

spaceship_cpp::problem1::Problem1SolutionBranch approximation_to_branch(
    const spaceship_cpp::problem1::Problem1RootApproximationResult& approximation
) {
    spaceship_cpp::problem1::Problem1SolutionBranch branch{};
    branch.valid = approximation.valid;
    branch.encounter_global_angle = approximation.predicted_encounter_global_angle;
    branch.target_arrival_true_anomaly = approximation.target_arrival_true_anomaly;
    branch.transfer_revolution = approximation.transfer_revolution;
    branch.target_revolution = approximation.target_revolution;
    branch.time_of_flight_seconds = approximation.transfer_time_seconds;
    branch.target_time_seconds = approximation.target_time_seconds;
    branch.residual_seconds = approximation.residual_seconds;
    branch.transfer_e = approximation.transfer_e;
    branch.transfer_p = approximation.transfer_p;
    branch.transfer_a = approximation.transfer_a;
    branch.theta_B = approximation.theta_B;
    return branch;
}

std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> route_a_baseline_query(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    const QuerySample& sample
) {
    namespace problem1 = spaceship_cpp::problem1;
    const auto nearest = problem1::find_nearest_problem1_root_table_node(
        table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
    std::vector<problem1::Problem1SolutionBranch> branches;
    if (!nearest.valid || nearest.cell == nullptr) {
        return branches;
    }
    for (const auto& source : nearest.cell->solutions_sorted_by_time_of_flight) {
        if (!source.valid) {
            continue;
        }
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            source,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (!attached.valid || !attached.derivatives_available) {
            continue;
        }
        const double dnu_A = normalize_angle_minus_pi_pi(sample.query_nu_A - nearest.node_nu_A);
        const double dnu_B = normalize_angle_minus_pi_pi(sample.query_nu_B - nearest.node_nu_B);
        const double dtheta_A = normalize_angle_minus_pi_pi(sample.query_theta_A - nearest.node_theta_A);
        const double alpha_linear = normalize_angle_0_2pi(
            attached.encounter_global_angle +
            attached.d_encounter_global_angle_d_nu_A * dnu_A +
            attached.d_encounter_global_angle_d_nu_B * dnu_B +
            attached.d_encounter_global_angle_d_theta_A * dtheta_A);
        const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
            table.config().departure_planet,
            table.config().target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            source.transfer_revolution,
            alpha_linear,
            source,
            table.config().max_target_revolution);
        const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            table.config().departure_planet,
            table.config().target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            source.transfer_revolution,
            q_selection.selected_q,
            alpha_linear,
            80,
            1e-2,
            1e-12,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (!refined.valid) {
            continue;
        }
        add_branch_if_not_duplicate(&branches, refined.branch);
    }
    std::sort(
        branches.begin(),
        branches.end(),
        [](const auto& lhs, const auto& rhs) {
            if (lhs.transfer_revolution != rhs.transfer_revolution) {
                return lhs.transfer_revolution < rhs.transfer_revolution;
            }
            if (lhs.time_of_flight_seconds != rhs.time_of_flight_seconds) {
                return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
            }
            return lhs.encounter_global_angle < rhs.encounter_global_angle;
        });
    return branches;
}

template <typename F>
auto time_call(F&& fn, double* elapsed_seconds) {
    const auto start = std::chrono::steady_clock::now();
    auto result = fn();
    const auto end = std::chrono::steady_clock::now();
    *elapsed_seconds = std::chrono::duration<double>(end - start).count();
    return result;
}

void record_branch_comparison(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& exact,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& approximate,
    MethodStats* stats
) {
    std::map<int, std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>> exact_by_k;
    std::map<int, std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>> approx_by_k;
    for (const auto& branch : exact) {
        if (branch.valid) {
            exact_by_k[branch.transfer_revolution].push_back(branch);
            stats->total_exact_branch_count += 1;
        }
    }
    for (const auto& branch : approximate) {
        if (branch.valid) {
            approx_by_k[branch.transfer_revolution].push_back(branch);
        }
    }
    for (auto& [k, branches] : exact_by_k) {
        std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    }
    for (auto& [k, branches] : approx_by_k) {
        std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    }
    for (const auto& [k, exact_group] : exact_by_k) {
        const auto approx_it = approx_by_k.find(k);
        if (approx_it == approx_by_k.end()) {
            continue;
        }
        const auto& approx_group = approx_it->second;
        const std::size_t pair_count = std::min(exact_group.size(), approx_group.size());
        for (std::size_t idx = 0; idx < pair_count; ++idx) {
            const auto& exact_branch = exact_group[idx];
            const auto& approx_branch = approx_group[idx];
            const double time_error = std::abs(approx_branch.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
            const double alpha_error = wrapped_alpha_distance(
                approx_branch.encounter_global_angle,
                exact_branch.encounter_global_angle);
            const double residual_abs = std::abs(approx_branch.residual_seconds);
            stats->paired_branch_count += 1;
            stats->sum_abs_time_error += time_error;
            stats->max_abs_time_error = std::max(stats->max_abs_time_error, time_error);
            stats->sum_alpha_wrapped_error += alpha_error;
            stats->max_alpha_wrapped_error = std::max(stats->max_alpha_wrapped_error, alpha_error);
            stats->sum_abs_residual_seconds += residual_abs;
            stats->max_abs_residual_seconds = std::max(stats->max_abs_residual_seconds, residual_abs);
            if (time_error <= kMatchTimeThresholdSeconds && alpha_error <= kMatchAlphaThreshold) {
                stats->matched_branch_count += 1;
            }
        }
    }
}

void print_method_summary(
    const std::string& group_name,
    const std::string& method_name,
    const MethodStats& stats,
    double route_a_avg_ms = 0.0
) {
    std::cout << "RouteBOneStepCompareSummary\n";
    std::cout << "group=" << group_name << '\n';
    std::cout << "method=" << method_name << '\n';
    std::cout << "total_queries=" << stats.total_queries << '\n';
    std::cout << "success_count=" << stats.success_count << '\n';
    std::cout << "fallback_required_count=" << stats.fallback_required_count << '\n';
    std::cout << "branch_count_complete_count=" << stats.branch_count_complete_count << '\n';
    std::cout << "no_valid_approximation_count=" << stats.no_valid_approximation_count << '\n';
    std::cout << "non_admissible_count=" << stats.non_admissible_count << '\n';
    for (const auto& [reason, count] : stats.invalid_reason_distribution) {
        std::cout << "invalid_reason[" << reason << "]=" << count << '\n';
    }
    std::cout << "matched_branch_count=" << stats.matched_branch_count << '\n';
    std::cout << "total_exact_branch_count=" << stats.total_exact_branch_count << '\n';
    std::cout << "coverage_ratio=" <<
        (stats.total_exact_branch_count > 0
            ? static_cast<double>(stats.matched_branch_count) / static_cast<double>(stats.total_exact_branch_count)
            : 0.0) << '\n';
    std::cout << "mean_abs_time_error=" <<
        (stats.paired_branch_count > 0 ? stats.sum_abs_time_error / stats.paired_branch_count : 0.0) << '\n';
    std::cout << "max_abs_time_error=" << stats.max_abs_time_error << '\n';
    std::cout << "mean_alpha_wrapped_error=" <<
        (stats.paired_branch_count > 0 ? stats.sum_alpha_wrapped_error / stats.paired_branch_count : 0.0) << '\n';
    std::cout << "max_alpha_wrapped_error=" << stats.max_alpha_wrapped_error << '\n';
    std::cout << "mean_abs_residual_seconds=" <<
        (stats.paired_branch_count > 0 ? stats.sum_abs_residual_seconds / stats.paired_branch_count : 0.0) << '\n';
    std::cout << "max_abs_residual_seconds=" << stats.max_abs_residual_seconds << '\n';
    const double avg_ms = stats.total_queries > 0
        ? 1000.0 * stats.total_time_seconds / static_cast<double>(stats.total_queries)
        : 0.0;
    std::cout << "avg_ms=" << avg_ms << '\n';
    if (route_a_avg_ms > 0.0) {
        std::cout << "speedup_vs_route_a=" << (avg_ms > 0.0 ? route_a_avg_ms / avg_ms : 0.0) << '\n';
    }
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
        return result;
    }
    return result;
}

spaceship_cpp::problem1::Problem1RootApproximationResult
evaluate_problem1_root_quadratic_one_step_corrected_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_target_revolution
) {
    namespace problem1 = spaceship_cpp::problem1;
    problem1::Problem1RootApproximationResult result{};
    result.method = "route_b_quadratic_one_step_corrected_test_only";
    result.transfer_revolution = node_branch.transfer_revolution;
    result.target_revolution = node_branch.target_revolution;
    result.diagnostics.source_target_revolution = node_branch.target_revolution;
    result.diagnostics.selected_target_revolution = node_branch.target_revolution;
    result.diagnostics.hessian_method = "debug_seconds_gate_tangent_finite_difference";
    result.diagnostics.hessian_step = kQuadraticOneStepHessianStep;

    if (!node_branch.valid) {
        result.invalid_reason = "node_branch_invalid";
        return result;
    }
    auto differentiated = node_branch;
    if (!differentiated.derivatives_available) {
        differentiated = problem1::attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated.derivatives_available) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }

    const auto hessian = build_test_only_hessian(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        differentiated,
        kQuadraticOneStepHessianStep);
    result.diagnostics.hessian_valid = hessian.valid;
    if (!hessian.valid) {
        result.invalid_reason = hessian.invalid_reason;
        return result;
    }

    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - node_theta_A);
    const double alpha_linear = normalize_angle_0_2pi(
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated.d_encounter_global_angle_d_theta_A * dx3);
    result.diagnostics.alpha_linear = alpha_linear;
    const double alpha_quadratic = normalize_angle_0_2pi(
        alpha_linear +
        0.5 * (
            hessian.H_aa * dx1 * dx1 +
            hessian.H_bb * dx2 * dx2 +
            hessian.H_cc * dx3 * dx3 +
            2.0 * hessian.H_ab * dx1 * dx2 +
            2.0 * hessian.H_ac * dx1 * dx3 +
            2.0 * hessian.H_bc * dx2 * dx3));
    result.diagnostics.alpha_quadratic = alpha_quadratic;

    const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        node_branch.transfer_revolution,
        alpha_quadratic,
        node_branch,
        max_target_revolution);
    result.diagnostics.q_sheet_selection_failed = q_selection.selection_failed;
    result.diagnostics.q_sheet_selection_changed = q_selection.q_changed;
    result.diagnostics.selected_q_continuity_error = q_selection.selected_continuity_error;
    result.diagnostics.source_q_continuity_error = q_selection.source_q_continuity_error;
    if (q_selection.selection_failed) {
        result.invalid_reason = "q_sheet_selection_failed";
        return result;
    }
    result.target_revolution = q_selection.selected_q;
    result.diagnostics.selected_target_revolution = q_selection.selected_q;

    const auto residual_before = problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        node_branch.transfer_revolution,
        q_selection.selected_q);
    if (!residual_before.valid || !std::isfinite(residual_before.residual_seconds)) {
        result.invalid_reason = "residual_before_correction_invalid";
        return result;
    }

    const auto derivatives_before = problem1::evaluate_problem1_root_residual_derivatives(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        node_branch.transfer_revolution,
        q_selection.selected_q);
    if (!derivatives_before.valid || !std::isfinite(derivatives_before.R_alpha)) {
        result.invalid_reason = "derivative_before_correction_invalid";
        return result;
    }
    if (std::abs(derivatives_before.R_alpha) <= 1e-12) {
        result.invalid_reason = "correction_derivative_too_small";
        return result;
    }

    const double residual_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    const double delta_alpha = -residual_scale_free / derivatives_before.R_alpha;
    const double alpha_corrected = normalize_angle_0_2pi(alpha_quadratic + delta_alpha);

    const auto residual_after = problem1::evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_corrected,
        node_branch.transfer_revolution,
        q_selection.selected_q);
    if (!residual_after.valid || !std::isfinite(residual_after.residual_seconds)) {
        result.invalid_reason = "residual_after_correction_invalid";
        return result;
    }
    if (std::abs(residual_after.residual_seconds) > kEngineeringResidualToleranceSeconds) {
        result.invalid_reason = "residual_after_correction_too_large";
        return result;
    }

    result.valid = true;
    result.predicted_encounter_global_angle = residual_after.encounter_global_angle;
    result.target_arrival_true_anomaly = residual_after.target_arrival_true_anomaly;
    result.residual_scale_free = residual_after.residual_scale_free;
    result.residual_seconds = residual_after.residual_seconds;
    result.transfer_time_seconds = residual_after.transfer_time_seconds;
    result.target_time_seconds = residual_after.target_time_seconds;
    result.transfer_e = residual_after.transfer_e;
    result.transfer_p = residual_after.transfer_p;
    result.transfer_a = residual_after.transfer_a;
    result.theta_B = residual_after.theta_B;
    result.diagnostics.raw_residual_scale_free = residual_after.residual_scale_free;
    result.diagnostics.raw_residual_seconds = residual_after.residual_seconds;
    result.diagnostics.admissible_for_fast_approximation = true;
    result.diagnostics.admissibility_reason = "admissible";
    return result;
}

TestOnlyRouteBQueryResult
query_problem1_root_table_route_b_quadratic_one_step_corrected_test_only(
    const spaceship_cpp::problem1::Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    TestOnlyRouteBQueryResult query_result{};
    auto& result = query_result.safe_result;
    result.method = "route_b_quadratic_one_step_corrected_test_only";

    const auto nearest = problem1::find_nearest_problem1_root_table_node(table, query_nu_A, query_nu_B, query_theta_A);
    if (!nearest.valid || nearest.cell == nullptr) {
        result.fallback_required = true;
        result.reason = "failed_to_locate_nearest_node";
        return query_result;
    }

    const int cell_i = lower_corner_index(nearest.i, normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A));
    const int cell_j = lower_corner_index(nearest.j, normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B));
    const int cell_k = lower_corner_index(nearest.k, normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A));
    result.cell_admissibility = problem1::evaluate_problem1_root_cell_admissibility(
        table,
        cell_i,
        cell_j,
        cell_k,
        table.config().max_transfer_revolution,
        table.config().max_target_revolution);
    if (!result.cell_admissibility.admissible) {
        result.fallback_required = true;
        result.reason = "cell_non_admissible_branch_count_inconsistent";
        return query_result;
    }

    std::vector<problem1::Problem1RootApproximationResult> approximations;
    for (const auto& node_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
        if (!node_branch.valid) {
            continue;
        }
        auto approximation = evaluate_problem1_root_quadratic_one_step_corrected_test_only(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            node_branch,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            table.config().max_target_revolution);
        if (!approximation.valid) {
            query_result.invalid_reason_distribution[approximation.invalid_reason] += 1;
            continue;
        }
        bool duplicate = false;
        for (auto& existing : approximations) {
            if (existing.transfer_revolution != approximation.transfer_revolution) {
                continue;
            }
            const double time_tol = branch_time_duplicate_tolerance_seconds(
                existing.transfer_time_seconds, approximation.transfer_time_seconds);
            const double time_diff = std::abs(existing.transfer_time_seconds - approximation.transfer_time_seconds);
            const double angle_diff = wrapped_alpha_distance(
                existing.predicted_encounter_global_angle, approximation.predicted_encounter_global_angle);
            if (time_diff <= time_tol && angle_diff <= 1e-4) {
                if (std::abs(approximation.residual_seconds) < std::abs(existing.residual_seconds)) {
                    existing = approximation;
                }
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            approximations.push_back(approximation);
        }
    }

    std::sort(
        approximations.begin(),
        approximations.end(),
        [](const auto& lhs, const auto& rhs) {
            if (lhs.transfer_revolution != rhs.transfer_revolution) {
                return lhs.transfer_revolution < rhs.transfer_revolution;
            }
            if (lhs.transfer_time_seconds != rhs.transfer_time_seconds) {
                return lhs.transfer_time_seconds < rhs.transfer_time_seconds;
            }
            return lhs.predicted_encounter_global_angle < rhs.predicted_encounter_global_angle;
        });

    result.expected_count_by_k = result.cell_admissibility.reference_root_count_by_k;
    for (const auto& approximation : approximations) {
        result.candidate_count_by_k[approximation.transfer_revolution] += 1;
    }
    result.branch_count_complete = true;
    for (const auto& [k, expected] : result.expected_count_by_k) {
        const int actual = result.candidate_count_by_k[k];
        if (actual != expected) {
            result.branch_count_complete = false;
            result.missing_count_by_k[k] = std::max(0, expected - actual);
            result.extra_count_by_k[k] = std::max(0, actual - expected);
            result.incomplete_by_k[k] = true;
        } else {
            result.incomplete_by_k[k] = false;
        }
    }

    if (approximations.empty()) {
        result.valid = false;
        result.fallback_required = true;
        result.reason = "route_b_quadratic_one_step_no_valid_approximations";
        return query_result;
    }

    result.approximations = approximations;
    result.valid = true;
    if (!result.branch_count_complete) {
        result.fallback_required = true;
        result.reason = "branch_count_incomplete";
        return query_result;
    }

        result.reason = "route_b_quadratic_one_step_ok";
    return query_result;
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
    const int samples_per_group = 12;
    const double quadratic_hessian_step = 5e-4;
    const double tangent_residual_tolerance = problem1::kProblem1RootExperimentalTangentResidualTolerance;

    const auto samples = build_samples(samples_per_group);
    std::map<std::string, MethodStats> route_a_stats;
    std::map<std::string, MethodStats> route_b_linear_stats;
    std::map<std::string, MethodStats> route_b_quadratic_stats;
    std::map<std::string, MethodStats> route_b_quadratic_one_step_stats;

    for (const auto& sample : samples) {
        auto local_table = build_local_cell_table_for_query(
            departure_planet,
            target_planet,
            sample,
            max_transfer_revolution,
            max_target_revolution);

        double exact_elapsed = 0.0;
        const auto exact_branches = time_call([&] {
            return problem1::solve_problem1_from_departure_anomalies(
                departure_planet,
                target_planet,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                max_transfer_revolution,
                max_target_revolution);
        }, &exact_elapsed);

        double route_a_elapsed = 0.0;
        const auto route_a_branches = time_call([&] {
            return route_a_baseline_query(local_table, sample);
        }, &route_a_elapsed);
        auto& route_a = route_a_stats[sample.group_name];
        route_a.total_queries += 1;
        route_a.total_time_seconds += route_a_elapsed;
        route_a.success_count += 1;
        route_a.branch_count_complete_count += 1;
        record_branch_comparison(exact_branches, route_a_branches, &route_a);

        double route_b_linear_elapsed = 0.0;
        const auto route_b_linear_result = time_call([&] {
            return problem1::query_problem1_root_table_route_b_linear_safe(
                local_table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
        }, &route_b_linear_elapsed);
        auto& route_b_linear = route_b_linear_stats[sample.group_name];
        route_b_linear.total_queries += 1;
        route_b_linear.total_time_seconds += route_b_linear_elapsed;
        if (route_b_linear_result.valid) {
            route_b_linear.success_count += 1;
        }
        if (route_b_linear_result.fallback_required) {
            route_b_linear.fallback_required_count += 1;
        }
        if (route_b_linear_result.branch_count_complete) {
            route_b_linear.branch_count_complete_count += 1;
        }
        if (route_b_linear_result.reason == "route_b_linear_no_valid_approximations") {
            route_b_linear.no_valid_approximation_count += 1;
        }
        if (route_b_linear_result.reason == "cell_non_admissible_branch_count_inconsistent") {
            route_b_linear.non_admissible_count += 1;
        }
        if (!route_b_linear_result.valid || route_b_linear_result.fallback_required) {
            route_b_linear.invalid_reason_distribution[route_b_linear_result.reason] += 1;
        }
        std::vector<problem1::Problem1SolutionBranch> route_b_linear_branches;
        for (const auto& approximation : route_b_linear_result.approximations) {
            route_b_linear_branches.push_back(approximation_to_branch(approximation));
        }
        record_branch_comparison(exact_branches, route_b_linear_branches, &route_b_linear);

        double route_b_quadratic_elapsed = 0.0;
        const auto route_b_quadratic_result = time_call([&] {
            return problem1::query_problem1_root_table_route_b_quadratic_safe(
                local_table,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                quadratic_hessian_step,
                problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                tangent_residual_tolerance);
        }, &route_b_quadratic_elapsed);
        auto& route_b_quadratic = route_b_quadratic_stats[sample.group_name];
        route_b_quadratic.total_queries += 1;
        route_b_quadratic.total_time_seconds += route_b_quadratic_elapsed;
        if (route_b_quadratic_result.valid) {
            route_b_quadratic.success_count += 1;
        }
        if (route_b_quadratic_result.fallback_required) {
            route_b_quadratic.fallback_required_count += 1;
        }
        if (route_b_quadratic_result.branch_count_complete) {
            route_b_quadratic.branch_count_complete_count += 1;
        }
        if (route_b_quadratic_result.reason == "route_b_quadratic_no_valid_approximations") {
            route_b_quadratic.no_valid_approximation_count += 1;
        }
        if (route_b_quadratic_result.reason == "cell_non_admissible_branch_count_inconsistent") {
            route_b_quadratic.non_admissible_count += 1;
        }
        if (!route_b_quadratic_result.valid || route_b_quadratic_result.fallback_required) {
            route_b_quadratic.invalid_reason_distribution[route_b_quadratic_result.reason] += 1;
        }
        std::vector<problem1::Problem1SolutionBranch> route_b_quadratic_branches;
        for (const auto& approximation : route_b_quadratic_result.approximations) {
            route_b_quadratic_branches.push_back(approximation_to_branch(approximation));
        }
        record_branch_comparison(exact_branches, route_b_quadratic_branches, &route_b_quadratic);

        double route_b_quadratic_one_step_elapsed = 0.0;
        const auto route_b_quadratic_one_step_query = time_call([&] {
            return query_problem1_root_table_route_b_quadratic_one_step_corrected_test_only(
                local_table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
        }, &route_b_quadratic_one_step_elapsed);
        const auto& route_b_quadratic_one_step_result = route_b_quadratic_one_step_query.safe_result;
        auto& route_b_quadratic_one_step = route_b_quadratic_one_step_stats[sample.group_name];
        route_b_quadratic_one_step.total_queries += 1;
        route_b_quadratic_one_step.total_time_seconds += route_b_quadratic_one_step_elapsed;
        if (route_b_quadratic_one_step_result.valid) {
            route_b_quadratic_one_step.success_count += 1;
        }
        if (route_b_quadratic_one_step_result.fallback_required) {
            route_b_quadratic_one_step.fallback_required_count += 1;
        }
        if (route_b_quadratic_one_step_result.branch_count_complete) {
            route_b_quadratic_one_step.branch_count_complete_count += 1;
        }
        if (route_b_quadratic_one_step_result.reason == "route_b_quadratic_one_step_no_valid_approximations") {
            route_b_quadratic_one_step.no_valid_approximation_count += 1;
        }
        if (route_b_quadratic_one_step_result.reason == "cell_non_admissible_branch_count_inconsistent") {
            route_b_quadratic_one_step.non_admissible_count += 1;
        }
        if (!route_b_quadratic_one_step_result.valid || route_b_quadratic_one_step_result.fallback_required) {
            route_b_quadratic_one_step.invalid_reason_distribution[route_b_quadratic_one_step_result.reason] += 1;
        }
        for (const auto& [reason, count] : route_b_quadratic_one_step_query.invalid_reason_distribution) {
            route_b_quadratic_one_step.invalid_reason_distribution[reason] += count;
        }
        std::vector<problem1::Problem1SolutionBranch> route_b_quadratic_one_step_branches;
        for (const auto& approximation : route_b_quadratic_one_step_result.approximations) {
            route_b_quadratic_one_step_branches.push_back(approximation_to_branch(approximation));
        }
        record_branch_comparison(exact_branches, route_b_quadratic_one_step_branches, &route_b_quadratic_one_step);
    }

    for (const std::string group_name : {"near_node", "mid_cell", "physical_launch"}) {
        const double route_a_avg_ms = route_a_stats[group_name].total_queries > 0
            ? 1000.0 * route_a_stats[group_name].total_time_seconds / route_a_stats[group_name].total_queries
            : 0.0;
        print_method_summary(group_name, "route_a", route_a_stats[group_name]);
        print_method_summary(group_name, "route_b_linear", route_b_linear_stats[group_name], route_a_avg_ms);
        print_method_summary(group_name, "route_b_quadratic_no_newton", route_b_quadratic_stats[group_name], route_a_avg_ms);
        print_method_summary(
            group_name,
            "route_b_quadratic_one_step_corrected_test_only",
            route_b_quadratic_one_step_stats[group_name],
            route_a_avg_ms);
    }

    return 0;
}
