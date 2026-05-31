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
#include <sstream>
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
constexpr std::array<double, 5> kLinearCoverageTimeThresholdsSeconds{{1.0, 1e2, 1e3, 1e4, 1e5}};
constexpr std::array<double, 3> kLinearCoverageAlphaThresholds{{1e-3, 1e-2, 1e-1}};

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
};

struct ThresholdCoverageStats {
    int total_exact_branch_count = 0;
    std::array<std::array<int, kLinearCoverageAlphaThresholds.size()>, kLinearCoverageTimeThresholdsSeconds.size()>
        matched_counts{};
};

struct QuadraticInvalidCase {
    std::string group_name;
    int sample_index = -1;
    int source_branch_index = -1;
    int transfer_revolution = 0;
    int source_target_revolution = 0;
    int selected_target_revolution = 0;
    double time_of_flight_seconds = std::numeric_limits<double>::quiet_NaN();
    std::string invalid_reason;
    std::string invalid_category;
    std::string hessian_method;
    bool hessian_valid = false;
    double tangent_residual_max_scale_free = std::numeric_limits<double>::quiet_NaN();
    double tangent_residual_max_seconds = std::numeric_limits<double>::quiet_NaN();
    double alpha_linear = std::numeric_limits<double>::quiet_NaN();
    double alpha_quadratic = std::numeric_limits<double>::quiet_NaN();
    bool q_sheet_selection_failed = false;
    bool residual_evaluated = false;
    bool residual_valid = false;
};

struct QuadraticDiagnosticStats {
    int total_quadratic_attempt_count = 0;
    int valid_quadratic_count = 0;
    int invalid_quadratic_count = 0;
    int query_non_admissible_count = 0;
    int query_branch_count_incomplete_count = 0;
    std::map<std::string, int> invalid_reason_distribution;
    std::map<std::string, int> invalid_category_distribution;
    std::vector<QuadraticInvalidCase> invalid_case_samples;
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

std::string threshold_key(double time_threshold_seconds, double alpha_threshold) {
    std::ostringstream oss;
    oss << "time<=" << std::scientific << time_threshold_seconds
        << ",alpha<=" << alpha_threshold;
    return oss.str();
}

std::string bool_text(bool value) {
    return value ? "true" : "false";
}

std::string quadratic_invalid_category(
    const spaceship_cpp::problem1::Problem1RootApproximationResult& approximation,
    bool residual_evaluated,
    bool residual_valid
) {
    if (!approximation.diagnostics.hessian_valid) {
        if (approximation.invalid_reason.find("tangent") != std::string::npos) {
            return "tangent_finite_difference_failed";
        }
        return "hessian_invalid";
    }
    if (!std::isfinite(approximation.diagnostics.alpha_quadratic)) {
        return "alpha_quadratic_became_invalid";
    }
    if (approximation.diagnostics.q_sheet_selection_failed) {
        return "q_sheet_selection_failed";
    }
    if (residual_evaluated && !residual_valid) {
        return "residual_evaluation_invalid";
    }
    if (approximation.invalid_reason.find("residual") != std::string::npos) {
        return "residual_evaluation_invalid";
    }
    return "other_invalid";
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
    MethodStats* stats,
    ThresholdCoverageStats* threshold_stats = nullptr
) {
    std::map<int, std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>> exact_by_k;
    std::map<int, std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>> approx_by_k;
    for (const auto& branch : exact) {
        if (branch.valid) {
            exact_by_k[branch.transfer_revolution].push_back(branch);
            stats->total_exact_branch_count += 1;
            if (threshold_stats != nullptr) {
                threshold_stats->total_exact_branch_count += 1;
            }
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
            if (threshold_stats != nullptr) {
                for (std::size_t time_i = 0; time_i < kLinearCoverageTimeThresholdsSeconds.size(); ++time_i) {
                    for (std::size_t alpha_i = 0; alpha_i < kLinearCoverageAlphaThresholds.size(); ++alpha_i) {
                        if (time_error <= kLinearCoverageTimeThresholdsSeconds[time_i] &&
                            alpha_error <= kLinearCoverageAlphaThresholds[alpha_i]) {
                            threshold_stats->matched_counts[time_i][alpha_i] += 1;
                        }
                    }
                }
            }
        }
    }
}

void record_quadratic_attempt_diagnostic(
    const QuerySample& sample,
    int sample_index,
    int source_branch_index,
    const spaceship_cpp::problem1::Problem1SolutionBranch& node_branch,
    const spaceship_cpp::problem1::Problem1RootApproximationResult& approximation,
    bool residual_evaluated,
    bool residual_valid,
    QuadraticDiagnosticStats* stats
) {
    stats->total_quadratic_attempt_count += 1;
    if (approximation.valid) {
        stats->valid_quadratic_count += 1;
        return;
    }

    stats->invalid_quadratic_count += 1;
    const std::string invalid_reason =
        approximation.invalid_reason.empty() ? "unknown_invalid_reason" : approximation.invalid_reason;
    const std::string invalid_category = quadratic_invalid_category(approximation, residual_evaluated, residual_valid);
    stats->invalid_reason_distribution[invalid_reason] += 1;
    stats->invalid_category_distribution[invalid_category] += 1;
    if (stats->invalid_case_samples.size() >= 3) {
        return;
    }

    QuadraticInvalidCase sample_case{};
    sample_case.group_name = sample.group_name;
    sample_case.sample_index = sample_index;
    sample_case.source_branch_index = source_branch_index;
    sample_case.transfer_revolution = node_branch.transfer_revolution;
    sample_case.source_target_revolution = node_branch.target_revolution;
    sample_case.selected_target_revolution = approximation.diagnostics.selected_target_revolution;
    sample_case.time_of_flight_seconds = node_branch.time_of_flight_seconds;
    sample_case.invalid_reason = invalid_reason;
    sample_case.invalid_category = invalid_category;
    sample_case.hessian_method = approximation.diagnostics.hessian_method;
    sample_case.hessian_valid = approximation.diagnostics.hessian_valid;
    sample_case.tangent_residual_max_scale_free = approximation.diagnostics.tangent_residual_max_scale_free;
    sample_case.tangent_residual_max_seconds = approximation.diagnostics.tangent_residual_max_seconds;
    sample_case.alpha_linear = approximation.diagnostics.alpha_linear;
    sample_case.alpha_quadratic = approximation.diagnostics.alpha_quadratic;
    sample_case.q_sheet_selection_failed = approximation.diagnostics.q_sheet_selection_failed;
    sample_case.residual_evaluated = residual_evaluated;
    sample_case.residual_valid = residual_valid;
    stats->invalid_case_samples.push_back(sample_case);
}

void print_method_summary(
    const std::string& group_name,
    const std::string& method_name,
    const MethodStats& stats,
    double route_a_avg_ms = 0.0
) {
    std::cout << "RouteBCompareSummary\n";
    std::cout << "group=" << group_name << '\n';
    std::cout << "method=" << method_name << '\n';
    std::cout << "total_queries=" << stats.total_queries << '\n';
    std::cout << "success_count=" << stats.success_count << '\n';
    std::cout << "fallback_required_count=" << stats.fallback_required_count << '\n';
    std::cout << "branch_count_complete_count=" << stats.branch_count_complete_count << '\n';
    std::cout << "no_valid_approximation_count=" << stats.no_valid_approximation_count << '\n';
    std::cout << "non_admissible_count=" << stats.non_admissible_count << '\n';
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

void print_linear_threshold_coverage(
    const std::string& group_name,
    const ThresholdCoverageStats& stats
) {
    std::cout << "RouteBLinearThresholdCoverage\n";
    std::cout << "group=" << group_name << '\n';
    std::cout << "total_exact_branch_count=" << stats.total_exact_branch_count << '\n';
    for (std::size_t time_i = 0; time_i < kLinearCoverageTimeThresholdsSeconds.size(); ++time_i) {
        for (std::size_t alpha_i = 0; alpha_i < kLinearCoverageAlphaThresholds.size(); ++alpha_i) {
            const double coverage = stats.total_exact_branch_count > 0
                ? static_cast<double>(stats.matched_counts[time_i][alpha_i]) / stats.total_exact_branch_count
                : 0.0;
            std::cout << threshold_key(
                kLinearCoverageTimeThresholdsSeconds[time_i],
                kLinearCoverageAlphaThresholds[alpha_i])
                << "_coverage=" << coverage << '\n';
        }
    }
}

void print_quadratic_diagnostic_summary(
    const std::string& group_name,
    const QuadraticDiagnosticStats& stats
) {
    std::cout << "RouteBQuadraticInvalidDiagnostic\n";
    std::cout << "group=" << group_name << '\n';
    std::cout << "total_quadratic_attempt_count=" << stats.total_quadratic_attempt_count << '\n';
    std::cout << "valid_quadratic_count=" << stats.valid_quadratic_count << '\n';
    std::cout << "invalid_quadratic_count=" << stats.invalid_quadratic_count << '\n';
    std::cout << "query_non_admissible_count=" << stats.query_non_admissible_count << '\n';
    std::cout << "query_branch_count_incomplete_count=" << stats.query_branch_count_incomplete_count << '\n';
    for (const auto& [reason, count] : stats.invalid_reason_distribution) {
        std::cout << "invalid_reason[" << reason << "]=" << count << '\n';
    }
    for (const auto& [category, count] : stats.invalid_category_distribution) {
        std::cout << "invalid_category[" << category << "]=" << count << '\n';
    }
    for (std::size_t i = 0; i < stats.invalid_case_samples.size(); ++i) {
        const auto& sample = stats.invalid_case_samples[i];
        std::cout << "invalid_case_rank=" << (i + 1) << '\n';
        std::cout << "invalid_case.group=" << sample.group_name << '\n';
        std::cout << "invalid_case.sample_index=" << sample.sample_index << '\n';
        std::cout << "invalid_case.source_branch_index=" << sample.source_branch_index << '\n';
        std::cout << "invalid_case.k=" << sample.transfer_revolution << '\n';
        std::cout << "invalid_case.source_q=" << sample.source_target_revolution << '\n';
        std::cout << "invalid_case.selected_q=" << sample.selected_target_revolution << '\n';
        std::cout << "invalid_case.time_of_flight=" << sample.time_of_flight_seconds << '\n';
        std::cout << "invalid_case.invalid_reason=" << sample.invalid_reason << '\n';
        std::cout << "invalid_case.invalid_category=" << sample.invalid_category << '\n';
        std::cout << "invalid_case.hessian_method=" << sample.hessian_method << '\n';
        std::cout << "invalid_case.hessian_valid=" << bool_text(sample.hessian_valid) << '\n';
        std::cout << "invalid_case.tangent_residual_max_scale_free=" << sample.tangent_residual_max_scale_free << '\n';
        std::cout << "invalid_case.tangent_residual_max_seconds=" << sample.tangent_residual_max_seconds << '\n';
        std::cout << "invalid_case.alpha_linear=" << sample.alpha_linear << '\n';
        std::cout << "invalid_case.alpha_quadratic=" << sample.alpha_quadratic << '\n';
        std::cout << "invalid_case.q_sheet_selection_failed=" << bool_text(sample.q_sheet_selection_failed) << '\n';
        std::cout << "invalid_case.residual.valid="
            << (sample.residual_evaluated ? bool_text(sample.residual_valid) : "not_evaluated") << '\n';
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
    const int samples_per_group = 12;
    const double hessian_step = 5e-4;
    const double tangent_residual_tolerance = problem1::kProblem1RootExperimentalTangentResidualTolerance;

    const auto samples = build_samples(samples_per_group);
    std::map<std::string, MethodStats> exact_stats;
    std::map<std::string, MethodStats> route_a_stats;
    std::map<std::string, MethodStats> route_b_linear_stats;
    std::map<std::string, MethodStats> route_b_quadratic_stats;
    std::map<std::string, ThresholdCoverageStats> route_b_linear_threshold_stats;
    std::map<std::string, QuadraticDiagnosticStats> route_b_quadratic_diagnostics;
    std::map<std::string, int> group_sample_indices;

    for (const auto& sample : samples) {
        const int sample_index = group_sample_indices[sample.group_name];
        group_sample_indices[sample.group_name] += 1;
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
        exact_stats[sample.group_name].total_queries += 1;
        exact_stats[sample.group_name].success_count += 1;
        exact_stats[sample.group_name].total_time_seconds += exact_elapsed;

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
                local_table,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A);
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
        std::vector<problem1::Problem1SolutionBranch> route_b_linear_branches;
        for (const auto& approximation : route_b_linear_result.approximations) {
            route_b_linear_branches.push_back(approximation_to_branch(approximation));
        }
        record_branch_comparison(
            exact_branches,
            route_b_linear_branches,
            &route_b_linear,
            &route_b_linear_threshold_stats[sample.group_name]);

        double route_b_quadratic_elapsed = 0.0;
        const auto route_b_quadratic_result = time_call([&] {
            return problem1::query_problem1_root_table_route_b_quadratic_safe(
                local_table,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                hessian_step,
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
        std::vector<problem1::Problem1SolutionBranch> route_b_quadratic_branches;
        for (const auto& approximation : route_b_quadratic_result.approximations) {
            route_b_quadratic_branches.push_back(approximation_to_branch(approximation));
        }
        record_branch_comparison(exact_branches, route_b_quadratic_branches, &route_b_quadratic);

        auto& quadratic_diagnostics = route_b_quadratic_diagnostics[sample.group_name];
        if (route_b_quadratic_result.reason == "cell_non_admissible_branch_count_inconsistent") {
            quadratic_diagnostics.query_non_admissible_count += 1;
        }
        if (route_b_quadratic_result.reason == "route_b_quadratic_incomplete_branch_count") {
            quadratic_diagnostics.query_branch_count_incomplete_count += 1;
        }

        const auto nearest = problem1::find_nearest_problem1_root_table_node(
            local_table,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A);
        if (!nearest.valid || nearest.cell == nullptr) {
            continue;
        }
        const int cell_i = lower_corner_index(
            nearest.i, normalize_angle_minus_pi_pi(sample.query_nu_A - nearest.node_nu_A));
        const int cell_j = lower_corner_index(
            nearest.j, normalize_angle_minus_pi_pi(sample.query_nu_B - nearest.node_nu_B));
        const int cell_k = lower_corner_index(
            nearest.k, normalize_angle_minus_pi_pi(sample.query_theta_A - nearest.node_theta_A));
        const auto admissibility = problem1::evaluate_problem1_root_cell_admissibility(
            local_table,
            cell_i,
            cell_j,
            cell_k,
            local_table.config().max_transfer_revolution,
            local_table.config().max_target_revolution);
        if (!admissibility.admissible) {
            continue;
        }
        for (std::size_t source_branch_index = 0;
             source_branch_index < nearest.cell->solutions_sorted_by_time_of_flight.size();
             ++source_branch_index) {
            const auto& node_branch = nearest.cell->solutions_sorted_by_time_of_flight[source_branch_index];
            if (!node_branch.valid) {
                continue;
            }
            const auto approximation = problem1::evaluate_problem1_root_quadratic_route_b_approximation_from_node(
                local_table.config().departure_planet,
                local_table.config().target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                local_table.config().max_target_revolution,
                hessian_step,
                problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                tangent_residual_tolerance);
            bool residual_evaluated = false;
            bool residual_valid = false;
            if (std::isfinite(approximation.diagnostics.alpha_quadratic)) {
                const int selected_q = approximation.diagnostics.q_sheet_selection_failed
                    ? node_branch.target_revolution
                    : approximation.diagnostics.selected_target_revolution;
                const auto residual = problem1::evaluate_problem1_root_residual(
                    local_table.config().departure_planet,
                    local_table.config().target_planet,
                    sample.query_nu_A,
                    sample.query_nu_B,
                    sample.query_theta_A,
                    approximation.diagnostics.alpha_quadratic,
                    node_branch.transfer_revolution,
                    selected_q);
                residual_evaluated = true;
                residual_valid = residual.valid;
            }
            record_quadratic_attempt_diagnostic(
                sample,
                sample_index,
                static_cast<int>(source_branch_index),
                node_branch,
                approximation,
                residual_evaluated,
                residual_valid,
                &quadratic_diagnostics);
        }
    }

    for (const std::string group_name : {"near_node", "mid_cell", "physical_launch"}) {
        const double route_a_avg_ms = route_a_stats[group_name].total_queries > 0
            ? 1000.0 * route_a_stats[group_name].total_time_seconds / route_a_stats[group_name].total_queries
            : 0.0;
        print_method_summary(group_name, "exact", exact_stats[group_name]);
        print_method_summary(group_name, "route_a", route_a_stats[group_name]);
        print_method_summary(group_name, "route_b_linear", route_b_linear_stats[group_name], route_a_avg_ms);
        print_method_summary(group_name, "route_b_quadratic", route_b_quadratic_stats[group_name], route_a_avg_ms);

        const double linear_coverage = route_b_linear_stats[group_name].total_exact_branch_count > 0
            ? static_cast<double>(route_b_linear_stats[group_name].matched_branch_count) /
                route_b_linear_stats[group_name].total_exact_branch_count
            : 0.0;
        const double quadratic_coverage = route_b_quadratic_stats[group_name].total_exact_branch_count > 0
            ? static_cast<double>(route_b_quadratic_stats[group_name].matched_branch_count) /
                route_b_quadratic_stats[group_name].total_exact_branch_count
            : 0.0;
        const double linear_avg_ms = route_b_linear_stats[group_name].total_queries > 0
            ? 1000.0 * route_b_linear_stats[group_name].total_time_seconds / route_b_linear_stats[group_name].total_queries
            : 0.0;
        const double quadratic_avg_ms = route_b_quadratic_stats[group_name].total_queries > 0
            ? 1000.0 * route_b_quadratic_stats[group_name].total_time_seconds / route_b_quadratic_stats[group_name].total_queries
            : 0.0;
        std::cout << "route_b_linear_usable_ratio=" <<
            (route_b_linear_stats[group_name].total_queries > 0
                ? static_cast<double>(route_b_linear_stats[group_name].success_count) /
                    route_b_linear_stats[group_name].total_queries
                : 0.0) << '\n';
        std::cout << "route_b_quadratic_usable_ratio=" <<
            (route_b_quadratic_stats[group_name].total_queries > 0
                ? static_cast<double>(route_b_quadratic_stats[group_name].success_count) /
                    route_b_quadratic_stats[group_name].total_queries
                : 0.0) << '\n';
        std::cout << "route_b_linear_speedup_vs_route_a=" <<
            (linear_avg_ms > 0.0 ? route_a_avg_ms / linear_avg_ms : 0.0) << '\n';
        std::cout << "route_b_quadratic_speedup_vs_route_a=" <<
            (quadratic_avg_ms > 0.0 ? route_a_avg_ms / quadratic_avg_ms : 0.0) << '\n';
        std::cout << "route_b_quadratic_accuracy_better_than_linear=" << (quadratic_coverage > linear_coverage) << '\n';
        print_quadratic_diagnostic_summary(group_name, route_b_quadratic_diagnostics[group_name]);
        print_linear_threshold_coverage(group_name, route_b_linear_threshold_stats[group_name]);
    }

    return 0;
}
