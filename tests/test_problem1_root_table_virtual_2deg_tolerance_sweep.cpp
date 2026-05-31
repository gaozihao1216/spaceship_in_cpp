#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kNearNodeOffsetDegrees = 0.1;
constexpr double kMidCellOffsetDegrees = 1.0;
constexpr double kStrictAlphaThreshold = 1e-5;
constexpr double kStrictTimeThresholdSeconds = 1e-1;
constexpr double kEngineeringResidualThresholdSeconds = 1e-2;
constexpr double kPhysicalTimeErrorThresholdSeconds = 1.0;
constexpr double kPhysicalAlphaErrorThresholdRadians = 1e-3;

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

struct QuerySample {
    std::string group_name;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

struct IndexedBranch {
    int original_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

enum class StopReason {
    ResidualToleranceSuccess,
    ResidualIncreaseStop,
    MaxIterationStop,
    DerivativeInvalidStop,
    ResidualInvalidStop,
    NonFiniteStepStop,
};

struct TestOnlyRefinementResult {
    bool valid = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
    StopReason stop_reason = StopReason::ResidualInvalidStop;
    double final_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
};

struct SourcePipelineRecord {
    int source_node_branch_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch source_node_branch;
    bool attach_derivative_ok = false;
    bool raw_prediction_finite = false;
    double raw_predicted_alpha = std::numeric_limits<double>::quiet_NaN();
    bool raw_residual_valid = false;
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    TestOnlyRefinementResult refined;
};

struct WrongRootCase {
    std::string group_name;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
    int k = 0;
    int q_exact = 0;
    int q_candidate = 0;
    double exact_time = 0.0;
    double candidate_time = 0.0;
    double time_error = 0.0;
    double exact_alpha = 0.0;
    double candidate_alpha = 0.0;
    double alpha_wrapped_error = 0.0;
    double residual_seconds = 0.0;
    StopReason stop_reason = StopReason::ResidualInvalidStop;
};

struct EngineeringStats {
    int exact_total_count = 0;
    int candidate_total_count = 0;
    int per_k_pair_count = 0;
    int engineering_coverage_count = 0;
    int strict_match_count = 0;
    int medium_match_count = 0;
    int loose_match_count = 0;
    int count_mismatch_count = 0;
    int missing_candidate_count = 0;
    int extra_candidate_count = 0;
    int physical_wrong_root_count = 0;
    int label_or_strict_mismatch_count = 0;
    int residual_increase_stop_count = 0;
    int residual_tolerance_success_count = 0;
    int max_iteration_stop_count = 0;
    int derivative_invalid_stop_count = 0;
    double max_residual_seconds = 0.0;
    double max_alpha_wrapped_error = 0.0;
    double max_time_error = 0.0;
    std::vector<WrongRootCase> wrong_root_cases;
};

struct ToleranceRunSummary {
    std::string group_name;
    double tolerance_seconds = 0.0;
    EngineeringStats engineering;
};

double wrapped_alpha_distance(double alpha1, double alpha2) {
    return std::abs(normalize_angle_minus_pi_pi(alpha1 - alpha2));
}

const char* stop_reason_name(StopReason reason) {
    switch (reason) {
        case StopReason::ResidualToleranceSuccess: return "residual_tolerance_success";
        case StopReason::ResidualIncreaseStop: return "residual_increase_stop";
        case StopReason::MaxIterationStop: return "max_iteration_stop";
        case StopReason::DerivativeInvalidStop: return "derivative_invalid_stop";
        case StopReason::ResidualInvalidStop: return "residual_invalid_stop";
        case StopReason::NonFiniteStepStop: return "nonfinite_step_stop";
    }
    return "unknown";
}

bool is_strict_match_1e6(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <=
            kStrictAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <=
            kStrictTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kEngineeringResidualThresholdSeconds;
}

bool is_engineering_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <=
            kPhysicalAlphaErrorThresholdRadians &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <=
            kPhysicalTimeErrorThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kEngineeringResidualThresholdSeconds;
}

bool is_medium_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <= 1e-3 &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= 1.0 &&
        std::abs(candidate.residual_seconds) <= 1e-3;
}

bool is_loose_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <= 1e-2 &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= 1e2 &&
        std::abs(candidate.residual_seconds) <= 1.0;
}

double branch_time_duplicate_tolerance_seconds(double t1, double t2) {
    return std::max(1e-4, 1e-10 * std::max(std::abs(t1), std::abs(t2)));
}

bool same_physical_branch_by_k_time(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return lhs.transfer_revolution == rhs.transfer_revolution &&
        std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds) <=
            branch_time_duplicate_tolerance_seconds(lhs.time_of_flight_seconds, rhs.time_of_flight_seconds);
}

bool safe_to_merge(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return same_physical_branch_by_k_time(lhs, rhs) &&
        wrapped_alpha_distance(lhs.encounter_global_angle, rhs.encounter_global_angle) <= 1e-4;
}

void add_if_not_duplicate(
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* branches,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    for (auto& existing : *branches) {
        if (!same_physical_branch_by_k_time(existing, branch)) {
            continue;
        }
        if (!safe_to_merge(existing, branch)) {
            continue;
        }
        if (std::abs(branch.residual_seconds) < std::abs(existing.residual_seconds)) {
            existing = branch;
        }
        return;
    }
    branches->push_back(branch);
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
    const double earth_period = spaceship_cpp::planet_params::planet_orbital_period(
        spaceship_cpp::planet_params::PlanetId::Earth);
    const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};
    for (int i = 0; i < samples_per_group; ++i) {
        const double launch_fraction = std::fmod(0.17 + 0.31 * static_cast<double>(i), 1.0);
        const double launch_time = launch_fraction * earth_period;
        const auto departure_state = spaceship_cpp::planet_params::planet_state_at_time(
            spaceship_cpp::planet_params::PlanetId::Earth, launch_time);
        const auto target_state = spaceship_cpp::planet_params::planet_state_at_time(
            spaceship_cpp::planet_params::PlanetId::Mars, launch_time);
        const double transfer_perihelion_angle =
            transfer_perihelion_angles[static_cast<std::size_t>(i) % transfer_perihelion_angles.size()];
        samples.push_back({"physical_launch",
                           departure_state.varphi,
                           target_state.varphi,
                           normalize_angle_0_2pi(departure_state.theta_global - transfer_perihelion_angle)});
    }
    return samples;
}

std::vector<IndexedBranch> collect_sorted_by_time(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::vector<IndexedBranch> result;
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (branches[i].valid) {
            result.push_back({static_cast<int>(i), branches[i]});
        }
    }
    std::sort(result.begin(), result.end(), [](const IndexedBranch& a, const IndexedBranch& b) {
        if (a.branch.time_of_flight_seconds < b.branch.time_of_flight_seconds) {
            return true;
        }
        if (a.branch.time_of_flight_seconds > b.branch.time_of_flight_seconds) {
            return false;
        }
        return a.branch.encounter_global_angle < b.branch.encounter_global_angle;
    });
    return result;
}

TestOnlyRefinementResult refine_problem1_root_branch_newton_residual_first_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_alpha,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance
) {
    namespace problem1 = spaceship_cpp::problem1;
    TestOnlyRefinementResult result{};
    double alpha_unwrapped = initial_alpha;
    double previous_abs_residual_scale_free = std::numeric_limits<double>::infinity();
    (void)alpha_step_tolerance;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const double alpha_wrapped = normalize_angle_0_2pi(alpha_unwrapped);
        const auto residual = problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            alpha_wrapped,
            transfer_revolution,
            target_revolution);
        if (!residual.valid || !std::isfinite(residual.residual_seconds) ||
            !std::isfinite(residual.residual_scale_free)) {
            result.stop_reason = StopReason::ResidualInvalidStop;
            result.final_residual_seconds = residual.residual_seconds;
            result.iterations = iteration + 1;
            return result;
        }
        result.final_residual_seconds = residual.residual_seconds;
        if (std::abs(residual.residual_seconds) <= residual_tolerance_seconds) {
            result.valid = true;
            result.stop_reason = StopReason::ResidualToleranceSuccess;
            result.iterations = iteration + 1;
            result.branch.valid = true;
            result.branch.encounter_global_angle = residual.encounter_global_angle;
            result.branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
            result.branch.transfer_revolution = transfer_revolution;
            result.branch.target_revolution = target_revolution;
            result.branch.time_of_flight_seconds = residual.transfer_time_seconds;
            result.branch.target_time_seconds = residual.target_time_seconds;
            result.branch.residual_seconds = residual.residual_seconds;
            result.branch.transfer_e = residual.transfer_e;
            result.branch.transfer_p = residual.transfer_p;
            result.branch.transfer_a = residual.transfer_a;
            result.branch.theta_B = residual.theta_B;
            return result;
        }
        const auto derivatives = problem1::evaluate_problem1_root_residual_derivatives_with_mode(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            alpha_wrapped,
            transfer_revolution,
            target_revolution,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (!derivatives.valid || !std::isfinite(derivatives.R_alpha) ||
            std::abs(derivatives.R_alpha) <= 1e-12) {
            result.stop_reason = StopReason::DerivativeInvalidStop;
            result.iterations = iteration + 1;
            return result;
        }
        const double delta_alpha = -residual.residual_scale_free / derivatives.R_alpha;
        if (!std::isfinite(delta_alpha)) {
            result.stop_reason = StopReason::NonFiniteStepStop;
            result.iterations = iteration + 1;
            return result;
        }
        const double trial_alpha_wrapped = normalize_angle_0_2pi(alpha_unwrapped + delta_alpha);
        const auto trial_residual = problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trial_alpha_wrapped,
            transfer_revolution,
            target_revolution);
        if (!trial_residual.valid || !std::isfinite(trial_residual.residual_scale_free)) {
            result.stop_reason = StopReason::ResidualInvalidStop;
            result.final_residual_seconds = trial_residual.residual_seconds;
            result.iterations = iteration + 1;
            return result;
        }
        const double current_abs = std::abs(residual.residual_scale_free);
        const double trial_abs = std::abs(trial_residual.residual_scale_free);
        if (trial_abs > current_abs || trial_abs > previous_abs_residual_scale_free) {
            result.stop_reason = StopReason::ResidualIncreaseStop;
            result.iterations = iteration + 1;
            result.final_residual_seconds = residual.residual_seconds;
            return result;
        }
        previous_abs_residual_scale_free = trial_abs;
        alpha_unwrapped += delta_alpha;
    }
    result.stop_reason = StopReason::MaxIterationStop;
    result.iterations = max_iterations;
    return result;
}

std::vector<SourcePipelineRecord> build_source_records(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& node_branches,
    double residual_tolerance_seconds
) {
    namespace problem1 = spaceship_cpp::problem1;
    std::vector<SourcePipelineRecord> records;
    for (std::size_t i = 0; i < node_branches.size(); ++i) {
        const auto& source = node_branches[i];
        if (!source.valid) {
            continue;
        }
        SourcePipelineRecord record{};
        record.source_node_branch_index = static_cast<int>(i);
        record.source_node_branch = source;
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            node.nu_A_node,
            node.nu_B_node,
            node.theta_A_node,
            source,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        record.attach_derivative_ok = attached.valid && attached.derivatives_available;
        if (!record.attach_derivative_ok) {
            records.push_back(record);
            continue;
        }
        record.raw_predicted_alpha = normalize_angle_0_2pi(
            attached.encounter_global_angle +
            attached.d_encounter_global_angle_d_nu_A * node.dnu_A +
            attached.d_encounter_global_angle_d_nu_B * node.dnu_B +
            attached.d_encounter_global_angle_d_theta_A * node.dtheta_A);
        record.raw_prediction_finite = std::isfinite(record.raw_predicted_alpha);
        if (!record.raw_prediction_finite) {
            records.push_back(record);
            continue;
        }
        const auto raw_residual = problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            record.raw_predicted_alpha,
            source.transfer_revolution,
            source.target_revolution);
        record.raw_residual_valid = raw_residual.valid;
        record.raw_residual_seconds = raw_residual.residual_seconds;
        record.refined = refine_problem1_root_branch_newton_residual_first_test_only(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            source.transfer_revolution,
            source.target_revolution,
            record.raw_predicted_alpha,
            80,
            residual_tolerance_seconds,
            1e-14);
        records.push_back(record);
    }
    return records;
}

std::map<int, std::vector<IndexedBranch>> group_by_k_sorted(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::map<int, std::vector<IndexedBranch>> grouped;
    const auto sorted = collect_sorted_by_time(branches);
    for (const auto& entry : sorted) {
        grouped[entry.branch.transfer_revolution].push_back(entry);
    }
    return grouped;
}

void collect_valid_candidates(
    const std::vector<SourcePipelineRecord>& records,
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* candidates
) {
    for (const auto& record : records) {
        if (!record.refined.valid) {
            continue;
        }
        add_if_not_duplicate(candidates, record.refined.branch);
    }
}

void accumulate_stop_counts(
    EngineeringStats* stats,
    const std::vector<SourcePipelineRecord>& records
) {
    for (const auto& record : records) {
        switch (record.refined.stop_reason) {
            case StopReason::ResidualToleranceSuccess:
                stats->residual_tolerance_success_count += 1;
                break;
            case StopReason::ResidualIncreaseStop:
                stats->residual_increase_stop_count += 1;
                break;
            case StopReason::MaxIterationStop:
                stats->max_iteration_stop_count += 1;
                break;
            case StopReason::DerivativeInvalidStop:
                stats->derivative_invalid_stop_count += 1;
                break;
            default:
                break;
        }
    }
}

void evaluate_engineering_policy(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& exact_branches,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& candidates,
    const std::vector<SourcePipelineRecord>& records,
    const QuerySample& sample,
    EngineeringStats* stats
) {
    EngineeringStats local{};
    for (const auto& b : exact_branches) {
        if (b.valid) {
            local.exact_total_count += 1;
        }
    }
    accumulate_stop_counts(&local, records);
    local.candidate_total_count = static_cast<int>(candidates.size());
    const auto exact_by_k = group_by_k_sorted(exact_branches);
    const auto candidate_by_k = group_by_k_sorted(candidates);
    for (const auto& [k, exact_group] : exact_by_k) {
        const auto it = candidate_by_k.find(k);
        const int candidate_count = it == candidate_by_k.end() ? 0 : static_cast<int>(it->second.size());
        const int exact_count = static_cast<int>(exact_group.size());
        if (candidate_count < exact_count) {
            local.missing_candidate_count += (exact_count - candidate_count);
        } else if (candidate_count > exact_count) {
            local.extra_candidate_count += (candidate_count - exact_count);
        }
        if (candidate_count != exact_count) {
            local.count_mismatch_count += 1;
        }
        const int pair_count = std::min(candidate_count, exact_count);
        for (int i = 0; i < pair_count; ++i) {
            const auto& exact = exact_group[static_cast<std::size_t>(i)].branch;
            const auto& candidate = it->second[static_cast<std::size_t>(i)].branch;
            local.per_k_pair_count += 1;
            const double alpha_err = wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle);
            const double time_err = std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds);
            const double residual_abs = std::abs(candidate.residual_seconds);
            local.max_residual_seconds = std::max(local.max_residual_seconds, residual_abs);
            local.max_alpha_wrapped_error = std::max(local.max_alpha_wrapped_error, alpha_err);
            local.max_time_error = std::max(local.max_time_error, time_err);
            if (is_strict_match_1e6(candidate, exact)) {
                local.strict_match_count += 1;
            }
            if (is_engineering_match(candidate, exact)) {
                local.engineering_coverage_count += 1;
            }
            if (is_medium_match(candidate, exact)) {
                local.medium_match_count += 1;
            }
            if (is_loose_match(candidate, exact)) {
                local.loose_match_count += 1;
            }
            if (alpha_err > kPhysicalAlphaErrorThresholdRadians || time_err > kPhysicalTimeErrorThresholdSeconds) {
                local.physical_wrong_root_count += 1;
                if (local.wrong_root_cases.size() < 5) {
                    int q_candidate = candidate.target_revolution;
                    StopReason stop_reason = StopReason::ResidualInvalidStop;
                    for (const auto& record : records) {
                        if (!record.refined.valid) {
                            continue;
                        }
                        if (record.refined.branch.transfer_revolution == candidate.transfer_revolution &&
                            std::abs(record.refined.branch.time_of_flight_seconds - candidate.time_of_flight_seconds) <= 1e-8 &&
                            wrapped_alpha_distance(record.refined.branch.encounter_global_angle, candidate.encounter_global_angle) <= 1e-10) {
                            stop_reason = record.refined.stop_reason;
                            break;
                        }
                    }
                    local.wrong_root_cases.push_back(
                        {sample.group_name,
                         sample.query_nu_A,
                         sample.query_nu_B,
                         sample.query_theta_A,
                         exact.transfer_revolution,
                         exact.target_revolution,
                         q_candidate,
                         exact.time_of_flight_seconds,
                         candidate.time_of_flight_seconds,
                         time_err,
                         exact.encounter_global_angle,
                         candidate.encounter_global_angle,
                         alpha_err,
                         candidate.residual_seconds,
                         stop_reason});
                }
            } else if (!is_strict_match_1e6(candidate, exact)) {
                local.label_or_strict_mismatch_count += 1;
            }
        }
    }

    stats->exact_total_count += local.exact_total_count;
    stats->candidate_total_count += local.candidate_total_count;
    stats->per_k_pair_count += local.per_k_pair_count;
    stats->engineering_coverage_count += local.engineering_coverage_count;
    stats->strict_match_count += local.strict_match_count;
    stats->medium_match_count += local.medium_match_count;
    stats->loose_match_count += local.loose_match_count;
    stats->count_mismatch_count += local.count_mismatch_count;
    stats->missing_candidate_count += local.missing_candidate_count;
    stats->extra_candidate_count += local.extra_candidate_count;
    stats->physical_wrong_root_count += local.physical_wrong_root_count;
    stats->label_or_strict_mismatch_count += local.label_or_strict_mismatch_count;
    stats->residual_increase_stop_count += local.residual_increase_stop_count;
    stats->residual_tolerance_success_count += local.residual_tolerance_success_count;
    stats->max_iteration_stop_count += local.max_iteration_stop_count;
    stats->derivative_invalid_stop_count += local.derivative_invalid_stop_count;
    stats->max_residual_seconds = std::max(stats->max_residual_seconds, local.max_residual_seconds);
    stats->max_alpha_wrapped_error = std::max(stats->max_alpha_wrapped_error, local.max_alpha_wrapped_error);
    stats->max_time_error = std::max(stats->max_time_error, local.max_time_error);
    for (const auto& item : local.wrong_root_cases) {
        if (stats->wrong_root_cases.size() < 5) {
            stats->wrong_root_cases.push_back(item);
        }
    }
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const auto samples = build_samples(12);
    const std::vector<double> tolerances{1e-5, 1e-3, 1e-2};

    std::map<double, std::map<std::string, ToleranceRunSummary>> summaries;

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "virtual_2deg_tolerance_sweep_test_only\n";

    for (double tolerance : tolerances) {
        for (const auto& sample : samples) {
            const auto node = find_nearest_virtual_root_table_node(
                sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
            auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
                departure_planet,
                target_planet,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                max_transfer_revolution,
                max_target_revolution);
            for (auto& branch : exact_branches) {
                if (!branch.valid) {
                    continue;
                }
                const auto polished = problem1::refine_problem1_root_branch_newton_seconds(
                    departure_planet,
                    target_planet,
                    sample.query_nu_A,
                    sample.query_nu_B,
                    sample.query_theta_A,
                    branch.transfer_revolution,
                    branch.target_revolution,
                    branch.encounter_global_angle,
                    30,
                    1e-6,
                    1e-14);
                if (polished.valid) {
                    branch = polished;
                }
            }
            const auto node_branches = problem1::solve_problem1_from_departure_anomalies(
                departure_planet,
                target_planet,
                node.nu_A_node,
                node.nu_B_node,
                node.theta_A_node,
                max_transfer_revolution,
                max_target_revolution);
            const auto records = build_source_records(
                departure_planet,
                target_planet,
                sample,
                node,
                node_branches,
                tolerance);
            std::vector<problem1::Problem1SolutionBranch> candidates;
            collect_valid_candidates(records, &candidates);

            auto& summary = summaries[tolerance][sample.group_name];
            summary.group_name = sample.group_name;
            summary.tolerance_seconds = tolerance;
            evaluate_engineering_policy(exact_branches, candidates, records, sample, &summary.engineering);

            if (std::abs(sample.query_nu_A - 4.555309e+00) < 1e-6 &&
                std::abs(sample.query_nu_B - 5.986479e+00) < 1e-6 &&
                std::abs(sample.query_theta_A - 2.949606e+00) < 1e-6 &&
                std::abs(tolerance - 1e-2) < 1e-12) {
                for (const auto& record : records) {
                    if (record.source_node_branch_index != 2 || !record.refined.valid) {
                        continue;
                    }
                    double best_time_error = std::numeric_limits<double>::infinity();
                    double best_alpha_error = std::numeric_limits<double>::infinity();
                    for (const auto& exact : exact_branches) {
                        if (!exact.valid || exact.transfer_revolution != record.refined.branch.transfer_revolution) {
                            continue;
                        }
                        const double terr = std::abs(exact.time_of_flight_seconds - record.refined.branch.time_of_flight_seconds);
                        if (terr < best_time_error) {
                            best_time_error = terr;
                            best_alpha_error = wrapped_alpha_distance(
                                exact.encounter_global_angle, record.refined.branch.encounter_global_angle);
                        }
                    }
                    std::cout << "fixed_branch2_case\n";
                    std::cout << "tolerance_seconds=" << tolerance << '\n';
                    std::cout << "source_branch_index=2\n";
                    std::cout << "final_residual_seconds=" << record.refined.final_residual_seconds << '\n';
                    std::cout << "matched_exact_time_error=" << best_time_error << '\n';
                    std::cout << "matched_exact_alpha_error=" << best_alpha_error << '\n';
                    std::cout << "stop_reason=" << stop_reason_name(record.refined.stop_reason) << '\n';
                }
            }
        }
    }

    // fixed branch 2 all tolerances
    const double query_nu_A = 4.555309e+00;
    const double query_nu_B = 5.986479e+00;
    const double query_theta_A = 2.949606e+00;
    const auto node = find_nearest_virtual_root_table_node(query_nu_A, query_nu_B, query_theta_A);
    const auto exact_fixed = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet, query_nu_A, query_nu_B, query_theta_A, max_transfer_revolution, max_target_revolution);
    const auto node_branches_fixed = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet, node.nu_A_node, node.nu_B_node, node.theta_A_node,
        max_transfer_revolution, max_target_revolution);
    std::vector<problem1::Problem1SolutionBranch> exact_k1;
    for (const auto& b : exact_fixed) {
        if (b.valid && b.transfer_revolution == 1) {
            exact_k1.push_back(b);
        }
    }
    std::sort(exact_k1.begin(), exact_k1.end(), [](const auto& a, const auto& b) {
        return a.time_of_flight_seconds < b.time_of_flight_seconds;
    });
    for (double tolerance : std::vector<double>{1e-5, 1e-3, 1e-2}) {
        const auto records = build_source_records(
            departure_planet, target_planet, {"mid_cell", query_nu_A, query_nu_B, query_theta_A}, node, node_branches_fixed, tolerance);
        for (const auto& record : records) {
            if (record.source_node_branch_index != 2) {
                continue;
            }
            double best_time_error = std::numeric_limits<double>::infinity();
            double best_alpha_error = std::numeric_limits<double>::infinity();
            if (record.refined.valid) {
                for (const auto& exact : exact_k1) {
                    const double terr = std::abs(exact.time_of_flight_seconds - record.refined.branch.time_of_flight_seconds);
                    if (terr < best_time_error) {
                        best_time_error = terr;
                        best_alpha_error = wrapped_alpha_distance(
                            exact.encounter_global_angle, record.refined.branch.encounter_global_angle);
                    }
                }
            }
            std::cout << "fixed_branch2_tolerance_run\n";
            std::cout << "tolerance_seconds=" << tolerance << '\n';
            std::cout << "valid=" << record.refined.valid << '\n';
            std::cout << "final_residual_seconds=" << record.refined.final_residual_seconds << '\n';
            std::cout << "final_time_error_to_exact=" << best_time_error << '\n';
            std::cout << "final_alpha_error_to_exact=" << best_alpha_error << '\n';
            std::cout << "stop_reason=" << stop_reason_name(record.refined.stop_reason) << '\n';
        }
    }

    for (const auto& [tolerance, groups] : summaries) {
        for (const auto& [group_name, summary] : groups) {
            std::cout << "tolerance_seconds=" << tolerance << '\n';
            std::cout << "group=" << group_name << '\n';
            const auto& s = summary.engineering;
            std::cout << "exact_total_count=" << s.exact_total_count << '\n';
            std::cout << "candidate_total_count=" << s.candidate_total_count << '\n';
            std::cout << "per_k_pair_count=" << s.per_k_pair_count << '\n';
            std::cout << "engineering_coverage_count=" << s.engineering_coverage_count << '\n';
            std::cout << "engineering_coverage_ratio=" <<
                (s.exact_total_count > 0 ? static_cast<double>(s.engineering_coverage_count) / s.exact_total_count : 0.0) << '\n';
            std::cout << "strict_match_count=" << s.strict_match_count << '\n';
            std::cout << "medium_match_count=" << s.medium_match_count << '\n';
            std::cout << "loose_match_count=" << s.loose_match_count << '\n';
            std::cout << "per_k_count_mismatch_count=" << s.count_mismatch_count << '\n';
            std::cout << "missing_candidate_count=" << s.missing_candidate_count << '\n';
            std::cout << "extra_candidate_count=" << s.extra_candidate_count << '\n';
            std::cout << "physical_wrong_root_count=" << s.physical_wrong_root_count << '\n';
            std::cout << "label_or_strict_mismatch_count=" << s.label_or_strict_mismatch_count << '\n';
            std::cout << "residual_increase_stop_count=" << s.residual_increase_stop_count << '\n';
            std::cout << "residual_tolerance_success_count=" << s.residual_tolerance_success_count << '\n';
            std::cout << "max_iteration_stop_count=" << s.max_iteration_stop_count << '\n';
            std::cout << "derivative_invalid_stop_count=" << s.derivative_invalid_stop_count << '\n';
            std::cout << "max_residual_seconds=" << s.max_residual_seconds << '\n';
            std::cout << "max_alpha_wrapped_error=" << s.max_alpha_wrapped_error << '\n';
            std::cout << "max_time_error_to_exact=" << s.max_time_error << '\n';
            for (std::size_t i = 0; i < s.wrong_root_cases.size(); ++i) {
                const auto& c = s.wrong_root_cases[i];
                std::cout << "wrong_root_case_rank=" << i << '\n';
                std::cout << "  group_name=" << c.group_name << '\n';
                std::cout << "  query_nu_A=" << c.query_nu_A << '\n';
                std::cout << "  query_nu_B=" << c.query_nu_B << '\n';
                std::cout << "  query_theta_A=" << c.query_theta_A << '\n';
                std::cout << "  k=" << c.k << '\n';
                std::cout << "  exact_time=" << c.exact_time << '\n';
                std::cout << "  candidate_time=" << c.candidate_time << '\n';
                std::cout << "  time_error=" << c.time_error << '\n';
                std::cout << "  exact_alpha=" << c.exact_alpha << '\n';
                std::cout << "  candidate_alpha=" << c.candidate_alpha << '\n';
                std::cout << "  alpha_wrapped_error=" << c.alpha_wrapped_error << '\n';
                std::cout << "  residual_seconds=" << c.residual_seconds << '\n';
                std::cout << "  q_exact=" << c.q_exact << '\n';
                std::cout << "  q_candidate=" << c.q_candidate << '\n';
                std::cout << "  stop_reason=" << stop_reason_name(c.stop_reason) << '\n';
            }
        }
    }

    return 0;
}
