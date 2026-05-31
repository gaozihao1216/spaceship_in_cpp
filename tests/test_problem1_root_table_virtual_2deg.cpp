#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepDegrees = 2.0;
constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kNearNodeOffsetDegrees = 0.1;
constexpr double kMidCellOffsetDegrees = 1.0;
constexpr double kRawFarAlphaThreshold = 0.05;
constexpr double kRawFarTimeThresholdSeconds = 86400.0;
constexpr double kStrictAlphaThreshold = 1e-5;
constexpr double kStrictTimeThresholdSeconds = 1e-2;
constexpr double kStrictResidualThresholdSeconds1e6 = 1e-6;
constexpr double kStrictResidualThresholdSeconds1e5 = 1e-5;

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

struct VirtualNodeKey {
    spaceship_cpp::planet_params::PlanetId departure_planet = spaceship_cpp::planet_params::PlanetId::Earth;
    spaceship_cpp::planet_params::PlanetId target_planet = spaceship_cpp::planet_params::PlanetId::Mars;
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    int max_transfer_revolution = 0;
    int max_target_revolution = 0;

    bool operator<(const VirtualNodeKey& other) const {
        return std::tie(
                   departure_planet,
                   target_planet,
                   nu_A_index,
                   nu_B_index,
                   theta_A_index,
                   max_transfer_revolution,
                   max_target_revolution) <
            std::tie(
                   other.departure_planet,
                   other.target_planet,
                   other.nu_A_index,
                   other.nu_B_index,
                   other.theta_A_index,
                   other.max_transfer_revolution,
                   other.max_target_revolution);
    }
};

struct QuerySample {
    std::string group_name;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

struct TimedNodeSolve {
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> branches;
    double solve_time_seconds = 0.0;
};

struct BranchStats {
    int refined_candidate_count_before_dedup = 0;
    int refined_candidate_count_after_dedup = 0;
    int same_group_exists_but_no_strict_match_count = 0;
    int no_group_candidate_count = 0;
    int strict_matched_count = 0;
    int medium_matched_count = 0;
    int loose_matched_count = 0;
    int best_error_count = 0;
    double best_alpha_error_sum = 0.0;
    double best_alpha_error_max = 0.0;
    double best_time_error_sum = 0.0;
    double best_time_error_max = 0.0;
    double best_residual_sum = 0.0;
    double best_residual_max = 0.0;
    int alpha_fail_count = 0;
    int time_fail_count = 0;
    int residual_fail_count = 0;
    int alpha_and_time_fail_count = 0;
    int strict_fail_but_medium_pass_count = 0;
    int strict_fail_but_loose_pass_count = 0;
};

struct NearNodeOutlierCase {
    int exact_index = -1;
    int exact_k = 0;
    int exact_q = 0;
    double exact_alpha = 0.0;
    double exact_time_of_flight = 0.0;
    std::size_t best_candidate_index = 0;
    int best_candidate_k = 0;
    int best_candidate_q = 0;
    double best_candidate_alpha = 0.0;
    double best_candidate_time_of_flight = 0.0;
    double alpha_error = 0.0;
    double time_error = 0.0;
    double residual_seconds = 0.0;
    bool strict_match = false;
    bool medium_match = false;
    bool loose_match = false;
};

struct IndexedBranch {
    int original_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct OrderMatchStats {
    int pair_count = 0;
    int strict_match_count = 0;
    int strict_match_count_relaxed = 0;
    int medium_match_count = 0;
    int loose_match_count = 0;
    int count_mismatch_count = 0;
    int error_count = 0;
    double alpha_error_sum = 0.0;
    double alpha_error_max = 0.0;
    double time_error_sum = 0.0;
    double time_error_max = 0.0;
    double residual_sum = 0.0;
    double residual_max = 0.0;
};

enum class CandidateQuality {
    Primary,
    Rescue,
    Rejected,
};

struct TwoTierPolicyStats {
    int pair_count = 0;
    int strict_1e6_match_count = 0;
    int strict_1e5_match_count = 0;
    int medium_match_count = 0;
    int loose_match_count = 0;
    int count_mismatch_count = 0;
    int wrong_root_count = 0;
    int extra_candidate_count = 0;
    int missing_candidate_count = 0;
};

struct SourcePipelineRecord {
    int source_node_branch_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch source_node_branch;
    bool attach_derivative_failed = false;
    bool alpha_linear_finite = false;
    double raw_predicted_alpha = std::numeric_limits<double>::quiet_NaN();
    bool raw_residual_valid = false;
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1SolutionBranch raw_branch;
    int selected_q = -1;
    bool q_selection_failed = false;
    bool q_changed = false;
    double source_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
    double selected_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
    bool newton_valid = false;
    spaceship_cpp::problem1::Problem1SolutionBranch final_branch;
    spaceship_cpp::problem1::Problem1RootRefinementResult refinement_result;
    std::string invalid_reason;
    bool dedup_merged = false;
    std::string pipeline_classification;
};

struct EngineeringCoverageStats {
    int exact_total_count = 0;
    int candidate_total_count = 0;
    int per_k_pair_count = 0;
    int per_k_count_mismatch_count = 0;
    int engineering_coverage_count = 0;
    int physical_wrong_root_count = 0;
    int label_or_strict_mismatch_count = 0;
    int missing_candidate_count = 0;
    int extra_candidate_count = 0;
    int residual_tolerance_success_count = 0;
    int residual_increase_stop_count = 0;
    int max_iteration_stop_count = 0;
    int derivative_invalid_stop_count = 0;
    int derivative_attach_failed_after_convergence_count = 0;
    int q_sheet_selection_changed_count = 0;
    int q_sheet_selection_failed_count = 0;
    int q_sheet_selection_source_q_kept_count = 0;
    int q_sheet_selection_successful_changed_count = 0;
    double max_residual_seconds = 0.0;
    double max_time_error_to_exact = 0.0;
    double max_alpha_wrapped_error = 0.0;
};

struct FailureClassification {
    int no_same_kq_seed_at_node = 0;
    int same_kq_seed_exists_but_raw_far = 0;
    int newton_failed = 0;
    int newton_converged_wrong_root = 0;
    int newton_converged_correct_root = 0;
    int derivative_failed_count = 0;
    int fallback_used_count = 0;
};

struct ModeGroupSummary {
    int sample_count = 0;
    int exact_total_branch_count = 0;
    int multi_root_same_kq_group_count = 0;
    int suspicious_time_duplicate_with_far_alpha_count = 0;
    BranchStats raw_same_kq_stats;
    BranchStats newton_same_kq_stats;
    BranchStats newton_same_k_time_stats;
    OrderMatchStats global_time_order_stats;
    OrderMatchStats per_k_time_order_stats;
    FailureClassification failures;
    std::vector<NearNodeOutlierCase> near_node_outliers;
    bool first_time_order_failure_recorded = false;
    std::string first_time_order_failure_report;
    std::string per_k_count_tables_report;
    bool first_per_k_count_mismatch_recorded = false;
    std::string first_per_k_count_mismatch_report;
    std::string first_missing_root_microscope_report;
    TwoTierPolicyStats primary_only_policy_stats;
    TwoTierPolicyStats primary_plus_rescue_policy_stats;
    TwoTierPolicyStats selective_rescue_policy_stats;
    std::string branch2_two_tier_report;
    EngineeringCoverageStats engineering_stats;
};

struct Problem1BranchMatchResult {
    bool matched = false;
    bool same_group_exists = false;
    std::size_t candidate_index = 0;
    double alpha_error = std::numeric_limits<double>::quiet_NaN();
    double time_error_seconds = std::numeric_limits<double>::quiet_NaN();
    double residual_seconds = std::numeric_limits<double>::quiet_NaN();
};

double average_or_zero(double sum, int count) {
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

double wrapped_alpha_distance(double alpha1, double alpha2) {
    return std::abs(normalize_angle_minus_pi_pi(alpha1 - alpha2));
}

struct QSheetSelectionResult {
    int selected_q = -1;
    bool selection_failed = false;
    bool q_changed = false;
    double selected_continuity_error = std::numeric_limits<double>::quiet_NaN();
    double source_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
};

QSheetSelectionResult select_q_by_target_time_sheet_continuity_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int transfer_revolution,
    double alpha_linear,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    int max_target_revolution
) {
    QSheetSelectionResult result{};
    const double source_time_reference =
        std::isfinite(source_branch.target_time_seconds)
            ? source_branch.target_time_seconds
            : source_branch.time_of_flight_seconds;

    double best_error = std::numeric_limits<double>::infinity();
    int best_q = source_branch.target_revolution;
    bool any_valid = false;
    for (int q = 0; q <= max_target_revolution; ++q) {
        const auto residual = spaceship_cpp::problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_linear,
            transfer_revolution,
            q);
        if (!residual.valid) {
            continue;
        }
        const double continuity_error = std::abs(residual.target_time_seconds - source_time_reference);
        if (q == source_branch.target_revolution) {
            result.source_q_continuity_error = continuity_error;
        }
        if (!any_valid || continuity_error < best_error) {
            any_valid = true;
            best_error = continuity_error;
            best_q = q;
        }
    }
    if (!any_valid) {
        result.selection_failed = true;
        result.selected_q = source_branch.target_revolution;
        result.selected_continuity_error = std::numeric_limits<double>::infinity();
        result.q_changed = false;
        return result;
    }
    result.selected_q = best_q;
    result.selected_continuity_error = best_error;
    result.q_changed = best_q != source_branch.target_revolution;
    return result;
}

void update_max(double value, double* current_max) {
    if (std::isfinite(value) && value > *current_max) {
        *current_max = value;
    }
}

void record_best_match_stats(
    BranchStats* stats,
    const Problem1BranchMatchResult& match,
    double strict_alpha_threshold,
    double strict_time_threshold_seconds,
    double strict_residual_threshold_seconds,
    double medium_alpha_threshold,
    double medium_time_threshold_seconds,
    double medium_residual_threshold_seconds,
    double loose_alpha_threshold,
    double loose_time_threshold_seconds,
    double loose_residual_threshold_seconds
) {
    if (!match.same_group_exists) {
        stats->no_group_candidate_count += 1;
        return;
    }
    stats->best_error_count += 1;
    stats->best_alpha_error_sum += match.alpha_error;
    stats->best_time_error_sum += match.time_error_seconds;
    stats->best_residual_sum += std::abs(match.residual_seconds);
    update_max(match.alpha_error, &stats->best_alpha_error_max);
    update_max(match.time_error_seconds, &stats->best_time_error_max);
    update_max(std::abs(match.residual_seconds), &stats->best_residual_max);

    const bool strict =
        match.alpha_error <= strict_alpha_threshold &&
        match.time_error_seconds <= strict_time_threshold_seconds &&
        std::abs(match.residual_seconds) <= strict_residual_threshold_seconds;
    const bool medium =
        match.alpha_error <= medium_alpha_threshold &&
        match.time_error_seconds <= medium_time_threshold_seconds &&
        std::abs(match.residual_seconds) <= medium_residual_threshold_seconds;
    const bool loose =
        match.alpha_error <= loose_alpha_threshold &&
        match.time_error_seconds <= loose_time_threshold_seconds &&
        std::abs(match.residual_seconds) <= loose_residual_threshold_seconds;

    if (strict) {
        stats->strict_matched_count += 1;
    } else {
        stats->same_group_exists_but_no_strict_match_count += 1;
        const bool alpha_fail = match.alpha_error > strict_alpha_threshold;
        const bool time_fail = match.time_error_seconds > strict_time_threshold_seconds;
        const bool residual_fail = std::abs(match.residual_seconds) > strict_residual_threshold_seconds;
        if (alpha_fail) {
            stats->alpha_fail_count += 1;
        }
        if (time_fail) {
            stats->time_fail_count += 1;
        }
        if (residual_fail) {
            stats->residual_fail_count += 1;
        }
        if (alpha_fail && time_fail) {
            stats->alpha_and_time_fail_count += 1;
        }
    }
    if (medium) {
        stats->medium_matched_count += 1;
        if (!strict) {
            stats->strict_fail_but_medium_pass_count += 1;
        }
    }
    if (loose) {
        stats->loose_matched_count += 1;
        if (!strict) {
            stats->strict_fail_but_loose_pass_count += 1;
        }
    }
}

int sample_count_from_environment() {
    const char* env = std::getenv("ROOT_TABLE_VIRTUAL_2DEG_SAMPLE_COUNT");
    if (env != nullptr) {
        const int parsed = std::atoi(env);
        if (parsed > 0) {
            return parsed;
        }
    }
    const char* expensive_env = std::getenv("RUN_EXPENSIVE_ROOT_TABLE_VIRTUAL_2DEG");
    if (expensive_env != nullptr && std::string(expensive_env) == "1") {
        return 12;
    }
    return 2;
}

Problem1RootVirtualGridNode find_nearest_virtual_root_table_node(
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double step_radians = kVirtualGridStepRadians
) {
    const int axis_count = static_cast<int>(std::llround(kTwoPi / step_radians));
    const auto nearest_index = [&](double angle) {
        long long index = std::llround(normalize_angle_0_2pi(angle) / step_radians);
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
    node.nu_A_node = normalize_angle_0_2pi(static_cast<double>(node.nu_A_index) * step_radians);
    node.nu_B_node = normalize_angle_0_2pi(static_cast<double>(node.nu_B_index) * step_radians);
    node.theta_A_node = normalize_angle_0_2pi(static_cast<double>(node.theta_A_index) * step_radians);
    node.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

Problem1BranchMatchResult find_best_refined_branch_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact_branch,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& refined_candidates,
    double alpha_threshold,
    double time_threshold_seconds,
    double residual_threshold_seconds
) {
    Problem1BranchMatchResult result{};
    double best_score = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < refined_candidates.size(); ++i) {
        const auto& candidate = refined_candidates[i];
        if (!candidate.valid ||
            candidate.transfer_revolution != exact_branch.transfer_revolution ||
            candidate.target_revolution != exact_branch.target_revolution) {
            continue;
        }
        result.same_group_exists = true;
        const double aerr = wrapped_alpha_distance(
            candidate.encounter_global_angle, exact_branch.encounter_global_angle);
        const double terr = std::abs(candidate.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
        const double rerr = std::abs(candidate.residual_seconds);
        const double score =
            aerr / alpha_threshold +
            terr / time_threshold_seconds +
            rerr / residual_threshold_seconds;
        if (score < best_score) {
            best_score = score;
            result.candidate_index = i;
            result.alpha_error = aerr;
            result.time_error_seconds = terr;
            result.residual_seconds = candidate.residual_seconds;
        }
    }
    if (!result.same_group_exists) {
        return result;
    }
    result.matched =
        result.alpha_error <= alpha_threshold &&
        result.time_error_seconds <= time_threshold_seconds &&
        std::abs(result.residual_seconds) <= residual_threshold_seconds;
    return result;
}

Problem1BranchMatchResult find_best_refined_branch_match_by_k_time(
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact_branch,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& refined_candidates,
    double alpha_threshold,
    double time_threshold_seconds,
    double residual_threshold_seconds
) {
    Problem1BranchMatchResult result{};
    double best_score = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < refined_candidates.size(); ++i) {
        const auto& candidate = refined_candidates[i];
        if (!candidate.valid ||
            candidate.transfer_revolution != exact_branch.transfer_revolution) {
            continue;
        }
        result.same_group_exists = true;
        const double aerr = wrapped_alpha_distance(
            candidate.encounter_global_angle, exact_branch.encounter_global_angle);
        const double terr = std::abs(candidate.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
        const double rerr = std::abs(candidate.residual_seconds);
        const double score =
            terr / time_threshold_seconds +
            aerr / alpha_threshold +
            rerr / residual_threshold_seconds;
        if (score < best_score) {
            best_score = score;
            result.candidate_index = i;
            result.alpha_error = aerr;
            result.time_error_seconds = terr;
            result.residual_seconds = candidate.residual_seconds;
        }
    }
    if (!result.same_group_exists) {
        return result;
    }
    result.matched =
        result.alpha_error <= alpha_threshold &&
        result.time_error_seconds <= time_threshold_seconds &&
        std::abs(result.residual_seconds) <= residual_threshold_seconds;
    return result;
}

double alpha_error(const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
                   const spaceship_cpp::problem1::Problem1SolutionBranch& rhs) {
    return wrapped_alpha_distance(lhs.encounter_global_angle, rhs.encounter_global_angle);
}

bool is_strict_match(const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
                     const spaceship_cpp::problem1::Problem1SolutionBranch& exact) {
    return alpha_error(candidate, exact) <= kStrictAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= kStrictTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kStrictResidualThresholdSeconds1e6;
}

bool is_strict_match_relaxed_residual(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return alpha_error(candidate, exact) <= kStrictAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= kStrictTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kStrictResidualThresholdSeconds1e5;
}

CandidateQuality classify_candidate_quality(const spaceship_cpp::problem1::Problem1SolutionBranch& branch) {
    const double abs_residual = std::abs(branch.residual_seconds);
    if (abs_residual <= kStrictResidualThresholdSeconds1e6) {
        return CandidateQuality::Primary;
    }
    if (abs_residual <= kStrictResidualThresholdSeconds1e5) {
        return CandidateQuality::Rescue;
    }
    return CandidateQuality::Rejected;
}

bool is_medium_match(const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
                     const spaceship_cpp::problem1::Problem1SolutionBranch& exact) {
    return alpha_error(candidate, exact) <= 1e-3 &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= 1e3 &&
        std::abs(candidate.residual_seconds) <= 1e-4;
}

bool is_loose_match(const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
                    const spaceship_cpp::problem1::Problem1SolutionBranch& exact) {
    return alpha_error(candidate, exact) <= 1e-2 &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= 1e5 &&
        std::abs(candidate.residual_seconds) <= 1e-2;
}

bool is_engineering_match(const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
                          const spaceship_cpp::problem1::Problem1SolutionBranch& exact) {
    return candidate.valid &&
        std::abs(candidate.residual_seconds) <= 1e-2 &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= 1.0 &&
        alpha_error(candidate, exact) <= 1e-3;
}

std::vector<IndexedBranch> collect_valid_sorted_by_time(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::vector<IndexedBranch> result;
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (!branches[i].valid) {
            continue;
        }
        result.push_back({static_cast<int>(i), branches[i]});
    }
    std::sort(
        result.begin(),
        result.end(),
        [](const IndexedBranch& lhs, const IndexedBranch& rhs) {
            if (lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds) {
                return true;
            }
            if (lhs.branch.time_of_flight_seconds > rhs.branch.time_of_flight_seconds) {
                return false;
            }
            return lhs.branch.encounter_global_angle < rhs.branch.encounter_global_angle;
        });
    return result;
}

void record_order_match(
    OrderMatchStats* stats,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact_branch,
    const spaceship_cpp::problem1::Problem1SolutionBranch& refined_branch
) {
    stats->pair_count += 1;
    const double aerr = alpha_error(refined_branch, exact_branch);
    const double terr = std::abs(refined_branch.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
    const double rerr = std::abs(refined_branch.residual_seconds);
    stats->error_count += 1;
    stats->alpha_error_sum += aerr;
    stats->time_error_sum += terr;
    stats->residual_sum += rerr;
    update_max(aerr, &stats->alpha_error_max);
    update_max(terr, &stats->time_error_max);
    update_max(rerr, &stats->residual_max);
    if (is_strict_match(refined_branch, exact_branch)) {
        stats->strict_match_count += 1;
    }
    if (is_strict_match_relaxed_residual(refined_branch, exact_branch)) {
        stats->strict_match_count_relaxed += 1;
    }
    if (is_medium_match(refined_branch, exact_branch)) {
        stats->medium_match_count += 1;
    }
    if (is_loose_match(refined_branch, exact_branch)) {
        stats->loose_match_count += 1;
    }
}

void record_two_tier_match(
    TwoTierPolicyStats* stats,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact_branch,
    const spaceship_cpp::problem1::Problem1SolutionBranch& refined_branch
) {
    stats->pair_count += 1;
    if (is_strict_match(refined_branch, exact_branch)) {
        stats->strict_1e6_match_count += 1;
    }
    if (is_strict_match_relaxed_residual(refined_branch, exact_branch)) {
        stats->strict_1e5_match_count += 1;
    }
    if (is_medium_match(refined_branch, exact_branch)) {
        stats->medium_match_count += 1;
    } else if (std::abs(refined_branch.residual_seconds) <= kStrictResidualThresholdSeconds1e5) {
        stats->wrong_root_count += 1;
    }
    if (is_loose_match(refined_branch, exact_branch)) {
        stats->loose_match_count += 1;
    }
}

void evaluate_per_k_two_tier_policy(
    const std::map<int, std::vector<IndexedBranch>>& exact_by_k,
    const std::map<int, std::vector<IndexedBranch>>& candidates_by_k,
    TwoTierPolicyStats* stats
) {
    for (const auto& [k, exact_group] : exact_by_k) {
        const auto it = candidates_by_k.find(k);
        const int candidate_count = it == candidates_by_k.end() ? 0 : static_cast<int>(it->second.size());
        const int exact_count = static_cast<int>(exact_group.size());
        if (candidate_count < exact_count) {
            stats->missing_candidate_count += (exact_count - candidate_count);
        } else if (candidate_count > exact_count) {
            stats->extra_candidate_count += (candidate_count - exact_count);
        }
        if (candidate_count != exact_count) {
            stats->count_mismatch_count += 1;
            continue;
        }
        for (std::size_t i = 0; i < exact_group.size(); ++i) {
            record_two_tier_match(stats, exact_group[i].branch, it->second[i].branch);
        }
    }
}

std::string build_time_order_failure_report(
    const std::string& group_name,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    const std::vector<IndexedBranch>& exact_sorted,
    const std::vector<IndexedBranch>& refined_sorted
) {
    std::ostringstream os;
    os << std::setprecision(6) << std::scientific;
    os << "time_order_failure_case\n";
    os << "group_name=" << group_name << '\n';
    os << "query_nu_A=" << sample.query_nu_A << '\n';
    os << "query_nu_B=" << sample.query_nu_B << '\n';
    os << "query_theta_A=" << sample.query_theta_A << '\n';
    os << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
    os << "node_offset_degrees=("
       << (node.dnu_A * 180.0 / kPi) << ","
       << (node.dnu_B * 180.0 / kPi) << ","
       << (node.dtheta_A * 180.0 / kPi) << ")\n";
    os << "exact_count=" << exact_sorted.size() << '\n';
    os << "refined_count=" << refined_sorted.size() << '\n';
    os << "exact_sorted_by_time\n";
    for (std::size_t i = 0; i < exact_sorted.size(); ++i) {
        const auto& entry = exact_sorted[i];
        os << "  exact_index_sorted=" << i
           << " original_index=" << entry.original_index
           << " k=" << entry.branch.transfer_revolution
           << " q=" << entry.branch.target_revolution
           << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
           << " alpha=" << entry.branch.encounter_global_angle
           << " residual_seconds=" << entry.branch.residual_seconds << '\n';
    }
    os << "refined_sorted_by_time\n";
    for (std::size_t i = 0; i < refined_sorted.size(); ++i) {
        const auto& entry = refined_sorted[i];
        os << "  refined_index_sorted=" << i
           << " original_index=" << entry.original_index
           << " k=" << entry.branch.transfer_revolution
           << " q=" << entry.branch.target_revolution
           << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
           << " alpha=" << entry.branch.encounter_global_angle
           << " residual_seconds=" << entry.branch.residual_seconds << '\n';
    }
    os << "paired_comparison\n";
    const std::size_t pair_count = std::min(exact_sorted.size(), refined_sorted.size());
    for (std::size_t i = 0; i < pair_count; ++i) {
        const auto& exact = exact_sorted[i].branch;
        const auto& refined = refined_sorted[i].branch;
        os << "  pair_index=" << i
           << " exact_k=" << exact.transfer_revolution
           << " exact_q=" << exact.target_revolution
           << " refined_k=" << refined.transfer_revolution
           << " refined_q=" << refined.target_revolution
           << " exact_time=" << exact.time_of_flight_seconds
           << " refined_time=" << refined.time_of_flight_seconds
           << " time_error=" << std::abs(refined.time_of_flight_seconds - exact.time_of_flight_seconds)
           << " alpha_wrapped_error=" << wrapped_alpha_distance(
               refined.encounter_global_angle, exact.encounter_global_angle)
           << " residual_seconds=" << refined.residual_seconds
           << " strict_match=" << is_strict_match(refined, exact)
           << " medium_match=" << is_medium_match(refined, exact)
           << " loose_match=" << is_loose_match(refined, exact) << '\n';
    }
    return os.str();
}

std::string build_per_k_count_table_report(
    const std::string& group_name,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    const std::map<int, std::vector<IndexedBranch>>& exact_by_k,
    const std::map<int, std::vector<IndexedBranch>>& refined_before_dedup_by_k,
    const std::map<int, std::vector<IndexedBranch>>& refined_after_dedup_by_k
) {
    std::ostringstream os;
    os << std::setprecision(6) << std::scientific;
    os << "per_k_count_table\n";
    os << "group_name=" << group_name << '\n';
    os << "query_nu_A=" << sample.query_nu_A << '\n';
    os << "query_nu_B=" << sample.query_nu_B << '\n';
    os << "query_theta_A=" << sample.query_theta_A << '\n';
    os << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
    os << "node_offset_degrees=("
       << (node.dnu_A * 180.0 / kPi) << ","
       << (node.dnu_B * 180.0 / kPi) << ","
       << (node.dtheta_A * 180.0 / kPi) << ")\n";
    std::map<int, bool> all_k;
    for (const auto& [k, _] : exact_by_k) {
        all_k[k] = true;
    }
    for (const auto& [k, _] : refined_before_dedup_by_k) {
        all_k[k] = true;
    }
    for (const auto& [k, _] : refined_after_dedup_by_k) {
        all_k[k] = true;
    }
    for (const auto& [k, _] : all_k) {
        const auto exact_it = exact_by_k.find(k);
        const auto before_it = refined_before_dedup_by_k.find(k);
        const auto after_it = refined_after_dedup_by_k.find(k);
        const int exact_count = exact_it == exact_by_k.end() ? 0 : static_cast<int>(exact_it->second.size());
        const int before_count = before_it == refined_before_dedup_by_k.end() ? 0 : static_cast<int>(before_it->second.size());
        const int after_count = after_it == refined_after_dedup_by_k.end() ? 0 : static_cast<int>(after_it->second.size());
        os << "  k=" << k
           << " exact_count_by_k=" << exact_count
           << " refined_count_before_dedup_by_k=" << before_count
           << " refined_count_after_dedup_by_k=" << after_count << '\n';
    }
    return os.str();
}

std::string build_per_k_count_mismatch_report(
    const std::string& group_name,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    int mismatch_k,
    const std::vector<IndexedBranch>& exact_in_k,
    const std::vector<IndexedBranch>& refined_before_dedup_in_k,
    const std::vector<IndexedBranch>& refined_after_dedup_in_k
) {
    std::ostringstream os;
    os << std::setprecision(6) << std::scientific;
    os << "per_k_count_mismatch_case\n";
    os << "group_name=" << group_name << '\n';
    os << "query_nu_A=" << sample.query_nu_A << '\n';
    os << "query_nu_B=" << sample.query_nu_B << '\n';
    os << "query_theta_A=" << sample.query_theta_A << '\n';
    os << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
    os << "node_offset_degrees=("
       << (node.dnu_A * 180.0 / kPi) << ","
       << (node.dnu_B * 180.0 / kPi) << ","
       << (node.dtheta_A * 180.0 / kPi) << ")\n";
    os << "mismatch_k=" << mismatch_k << '\n';
    os << "exact_count_in_k=" << exact_in_k.size() << '\n';
    os << "refined_count_before_dedup_in_k=" << refined_before_dedup_in_k.size() << '\n';
    os << "refined_count_after_dedup_in_k=" << refined_after_dedup_in_k.size() << '\n';
    os << "exact_branches_in_k_sorted_by_time\n";
    for (const auto& entry : exact_in_k) {
        os << "  original_index=" << entry.original_index
           << " k=" << entry.branch.transfer_revolution
           << " q=" << entry.branch.target_revolution
           << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
           << " alpha=" << entry.branch.encounter_global_angle
           << " residual_seconds=" << entry.branch.residual_seconds << '\n';
    }
    os << "exact_adjacent_time_gaps\n";
    for (std::size_t i = 1; i < exact_in_k.size(); ++i) {
        os << "  gap_index=" << (i - 1)
           << " gap_seconds=" << (exact_in_k[i].branch.time_of_flight_seconds -
                                  exact_in_k[i - 1].branch.time_of_flight_seconds) << '\n';
    }
    os << "refined_branches_before_dedup_in_k_sorted_by_time\n";
    for (const auto& entry : refined_before_dedup_in_k) {
        os << "  original_index=" << entry.original_index
           << " k=" << entry.branch.transfer_revolution
           << " q=" << entry.branch.target_revolution
           << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
           << " alpha=" << entry.branch.encounter_global_angle
           << " residual_seconds=" << entry.branch.residual_seconds << '\n';
    }
    os << "refined_before_dedup_adjacent_time_gaps\n";
    for (std::size_t i = 1; i < refined_before_dedup_in_k.size(); ++i) {
        os << "  gap_index=" << (i - 1)
           << " gap_seconds=" << (refined_before_dedup_in_k[i].branch.time_of_flight_seconds -
                                  refined_before_dedup_in_k[i - 1].branch.time_of_flight_seconds) << '\n';
    }
    os << "refined_branches_after_dedup_in_k_sorted_by_time\n";
    for (const auto& entry : refined_after_dedup_in_k) {
        os << "  original_index=" << entry.original_index
           << " k=" << entry.branch.transfer_revolution
           << " q=" << entry.branch.target_revolution
           << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
           << " alpha=" << entry.branch.encounter_global_angle
           << " residual_seconds=" << entry.branch.residual_seconds << '\n';
    }
    os << "refined_after_dedup_adjacent_time_gaps\n";
    for (std::size_t i = 1; i < refined_after_dedup_in_k.size(); ++i) {
        os << "  gap_index=" << (i - 1)
           << " gap_seconds=" << (refined_after_dedup_in_k[i].branch.time_of_flight_seconds -
                                  refined_after_dedup_in_k[i - 1].branch.time_of_flight_seconds) << '\n';
    }
    if (refined_before_dedup_in_k.size() == exact_in_k.size() &&
        refined_after_dedup_in_k.size() != exact_in_k.size()) {
        os << "mismatch_source=dedup_removed_a_candidate\n";
    } else if (refined_before_dedup_in_k.size() < exact_in_k.size()) {
        os << "mismatch_source=route_a_failed_to_generate_some_root\n";
    } else if (refined_before_dedup_in_k.size() > exact_in_k.size()) {
        os << "mismatch_source=route_a_generated_duplicate_or_extra_roots\n";
    } else {
        os << "mismatch_source=ordering_or_threshold_issue\n";
    }
    return os.str();
}

std::vector<IndexedBranch> collect_sorted_k_group(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    int k
) {
    std::vector<IndexedBranch> result;
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (!branches[i].valid || branches[i].transfer_revolution != k) {
            continue;
        }
        result.push_back({static_cast<int>(i), branches[i]});
    }
    std::sort(
        result.begin(),
        result.end(),
        [](const IndexedBranch& lhs, const IndexedBranch& rhs) {
            if (lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds) {
                return true;
            }
            if (lhs.branch.time_of_flight_seconds > rhs.branch.time_of_flight_seconds) {
                return false;
            }
            return lhs.branch.encounter_global_angle < rhs.branch.encounter_global_angle;
        });
    return result;
}

double branch_time_duplicate_tolerance_seconds(
    double time_1,
    double time_2
) {
    return std::max(1e-4, 1e-10 * std::max(std::abs(time_1), std::abs(time_2)));
}

bool are_same_physical_branch_by_k_time(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return lhs.transfer_revolution == rhs.transfer_revolution &&
        std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds) <=
        branch_time_duplicate_tolerance_seconds(lhs.time_of_flight_seconds, rhs.time_of_flight_seconds);
}

bool safe_to_merge_branch(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs,
    double alpha_merge_tolerance = 1e-4
) {
    return are_same_physical_branch_by_k_time(lhs, rhs) &&
        wrapped_alpha_distance(lhs.encounter_global_angle, rhs.encounter_global_angle) <= alpha_merge_tolerance;
}

bool add_refined_branch_if_not_duplicate(
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* branches,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch,
    int* suspicious_time_duplicate_with_far_alpha_count
) {
    constexpr double kAlphaWarningTolerance = 1e-4;
    for (auto& existing : *branches) {
        if (!are_same_physical_branch_by_k_time(existing, branch)) {
            continue;
        }
        if (!safe_to_merge_branch(existing, branch, kAlphaWarningTolerance)) {
            if (suspicious_time_duplicate_with_far_alpha_count != nullptr) {
                *suspicious_time_duplicate_with_far_alpha_count += 1;
            }
            continue;
        }
        if (std::abs(branch.residual_seconds) < std::abs(existing.residual_seconds)) {
            existing = branch;
        }
        return true;
    }
    branches->push_back(branch);
    return false;
}

std::vector<SourcePipelineRecord> build_route_a_source_pipeline_records(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& node_branches,
    spaceship_cpp::problem1::Problem1RootDerivativeMode derivative_mode,
    int max_target_revolution = 1
) {
    namespace problem1 = spaceship_cpp::problem1;
    std::vector<SourcePipelineRecord> records;
    records.reserve(node_branches.size());
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
            derivative_mode,
            1e-6);
        if (!attached.valid || !attached.derivatives_available) {
            record.attach_derivative_failed = true;
            record.pipeline_classification = "attach_derivative_failed";
            records.push_back(record);
            continue;
        }
        const double delta_alpha =
            attached.d_encounter_global_angle_d_nu_A * node.dnu_A +
            attached.d_encounter_global_angle_d_nu_B * node.dnu_B +
            attached.d_encounter_global_angle_d_theta_A * node.dtheta_A;
        const double alpha_linear = normalize_angle_0_2pi(attached.encounter_global_angle + delta_alpha);
        const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            source.transfer_revolution,
            alpha_linear,
            source,
            max_target_revolution);
        record.selected_q = q_selection.selected_q;
        record.q_selection_failed = q_selection.selection_failed;
        record.q_changed = q_selection.q_changed;
        record.source_q_continuity_error = q_selection.source_q_continuity_error;
        record.selected_q_continuity_error = q_selection.selected_continuity_error;
        record.raw_predicted_alpha = alpha_linear;
        record.alpha_linear_finite = std::isfinite(alpha_linear);
        if (!record.alpha_linear_finite) {
            record.pipeline_classification = "raw_prediction_invalid";
            records.push_back(record);
            continue;
        }
        const auto residual = problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            alpha_linear,
            source.transfer_revolution,
            q_selection.selected_q);
        record.raw_residual_valid = residual.valid;
        record.raw_residual_seconds = residual.residual_seconds;
        if (!residual.valid) {
            record.pipeline_classification = "raw_residual_invalid";
            records.push_back(record);
            continue;
        }
        record.raw_branch.valid = true;
        record.raw_branch.encounter_global_angle = residual.encounter_global_angle;
        record.raw_branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
        record.raw_branch.transfer_revolution = source.transfer_revolution;
        record.raw_branch.target_revolution = source.target_revolution;
        record.raw_branch.time_of_flight_seconds = residual.transfer_time_seconds;
        record.raw_branch.target_time_seconds = residual.target_time_seconds;
        record.raw_branch.residual_seconds = residual.residual_seconds;
        record.raw_branch.transfer_e = residual.transfer_e;
        record.raw_branch.transfer_p = residual.transfer_p;
        record.raw_branch.transfer_a = residual.transfer_a;
        record.raw_branch.theta_B = residual.theta_B;
        const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            source.transfer_revolution,
            q_selection.selected_q,
            alpha_linear,
            80,
            1e-2,
            1e-12,
            derivative_mode,
            1e-6);
        record.newton_valid = refined.valid;
        record.refinement_result = refined;
        record.invalid_reason = refined.valid ? std::string{} : "target_time_q_sheet_refinement_invalid";
        if (!record.newton_valid) {
            record.pipeline_classification = "newton_invalid";
            records.push_back(record);
            continue;
        }
        record.final_branch = refined.branch;
        records.push_back(record);
    }
    return records;
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

struct DerivativeModeConfig {
    std::string name;
    spaceship_cpp::problem1::Problem1RootDerivativeMode mode;
};

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;
    using clock = std::chrono::steady_clock;

    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const int samples_per_group = sample_count_from_environment();
    const std::vector<QuerySample> samples = build_samples(samples_per_group);
    const std::vector<DerivativeModeConfig> derivative_modes{
        {"AnalyticOnly", problem1::Problem1RootDerivativeMode::AnalyticOnly},
        {"AnalyticWithFiniteDifferenceFallback", problem1::Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback},
    };

    std::map<VirtualNodeKey, TimedNodeSolve> solve_cache;
    double node_solve_time_total = 0.0;
    double node_solve_time_min = std::numeric_limits<double>::infinity();
    double node_solve_time_max = 0.0;
    double query_exact_solve_time_total = 0.0;
    int node_total_branch_count = 0;
    int exact_total_branch_count = 0;

    std::map<std::string, std::map<std::string, ModeGroupSummary>> summaries;

    for (const QuerySample& sample : samples) {
        const Problem1RootVirtualGridNode node =
            find_nearest_virtual_root_table_node(sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
        const VirtualNodeKey key{
            departure_planet, target_planet, node.nu_A_index, node.nu_B_index, node.theta_A_index,
            max_transfer_revolution, max_target_revolution};
        auto solve_it = solve_cache.find(key);
        if (solve_it == solve_cache.end()) {
            const auto start = clock::now();
            TimedNodeSolve timed{};
            timed.branches = problem1::solve_problem1_from_departure_anomalies(
                departure_planet, target_planet, node.nu_A_node, node.nu_B_node, node.theta_A_node,
                max_transfer_revolution, max_target_revolution);
            const auto end = clock::now();
            timed.solve_time_seconds = std::chrono::duration<double>(end - start).count();
            node_solve_time_total += timed.solve_time_seconds;
            node_solve_time_min = std::min(node_solve_time_min, timed.solve_time_seconds);
            node_solve_time_max = std::max(node_solve_time_max, timed.solve_time_seconds);
            node_total_branch_count += static_cast<int>(timed.branches.size());
            solve_it = solve_cache.emplace(key, std::move(timed)).first;
        }

        const auto exact_start = clock::now();
        std::vector<problem1::Problem1SolutionBranch> exact_branches =
            problem1::solve_problem1_from_departure_anomalies(
                departure_planet, target_planet, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                max_transfer_revolution, max_target_revolution);
        for (auto& exact_branch : exact_branches) {
            if (!exact_branch.valid) {
                continue;
            }
            const auto polished = problem1::refine_problem1_root_branch_newton_seconds(
                departure_planet,
                target_planet,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                exact_branch.transfer_revolution,
                exact_branch.target_revolution,
                exact_branch.encounter_global_angle,
                30,
                1e-6,
                1e-14);
            if (polished.valid) {
                exact_branch = polished;
            }
        }
        const auto exact_end = clock::now();
        query_exact_solve_time_total += std::chrono::duration<double>(exact_end - exact_start).count();
        exact_total_branch_count += static_cast<int>(exact_branches.size());

        for (const auto& mode_config : derivative_modes) {
            ModeGroupSummary& summary = summaries[mode_config.name][sample.group_name];
            summary.sample_count += 1;
            summary.exact_total_branch_count += static_cast<int>(exact_branches.size());
            {
                std::map<std::pair<int, int>, int> exact_kq_counts;
                std::map<std::pair<int, int>, int> node_kq_counts;
                for (const auto& exact_branch : exact_branches) {
                    if (exact_branch.valid) {
                        exact_kq_counts[{exact_branch.transfer_revolution, exact_branch.target_revolution}] += 1;
                    }
                }
                for (const auto& node_branch : solve_it->second.branches) {
                    if (node_branch.valid) {
                        node_kq_counts[{node_branch.transfer_revolution, node_branch.target_revolution}] += 1;
                    }
                }
                bool has_multi_root_same_kq_group = false;
                for (const auto& [kq, count] : exact_kq_counts) {
                    if (count > 1) {
                        has_multi_root_same_kq_group = true;
                        break;
                    }
                }
                if (!has_multi_root_same_kq_group) {
                    for (const auto& [kq, count] : node_kq_counts) {
                        if (count > 1) {
                            has_multi_root_same_kq_group = true;
                            break;
                        }
                    }
                }
                if (has_multi_root_same_kq_group) {
                    summary.multi_root_same_kq_group_count += 1;
                }
            }

            std::vector<problem1::Problem1SolutionBranch> raw_predictions;
            struct RefinedPrediction {
                problem1::Problem1SolutionBranch branch;
                problem1::Problem1RootRefinementResult refinement_result;
                int source_node_branch_index = -1;
                problem1::Problem1SolutionBranch source_node_branch;
                double raw_predicted_alpha = std::numeric_limits<double>::quiet_NaN();
                double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
                std::string invalid_reason;
            };
            std::vector<RefinedPrediction> refined_predictions;
            std::vector<SourcePipelineRecord> source_pipeline_records =
                build_route_a_source_pipeline_records(
                    departure_planet,
                    target_planet,
                    sample,
                    node,
                    solve_it->second.branches,
                    mode_config.mode);
            for (const auto& record : source_pipeline_records) {
                if (record.raw_residual_valid) {
                    raw_predictions.push_back(record.raw_branch);
                }
                if (!record.newton_valid) {
                    continue;
                }
                refined_predictions.push_back({
                    record.final_branch,
                    record.refinement_result,
                    record.source_node_branch_index,
                    record.source_node_branch,
                    record.raw_predicted_alpha,
                    record.raw_residual_seconds,
                    record.invalid_reason});
            }

            summary.raw_same_kq_stats.refined_candidate_count_before_dedup +=
                static_cast<int>(raw_predictions.size());
            summary.newton_same_kq_stats.refined_candidate_count_before_dedup +=
                static_cast<int>(refined_predictions.size());
            summary.newton_same_k_time_stats.refined_candidate_count_before_dedup +=
                static_cast<int>(refined_predictions.size());

            std::vector<problem1::Problem1SolutionBranch> raw_predictions_deduped;
            for (const auto& raw : raw_predictions) {
                add_refined_branch_if_not_duplicate(
                    &raw_predictions_deduped,
                    raw,
                    &summary.suspicious_time_duplicate_with_far_alpha_count);
            }
            summary.raw_same_kq_stats.refined_candidate_count_after_dedup +=
                static_cast<int>(raw_predictions_deduped.size());

            std::vector<problem1::Problem1SolutionBranch> refined_only;
            for (auto& item : refined_predictions) {
                const bool merged = add_refined_branch_if_not_duplicate(
                    &refined_only,
                    item.branch,
                    &summary.suspicious_time_duplicate_with_far_alpha_count);
                for (auto& record : source_pipeline_records) {
                    if (record.source_node_branch_index == item.source_node_branch_index &&
                        record.newton_valid &&
                        record.final_branch.time_of_flight_seconds == item.branch.time_of_flight_seconds &&
                        record.final_branch.transfer_revolution == item.branch.transfer_revolution &&
                        record.final_branch.target_revolution == item.branch.target_revolution) {
                        record.dedup_merged = merged;
                        break;
                    }
                }
            }
            for (auto& record : source_pipeline_records) {
                if (record.newton_valid) {
                    if (record.dedup_merged) {
                        record.pipeline_classification = "newton_converged_duplicate";
                    } else {
                        double best_time_error = std::numeric_limits<double>::infinity();
                        for (const auto& exact_branch : exact_branches) {
                            if (!exact_branch.valid ||
                                exact_branch.transfer_revolution != record.final_branch.transfer_revolution) {
                                continue;
                            }
                            best_time_error = std::min(
                                best_time_error,
                                std::abs(exact_branch.time_of_flight_seconds - record.final_branch.time_of_flight_seconds));
                        }
                        if (std::isfinite(best_time_error) && best_time_error > 1e5) {
                            record.pipeline_classification = "newton_converged_wrong_time_region";
                        } else {
                            record.pipeline_classification = "newton_converged_unique";
                        }
                    }
                }
            }
            for (const auto& record : source_pipeline_records) {
                if (record.q_selection_failed) {
                    summary.engineering_stats.q_sheet_selection_failed_count += 1;
                } else if (record.q_changed) {
                    summary.engineering_stats.q_sheet_selection_changed_count += 1;
                    if (record.newton_valid) {
                        summary.engineering_stats.q_sheet_selection_successful_changed_count += 1;
                    }
                } else {
                    summary.engineering_stats.q_sheet_selection_source_q_kept_count += 1;
                }
                if (record.newton_valid) {
                    summary.engineering_stats.residual_tolerance_success_count += 1;
                    if (record.refinement_result.diagnostic.derivative_attach_failed_after_convergence) {
                        summary.engineering_stats.derivative_attach_failed_after_convergence_count += 1;
                    }
                    continue;
                }
                const std::string& reason = record.invalid_reason;
                if (reason == "residual_increase_stop") {
                    summary.engineering_stats.residual_increase_stop_count += 1;
                } else if (reason == "Newton refinement exceeded maximum iterations") {
                    summary.engineering_stats.max_iteration_stop_count += 1;
                } else if (record.refinement_result.diagnostic.derivative_failed) {
                    summary.engineering_stats.derivative_invalid_stop_count += 1;
                }
            }
            std::vector<problem1::Problem1SolutionBranch> refined_primary_only;
            std::vector<problem1::Problem1SolutionBranch> refined_primary_plus_rescue;
            for (const auto& branch : refined_only) {
                const CandidateQuality quality = classify_candidate_quality(branch);
                if (quality == CandidateQuality::Primary) {
                    refined_primary_only.push_back(branch);
                    refined_primary_plus_rescue.push_back(branch);
                } else if (quality == CandidateQuality::Rescue) {
                    refined_primary_plus_rescue.push_back(branch);
                }
            }
            summary.newton_same_kq_stats.refined_candidate_count_after_dedup +=
                static_cast<int>(refined_only.size());
            summary.newton_same_k_time_stats.refined_candidate_count_after_dedup +=
                static_cast<int>(refined_only.size());

            const auto exact_sorted = collect_valid_sorted_by_time(exact_branches);
            std::map<int, std::vector<IndexedBranch>> exact_by_k;
            for (const auto& entry : exact_sorted) {
                exact_by_k[entry.branch.transfer_revolution].push_back(entry);
            }

            std::vector<IndexedBranch> refined_before_dedup_sorted;
            refined_before_dedup_sorted.reserve(refined_predictions.size());
            for (std::size_t i = 0; i < refined_predictions.size(); ++i) {
                if (!refined_predictions[i].branch.valid) {
                    continue;
                }
                refined_before_dedup_sorted.push_back(
                    {static_cast<int>(i), refined_predictions[i].branch});
            }
            std::sort(
                refined_before_dedup_sorted.begin(),
                refined_before_dedup_sorted.end(),
                [](const IndexedBranch& lhs, const IndexedBranch& rhs) {
                    if (lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds) {
                        return true;
                    }
                    if (lhs.branch.time_of_flight_seconds > rhs.branch.time_of_flight_seconds) {
                        return false;
                    }
                    return lhs.branch.encounter_global_angle < rhs.branch.encounter_global_angle;
                });
            std::map<int, std::vector<IndexedBranch>> refined_before_dedup_by_k;
            for (const auto& entry : refined_before_dedup_sorted) {
                refined_before_dedup_by_k[entry.branch.transfer_revolution].push_back(entry);
            }

            const auto refined_sorted = collect_valid_sorted_by_time(refined_only);
            std::map<int, std::vector<IndexedBranch>> refined_by_k;
            for (const auto& entry : refined_sorted) {
                refined_by_k[entry.branch.transfer_revolution].push_back(entry);
            }
            summary.engineering_stats.exact_total_count += static_cast<int>(exact_sorted.size());
            summary.engineering_stats.candidate_total_count += static_cast<int>(refined_sorted.size());
            const auto refined_primary_only_sorted = collect_valid_sorted_by_time(refined_primary_only);
            std::map<int, std::vector<IndexedBranch>> refined_primary_only_by_k;
            for (const auto& entry : refined_primary_only_sorted) {
                refined_primary_only_by_k[entry.branch.transfer_revolution].push_back(entry);
            }
            const auto refined_primary_plus_rescue_sorted = collect_valid_sorted_by_time(refined_primary_plus_rescue);
            std::map<int, std::vector<IndexedBranch>> refined_primary_plus_rescue_by_k;
            for (const auto& entry : refined_primary_plus_rescue_sorted) {
                refined_primary_plus_rescue_by_k[entry.branch.transfer_revolution].push_back(entry);
            }
            std::map<int, std::vector<IndexedBranch>> refined_selective_rescue_by_k = refined_primary_only_by_k;
            for (const auto& [k, exact_group] : exact_by_k) {
                auto& selected = refined_selective_rescue_by_k[k];
                const int need = static_cast<int>(exact_group.size()) - static_cast<int>(selected.size());
                if (need <= 0) {
                    continue;
                }
                const auto pool_it = refined_primary_plus_rescue_by_k.find(k);
                if (pool_it == refined_primary_plus_rescue_by_k.end()) {
                    continue;
                }
                for (const auto& candidate : pool_it->second) {
                    bool already_selected = false;
                    for (const auto& existing : selected) {
                        if (safe_to_merge_branch(existing.branch, candidate.branch)) {
                            already_selected = true;
                            break;
                        }
                    }
                    if (already_selected) {
                        continue;
                    }
                    selected.push_back(candidate);
                    if (static_cast<int>(selected.size()) >= static_cast<int>(exact_group.size())) {
                        break;
                    }
                }
                std::sort(
                    selected.begin(),
                    selected.end(),
                    [](const IndexedBranch& lhs, const IndexedBranch& rhs) {
                        if (lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds) {
                            return true;
                        }
                        if (lhs.branch.time_of_flight_seconds > rhs.branch.time_of_flight_seconds) {
                            return false;
                        }
                        return lhs.branch.encounter_global_angle < rhs.branch.encounter_global_angle;
                    });
            }
            summary.per_k_count_tables_report += build_per_k_count_table_report(
                sample.group_name,
                sample,
                node,
                exact_by_k,
                refined_before_dedup_by_k,
                refined_by_k);
            evaluate_per_k_two_tier_policy(
                exact_by_k,
                refined_primary_only_by_k,
                &summary.primary_only_policy_stats);
            evaluate_per_k_two_tier_policy(
                exact_by_k,
                refined_primary_plus_rescue_by_k,
                &summary.primary_plus_rescue_policy_stats);
            evaluate_per_k_two_tier_policy(
                exact_by_k,
                refined_selective_rescue_by_k,
                &summary.selective_rescue_policy_stats);
            if (exact_sorted.size() != refined_sorted.size()) {
                summary.global_time_order_stats.count_mismatch_count += 1;
                summary.per_k_time_order_stats.count_mismatch_count += 1;
            }
            for (std::size_t i = 0; i < std::min(exact_sorted.size(), refined_sorted.size()); ++i) {
                record_order_match(
                    &summary.global_time_order_stats,
                    exact_sorted[i].branch,
                    refined_sorted[i].branch);
            }
            for (const auto& [k, exact_group] : exact_by_k) {
                const auto refined_it = refined_by_k.find(k);
                const int refined_count_in_k =
                    refined_it == refined_by_k.end() ? 0 : static_cast<int>(refined_it->second.size());
                if (refined_count_in_k < static_cast<int>(exact_group.size())) {
                    summary.engineering_stats.missing_candidate_count +=
                        static_cast<int>(exact_group.size()) - refined_count_in_k;
                } else if (refined_count_in_k > static_cast<int>(exact_group.size())) {
                    summary.engineering_stats.extra_candidate_count +=
                        refined_count_in_k - static_cast<int>(exact_group.size());
                }
                if (refined_it == refined_by_k.end() ||
                    refined_it->second.size() != exact_group.size()) {
                    summary.per_k_time_order_stats.count_mismatch_count += 1;
                    summary.engineering_stats.per_k_count_mismatch_count += 1;
                    if (!summary.first_per_k_count_mismatch_recorded) {
                        summary.first_per_k_count_mismatch_recorded = true;
                        const auto before_it = refined_before_dedup_by_k.find(k);
                        const std::vector<IndexedBranch> empty_group;
                        summary.first_per_k_count_mismatch_report = build_per_k_count_mismatch_report(
                            sample.group_name,
                            sample,
                            node,
                            k,
                            exact_group,
                            before_it == refined_before_dedup_by_k.end() ? empty_group : before_it->second,
                            refined_it == refined_by_k.end() ? empty_group : refined_it->second);
                        if (mode_config.name == "AnalyticOnly" &&
                            summary.first_missing_root_microscope_report.empty()) {
                            const auto node_group = collect_sorted_k_group(solve_it->second.branches, k);
                            std::vector<const SourcePipelineRecord*> pipeline_records_in_k;
                            for (const auto& record : source_pipeline_records) {
                                if (record.source_node_branch.transfer_revolution == k) {
                                    pipeline_records_in_k.push_back(&record);
                                }
                            }
                            std::sort(
                                pipeline_records_in_k.begin(),
                                pipeline_records_in_k.end(),
                                [](const SourcePipelineRecord* lhs, const SourcePipelineRecord* rhs) {
                                    const double lhs_time = lhs->newton_valid
                                        ? lhs->final_branch.time_of_flight_seconds
                                        : lhs->source_node_branch.time_of_flight_seconds;
                                    const double rhs_time = rhs->newton_valid
                                        ? rhs->final_branch.time_of_flight_seconds
                                        : rhs->source_node_branch.time_of_flight_seconds;
                                    return lhs_time < rhs_time;
                                });
                            std::ostringstream os;
                            os << std::setprecision(6) << std::scientific;
                            os << "missing_root_microscope_case\n";
                            os << "group_name=" << sample.group_name << '\n';
                            os << "query_nu_A=" << sample.query_nu_A << '\n';
                            os << "query_nu_B=" << sample.query_nu_B << '\n';
                            os << "query_theta_A=" << sample.query_theta_A << '\n';
                            os << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
                            os << "node_offset_degrees=("
                               << (node.dnu_A * 180.0 / kPi) << ","
                               << (node.dnu_B * 180.0 / kPi) << ","
                               << (node.dtheta_A * 180.0 / kPi) << ")\n";
                            os << "mismatch_k=" << k << '\n';
                            os << "exact_count_in_k=" << exact_group.size() << '\n';
                            os << "nearest_node_solve_count_in_k=" << node_group.size() << '\n';
                            os << "route_a_refined_count_before_dedup_in_k="
                               << (before_it == refined_before_dedup_by_k.end() ? 0 : before_it->second.size()) << '\n';
                            os << "exact_branches_in_mismatch_k\n";
                            for (const auto& entry : exact_group) {
                                os << "  original_index=" << entry.original_index
                                   << " k=" << entry.branch.transfer_revolution
                                   << " q=" << entry.branch.target_revolution
                                   << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
                                   << " alpha=" << entry.branch.encounter_global_angle
                                   << " residual_seconds=" << entry.branch.residual_seconds << '\n';
                            }
                            os << "nearest_node_solve_branches_in_mismatch_k\n";
                            for (const auto& entry : node_group) {
                                os << "  node_original_index=" << entry.original_index
                                   << " k=" << entry.branch.transfer_revolution
                                   << " q=" << entry.branch.target_revolution
                                   << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
                                   << " alpha=" << entry.branch.encounter_global_angle
                                   << " residual_seconds=" << entry.branch.residual_seconds
                                   << " derivatives_available=" << entry.branch.derivatives_available << '\n';
                            }
                            os << "route_a_refined_candidates_before_dedup_in_mismatch_k\n";
                            for (const auto* record : pipeline_records_in_k) {
                                if (!record->newton_valid) {
                                    continue;
                                }
                                os << "  candidate_original_index=" << record->source_node_branch_index
                                   << " source_node_branch_index=" << record->source_node_branch_index
                                   << " source_node_k=" << record->source_node_branch.transfer_revolution
                                   << " source_node_q=" << record->source_node_branch.target_revolution
                                   << " source_node_time=" << record->source_node_branch.time_of_flight_seconds
                                   << " source_node_alpha=" << record->source_node_branch.encounter_global_angle
                                   << " raw_predicted_alpha=" << record->raw_predicted_alpha
                                   << " raw_residual_seconds=" << record->raw_residual_seconds
                                   << " newton_valid=" << record->newton_valid
                                   << " final_k=" << record->final_branch.transfer_revolution
                                   << " final_q=" << record->final_branch.target_revolution
                                   << " final_time_of_flight_seconds=" << record->final_branch.time_of_flight_seconds
                                   << " final_alpha=" << record->final_branch.encounter_global_angle
                                   << " final_residual_seconds=" << record->final_branch.residual_seconds
                                   << " invalid_reason=" << record->invalid_reason << '\n';
                            }
                            os << "source_node_branch_pipeline_status_in_mismatch_k\n";
                            for (const auto* record : pipeline_records_in_k) {
                                os << "  source_node_branch_index=" << record->source_node_branch_index
                                   << " source_node_k=" << record->source_node_branch.transfer_revolution
                                   << " source_node_q=" << record->source_node_branch.target_revolution
                                   << " source_node_time=" << record->source_node_branch.time_of_flight_seconds
                                   << " source_node_alpha=" << record->source_node_branch.encounter_global_angle
                                   << " attach_derivative_failed=" << record->attach_derivative_failed
                                   << " alpha_linear_finite=" << record->alpha_linear_finite
                                   << " raw_residual_valid=" << record->raw_residual_valid
                                   << " raw_residual_seconds=" << record->raw_residual_seconds
                                   << " newton_valid=" << record->newton_valid
                                   << " final_time_of_flight_seconds="
                                   << (record->newton_valid ? record->final_branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN())
                                   << " final_alpha="
                                   << (record->newton_valid ? record->final_branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN())
                                   << " dedup_merged=" << record->dedup_merged
                                   << " pipeline_classification=" << record->pipeline_classification
                                   << " invalid_reason=" << record->invalid_reason << '\n';
                            }

                            os << "missing_exact_root_comparison_to_nearest_node\n";
                            const std::vector<IndexedBranch>& refined_after_group =
                                (refined_it == refined_by_k.end() ? empty_group : refined_it->second);
                            for (const auto& exact_entry : exact_group) {
                                bool strict_found = false;
                                for (const auto& refined_entry : refined_after_group) {
                                    if (is_strict_match(refined_entry.branch, exact_entry.branch)) {
                                        strict_found = true;
                                        break;
                                    }
                                }
                                if (strict_found) {
                                    continue;
                                }
                                double best_time_gap = std::numeric_limits<double>::infinity();
                                double nearest_node_time = std::numeric_limits<double>::quiet_NaN();
                                double best_alpha_gap = std::numeric_limits<double>::infinity();
                                double nearest_node_alpha = std::numeric_limits<double>::quiet_NaN();
                                bool same_q_exists = false;
                                for (const auto& node_entry : node_group) {
                                    const double tgap = std::abs(
                                        node_entry.branch.time_of_flight_seconds - exact_entry.branch.time_of_flight_seconds);
                                    if (tgap < best_time_gap) {
                                        best_time_gap = tgap;
                                        nearest_node_time = node_entry.branch.time_of_flight_seconds;
                                    }
                                    const double agap = wrapped_alpha_distance(
                                        node_entry.branch.encounter_global_angle,
                                        exact_entry.branch.encounter_global_angle);
                                    if (agap < best_alpha_gap) {
                                        best_alpha_gap = agap;
                                        nearest_node_alpha = node_entry.branch.encounter_global_angle;
                                    }
                                    if (node_entry.branch.target_revolution == exact_entry.branch.target_revolution) {
                                        same_q_exists = true;
                                    }
                                }
                                os << "  exact_missing_root_time=" << exact_entry.branch.time_of_flight_seconds
                                   << " exact_missing_root_alpha=" << exact_entry.branch.encounter_global_angle
                                   << " exact_missing_root_q=" << exact_entry.branch.target_revolution
                                   << " nearest_node_by_time_time=" << nearest_node_time
                                   << " time_gap=" << best_time_gap
                                   << " nearest_node_by_alpha_alpha_error=" << best_alpha_gap
                                   << " nearest_node_by_alpha_alpha=" << nearest_node_alpha
                                   << " same_q_exists=" << same_q_exists << '\n';
                            }

                            auto run_node_pipeline = [&](const Problem1RootVirtualGridNode& pool_node) {
                                const auto pool_branches = problem1::solve_problem1_from_departure_anomalies(
                                    departure_planet,
                                    target_planet,
                                    pool_node.nu_A_node,
                                    pool_node.nu_B_node,
                                    pool_node.theta_A_node,
                                    max_transfer_revolution,
                                    max_target_revolution);
                                auto pool_records = build_route_a_source_pipeline_records(
                                    departure_planet,
                                    target_planet,
                                    sample,
                                    pool_node,
                                    pool_branches,
                                    mode_config.mode);
                                std::vector<problem1::Problem1SolutionBranch> pool_refined;
                                int ignore_suspicious = 0;
                                for (const auto& record : pool_records) {
                                    if (!record.newton_valid) {
                                        continue;
                                    }
                                    add_refined_branch_if_not_duplicate(
                                        &pool_refined,
                                        record.final_branch,
                                        &ignore_suspicious);
                                }
                                return pool_refined;
                            };
                            const auto one_node_pool = refined_only;
                            std::vector<problem1::Problem1SolutionBranch> neighbor_pool;
                            const int axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));
                            for (int dA = -1; dA <= 1; ++dA) {
                                for (int dB = -1; dB <= 1; ++dB) {
                                    for (int dT = -1; dT <= 1; ++dT) {
                                        Problem1RootVirtualGridNode neighbor = node;
                                        auto wrap_index = [&](int base, int delta) {
                                            int idx = (base + delta) % axis_count;
                                            if (idx < 0) {
                                                idx += axis_count;
                                            }
                                            return idx;
                                        };
                                        neighbor.nu_A_index = wrap_index(node.nu_A_index, dA);
                                        neighbor.nu_B_index = wrap_index(node.nu_B_index, dB);
                                        neighbor.theta_A_index = wrap_index(node.theta_A_index, dT);
                                        neighbor.nu_A_node = normalize_angle_0_2pi(neighbor.nu_A_index * kVirtualGridStepRadians);
                                        neighbor.nu_B_node = normalize_angle_0_2pi(neighbor.nu_B_index * kVirtualGridStepRadians);
                                        neighbor.theta_A_node = normalize_angle_0_2pi(neighbor.theta_A_index * kVirtualGridStepRadians);
                                        neighbor.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(sample.query_nu_A) - neighbor.nu_A_node);
                                        neighbor.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(sample.query_nu_B) - neighbor.nu_B_node);
                                        neighbor.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(sample.query_theta_A) - neighbor.theta_A_node);
                                        auto pool_refined = run_node_pipeline(neighbor);
                                        int ignore_suspicious = 0;
                                        for (const auto& branch : pool_refined) {
                                            add_refined_branch_if_not_duplicate(&neighbor_pool, branch, &ignore_suspicious);
                                        }
                                    }
                                }
                            }
                            const auto one_node_sorted_k = collect_sorted_k_group(one_node_pool, k);
                            const auto neighbor_pool_sorted_k = collect_sorted_k_group(neighbor_pool, k);
                            int one_node_strict = 0;
                            for (std::size_t i = 0; i < std::min(exact_group.size(), one_node_sorted_k.size()); ++i) {
                                if (is_strict_match(one_node_sorted_k[i].branch, exact_group[i].branch)) {
                                    one_node_strict += 1;
                                }
                            }
                            int neighbor_pool_strict = 0;
                            for (std::size_t i = 0; i < std::min(exact_group.size(), neighbor_pool_sorted_k.size()); ++i) {
                                if (is_strict_match(neighbor_pool_sorted_k[i].branch, exact_group[i].branch)) {
                                    neighbor_pool_strict += 1;
                                }
                            }
                            int missing_roots_recovered = 0;
                            for (const auto& exact_entry : exact_group) {
                                bool one_node_has = false;
                                for (const auto& refined_entry : one_node_sorted_k) {
                                    if (is_strict_match(refined_entry.branch, exact_entry.branch)) {
                                        one_node_has = true;
                                        break;
                                    }
                                }
                                if (one_node_has) {
                                    continue;
                                }
                                for (const auto& neighbor_entry : neighbor_pool_sorted_k) {
                                    if (is_strict_match(neighbor_entry.branch, exact_entry.branch)) {
                                        missing_roots_recovered += 1;
                                        break;
                                    }
                                }
                            }
                            os << "neighbor_pool_diagnostic\n";
                            os << "  one_node_refined_count_in_k=" << one_node_sorted_k.size() << '\n';
                            os << "  neighbor_pool_refined_count_in_k=" << neighbor_pool_sorted_k.size() << '\n';
                            os << "  exact_count_in_k=" << exact_group.size() << '\n';
                            os << "  one_node_coverage_in_k=" << one_node_strict << "/" << exact_group.size() << '\n';
                            os << "  neighbor_pool_coverage_in_k=" << neighbor_pool_strict << "/" << exact_group.size() << '\n';
                            os << "  missing_roots_recovered_by_neighbor_pool=" << missing_roots_recovered << '\n';
                            summary.first_missing_root_microscope_report = os.str();
                        }
                    }
                }
                const std::size_t pair_count = refined_it == refined_by_k.end()
                    ? 0
                    : std::min(exact_group.size(), refined_it->second.size());
                for (std::size_t i = 0; i < pair_count; ++i) {
                    record_order_match(
                        &summary.per_k_time_order_stats,
                        exact_group[i].branch,
                        refined_it->second[i].branch);
                    summary.engineering_stats.per_k_pair_count += 1;
                    const auto& exact_branch = exact_group[i].branch;
                    const auto& candidate_branch = refined_it->second[i].branch;
                    const double residual_abs = std::abs(candidate_branch.residual_seconds);
                    const double time_error =
                        std::abs(candidate_branch.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
                    const double alpha_wrapped_error = alpha_error(candidate_branch, exact_branch);
                    update_max(residual_abs, &summary.engineering_stats.max_residual_seconds);
                    update_max(time_error, &summary.engineering_stats.max_time_error_to_exact);
                    update_max(alpha_wrapped_error, &summary.engineering_stats.max_alpha_wrapped_error);
                    if (is_engineering_match(candidate_branch, exact_branch)) {
                        summary.engineering_stats.engineering_coverage_count += 1;
                    } else if (residual_abs <= 1e-2 &&
                               time_error <= 1.0 &&
                               alpha_wrapped_error <= 1e-3) {
                        summary.engineering_stats.label_or_strict_mismatch_count += 1;
                    } else {
                        summary.engineering_stats.physical_wrong_root_count += 1;
                    }
                }
            }
            if (!summary.first_time_order_failure_recorded) {
                bool should_record = exact_sorted.size() != refined_sorted.size();
                if (!should_record) {
                    for (std::size_t i = 0; i < std::min(exact_sorted.size(), refined_sorted.size()); ++i) {
                        if (!is_strict_match(refined_sorted[i].branch, exact_sorted[i].branch)) {
                            should_record = true;
                            break;
                        }
                    }
                }
                if (should_record) {
                    summary.first_time_order_failure_recorded = true;
                    summary.first_time_order_failure_report = build_time_order_failure_report(
                        sample.group_name,
                        sample,
                        node,
                        exact_sorted,
                        refined_sorted);
                }
            }

            for (std::size_t exact_index = 0; exact_index < exact_branches.size(); ++exact_index) {
                const auto& exact_branch = exact_branches[exact_index];
                if (!exact_branch.valid) {
                    continue;
                }
                bool same_kq_exists_at_node = false;
                for (const auto& node_branch : solve_it->second.branches) {
                    if (node_branch.valid &&
                        node_branch.transfer_revolution == exact_branch.transfer_revolution &&
                        node_branch.target_revolution == exact_branch.target_revolution) {
                        same_kq_exists_at_node = true;
                        break;
                    }
                }
                const auto raw_match = find_best_refined_branch_match(
                    exact_branch,
                    raw_predictions_deduped,
                    kStrictAlphaThreshold,
                    kStrictTimeThresholdSeconds,
                    kStrictResidualThresholdSeconds1e6);
                record_best_match_stats(
                    &summary.raw_same_kq_stats,
                    raw_match,
                    kStrictAlphaThreshold,
                    kStrictTimeThresholdSeconds,
                    kStrictResidualThresholdSeconds1e6,
                    1e-3,
                    1e3,
                    1e-4,
                    1e-2,
                    1e5,
                    1e-2);
                const auto newton_match = find_best_refined_branch_match(
                    exact_branch,
                    refined_only,
                    kStrictAlphaThreshold,
                    kStrictTimeThresholdSeconds,
                    kStrictResidualThresholdSeconds1e6);
                record_best_match_stats(
                    &summary.newton_same_kq_stats,
                    newton_match,
                    kStrictAlphaThreshold,
                    kStrictTimeThresholdSeconds,
                    kStrictResidualThresholdSeconds1e6,
                    1e-3,
                    1e3,
                    1e-4,
                    1e-2,
                    1e5,
                    1e-2);
                const auto newton_time_match = find_best_refined_branch_match_by_k_time(
                    exact_branch,
                    refined_only,
                    kStrictAlphaThreshold,
                    kStrictTimeThresholdSeconds,
                    kStrictResidualThresholdSeconds1e6);
                record_best_match_stats(
                    &summary.newton_same_k_time_stats,
                    newton_time_match,
                    kStrictAlphaThreshold,
                    kStrictTimeThresholdSeconds,
                    kStrictResidualThresholdSeconds1e6,
                    1e-3,
                    1e3,
                    1e-4,
                    1e-2,
                    1e5,
                    1e-2);
                if (sample.group_name == "near_node" && newton_time_match.same_group_exists) {
                    const auto& best_candidate = refined_only[newton_time_match.candidate_index];
                    summary.near_node_outliers.push_back({
                        static_cast<int>(exact_index),
                        exact_branch.transfer_revolution,
                        exact_branch.target_revolution,
                        exact_branch.encounter_global_angle,
                        exact_branch.time_of_flight_seconds,
                        newton_time_match.candidate_index,
                        best_candidate.transfer_revolution,
                        best_candidate.target_revolution,
                        best_candidate.encounter_global_angle,
                        best_candidate.time_of_flight_seconds,
                        newton_time_match.alpha_error,
                        newton_time_match.time_error_seconds,
                        newton_time_match.residual_seconds,
                        newton_time_match.matched,
                        newton_time_match.alpha_error <= 1e-3 &&
                            newton_time_match.time_error_seconds <= 1e3 &&
                            std::abs(newton_time_match.residual_seconds) <= 1e-4,
                        newton_time_match.alpha_error <= 1e-2 &&
                            newton_time_match.time_error_seconds <= 1e5 &&
                            std::abs(newton_time_match.residual_seconds) <= 1e-2});
                }

                if (mode_config.name == "AnalyticOnly" &&
                    summary.branch2_two_tier_report.empty() &&
                    sample.group_name == "mid_cell" &&
                    std::abs(sample.query_nu_A - 4.555309e+00) < 1e-6 &&
                    std::abs(sample.query_nu_B - 5.986479e+00) < 1e-6 &&
                    std::abs(sample.query_theta_A - 2.949606e+00) < 1e-6) {
                    for (const auto& record : source_pipeline_records) {
                        if (record.source_node_branch_index != 2 || !record.newton_valid) {
                            continue;
                        }
                        double best_time_error = std::numeric_limits<double>::infinity();
                        double best_alpha_error = std::numeric_limits<double>::infinity();
                        for (const auto& exact_branch : exact_branches) {
                            if (!exact_branch.valid ||
                                exact_branch.transfer_revolution != record.final_branch.transfer_revolution) {
                                continue;
                            }
                            const double time_error =
                                std::abs(exact_branch.time_of_flight_seconds - record.final_branch.time_of_flight_seconds);
                            if (time_error < best_time_error) {
                                best_time_error = time_error;
                                best_alpha_error = wrapped_alpha_distance(
                                    exact_branch.encounter_global_angle,
                                    record.final_branch.encounter_global_angle);
                            }
                        }
                        std::ostringstream os;
                        os << std::setprecision(6) << std::scientific;
                        os << "fixed_branch2_two_tier_case\n";
                        os << "source_branch_index=2\n";
                        os << "final_residual_seconds=" << record.final_branch.residual_seconds << '\n';
                        const CandidateQuality quality = classify_candidate_quality(record.final_branch);
                        os << "engineering_candidate_accepted="
                           << is_engineering_match(record.final_branch, record.final_branch) << '\n';
                        os << "engineering_residual_threshold_seconds=1.000000e-02\n";
                        os << "legacy_two_tier_quality="
                           << (quality == CandidateQuality::Primary ? "Primary"
                               : quality == CandidateQuality::Rescue ? "Rescue" : "Rejected") << '\n';
                        os << "matched_exact_time_error=" << best_time_error << '\n';
                        os << "matched_exact_alpha_error=" << best_alpha_error << '\n';
                        summary.branch2_two_tier_report = os.str();
                        break;
                    }
                }

                if (!same_kq_exists_at_node) {
                    summary.failures.no_same_kq_seed_at_node += 1;
                    continue;
                }
                if (!raw_match.same_group_exists) {
                    summary.failures.same_kq_seed_exists_but_raw_far += 1;
                    continue;
                }
                const double raw_alpha_error = raw_match.alpha_error;
                const double raw_time_error = raw_match.time_error_seconds;
                if (raw_alpha_error > kRawFarAlphaThreshold || raw_time_error > kRawFarTimeThresholdSeconds) {
                    summary.failures.same_kq_seed_exists_but_raw_far += 1;
                    continue;
                }
                if (!newton_time_match.same_group_exists) {
                    summary.failures.newton_failed += 1;
                } else if (newton_time_match.matched) {
                    summary.failures.newton_converged_correct_root += 1;
                } else {
                    summary.failures.newton_converged_wrong_root += 1;
                }
            }
        }
    }

    const int virtual_grid_axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));
    const long long virtual_grid_total_nodes =
        static_cast<long long>(virtual_grid_axis_count) *
        static_cast<long long>(virtual_grid_axis_count) *
        static_cast<long long>(virtual_grid_axis_count);
    const int unique_virtual_node_count = static_cast<int>(solve_cache.size());
    const double node_solve_time_avg_per_unique_node = average_or_zero(node_solve_time_total, unique_virtual_node_count);
    const double estimated_total_build_seconds =
        node_solve_time_avg_per_unique_node * static_cast<double>(virtual_grid_total_nodes);
    const double estimated_total_build_hours = estimated_total_build_seconds / 3600.0;
    const unsigned int hardware_threads = std::thread::hardware_concurrency();

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Virtual 2-degree nearest-node experiment with derivative mode\n";
    std::cout << "virtual_grid_step_degrees=" << kVirtualGridStepDegrees << '\n';
    std::cout << "virtual_grid_axis_count=" << virtual_grid_axis_count << '\n';
    std::cout << "virtual_grid_total_nodes=" << virtual_grid_total_nodes << '\n';
    std::cout << "samples_per_group=" << samples_per_group << '\n';
    std::cout << "sample_count_total=" << samples.size() << '\n';
    std::cout << "unique_virtual_node_count=" << unique_virtual_node_count << '\n';
    std::cout << "node_solve_time_avg_per_unique_node=" << node_solve_time_avg_per_unique_node << '\n';
    std::cout << "node_solve_time_min=" << node_solve_time_min << '\n';
    std::cout << "node_solve_time_max=" << node_solve_time_max << '\n';
    std::cout << "query_exact_solve_time_total=" << query_exact_solve_time_total << '\n';
    std::cout << "estimated_total_build_seconds=" << estimated_total_build_seconds << '\n';
    std::cout << "estimated_total_build_hours=" << estimated_total_build_hours << '\n';
    if (hardware_threads > 0U) {
        std::cout << "estimated_build_hours_hardware_threads=" <<
            (estimated_total_build_hours / static_cast<double>(hardware_threads)) << '\n';
    }
    std::cout << "exact_total_branch_count=" << exact_total_branch_count << '\n';
    std::cout << "node_total_branch_count=" << node_total_branch_count << '\n';

    for (const auto& [mode_name, groups] : summaries) {
        for (const auto& [group_name, summary] : groups) {
            std::cout << "mode=" << mode_name << '\n';
            std::cout << "group=" << group_name << '\n';
            std::cout << "group_sample_count=" << summary.sample_count << '\n';
            std::cout << "group_exact_total_branch_count=" << summary.exact_total_branch_count << '\n';
            std::cout << "multi_root_same_kq_group_count=" << summary.multi_root_same_kq_group_count << '\n';
            std::cout << "suspicious_time_duplicate_with_far_alpha_count=" <<
                summary.suspicious_time_duplicate_with_far_alpha_count << '\n';
            std::cout << "time_based_dedup_before_count=" <<
                summary.newton_same_k_time_stats.refined_candidate_count_before_dedup << '\n';
            std::cout << "time_based_dedup_after_count=" <<
                summary.newton_same_k_time_stats.refined_candidate_count_after_dedup << '\n';
            std::cout << "raw_refined_candidate_count_before_dedup=" <<
                summary.raw_same_kq_stats.refined_candidate_count_before_dedup << '\n';
            std::cout << "raw_refined_candidate_count_after_dedup=" <<
                summary.raw_same_kq_stats.refined_candidate_count_after_dedup << '\n';
            std::cout << "newton_refined_candidate_count_before_dedup=" <<
                summary.newton_same_kq_stats.refined_candidate_count_before_dedup << '\n';
            std::cout << "newton_refined_candidate_count_after_dedup=" <<
                summary.newton_same_kq_stats.refined_candidate_count_after_dedup << '\n';
            std::cout << "coverage_ratio_raw=" <<
                average_or_zero(static_cast<double>(summary.raw_same_kq_stats.strict_matched_count), summary.exact_total_branch_count) << '\n';
            std::cout << "same_kq_strict_coverage=" <<
                average_or_zero(static_cast<double>(summary.newton_same_kq_stats.strict_matched_count), summary.exact_total_branch_count) << '\n';
            std::cout << "same_k_time_strict_coverage=" <<
                average_or_zero(static_cast<double>(summary.newton_same_k_time_stats.strict_matched_count), summary.exact_total_branch_count) << '\n';
            std::cout << "same_k_time_medium_coverage=" <<
                average_or_zero(static_cast<double>(summary.newton_same_k_time_stats.medium_matched_count), summary.exact_total_branch_count) << '\n';
            std::cout << "same_k_time_loose_coverage=" <<
                average_or_zero(static_cast<double>(summary.newton_same_k_time_stats.loose_matched_count), summary.exact_total_branch_count) << '\n';
            std::cout << "global_time_order_pair_count=" << summary.global_time_order_stats.pair_count << '\n';
            std::cout << "global_time_order_strict_1e_6_match_count=" <<
                summary.global_time_order_stats.strict_match_count << '\n';
            std::cout << "global_time_order_strict_1e_5_match_count=" <<
                summary.global_time_order_stats.strict_match_count_relaxed << '\n';
            std::cout << "global_time_order_medium_match_count=" <<
                summary.global_time_order_stats.medium_match_count << '\n';
            std::cout << "global_time_order_loose_match_count=" <<
                summary.global_time_order_stats.loose_match_count << '\n';
            std::cout << "global_time_order_alpha_error_avg=" <<
                average_or_zero(summary.global_time_order_stats.alpha_error_sum, summary.global_time_order_stats.error_count) << '\n';
            std::cout << "global_time_order_alpha_error_max=" << summary.global_time_order_stats.alpha_error_max << '\n';
            std::cout << "global_time_order_time_error_avg=" <<
                average_or_zero(summary.global_time_order_stats.time_error_sum, summary.global_time_order_stats.error_count) << '\n';
            std::cout << "global_time_order_time_error_max=" << summary.global_time_order_stats.time_error_max << '\n';
            std::cout << "global_time_order_residual_avg=" <<
                average_or_zero(summary.global_time_order_stats.residual_sum, summary.global_time_order_stats.error_count) << '\n';
            std::cout << "global_time_order_residual_max=" << summary.global_time_order_stats.residual_max << '\n';
            std::cout << "per_k_time_order_pair_count=" << summary.per_k_time_order_stats.pair_count << '\n';
            std::cout << "per_k_time_order_strict_1e_6_match_count=" <<
                summary.per_k_time_order_stats.strict_match_count << '\n';
            std::cout << "per_k_time_order_strict_1e_5_match_count=" <<
                summary.per_k_time_order_stats.strict_match_count_relaxed << '\n';
            std::cout << "per_k_time_order_medium_match_count=" <<
                summary.per_k_time_order_stats.medium_match_count << '\n';
            std::cout << "per_k_time_order_loose_match_count=" <<
                summary.per_k_time_order_stats.loose_match_count << '\n';
            std::cout << "per_k_time_order_alpha_error_avg=" <<
                average_or_zero(summary.per_k_time_order_stats.alpha_error_sum, summary.per_k_time_order_stats.error_count) << '\n';
            std::cout << "per_k_time_order_alpha_error_max=" << summary.per_k_time_order_stats.alpha_error_max << '\n';
            std::cout << "per_k_time_order_time_error_avg=" <<
                average_or_zero(summary.per_k_time_order_stats.time_error_sum, summary.per_k_time_order_stats.error_count) << '\n';
            std::cout << "per_k_time_order_time_error_max=" << summary.per_k_time_order_stats.time_error_max << '\n';
            std::cout << "per_k_time_order_residual_avg=" <<
                average_or_zero(summary.per_k_time_order_stats.residual_sum, summary.per_k_time_order_stats.error_count) << '\n';
            std::cout << "per_k_time_order_residual_max=" << summary.per_k_time_order_stats.residual_max << '\n';
            std::cout << "per_k_time_order_count_mismatch_count=" <<
                summary.per_k_time_order_stats.count_mismatch_count << '\n';
            std::cout << "no_same_kq_candidate_count=" << summary.newton_same_kq_stats.no_group_candidate_count << '\n';
            std::cout << "same_kq_exists_but_no_strict_match_count=" <<
                summary.newton_same_kq_stats.same_group_exists_but_no_strict_match_count << '\n';
            std::cout << "no_same_k_time_candidate_count=" <<
                summary.newton_same_k_time_stats.no_group_candidate_count << '\n';
            std::cout << "same_k_time_exists_but_no_strict_match_count=" <<
                summary.newton_same_k_time_stats.same_group_exists_but_no_strict_match_count << '\n';
            std::cout << "alpha_fail_count=" << summary.newton_same_k_time_stats.alpha_fail_count << '\n';
            std::cout << "time_fail_count=" << summary.newton_same_k_time_stats.time_fail_count << '\n';
            std::cout << "residual_fail_count=" << summary.newton_same_k_time_stats.residual_fail_count << '\n';
            std::cout << "alpha_and_time_fail_count=" <<
                summary.newton_same_k_time_stats.alpha_and_time_fail_count << '\n';
            std::cout << "strict_fail_but_medium_pass_count=" <<
                summary.newton_same_k_time_stats.strict_fail_but_medium_pass_count << '\n';
            std::cout << "strict_fail_but_loose_pass_count=" <<
                summary.newton_same_k_time_stats.strict_fail_but_loose_pass_count << '\n';
            std::cout << "best_alpha_error_avg=" <<
                average_or_zero(
                    summary.newton_same_k_time_stats.best_alpha_error_sum,
                    summary.newton_same_k_time_stats.best_error_count) << '\n';
            std::cout << "best_alpha_error_max=" << summary.newton_same_k_time_stats.best_alpha_error_max << '\n';
            std::cout << "best_time_error_avg=" <<
                average_or_zero(
                    summary.newton_same_k_time_stats.best_time_error_sum,
                    summary.newton_same_k_time_stats.best_error_count) << '\n';
            std::cout << "best_time_error_max=" << summary.newton_same_k_time_stats.best_time_error_max << '\n';
            std::cout << "best_residual_avg=" <<
                average_or_zero(
                    summary.newton_same_k_time_stats.best_residual_sum,
                    summary.newton_same_k_time_stats.best_error_count) << '\n';
            std::cout << "best_residual_max=" << summary.newton_same_k_time_stats.best_residual_max << '\n';
            std::cout << "newton_failed=" << summary.failures.newton_failed << '\n';
            std::cout << "newton_converged_wrong_root=" << summary.failures.newton_converged_wrong_root << '\n';
            std::cout << "newton_converged_correct_root=" << summary.failures.newton_converged_correct_root << '\n';
            std::cout << "same_kq_seed_exists_but_raw_far=" << summary.failures.same_kq_seed_exists_but_raw_far << '\n';
            std::cout << "derivative_failed_count=" << summary.failures.derivative_failed_count << '\n';
            std::cout << "fallback_used_count=" << summary.failures.fallback_used_count << '\n';
            std::cout << "primary_only_per_k_pair_count=" << summary.primary_only_policy_stats.pair_count << '\n';
            std::cout << "primary_only_per_k_strict_1e_6_match_count=" << summary.primary_only_policy_stats.strict_1e6_match_count << '\n';
            std::cout << "primary_only_per_k_strict_1e_5_match_count=" << summary.primary_only_policy_stats.strict_1e5_match_count << '\n';
            std::cout << "primary_only_per_k_medium_match_count=" << summary.primary_only_policy_stats.medium_match_count << '\n';
            std::cout << "primary_only_per_k_loose_match_count=" << summary.primary_only_policy_stats.loose_match_count << '\n';
            std::cout << "primary_only_per_k_count_mismatch_count=" << summary.primary_only_policy_stats.count_mismatch_count << '\n';
            std::cout << "primary_only_wrong_root_count=" << summary.primary_only_policy_stats.wrong_root_count << '\n';
            std::cout << "primary_only_extra_candidate_count=" << summary.primary_only_policy_stats.extra_candidate_count << '\n';
            std::cout << "primary_only_missing_candidate_count=" << summary.primary_only_policy_stats.missing_candidate_count << '\n';
            std::cout << "primary_plus_rescue_per_k_pair_count=" << summary.primary_plus_rescue_policy_stats.pair_count << '\n';
            std::cout << "primary_plus_rescue_per_k_strict_1e_6_match_count=" << summary.primary_plus_rescue_policy_stats.strict_1e6_match_count << '\n';
            std::cout << "primary_plus_rescue_per_k_strict_1e_5_match_count=" << summary.primary_plus_rescue_policy_stats.strict_1e5_match_count << '\n';
            std::cout << "primary_plus_rescue_per_k_medium_match_count=" << summary.primary_plus_rescue_policy_stats.medium_match_count << '\n';
            std::cout << "primary_plus_rescue_per_k_loose_match_count=" << summary.primary_plus_rescue_policy_stats.loose_match_count << '\n';
            std::cout << "primary_plus_rescue_per_k_count_mismatch_count=" << summary.primary_plus_rescue_policy_stats.count_mismatch_count << '\n';
            std::cout << "primary_plus_rescue_wrong_root_count=" << summary.primary_plus_rescue_policy_stats.wrong_root_count << '\n';
            std::cout << "primary_plus_rescue_extra_candidate_count=" << summary.primary_plus_rescue_policy_stats.extra_candidate_count << '\n';
            std::cout << "primary_plus_rescue_missing_candidate_count=" << summary.primary_plus_rescue_policy_stats.missing_candidate_count << '\n';
            std::cout << "selective_rescue_per_k_pair_count=" << summary.selective_rescue_policy_stats.pair_count << '\n';
            std::cout << "selective_rescue_per_k_strict_1e_6_match_count=" << summary.selective_rescue_policy_stats.strict_1e6_match_count << '\n';
            std::cout << "selective_rescue_per_k_strict_1e_5_match_count=" << summary.selective_rescue_policy_stats.strict_1e5_match_count << '\n';
            std::cout << "selective_rescue_per_k_medium_match_count=" << summary.selective_rescue_policy_stats.medium_match_count << '\n';
            std::cout << "selective_rescue_per_k_loose_match_count=" << summary.selective_rescue_policy_stats.loose_match_count << '\n';
            std::cout << "selective_rescue_per_k_count_mismatch_count=" << summary.selective_rescue_policy_stats.count_mismatch_count << '\n';
            std::cout << "selective_rescue_wrong_root_count=" << summary.selective_rescue_policy_stats.wrong_root_count << '\n';
            std::cout << "selective_rescue_extra_candidate_count=" << summary.selective_rescue_policy_stats.extra_candidate_count << '\n';
            std::cout << "selective_rescue_missing_candidate_count=" << summary.selective_rescue_policy_stats.missing_candidate_count << '\n';
            std::cout << "engineering_exact_total_count=" << summary.engineering_stats.exact_total_count << '\n';
            std::cout << "engineering_candidate_total_count=" << summary.engineering_stats.candidate_total_count << '\n';
            std::cout << "engineering_per_k_pair_count=" << summary.engineering_stats.per_k_pair_count << '\n';
            std::cout << "engineering_per_k_count_mismatch_count=" << summary.engineering_stats.per_k_count_mismatch_count << '\n';
            std::cout << "engineering_coverage_count=" << summary.engineering_stats.engineering_coverage_count << '\n';
            std::cout << "engineering_coverage_ratio=" <<
                average_or_zero(
                    static_cast<double>(summary.engineering_stats.engineering_coverage_count),
                    summary.engineering_stats.exact_total_count) << '\n';
            std::cout << "engineering_physical_wrong_root_count=" << summary.engineering_stats.physical_wrong_root_count << '\n';
            std::cout << "engineering_label_or_strict_mismatch_count=" << summary.engineering_stats.label_or_strict_mismatch_count << '\n';
            std::cout << "engineering_missing_candidate_count=" << summary.engineering_stats.missing_candidate_count << '\n';
            std::cout << "engineering_extra_candidate_count=" << summary.engineering_stats.extra_candidate_count << '\n';
            std::cout << "engineering_residual_tolerance_success_count=" << summary.engineering_stats.residual_tolerance_success_count << '\n';
            std::cout << "engineering_residual_increase_stop_count=" << summary.engineering_stats.residual_increase_stop_count << '\n';
            std::cout << "engineering_max_iteration_stop_count=" << summary.engineering_stats.max_iteration_stop_count << '\n';
            std::cout << "engineering_derivative_invalid_stop_count=" << summary.engineering_stats.derivative_invalid_stop_count << '\n';
            std::cout << "engineering_derivative_attach_failed_after_convergence_count=" <<
                summary.engineering_stats.derivative_attach_failed_after_convergence_count << '\n';
            std::cout << "engineering_q_sheet_selection_changed_count=" <<
                summary.engineering_stats.q_sheet_selection_changed_count << '\n';
            std::cout << "engineering_q_sheet_selection_failed_count=" <<
                summary.engineering_stats.q_sheet_selection_failed_count << '\n';
            std::cout << "engineering_q_sheet_selection_source_q_kept_count=" <<
                summary.engineering_stats.q_sheet_selection_source_q_kept_count << '\n';
            std::cout << "engineering_q_sheet_selection_successful_changed_count=" <<
                summary.engineering_stats.q_sheet_selection_successful_changed_count << '\n';
            std::cout << "engineering_max_residual_seconds=" << summary.engineering_stats.max_residual_seconds << '\n';
            std::cout << "engineering_max_time_error_to_exact=" << summary.engineering_stats.max_time_error_to_exact << '\n';
            std::cout << "engineering_max_alpha_wrapped_error=" << summary.engineering_stats.max_alpha_wrapped_error << '\n';
            if (group_name == "near_node" && !summary.near_node_outliers.empty()) {
                std::vector<NearNodeOutlierCase> outliers = summary.near_node_outliers;
                std::sort(
                    outliers.begin(),
                    outliers.end(),
                    [](const NearNodeOutlierCase& lhs, const NearNodeOutlierCase& rhs) {
                        return lhs.alpha_error > rhs.alpha_error;
                    });
                const std::size_t limit = std::min<std::size_t>(5, outliers.size());
                for (std::size_t i = 0; i < limit; ++i) {
                    const auto& outlier = outliers[i];
                    std::cout << "near_node_top_alpha_outlier_rank=" << i << '\n';
                    std::cout << "  exact_index=" << outlier.exact_index << '\n';
                    std::cout << "  exact_k=" << outlier.exact_k << '\n';
                    std::cout << "  exact_q=" << outlier.exact_q << '\n';
                    std::cout << "  exact_alpha=" << outlier.exact_alpha << '\n';
                    std::cout << "  exact_time_of_flight=" << outlier.exact_time_of_flight << '\n';
                    std::cout << "  best_candidate_index=" << outlier.best_candidate_index << '\n';
                    std::cout << "  best_candidate_k=" << outlier.best_candidate_k << '\n';
                    std::cout << "  best_candidate_q=" << outlier.best_candidate_q << '\n';
                    std::cout << "  best_candidate_alpha=" << outlier.best_candidate_alpha << '\n';
                    std::cout << "  best_candidate_time_of_flight=" << outlier.best_candidate_time_of_flight << '\n';
                    std::cout << "  alpha_error=" << outlier.alpha_error << '\n';
                    std::cout << "  time_error=" << outlier.time_error << '\n';
                    std::cout << "  residual_seconds=" << outlier.residual_seconds << '\n';
                    std::cout << "  strict_match=" << outlier.strict_match << '\n';
                    std::cout << "  medium_match=" << outlier.medium_match << '\n';
                    std::cout << "  loose_match=" << outlier.loose_match << '\n';
                }
            }
            if (!summary.first_time_order_failure_report.empty()) {
                std::cout << summary.first_time_order_failure_report;
            }
            if (!summary.per_k_count_tables_report.empty()) {
                std::cout << summary.per_k_count_tables_report;
            }
            if (!summary.first_per_k_count_mismatch_report.empty()) {
                std::cout << summary.first_per_k_count_mismatch_report;
            }
            if (!summary.first_missing_root_microscope_report.empty()) {
                std::cout << summary.first_missing_root_microscope_report;
            }
            if (!summary.branch2_two_tier_report.empty()) {
                std::cout << summary.branch2_two_tier_report;
            }
        }
    }

    assert(unique_virtual_node_count > 0);
    assert(node_solve_time_avg_per_unique_node > 0.0);
    assert(std::isfinite(estimated_total_build_seconds));
    assert(exact_total_branch_count > 0);
    return 0;
}
