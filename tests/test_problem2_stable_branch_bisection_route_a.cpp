#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kThetaSampleCount = 32;
constexpr int kMaxRecursionDepth = 8;
constexpr double kTopologyEpsilon = 1e-5;
constexpr int kBisectionMaxIterations = 60;
constexpr double kThetaTolerance = 1e-10;
constexpr double kResidualTolerance = 1e-8;

struct Problem2SlingshotResidualBranchTestOnly {
    bool valid = false;
    std::string invalid_reason;

    double theta_prime = 0.0;
    double alpha = 0.0;
    double encounter_global_angle = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double problem1_residual_seconds = 0.0;

    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;

    double incoming_invariant = 0.0;
    double outgoing_invariant = 0.0;
    double slingshot_residual = 0.0;
    bool derivatives_available = false;
    double d_alpha_d_theta_A = 0.0;
};

struct ThetaPrimeSampleSummary {
    double theta_prime = 0.0;
    int total_problem1_branch_count = 0;
    int valid_slingshot_branch_count = 0;
    std::map<int, int> valid_count_by_k;
    std::map<int, int> problem1_count_by_k;
    std::vector<Problem2SlingshotResidualBranchTestOnly> branches;
};

struct Problem2StableSignChangePair {
    double theta_left = 0.0;
    double theta_right = 0.0;
    int k = 0;
    int rank_in_k = 0;
    Problem2SlingshotResidualBranchTestOnly left_branch;
    Problem2SlingshotResidualBranchTestOnly right_branch;
};

struct Problem2StableBranchBisectionSolution {
    bool valid = false;
    std::string invalid_reason;
    int k = 0;
    int rank_in_k = 0;
    double theta_left_initial = 0.0;
    double theta_right_initial = 0.0;
    double theta_root = 0.0;
    double alpha_root = 0.0;
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double slingshot_residual = 0.0;
    int iterations = 0;
    double final_width = 0.0;
    double time_of_flight_seconds = 0.0;
    int target_revolution = 0;
    double problem1_residual_seconds = 0.0;
    double last_theta = 0.0;
    double last_residual = 0.0;
};

struct PairCollectionStats {
    int raw_sign_change_pair_count = 0;
    int continuation_stable_sign_change_pair_count = 0;
    int unstable_candidate_pair_count = 0;
    int midpoint_left_continuation_invalid_count = 0;
    int midpoint_right_continuation_invalid_count = 0;
    int midpoint_both_continuation_invalid_count = 0;
    int midpoint_one_sided_only_count = 0;
    int midpoint_two_sided_mismatch_count = 0;
    int midpoint_slingshot_geometry_invalid_count = 0;
    int midpoint_endpoint_derivatives_unavailable_count = 0;
    int midpoint_continuation_refine_failed_count = 0;
    int midpoint_continuation_problem1_residual_too_large_count = 0;
    int midpoint_continuation_derivatives_unavailable_count = 0;
    int q_diagnostic_endpoint_q_best_count = 0;
    int q_diagnostic_changed_q_best_count = 0;
};

struct ContinuationAttemptDiagnostic {
    bool valid = false;
    std::string invalid_reason;
    bool refined_valid = false;
    double refined_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double endpoint_theta_prime = 0.0;
    double theta_prime_new = 0.0;
    double delta_theta_prime = 0.0;
    double endpoint_alpha = 0.0;
    double endpoint_encounter_global_angle = 0.0;
    double d_alpha_d_theta_A = 0.0;
    double predicted_encounter_alpha_seed = 0.0;
    double refined_encounter_global_angle = std::numeric_limits<double>::quiet_NaN();
    double refined_target_arrival_true_anomaly = std::numeric_limits<double>::quiet_NaN();
    int endpoint_q = 0;
    int best_q = -1;
    bool best_q_valid = false;
    double best_q_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double best_q_alpha = std::numeric_limits<double>::quiet_NaN();
    double best_q_time_of_flight_seconds = std::numeric_limits<double>::quiet_NaN();
    bool best_q_same_as_endpoint = false;
};

struct TwoSidedContinuationDiagnostic {
    bool valid = false;
    std::string invalid_reason;
    bool left_valid = false;
    std::string left_invalid_reason;
    bool right_valid = false;
    std::string right_invalid_reason;
    double left_alpha = std::numeric_limits<double>::quiet_NaN();
    double right_alpha = std::numeric_limits<double>::quiet_NaN();
    double alpha_diff = std::numeric_limits<double>::quiet_NaN();
    double left_time = std::numeric_limits<double>::quiet_NaN();
    double right_time = std::numeric_limits<double>::quiet_NaN();
    double time_diff = std::numeric_limits<double>::quiet_NaN();
    double left_e_prime = std::numeric_limits<double>::quiet_NaN();
    double right_e_prime = std::numeric_limits<double>::quiet_NaN();
    double e_prime_diff = std::numeric_limits<double>::quiet_NaN();
    double left_slingshot_residual = std::numeric_limits<double>::quiet_NaN();
    double right_slingshot_residual = std::numeric_limits<double>::quiet_NaN();
    double slingshot_residual_diff = std::numeric_limits<double>::quiet_NaN();
    ContinuationAttemptDiagnostic left_continuation;
    ContinuationAttemptDiagnostic right_continuation;
    Problem2SlingshotResidualBranchTestOnly selected_branch;
};

constexpr int kMaxRejectedCasePrintCount = 20;
constexpr int kMaxContinuationFailurePrintCount = 20;
int g_rejected_case_print_count = 0;
int g_continuation_failure_print_count = 0;

Problem2SlingshotResidualBranchTestOnly continue_problem2_branch_to_theta_prime_route_a_test_only(
    const Problem2SlingshotResidualBranchTestOnly& endpoint_branch,
    double theta_prime_new,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution,
    ContinuationAttemptDiagnostic* diagnostic = nullptr
);

TwoSidedContinuationDiagnostic evaluate_midpoint_branch_by_two_sided_continuation_test_only(
    double theta_mid,
    const Problem2SlingshotResidualBranchTestOnly& left_branch,
    const Problem2SlingshotResidualBranchTestOnly& right_branch,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution
);

double wrapped_angle_distance(double lhs, double rhs) {
    const double delta = std::fmod(std::abs(lhs - rhs), kTwoPi);
    return std::min(delta, kTwoPi - delta);
}

double wrapped_angle_midpoint(double lhs, double rhs) {
    const double forward = std::fmod((rhs - lhs) + kTwoPi, kTwoPi);
    return normalize_angle_0_2pi(lhs + 0.5 * forward);
}

bool same_root_for_dedup(
    const Problem2StableBranchBisectionSolution& lhs,
    const Problem2StableBranchBisectionSolution& rhs
) {
    return lhs.k == rhs.k &&
           std::abs(lhs.theta_root - rhs.theta_root) <= 1e-6 &&
           wrapped_angle_distance(lhs.alpha_root, rhs.alpha_root) <= 1e-6;
}

std::vector<Problem2SlingshotResidualBranchTestOnly>
evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    int max_transfer_revolution,
    int max_target_revolution
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);

    const std::vector<problem1::Problem1SolutionBranch> branches =
        problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            encounter_anomaly_phi,
            target_current_anomaly_beta,
            theta_A,
            max_transfer_revolution,
            max_target_revolution);

    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);

    std::vector<Problem2SlingshotResidualBranchTestOnly> results;
    results.reserve(branches.size());
    for (const auto& branch : branches) {
        Problem2SlingshotResidualBranchTestOnly residual_branch{};
        residual_branch.theta_prime = theta_prime;
        residual_branch.alpha = branch.target_arrival_true_anomaly;
        residual_branch.encounter_global_angle = branch.encounter_global_angle;
        residual_branch.transfer_revolution = branch.transfer_revolution;
        residual_branch.target_revolution = branch.target_revolution;
        residual_branch.time_of_flight_seconds = branch.time_of_flight_seconds;
        residual_branch.target_time_seconds = branch.target_time_seconds;
        residual_branch.problem1_residual_seconds = branch.residual_seconds;

        if (!branch.valid) {
            residual_branch.invalid_reason =
                branch.invalid_reason.empty() ? "problem1_branch_invalid" : branch.invalid_reason;
            results.push_back(residual_branch);
            continue;
        }

        const auto differentiated = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            encounter_anomaly_phi,
            target_current_anomaly_beta,
            theta_A,
            branch,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (differentiated.valid && differentiated.derivatives_available) {
            residual_branch.derivatives_available = true;
            // target_arrival_true_anomaly differs from encounter_global_angle by the target-orbit
            // reference angle, so their theta_A derivative is the same for a fixed target state.
            residual_branch.d_alpha_d_theta_A = differentiated.d_encounter_global_angle_d_theta_A;
            residual_branch.encounter_global_angle = differentiated.encounter_global_angle;
        }

        const auto theta_alpha_residual =
            problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
                departure_params.orbit.p,
                departure_params.orbit.e,
                target_params.orbit.p,
                target_params.orbit.e,
                encounter_anomaly_phi,
                branch.target_arrival_true_anomaly,
                incoming_e,
                incoming_theta,
                theta_prime);
        if (!theta_alpha_residual.valid) {
            residual_branch.invalid_reason = theta_alpha_residual.invalid_reason;
            results.push_back(residual_branch);
            continue;
        }

        residual_branch.valid = true;
        residual_branch.outgoing_eccentricity = theta_alpha_residual.outgoing_eccentricity;
        residual_branch.outgoing_semi_latus_rectum = theta_alpha_residual.outgoing_semi_latus_rectum;
        residual_branch.incoming_invariant = theta_alpha_residual.incoming_invariant;
        residual_branch.outgoing_invariant = theta_alpha_residual.outgoing_invariant;
        residual_branch.slingshot_residual = theta_alpha_residual.slingshot_residual;
        results.push_back(residual_branch);
    }
    return results;
}

ThetaPrimeSampleSummary build_theta_summary(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    int max_transfer_revolution,
    int max_target_revolution
) {
    ThetaPrimeSampleSummary summary{};
    summary.theta_prime = theta_prime;
    summary.branches = evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
        departure_planet,
        target_planet,
        encounter_time,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        incoming_e,
        incoming_theta,
        theta_prime,
        max_transfer_revolution,
        max_target_revolution);
    summary.total_problem1_branch_count = static_cast<int>(summary.branches.size());
    for (const auto& branch : summary.branches) {
        summary.problem1_count_by_k[branch.transfer_revolution] += 1;
        if (branch.valid) {
            summary.valid_slingshot_branch_count += 1;
            summary.valid_count_by_k[branch.transfer_revolution] += 1;
        }
    }
    return summary;
}

bool same_count_by_k(const ThetaPrimeSampleSummary& lhs, const ThetaPrimeSampleSummary& rhs) {
    return lhs.valid_count_by_k == rhs.valid_count_by_k;
}

std::map<int, std::vector<Problem2SlingshotResidualBranchTestOnly>>
group_valid_branches_by_k_sorted(const std::vector<Problem2SlingshotResidualBranchTestOnly>& branches) {
    std::map<int, std::vector<Problem2SlingshotResidualBranchTestOnly>> by_k;
    for (const auto& branch : branches) {
        if (branch.valid) {
            by_k[branch.transfer_revolution].push_back(branch);
        }
    }
    for (auto& [k, group] : by_k) {
        (void)k;
        std::sort(
            group.begin(),
            group.end(),
            [](const auto& a, const auto& b) { return a.time_of_flight_seconds < b.time_of_flight_seconds; });
    }
    return by_k;
}

template <class SummaryBuilder>
ThetaPrimeSampleSummary evaluate_cached_or_build(
    std::map<double, ThetaPrimeSampleSummary>* cache,
    double theta_prime,
    SummaryBuilder&& builder
) {
    const auto it = cache->find(theta_prime);
    if (it != cache->end()) {
        return it->second;
    }
    ThetaPrimeSampleSummary summary = builder(theta_prime);
    cache->emplace(theta_prime, summary);
    return summary;
}

void collect_stable_sign_change_pairs(
    const ThetaPrimeSampleSummary& left,
    const ThetaPrimeSampleSummary& right,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution,
    std::vector<Problem2StableSignChangePair>* pairs,
    PairCollectionStats* stats
) {
    const auto left_by_k = group_valid_branches_by_k_sorted(left.branches);
    const auto right_by_k = group_valid_branches_by_k_sorted(right.branches);
    for (const auto& [k, left_group] : left_by_k) {
        const auto right_it = right_by_k.find(k);
        if (right_it == right_by_k.end()) {
            continue;
        }
        const auto& right_group = right_it->second;
        const int pair_count = std::min<int>(left_group.size(), right_group.size());
        for (int rank = 0; rank < pair_count; ++rank) {
            const auto& l = left_group[static_cast<std::size_t>(rank)];
            const auto& r = right_group[static_cast<std::size_t>(rank)];
            if (!std::isfinite(l.slingshot_residual) || !std::isfinite(r.slingshot_residual)) {
                continue;
            }
            const bool sign_change =
                (l.slingshot_residual < 0.0 && r.slingshot_residual > 0.0) ||
                (l.slingshot_residual > 0.0 && r.slingshot_residual < 0.0);
            if (!sign_change) {
                continue;
            }
            stats->raw_sign_change_pair_count += 1;
            const double theta_mid = 0.5 * (left.theta_prime + right.theta_prime);
            const auto mid = evaluate_midpoint_branch_by_two_sided_continuation_test_only(
                theta_mid,
                l,
                r,
                departure_planet,
                target_planet,
                encounter_time,
                encounter_anomaly_phi,
                target_current_anomaly_beta,
                incoming_e,
                incoming_theta,
                max_transfer_revolution,
                max_target_revolution);
            if (!mid.valid) {
                stats->unstable_candidate_pair_count += 1;
                const bool left_invalid = !mid.left_valid;
                const bool right_invalid = !mid.right_valid;
                if (left_invalid && right_invalid) {
                    stats->midpoint_both_continuation_invalid_count += 1;
                } else if (left_invalid) {
                    stats->midpoint_left_continuation_invalid_count += 1;
                } else if (right_invalid) {
                    stats->midpoint_right_continuation_invalid_count += 1;
                }
                if (mid.invalid_reason == "one_sided_continuation_only") {
                    stats->midpoint_one_sided_only_count += 1;
                } else if (mid.invalid_reason == "two_sided_continuation_mismatch") {
                    stats->midpoint_two_sided_mismatch_count += 1;
                }
                const auto bump_reason = [&](const std::string& reason) {
                    if (reason == "endpoint_derivatives_unavailable") {
                        stats->midpoint_endpoint_derivatives_unavailable_count += 1;
                    } else if (reason == "continuation_refine_failed") {
                        stats->midpoint_continuation_refine_failed_count += 1;
                    } else if (reason == "continuation_problem1_residual_too_large") {
                        stats->midpoint_continuation_problem1_residual_too_large_count += 1;
                    } else if (reason == "continuation_derivatives_unavailable") {
                        stats->midpoint_continuation_derivatives_unavailable_count += 1;
                    } else if (reason == "continuation_slingshot_geometry_invalid") {
                        stats->midpoint_slingshot_geometry_invalid_count += 1;
                    }
                };
                bump_reason(mid.left_invalid_reason);
                bump_reason(mid.right_invalid_reason);
                if (mid.left_continuation.best_q_valid) {
                    if (mid.left_continuation.best_q_same_as_endpoint) {
                        stats->q_diagnostic_endpoint_q_best_count += 1;
                    } else {
                        stats->q_diagnostic_changed_q_best_count += 1;
                    }
                }
                if (mid.right_continuation.best_q_valid) {
                    if (mid.right_continuation.best_q_same_as_endpoint) {
                        stats->q_diagnostic_endpoint_q_best_count += 1;
                    } else {
                        stats->q_diagnostic_changed_q_best_count += 1;
                    }
                }
                if (g_rejected_case_print_count < kMaxRejectedCasePrintCount) {
                    std::cout << "Problem2MidpointContinuationRejectedCase\n";
                    std::cout << "theta_left=" << left.theta_prime << '\n';
                    std::cout << "theta_mid=" << theta_mid << '\n';
                    std::cout << "theta_right=" << right.theta_prime << '\n';
                    std::cout << "k=" << k << '\n';
                    std::cout << "rank_in_k=" << rank << '\n';
                    std::cout << "left_endpoint_alpha=" << l.alpha << '\n';
                    std::cout << "right_endpoint_alpha=" << r.alpha << '\n';
                    std::cout << "left_endpoint_residual=" << l.slingshot_residual << '\n';
                    std::cout << "right_endpoint_residual=" << r.slingshot_residual << '\n';
                    std::cout << "left_valid=" << (mid.left_valid ? 1 : 0) << '\n';
                    std::cout << "left_invalid_reason=" << mid.left_invalid_reason << '\n';
                    std::cout << "right_valid=" << (mid.right_valid ? 1 : 0) << '\n';
                    std::cout << "right_invalid_reason=" << mid.right_invalid_reason << '\n';
                    std::cout << "two_sided_invalid_reason=" << mid.invalid_reason << '\n';
                    std::cout << "left_mid_alpha=" << mid.left_alpha << '\n';
                    std::cout << "right_mid_alpha=" << mid.right_alpha << '\n';
                    std::cout << "alpha_diff=" << mid.alpha_diff << '\n';
                    std::cout << "left_mid_time=" << mid.left_time << '\n';
                    std::cout << "right_mid_time=" << mid.right_time << '\n';
                    std::cout << "time_diff=" << mid.time_diff << '\n';
                    std::cout << "left_mid_e_prime=" << mid.left_e_prime << '\n';
                    std::cout << "right_mid_e_prime=" << mid.right_e_prime << '\n';
                    std::cout << "e_prime_diff=" << mid.e_prime_diff << '\n';
                    std::cout << "left_mid_slingshot_residual=" << mid.left_slingshot_residual << '\n';
                    std::cout << "right_mid_slingshot_residual=" << mid.right_slingshot_residual << '\n';
                    std::cout << "slingshot_residual_diff=" << mid.slingshot_residual_diff << '\n';
                    ++g_rejected_case_print_count;
                }
                const auto print_continuation_failure = [&](const char* side, const ContinuationAttemptDiagnostic& d) {
                    if (d.valid || g_continuation_failure_print_count >= kMaxContinuationFailurePrintCount) {
                        return;
                    }
                    std::cout << "Problem2EndpointContinuationFailureCase\n";
                    std::cout << "side=" << side << '\n';
                    std::cout << "endpoint_theta_prime=" << d.endpoint_theta_prime << '\n';
                    std::cout << "theta_prime_new=" << d.theta_prime_new << '\n';
                    std::cout << "delta_theta_prime=" << d.delta_theta_prime << '\n';
                    std::cout << "endpoint_alpha=" << d.endpoint_alpha << '\n';
                    std::cout << "endpoint_encounter_global_angle=" << d.endpoint_encounter_global_angle << '\n';
                    std::cout << "d_alpha_d_theta_A=" << d.d_alpha_d_theta_A << '\n';
                    std::cout << "predicted_encounter_alpha_seed=" << d.predicted_encounter_alpha_seed << '\n';
                    std::cout << "invalid_reason=" << d.invalid_reason << '\n';
                    std::cout << "refined_valid=" << (d.refined_valid ? 1 : 0) << '\n';
                    std::cout << "refined_residual_seconds=" << d.refined_residual_seconds << '\n';
                    std::cout << "Problem2ContinuationQDiagnostic\n";
                    std::cout << "side=" << side << '\n';
                    std::cout << "endpoint_q=" << d.endpoint_q << '\n';
                    std::cout << "best_q=" << d.best_q << '\n';
                    std::cout << "best_valid=" << (d.best_q_valid ? 1 : 0) << '\n';
                    std::cout << "best_residual_seconds=" << d.best_q_residual_seconds << '\n';
                    std::cout << "best_alpha=" << d.best_q_alpha << '\n';
                    std::cout << "best_time_of_flight_seconds=" << d.best_q_time_of_flight_seconds << '\n';
                    std::cout << "best_q_same_as_endpoint=" << (d.best_q_same_as_endpoint ? 1 : 0) << '\n';
                    ++g_continuation_failure_print_count;
                };
                print_continuation_failure("left", mid.left_continuation);
                print_continuation_failure("right", mid.right_continuation);
                continue;
            }
            stats->continuation_stable_sign_change_pair_count += 1;
            Problem2StableSignChangePair pair{};
            pair.theta_left = left.theta_prime;
            pair.theta_right = right.theta_prime;
            pair.k = k;
            pair.rank_in_k = rank;
            pair.left_branch = l;
            pair.right_branch = r;
            pairs->push_back(pair);
        }
    }
}

template <class SummaryBuilder>
void collect_pairs_from_interval(
    const ThetaPrimeSampleSummary& left,
    const ThetaPrimeSampleSummary& right,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution,
    int depth,
    std::map<double, ThetaPrimeSampleSummary>* cache,
    SummaryBuilder&& builder,
    std::vector<Problem2StableSignChangePair>* pairs,
    PairCollectionStats* stats
) {
    const double width = right.theta_prime - left.theta_prime;
    if (same_count_by_k(left, right)) {
        if (width >= kTopologyEpsilon && depth < kMaxRecursionDepth) {
            const double theta_mid = 0.5 * (left.theta_prime + right.theta_prime);
            const ThetaPrimeSampleSummary mid = evaluate_cached_or_build(cache, theta_mid, builder);
            if (!same_count_by_k(left, mid) || !same_count_by_k(mid, right)) {
                collect_pairs_from_interval(
                    left, mid,
                    departure_planet, target_planet, encounter_time, encounter_anomaly_phi,
                    target_current_anomaly_beta, incoming_e, incoming_theta,
                    max_transfer_revolution, max_target_revolution,
                    depth + 1, cache, builder, pairs, stats);
                collect_pairs_from_interval(
                    mid, right,
                    departure_planet, target_planet, encounter_time, encounter_anomaly_phi,
                    target_current_anomaly_beta, incoming_e, incoming_theta,
                    max_transfer_revolution, max_target_revolution,
                    depth + 1, cache, builder, pairs, stats);
                return;
            }
        }
        collect_stable_sign_change_pairs(
            left,
            right,
            departure_planet,
            target_planet,
            encounter_time,
            encounter_anomaly_phi,
            target_current_anomaly_beta,
            incoming_e,
            incoming_theta,
            max_transfer_revolution,
            max_target_revolution,
            pairs,
            stats);
        return;
    }
    if (width < kTopologyEpsilon || depth >= kMaxRecursionDepth) {
        return;
    }

    const double theta_mid = 0.5 * (left.theta_prime + right.theta_prime);
    const ThetaPrimeSampleSummary mid = evaluate_cached_or_build(cache, theta_mid, builder);
    collect_pairs_from_interval(
        left, mid,
        departure_planet, target_planet, encounter_time, encounter_anomaly_phi,
        target_current_anomaly_beta, incoming_e, incoming_theta,
        max_transfer_revolution, max_target_revolution,
        depth + 1, cache, builder, pairs, stats);
    collect_pairs_from_interval(
        mid, right,
        departure_planet, target_planet, encounter_time, encounter_anomaly_phi,
        target_current_anomaly_beta, incoming_e, incoming_theta,
        max_transfer_revolution, max_target_revolution,
        depth + 1, cache, builder, pairs, stats);
}

Problem2SlingshotResidualBranchTestOnly
evaluate_problem2_slingshot_residual_with_branch_hint_route_a_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    int max_transfer_revolution,
    int max_target_revolution,
    int k,
    double branch_hint_alpha,
    double branch_hint_time_of_flight,
    double branch_hint_e_prime,
    double branch_hint_residual
) {
    Problem2SlingshotResidualBranchTestOnly result{};
    result.theta_prime = theta_prime;
    result.transfer_revolution = k;
    result.invalid_reason = "branch_hint_no_same_k_candidate";

    const auto branches = evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
        departure_planet,
        target_planet,
        encounter_time,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        incoming_e,
        incoming_theta,
        theta_prime,
        max_transfer_revolution,
        max_target_revolution);
    const auto by_k = group_valid_branches_by_k_sorted(branches);
    const auto it = by_k.find(k);
    if (it == by_k.end()) {
        return result;
    }

    constexpr double kTimeScale = 1e7;
    constexpr double kEScale = 1.0;
    constexpr double kResidualScale = 1.0;

    double best_score = std::numeric_limits<double>::infinity();
    for (const auto& candidate : it->second) {
        const double alpha_score = wrapped_angle_distance(candidate.alpha, branch_hint_alpha);
        const double time_score =
            std::abs(candidate.time_of_flight_seconds - branch_hint_time_of_flight) / kTimeScale;
        const double e_score = std::abs(candidate.outgoing_eccentricity - branch_hint_e_prime) / kEScale;
        const double residual_score = std::abs(candidate.slingshot_residual - branch_hint_residual) / kResidualScale;
        const double score = alpha_score + time_score + e_score + residual_score;
        if (score < best_score) {
            best_score = score;
            result = candidate;
        }
    }
    return result;
}

Problem2SlingshotResidualBranchTestOnly
continue_problem2_branch_to_theta_prime_route_a_test_only(
    const Problem2SlingshotResidualBranchTestOnly& endpoint_branch,
    double theta_prime_new,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution,
    ContinuationAttemptDiagnostic* diagnostic
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    Problem2SlingshotResidualBranchTestOnly result{};
    result.theta_prime = theta_prime_new;
    result.transfer_revolution = endpoint_branch.transfer_revolution;
    result.target_revolution = endpoint_branch.target_revolution;
    result.alpha = endpoint_branch.alpha;
    if (diagnostic != nullptr) {
        diagnostic->endpoint_theta_prime = endpoint_branch.theta_prime;
        diagnostic->theta_prime_new = theta_prime_new;
        diagnostic->endpoint_alpha = endpoint_branch.alpha;
        diagnostic->endpoint_encounter_global_angle = endpoint_branch.encounter_global_angle;
        diagnostic->d_alpha_d_theta_A = endpoint_branch.d_alpha_d_theta_A;
        diagnostic->endpoint_q = endpoint_branch.target_revolution;
    }
    if (!endpoint_branch.derivatives_available) {
        result.invalid_reason = "endpoint_derivatives_unavailable";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = result.invalid_reason;
        }
        return result;
    }

    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A_new = normalize_angle_0_2pi(departure_state.theta_global - theta_prime_new);
    const double delta_theta_prime = normalize_angle_minus_pi_pi(theta_prime_new - endpoint_branch.theta_prime);
    const double encounter_alpha_seed = normalize_angle_0_2pi(
        endpoint_branch.encounter_global_angle - endpoint_branch.d_alpha_d_theta_A * delta_theta_prime);
    if (diagnostic != nullptr) {
        diagnostic->delta_theta_prime = delta_theta_prime;
        diagnostic->predicted_encounter_alpha_seed = encounter_alpha_seed;
    }

    const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
        departure_planet,
        target_planet,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        theta_A_new,
        endpoint_branch.transfer_revolution,
        endpoint_branch.target_revolution,
        encounter_alpha_seed,
        80,
        1e-2,
        1e-12,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    if (diagnostic != nullptr) {
        diagnostic->refined_valid = refined.valid && refined.branch.valid;
        diagnostic->refined_residual_seconds = refined.branch.residual_seconds;
        diagnostic->refined_encounter_global_angle = refined.branch.encounter_global_angle;
        diagnostic->refined_target_arrival_true_anomaly = refined.branch.target_arrival_true_anomaly;
    }
    if (!refined.valid || !refined.branch.valid) {
        result.invalid_reason = "continuation_refine_failed";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = result.invalid_reason;
            problem1::Problem1RootRefinementResult best_refined{};
            double best_abs_residual = std::numeric_limits<double>::infinity();
            for (int q = 0; q <= max_target_revolution; ++q) {
                const auto q_refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
                    departure_planet,
                    target_planet,
                    encounter_anomaly_phi,
                    target_current_anomaly_beta,
                    theta_A_new,
                    endpoint_branch.transfer_revolution,
                    q,
                    encounter_alpha_seed,
                    80,
                    1e-2,
                    1e-12,
                    problem1::Problem1RootDerivativeMode::AnalyticOnly,
                    1e-6);
                if (q_refined.valid && q_refined.branch.valid) {
                    const double abs_res = std::abs(q_refined.branch.residual_seconds);
                    if (abs_res < best_abs_residual) {
                        best_abs_residual = abs_res;
                        best_refined = q_refined;
                        diagnostic->best_q = q;
                        diagnostic->best_q_valid = true;
                        diagnostic->best_q_residual_seconds = q_refined.branch.residual_seconds;
                        diagnostic->best_q_alpha = q_refined.branch.target_arrival_true_anomaly;
                        diagnostic->best_q_time_of_flight_seconds = q_refined.branch.time_of_flight_seconds;
                    }
                }
            }
            diagnostic->best_q_same_as_endpoint = diagnostic->best_q == diagnostic->endpoint_q;
        }
        return result;
    }
    if (std::abs(refined.branch.residual_seconds) > 1e-2) {
        result.invalid_reason = "continuation_problem1_residual_too_large";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = result.invalid_reason;
        }
        return result;
    }

    const auto differentiated = problem1::attach_problem1_root_derivatives_with_mode(
        departure_planet,
        target_planet,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        theta_A_new,
        refined.branch,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    if (!differentiated.valid || !differentiated.derivatives_available) {
        result.invalid_reason = "continuation_derivatives_unavailable";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = result.invalid_reason;
        }
        return result;
    }

    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);
    const auto theta_alpha_residual =
        problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            encounter_anomaly_phi,
            differentiated.target_arrival_true_anomaly,
            incoming_e,
            incoming_theta,
            theta_prime_new);
    if (!theta_alpha_residual.valid) {
        result.invalid_reason = "continuation_slingshot_geometry_invalid";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = result.invalid_reason;
        }
        return result;
    }

    result.valid = true;
    result.alpha = differentiated.target_arrival_true_anomaly;
    result.encounter_global_angle = differentiated.encounter_global_angle;
    result.transfer_revolution = differentiated.transfer_revolution;
    result.target_revolution = differentiated.target_revolution;
    result.time_of_flight_seconds = differentiated.time_of_flight_seconds;
    result.target_time_seconds = differentiated.target_time_seconds;
    result.problem1_residual_seconds = differentiated.residual_seconds;
    result.outgoing_eccentricity = theta_alpha_residual.outgoing_eccentricity;
    result.outgoing_semi_latus_rectum = theta_alpha_residual.outgoing_semi_latus_rectum;
    result.incoming_invariant = theta_alpha_residual.incoming_invariant;
    result.outgoing_invariant = theta_alpha_residual.outgoing_invariant;
    result.slingshot_residual = theta_alpha_residual.slingshot_residual;
    result.derivatives_available = true;
    result.d_alpha_d_theta_A = differentiated.d_encounter_global_angle_d_theta_A;
    if (diagnostic != nullptr) {
        diagnostic->valid = true;
    }
    return result;
}

TwoSidedContinuationDiagnostic
evaluate_midpoint_branch_by_two_sided_continuation_test_only(
    double theta_mid,
    const Problem2SlingshotResidualBranchTestOnly& left_branch,
    const Problem2SlingshotResidualBranchTestOnly& right_branch,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution
) {
    TwoSidedContinuationDiagnostic diagnostic{};
    auto left_mid = continue_problem2_branch_to_theta_prime_route_a_test_only(
        left_branch,
        theta_mid,
        departure_planet,
        target_planet,
        encounter_time,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        incoming_e,
        incoming_theta,
        max_transfer_revolution,
        max_target_revolution,
        &diagnostic.left_continuation);
    auto right_mid = continue_problem2_branch_to_theta_prime_route_a_test_only(
        right_branch,
        theta_mid,
        departure_planet,
        target_planet,
        encounter_time,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        incoming_e,
        incoming_theta,
        max_transfer_revolution,
        max_target_revolution,
        &diagnostic.right_continuation);

    diagnostic.left_valid = left_mid.valid;
    diagnostic.left_invalid_reason = left_mid.invalid_reason;
    diagnostic.right_valid = right_mid.valid;
    diagnostic.right_invalid_reason = right_mid.invalid_reason;
    if (!left_mid.valid && !right_mid.valid) {
        diagnostic.invalid_reason = "two_sided_continuation_failed";
        return diagnostic;
    }
    if (!left_mid.valid || !right_mid.valid) {
        diagnostic.invalid_reason = "one_sided_continuation_only";
        return diagnostic;
    }

    diagnostic.left_alpha = left_mid.alpha;
    diagnostic.right_alpha = right_mid.alpha;
    diagnostic.alpha_diff = wrapped_angle_distance(left_mid.alpha, right_mid.alpha);
    diagnostic.left_time = left_mid.time_of_flight_seconds;
    diagnostic.right_time = right_mid.time_of_flight_seconds;
    diagnostic.time_diff = std::abs(left_mid.time_of_flight_seconds - right_mid.time_of_flight_seconds);
    diagnostic.left_e_prime = left_mid.outgoing_eccentricity;
    diagnostic.right_e_prime = right_mid.outgoing_eccentricity;
    diagnostic.e_prime_diff = std::abs(left_mid.outgoing_eccentricity - right_mid.outgoing_eccentricity);
    diagnostic.left_slingshot_residual = left_mid.slingshot_residual;
    diagnostic.right_slingshot_residual = right_mid.slingshot_residual;
    diagnostic.slingshot_residual_diff = std::abs(left_mid.slingshot_residual - right_mid.slingshot_residual);
    const bool same_branch =
        left_mid.transfer_revolution == right_mid.transfer_revolution &&
        diagnostic.alpha_diff <= 1e-5 &&
        diagnostic.time_diff <= 1e3 &&
        diagnostic.e_prime_diff <= 1e-4 &&
        diagnostic.slingshot_residual_diff <= 1e-4;
    if (!same_branch) {
        diagnostic.invalid_reason = "two_sided_continuation_mismatch";
        return diagnostic;
    }
    diagnostic.valid = true;
    diagnostic.selected_branch =
        std::abs(left_mid.problem1_residual_seconds) <= std::abs(right_mid.problem1_residual_seconds) ? left_mid
                                                                                                       : right_mid;
    return diagnostic;
}

Problem2StableBranchBisectionSolution bisect_sign_change_pair(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    int max_transfer_revolution,
    int max_target_revolution,
    const Problem2StableSignChangePair& pair
) {
    Problem2StableBranchBisectionSolution solution{};
    solution.k = pair.k;
    solution.rank_in_k = pair.rank_in_k;
    solution.theta_left_initial = pair.theta_left;
    solution.theta_right_initial = pair.theta_right;

    double theta_left = pair.theta_left;
    double theta_right = pair.theta_right;
    Problem2SlingshotResidualBranchTestOnly left_branch = pair.left_branch;
    Problem2SlingshotResidualBranchTestOnly right_branch = pair.right_branch;
    double residual_left = left_branch.slingshot_residual;
    for (int iteration = 0; iteration < kBisectionMaxIterations; ++iteration) {
        const double theta_mid = 0.5 * (theta_left + theta_right);
        const auto mid = evaluate_midpoint_branch_by_two_sided_continuation_test_only(
            theta_mid,
            left_branch,
            right_branch,
            departure_planet,
            target_planet,
            encounter_time,
            encounter_anomaly_phi,
            target_current_anomaly_beta,
            incoming_e,
            incoming_theta,
            max_transfer_revolution,
            max_target_revolution);
        solution.iterations = iteration + 1;
        solution.last_theta = theta_mid;
        if (!mid.valid) {
            solution.invalid_reason = mid.invalid_reason;
            return solution;
        }
        const auto& mid_branch = mid.selected_branch;

        solution.last_residual = mid_branch.slingshot_residual;
        if (std::abs(mid_branch.slingshot_residual) <= kResidualTolerance) {
            solution.valid = true;
            solution.theta_root = theta_mid;
            solution.alpha_root = mid_branch.alpha;
            solution.outgoing_eccentricity = mid_branch.outgoing_eccentricity;
            solution.outgoing_semi_latus_rectum = mid_branch.outgoing_semi_latus_rectum;
            solution.slingshot_residual = mid_branch.slingshot_residual;
            solution.final_width = theta_right - theta_left;
            solution.time_of_flight_seconds = mid_branch.time_of_flight_seconds;
            solution.target_revolution = mid_branch.target_revolution;
            solution.problem1_residual_seconds = mid_branch.problem1_residual_seconds;
            return solution;
        }
        if ((theta_right - theta_left) <= kThetaTolerance) {
            solution.invalid_reason = "bisection_width_converged_but_residual_too_large";
            solution.final_width = theta_right - theta_left;
            return solution;
        }

        const bool left_mid_sign_change =
            (residual_left < 0.0 && mid_branch.slingshot_residual > 0.0) ||
            (residual_left > 0.0 && mid_branch.slingshot_residual < 0.0);
        if (left_mid_sign_change) {
            theta_right = theta_mid;
            right_branch = mid_branch;
        } else {
            theta_left = theta_mid;
            residual_left = mid_branch.slingshot_residual;
            left_branch = mid_branch;
        }
    }

    solution.invalid_reason = "bisection_max_iterations_residual_too_large";
    solution.final_width = theta_right - theta_left;
    return solution;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;

    std::cout << std::setprecision(6) << std::scientific;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time =
        0.17 * planet_params::planet_orbital_period(departure_planet);
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double incoming_e = 0.3;
    const double incoming_theta = 0.4;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;

    auto builder = [&](double theta_prime) {
        return build_theta_summary(
            departure_planet,
            target_planet,
            encounter_time,
            phi,
            beta,
            incoming_e,
            incoming_theta,
            theta_prime,
            max_transfer_revolution,
            max_target_revolution);
    };

    std::map<double, ThetaPrimeSampleSummary> cache;
    std::vector<ThetaPrimeSampleSummary> samples;
    samples.reserve(kThetaSampleCount);
    for (int i = 0; i < kThetaSampleCount; ++i) {
        const double theta_prime = kTwoPi * static_cast<double>(i) / static_cast<double>(kThetaSampleCount);
        samples.push_back(evaluate_cached_or_build(&cache, theta_prime, builder));
    }

    std::vector<Problem2StableSignChangePair> sign_change_pairs;
    PairCollectionStats pair_stats{};
    for (int i = 0; i + 1 < kThetaSampleCount; ++i) {
        collect_pairs_from_interval(
            samples[static_cast<std::size_t>(i)],
            samples[static_cast<std::size_t>(i + 1)],
            departure_planet,
            target_planet,
            encounter_time,
            phi,
            beta,
            incoming_e,
            incoming_theta,
            max_transfer_revolution,
            max_target_revolution,
            0,
            &cache,
            builder,
            &sign_change_pairs,
            &pair_stats);
    }

    int strict_success_count = 0;
    int failure_count = 0;
    int endpoint_derivatives_unavailable_count = 0;
    int continuation_refine_failed_count = 0;
    int continuation_derivatives_unavailable_count = 0;
    int continuation_slingshot_geometry_invalid_count = 0;
    int two_sided_continuation_failed_count = 0;
    int one_sided_continuation_only_count = 0;
    int two_sided_continuation_mismatch_count = 0;
    int branch_hint_missing_during_bisection_count = 0;
    int width_converged_residual_too_large_count = 0;
    int max_iterations_residual_too_large_count = 0;
    double max_abs_slingshot_residual_at_root = 0.0;
    double max_abs_problem1_residual_seconds_at_root = 0.0;
    double min_theta_root = std::numeric_limits<double>::infinity();
    double max_theta_root = -std::numeric_limits<double>::infinity();
    std::vector<Problem2StableBranchBisectionSolution> strict_successes;
    for (const auto& pair : sign_change_pairs) {
        const auto solution = bisect_sign_change_pair(
            departure_planet,
            target_planet,
            encounter_time,
            phi,
            beta,
            incoming_e,
            incoming_theta,
            max_transfer_revolution,
            max_target_revolution,
            pair);
        if (solution.valid) {
            strict_success_count += 1;
            strict_successes.push_back(solution);
            max_abs_slingshot_residual_at_root =
                std::max(max_abs_slingshot_residual_at_root, std::abs(solution.slingshot_residual));
            max_abs_problem1_residual_seconds_at_root =
                std::max(max_abs_problem1_residual_seconds_at_root, std::abs(solution.problem1_residual_seconds));
            min_theta_root = std::min(min_theta_root, solution.theta_root);
            max_theta_root = std::max(max_theta_root, solution.theta_root);

            std::cout << "Problem2StableBranchBisectionSolution\n";
            std::cout << "valid=1\n";
            std::cout << "k=" << solution.k << '\n';
            std::cout << "rank_in_k=" << solution.rank_in_k << '\n';
            std::cout << "theta_left_initial=" << solution.theta_left_initial << '\n';
            std::cout << "theta_right_initial=" << solution.theta_right_initial << '\n';
            std::cout << "theta_root=" << solution.theta_root << '\n';
            std::cout << "alpha_root=" << solution.alpha_root << '\n';
            std::cout << "target_revolution=" << solution.target_revolution << '\n';
            std::cout << "time_of_flight_seconds=" << solution.time_of_flight_seconds << '\n';
            std::cout << "outgoing_eccentricity=" << solution.outgoing_eccentricity << '\n';
            std::cout << "outgoing_semi_latus_rectum=" << solution.outgoing_semi_latus_rectum << '\n';
            std::cout << "slingshot_residual=" << solution.slingshot_residual << '\n';
            std::cout << "problem1_residual_seconds=" << solution.problem1_residual_seconds << '\n';
            std::cout << "iterations=" << solution.iterations << '\n';
            std::cout << "final_width=" << solution.final_width << '\n';
        } else {
            failure_count += 1;
            if (solution.invalid_reason == "endpoint_derivatives_unavailable") {
                endpoint_derivatives_unavailable_count += 1;
            } else if (solution.invalid_reason == "continuation_refine_failed") {
                continuation_refine_failed_count += 1;
            } else if (solution.invalid_reason == "continuation_derivatives_unavailable") {
                continuation_derivatives_unavailable_count += 1;
            } else if (solution.invalid_reason == "continuation_slingshot_geometry_invalid") {
                continuation_slingshot_geometry_invalid_count += 1;
            } else if (solution.invalid_reason == "two_sided_continuation_failed") {
                two_sided_continuation_failed_count += 1;
            } else if (solution.invalid_reason == "one_sided_continuation_only") {
                one_sided_continuation_only_count += 1;
            } else if (solution.invalid_reason == "two_sided_continuation_mismatch") {
                two_sided_continuation_mismatch_count += 1;
            } else if (solution.invalid_reason == "branch_hint_missing_during_bisection") {
                branch_hint_missing_during_bisection_count += 1;
            } else if (solution.invalid_reason == "bisection_width_converged_but_residual_too_large") {
                width_converged_residual_too_large_count += 1;
            } else if (solution.invalid_reason == "bisection_max_iterations_residual_too_large") {
                max_iterations_residual_too_large_count += 1;
            }
            std::cout << "Problem2StableBranchBisectionFailure\n";
            std::cout << "k=" << solution.k << '\n';
            std::cout << "rank_in_k=" << solution.rank_in_k << '\n';
            std::cout << "theta_left_initial=" << solution.theta_left_initial << '\n';
            std::cout << "theta_right_initial=" << solution.theta_right_initial << '\n';
            std::cout << "invalid_reason=" << solution.invalid_reason << '\n';
            std::cout << "last_theta=" << solution.last_theta << '\n';
            std::cout << "last_residual=" << solution.last_residual << '\n';
        }
    }

    std::vector<Problem2StableBranchBisectionSolution> dedup_successes;
    for (const auto& solution : strict_successes) {
        bool duplicate = false;
        for (const auto& existing : dedup_successes) {
            if (same_root_for_dedup(solution, existing)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            dedup_successes.push_back(solution);
        }
    }

    std::cout << "Problem2StableBranchBisectionSummary\n";
    std::cout << "theta_sample_count=" << kThetaSampleCount << '\n';
    std::cout << "raw_sign_change_pair_count=" << pair_stats.raw_sign_change_pair_count << '\n';
    std::cout << "continuation_stable_sign_change_pair_count="
              << pair_stats.continuation_stable_sign_change_pair_count << '\n';
    std::cout << "unstable_candidate_pair_count=" << pair_stats.unstable_candidate_pair_count << '\n';
    std::cout << "strict_success_count=" << strict_success_count << '\n';
    std::cout << "raw_success_count=" << strict_success_count << '\n';
    std::cout << "dedup_success_count=" << dedup_successes.size() << '\n';
    std::cout << "bisection_failure_count=" << failure_count << '\n';
    std::cout << "endpoint_derivatives_unavailable_count=" << endpoint_derivatives_unavailable_count << '\n';
    std::cout << "continuation_refine_failed_count=" << continuation_refine_failed_count << '\n';
    std::cout << "continuation_derivatives_unavailable_count=" << continuation_derivatives_unavailable_count << '\n';
    std::cout << "continuation_slingshot_geometry_invalid_count="
              << continuation_slingshot_geometry_invalid_count << '\n';
    std::cout << "midpoint_left_continuation_invalid_count="
              << pair_stats.midpoint_left_continuation_invalid_count << '\n';
    std::cout << "midpoint_right_continuation_invalid_count="
              << pair_stats.midpoint_right_continuation_invalid_count << '\n';
    std::cout << "midpoint_both_continuation_invalid_count="
              << pair_stats.midpoint_both_continuation_invalid_count << '\n';
    std::cout << "midpoint_one_sided_only_count=" << pair_stats.midpoint_one_sided_only_count << '\n';
    std::cout << "midpoint_two_sided_mismatch_count=" << pair_stats.midpoint_two_sided_mismatch_count << '\n';
    std::cout << "midpoint_continuation_refine_failed_count="
              << pair_stats.midpoint_continuation_refine_failed_count << '\n';
    std::cout << "midpoint_continuation_problem1_residual_too_large_count="
              << pair_stats.midpoint_continuation_problem1_residual_too_large_count << '\n';
    std::cout << "midpoint_continuation_derivatives_unavailable_count="
              << pair_stats.midpoint_continuation_derivatives_unavailable_count << '\n';
    std::cout << "midpoint_slingshot_geometry_invalid_count="
              << pair_stats.midpoint_slingshot_geometry_invalid_count << '\n';
    std::cout << "midpoint_endpoint_derivatives_unavailable_count="
              << pair_stats.midpoint_endpoint_derivatives_unavailable_count << '\n';
    std::cout << "q_diagnostic_endpoint_q_best_count="
              << pair_stats.q_diagnostic_endpoint_q_best_count << '\n';
    std::cout << "q_diagnostic_changed_q_best_count="
              << pair_stats.q_diagnostic_changed_q_best_count << '\n';
    std::cout << "two_sided_continuation_failed_count=" << two_sided_continuation_failed_count << '\n';
    std::cout << "one_sided_continuation_only_count=" << one_sided_continuation_only_count << '\n';
    std::cout << "two_sided_continuation_mismatch_count=" << two_sided_continuation_mismatch_count << '\n';
    std::cout << "branch_hint_missing_during_bisection_count="
              << branch_hint_missing_during_bisection_count << '\n';
    std::cout << "width_converged_residual_too_large_count="
              << width_converged_residual_too_large_count << '\n';
    std::cout << "max_iterations_residual_too_large_count="
              << max_iterations_residual_too_large_count << '\n';
    std::cout << "max_abs_slingshot_residual_at_strict_root=" << max_abs_slingshot_residual_at_root << '\n';
    std::cout << "max_abs_problem1_residual_seconds_at_strict_root="
              << max_abs_problem1_residual_seconds_at_root << '\n';
    std::cout << "min_theta_root=" << min_theta_root << '\n';
    std::cout << "max_theta_root=" << max_theta_root << '\n';

    assert(pair_stats.raw_sign_change_pair_count > 0);
    if (strict_success_count > 0) {
        assert(std::isfinite(max_abs_slingshot_residual_at_root));
        assert(std::isfinite(max_abs_problem1_residual_seconds_at_root));
    }

    std::cout << "problem2_stable_branch_bisection_route_a_ok\n";
    return 0;
}
