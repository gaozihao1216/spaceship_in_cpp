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

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr int kThetaSampleCount = 32;
constexpr int kMaxRecursionDepth = 4;
constexpr double kTopologyEpsilon = 1e-5;

struct Problem2SlingshotResidualBranchTestOnly {
    bool valid = false;
    std::string invalid_reason;

    double theta_prime = 0.0;
    double alpha = 0.0;

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
};

struct ThetaPrimeSampleSummary {
    double theta_prime = 0.0;
    int total_problem1_branch_count = 0;
    int valid_slingshot_branch_count = 0;
    std::map<int, int> valid_count_by_k;
    std::map<int, int> problem1_count_by_k;
    std::vector<Problem2SlingshotResidualBranchTestOnly> branches;
};

struct TopologyStats {
    int stable_interval_count = 0;
    int topology_change_interval_count = 0;
    int topology_skipped_count = 0;
    int topology_unresolved_count = 0;
    int stable_branch_pair_count = 0;
    int sign_change_pair_count = 0;
};

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

std::string map_to_string(const std::map<int, int>& counts) {
    std::string out = "{";
    bool first = true;
    for (const auto& [k, count] : counts) {
        if (!first) {
            out += ",";
        }
        first = false;
        out += "k" + std::to_string(k) + "=" + std::to_string(count);
    }
    out += "}";
    return out;
}

void print_stable_branch_pair_preview(
    const ThetaPrimeSampleSummary& left,
    const ThetaPrimeSampleSummary& right,
    TopologyStats* stats
) {
    std::map<int, std::vector<Problem2SlingshotResidualBranchTestOnly>> left_by_k;
    std::map<int, std::vector<Problem2SlingshotResidualBranchTestOnly>> right_by_k;
    for (const auto& branch : left.branches) {
        if (branch.valid) {
            left_by_k[branch.transfer_revolution].push_back(branch);
        }
    }
    for (const auto& branch : right.branches) {
        if (branch.valid) {
            right_by_k[branch.transfer_revolution].push_back(branch);
        }
    }
    for (auto& [k, branches] : left_by_k) {
        auto& other = right_by_k[k];
        std::sort(
            branches.begin(),
            branches.end(),
            [](const auto& a, const auto& b) { return a.time_of_flight_seconds < b.time_of_flight_seconds; });
        std::sort(
            other.begin(),
            other.end(),
            [](const auto& a, const auto& b) { return a.time_of_flight_seconds < b.time_of_flight_seconds; });
        const int pair_count = std::min<int>(branches.size(), other.size());
        stats->stable_branch_pair_count += pair_count;
        const int preview_count = std::min(3, pair_count);
        for (int rank = 0; rank < pair_count; ++rank) {
            const auto& l = branches[static_cast<std::size_t>(rank)];
            const auto& r = other[static_cast<std::size_t>(rank)];
            const bool sign_change =
                std::isfinite(l.slingshot_residual) && std::isfinite(r.slingshot_residual) &&
                ((l.slingshot_residual < 0.0 && r.slingshot_residual > 0.0) ||
                 (l.slingshot_residual > 0.0 && r.slingshot_residual < 0.0));
            if (sign_change) {
                stats->sign_change_pair_count += 1;
            }
            if (rank >= preview_count) {
                continue;
            }
            std::cout << "Problem2StableBranchPairPreview\n";
            std::cout << "theta_left=" << left.theta_prime << '\n';
            std::cout << "theta_right=" << right.theta_prime << '\n';
            std::cout << "k=" << k << '\n';
            std::cout << "rank_in_k=" << rank << '\n';
            std::cout << "left_alpha=" << l.alpha << '\n';
            std::cout << "right_alpha=" << r.alpha << '\n';
            std::cout << "left_time_of_flight=" << l.time_of_flight_seconds << '\n';
            std::cout << "right_time_of_flight=" << r.time_of_flight_seconds << '\n';
            std::cout << "left_slingshot_residual=" << l.slingshot_residual << '\n';
            std::cout << "right_slingshot_residual=" << r.slingshot_residual << '\n';
            std::cout << "residual_sign_change=" << (sign_change ? 1 : 0) << '\n';
            std::cout << "left_e_prime=" << l.outgoing_eccentricity << '\n';
            std::cout << "right_e_prime=" << r.outgoing_eccentricity << '\n';
        }
    }
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

template <class SummaryBuilder>
void process_interval(
    const ThetaPrimeSampleSummary& left,
    const ThetaPrimeSampleSummary& right,
    int depth,
    std::map<double, ThetaPrimeSampleSummary>* cache,
    SummaryBuilder&& builder,
    TopologyStats* stats
) {
    const double width = right.theta_prime - left.theta_prime;
    if (same_count_by_k(left, right)) {
        stats->stable_interval_count += 1;
        std::cout << "Problem2StableThetaInterval\n";
        std::cout << "theta_left=" << left.theta_prime << '\n';
        std::cout << "theta_right=" << right.theta_prime << '\n';
        std::cout << "width=" << width << '\n';
        std::cout << "count_by_k=" << map_to_string(left.valid_count_by_k) << '\n';
        std::cout << "total_valid_count=" << left.valid_slingshot_branch_count << '\n';
        print_stable_branch_pair_preview(left, right, stats);
        return;
    }

    stats->topology_change_interval_count += 1;
    int max_count_difference_by_k = 0;
    std::map<int, int> merged = left.valid_count_by_k;
    for (const auto& [k, count] : right.valid_count_by_k) {
        merged.emplace(k, count);
    }
    for (const auto& [k, _] : merged) {
        const int left_count = left.valid_count_by_k.count(k) ? left.valid_count_by_k.at(k) : 0;
        const int right_count = right.valid_count_by_k.count(k) ? right.valid_count_by_k.at(k) : 0;
        max_count_difference_by_k = std::max(max_count_difference_by_k, std::abs(left_count - right_count));
    }
    std::cout << "Problem2TopologyChangeThetaInterval\n";
    std::cout << "theta_left=" << left.theta_prime << '\n';
    std::cout << "theta_right=" << right.theta_prime << '\n';
    std::cout << "width=" << width << '\n';
    std::cout << "left_count_by_k=" << map_to_string(left.valid_count_by_k) << '\n';
    std::cout << "right_count_by_k=" << map_to_string(right.valid_count_by_k) << '\n';
    std::cout << "left_total_valid_count=" << left.valid_slingshot_branch_count << '\n';
    std::cout << "right_total_valid_count=" << right.valid_slingshot_branch_count << '\n';
    std::cout << "count_difference_total="
              << std::abs(left.valid_slingshot_branch_count - right.valid_slingshot_branch_count) << '\n';
    std::cout << "max_count_difference_by_k=" << max_count_difference_by_k << '\n';

    if (width < kTopologyEpsilon) {
        stats->topology_skipped_count += 1;
        std::cout << "Problem2TopologyTransitionSkipped\n";
        std::cout << "theta_left=" << left.theta_prime << '\n';
        std::cout << "theta_right=" << right.theta_prime << '\n';
        std::cout << "width=" << width << '\n';
        std::cout << "left_count_by_k=" << map_to_string(left.valid_count_by_k) << '\n';
        std::cout << "right_count_by_k=" << map_to_string(right.valid_count_by_k) << '\n';
        std::cout << "reason=topology_interval_too_small\n";
        return;
    }

    const double theta_mid = 0.5 * (left.theta_prime + right.theta_prime);
    const ThetaPrimeSampleSummary mid =
        evaluate_cached_or_build(cache, theta_mid, builder);
    std::string likely_change_side = "both_or_unclear";
    if (same_count_by_k(left, mid) && !same_count_by_k(mid, right)) {
        likely_change_side = "right_half";
    } else if (!same_count_by_k(left, mid) && same_count_by_k(mid, right)) {
        likely_change_side = "left_half";
    }
    std::cout << "Problem2TopologyMidpointDiagnostic\n";
    std::cout << "theta_left=" << left.theta_prime << '\n';
    std::cout << "theta_mid=" << mid.theta_prime << '\n';
    std::cout << "theta_right=" << right.theta_prime << '\n';
    std::cout << "left_count_by_k=" << map_to_string(left.valid_count_by_k) << '\n';
    std::cout << "mid_count_by_k=" << map_to_string(mid.valid_count_by_k) << '\n';
    std::cout << "right_count_by_k=" << map_to_string(right.valid_count_by_k) << '\n';
    std::cout << "left_total=" << left.valid_slingshot_branch_count << '\n';
    std::cout << "mid_total=" << mid.valid_slingshot_branch_count << '\n';
    std::cout << "right_total=" << right.valid_slingshot_branch_count << '\n';
    std::cout << "likely_change_side=" << likely_change_side << '\n';

    if (depth >= kMaxRecursionDepth) {
        stats->topology_unresolved_count += 1;
        std::cout << "Problem2TopologyTransitionUnresolved\n";
        std::cout << "theta_left=" << left.theta_prime << '\n';
        std::cout << "theta_right=" << right.theta_prime << '\n';
        std::cout << "width=" << width << '\n';
        std::cout << "depth=" << depth << '\n';
        std::cout << "left_count_by_k=" << map_to_string(left.valid_count_by_k) << '\n';
        std::cout << "right_count_by_k=" << map_to_string(right.valid_count_by_k) << '\n';
        return;
    }

    process_interval(left, mid, depth + 1, cache, builder, stats);
    process_interval(mid, right, depth + 1, cache, builder, stats);
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
        ThetaPrimeSampleSummary summary = evaluate_cached_or_build(&cache, theta_prime, builder);
        samples.push_back(summary);
    }

    TopologyStats stats{};
    for (int i = 0; i + 1 < kThetaSampleCount; ++i) {
        process_interval(samples[static_cast<std::size_t>(i)],
                         samples[static_cast<std::size_t>(i + 1)],
                         0,
                         &cache,
                         builder,
                         &stats);
    }

    std::cout << "Problem2ThetaPrimeTopologySummary\n";
    std::cout << "theta_sample_count=" << kThetaSampleCount << '\n';
    std::cout << "stable_interval_count=" << stats.stable_interval_count << '\n';
    std::cout << "topology_change_interval_count=" << stats.topology_change_interval_count << '\n';
    std::cout << "topology_skipped_count=" << stats.topology_skipped_count << '\n';
    std::cout << "topology_unresolved_count=" << stats.topology_unresolved_count << '\n';
    std::cout << "stable_branch_pair_count=" << stats.stable_branch_pair_count << '\n';
    std::cout << "sign_change_pair_count=" << stats.sign_change_pair_count << '\n';

    std::cout << "problem2_theta_prime_branch_count_diagnostic_ok\n";
    return 0;
}
