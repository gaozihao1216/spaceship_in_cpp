#pragma once

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace problem1_5deg_test {

using Clock = std::chrono::steady_clock;

constexpr double kFiveDegRadians = 5.0 * spaceship_cpp::common::kPi / 180.0;
constexpr double kMatchAngleTolerance = 1e-6;
constexpr double kMatchTimeToleranceSeconds = 1e3;
constexpr int kMaxNewtonIterations = 80;
constexpr double kResidualToleranceSeconds = 1e-2;
constexpr double kAlphaStepTolerance = 1e-12;
constexpr double kDerivativeFiniteDifferenceStep = 1e-6;

struct PairCase {
    spaceship_cpp::planet_params::PlanetId from;
    spaceship_cpp::planet_params::PlanetId to;
};

struct QueryPoint {
    spaceship_cpp::planet_params::PlanetId from;
    spaceship_cpp::planet_params::PlanetId to;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
};

struct BranchMatchStats {
    int direct_branch_count = 0;
    int candidate_branch_count = 0;
    int matched_branch_count = 0;
    int missing_count = 0;
    int new_count = 0;
    int branch_order_mismatch_count = 0;
    double theta_prime_error_sum = 0.0;
    double theta_prime_error_max = 0.0;
    double alpha_error_sum = 0.0;
    double alpha_error_max = 0.0;
    double arrival_time_error_sum = 0.0;
    double arrival_time_error_max = 0.0;
};

struct SimulatedFiveDegResult {
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> seed_branches;
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> refined_branches;
    int refine_attempt_count = 0;
    int refine_success_count = 0;
    int refine_fail_count = 0;
    int derivative_attach_attempt_count = 0;
    int derivative_attach_success_count = 0;
    double seed_ms = 0.0;
    double refine_ms = 0.0;
    double derivative_attach_ms = 0.0;
    double dedup_ms = 0.0;
};

inline int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    const int value = std::atoi(raw);
    return value > 0 ? value : fallback;
}

inline unsigned env_uint(const char* name, unsigned fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    const unsigned long value = std::strtoul(raw, nullptr, 10);
    return value > 0 ? static_cast<unsigned>(value) : fallback;
}

inline double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

inline double angle_distance(double lhs, double rhs) {
    double diff = std::fmod(lhs - rhs + spaceship_cpp::common::kPi, spaceship_cpp::common::kTwoPi);
    if (diff < 0.0) diff += spaceship_cpp::common::kTwoPi;
    return std::abs(diff - spaceship_cpp::common::kPi);
}

inline double round_to_5deg(double angle) {
    const double wrapped = spaceship_cpp::common::normalize_angle_0_2pi(angle);
    const double index = std::round(wrapped / kFiveDegRadians);
    return spaceship_cpp::common::normalize_angle_0_2pi(index * kFiveDegRadians);
}

inline QueryPoint nearest_5deg_query(const QueryPoint& query) {
    QueryPoint out = query;
    out.nu_A = round_to_5deg(query.nu_A);
    out.nu_B = round_to_5deg(query.nu_B);
    out.theta_A = round_to_5deg(query.theta_A);
    return out;
}

inline std::vector<PairCase> default_pairs() {
    namespace planet = spaceship_cpp::planet_params;
    return {
        {planet::PlanetId::Earth, planet::PlanetId::Mars},
        {planet::PlanetId::Mars, planet::PlanetId::Earth},
        {planet::PlanetId::Earth, planet::PlanetId::Mercury},
        {planet::PlanetId::Venus, planet::PlanetId::Mars},
        {planet::PlanetId::Venus, planet::PlanetId::Mercury},
        {planet::PlanetId::Mars, planet::PlanetId::Mercury},
    };
}

inline std::vector<PairCase> selected_pairs_from_env() {
    auto pairs = default_pairs();
    const int pair_limit = env_int("ROOT_TABLE_5DEG_TEST_PAIR_LIMIT", static_cast<int>(pairs.size()));
    if (pair_limit < static_cast<int>(pairs.size())) {
        pairs.resize(static_cast<std::size_t>(std::max(1, pair_limit)));
    }
    return pairs;
}

inline std::vector<QueryPoint> make_random_queries_for_pair(
    const PairCase& pair,
    int count,
    std::mt19937_64* rng
) {
    std::uniform_real_distribution<double> angle_dist(0.0, spaceship_cpp::common::kTwoPi);
    std::vector<QueryPoint> queries;
    queries.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        queries.push_back(QueryPoint{
            pair.from,
            pair.to,
            angle_dist(*rng),
            angle_dist(*rng),
            angle_dist(*rng)});
    }
    return queries;
}

inline std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> solve_direct(const QueryPoint& query) {
    return spaceship_cpp::problem1::solve_problem1_from_departure_anomalies(
        query.from, query.to, query.nu_A, query.nu_B, query.theta_A, 1, 1);
}

inline void sort_branches(std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* branches) {
    std::sort(branches->begin(), branches->end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.target_time_seconds,
                        lhs.encounter_global_angle) <
               std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.target_time_seconds,
                        rhs.encounter_global_angle);
    });
}

inline void dedup_branches(std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* branches) {
    sort_branches(branches);
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> deduped;
    for (const auto& branch : *branches) {
        const auto found = std::find_if(deduped.begin(), deduped.end(), [&](const auto& existing) {
            return existing.transfer_revolution == branch.transfer_revolution &&
                   existing.target_revolution == branch.target_revolution &&
                   angle_distance(existing.encounter_global_angle, branch.encounter_global_angle) <= 1e-8 &&
                   std::abs(existing.target_time_seconds - branch.target_time_seconds) <= 1e-3;
        });
        if (found == deduped.end()) {
            deduped.push_back(branch);
        }
    }
    *branches = std::move(deduped);
}

inline SimulatedFiveDegResult simulate_5deg_nearest_refine(const QueryPoint& query) {
    namespace problem1 = spaceship_cpp::problem1;
    SimulatedFiveDegResult result{};
    const QueryPoint seed_query = nearest_5deg_query(query);

    const auto seed_start = Clock::now();
    result.seed_branches = solve_direct(seed_query);
    result.seed_ms = elapsed_ms(seed_start, Clock::now());

    const auto refine_start = Clock::now();
    for (const auto& seed_branch : result.seed_branches) {
        if (!seed_branch.valid) continue;
        result.refine_attempt_count += 1;
        const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            query.from,
            query.to,
            query.nu_A,
            query.nu_B,
            query.theta_A,
            seed_branch.transfer_revolution,
            seed_branch.target_revolution,
            seed_branch.encounter_global_angle,
            kMaxNewtonIterations,
            kResidualToleranceSeconds,
            kAlphaStepTolerance,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            kDerivativeFiniteDifferenceStep);
        if (!refined.valid || !refined.branch.valid) {
            result.refine_fail_count += 1;
            continue;
        }
        result.refine_success_count += 1;
        result.refined_branches.push_back(refined.branch);
    }
    result.refine_ms = elapsed_ms(refine_start, Clock::now());

    const auto attach_start = Clock::now();
    for (auto& branch : result.refined_branches) {
        result.derivative_attach_attempt_count += 1;
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            query.from,
            query.to,
            query.nu_A,
            query.nu_B,
            query.theta_A,
            branch,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            kDerivativeFiniteDifferenceStep);
        if (attached.valid && attached.derivatives_available) {
            branch = attached;
            result.derivative_attach_success_count += 1;
        }
    }
    result.derivative_attach_ms = elapsed_ms(attach_start, Clock::now());

    const auto dedup_start = Clock::now();
    dedup_branches(&result.refined_branches);
    result.dedup_ms = elapsed_ms(dedup_start, Clock::now());
    return result;
}

inline BranchMatchStats compare_branches(
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> direct,
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> candidate
) {
    sort_branches(&direct);
    sort_branches(&candidate);
    BranchMatchStats stats{};
    stats.direct_branch_count = static_cast<int>(direct.size());
    stats.candidate_branch_count = static_cast<int>(candidate.size());
    const int common_count = std::min<int>(direct.size(), candidate.size());
    for (int i = 0; i < common_count; ++i) {
        if (direct[static_cast<std::size_t>(i)].transfer_revolution != candidate[static_cast<std::size_t>(i)].transfer_revolution ||
            direct[static_cast<std::size_t>(i)].target_revolution != candidate[static_cast<std::size_t>(i)].target_revolution) {
            stats.branch_order_mismatch_count += 1;
        }
    }

    std::vector<bool> used(candidate.size(), false);
    for (const auto& lhs : direct) {
        int best = -1;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < candidate.size(); ++i) {
            if (used[i]) continue;
            const auto& rhs = candidate[i];
            if (lhs.transfer_revolution != rhs.transfer_revolution ||
                lhs.target_revolution != rhs.target_revolution) {
                continue;
            }
            const double theta_error = angle_distance(lhs.target_arrival_true_anomaly, rhs.target_arrival_true_anomaly);
            const double alpha_error = angle_distance(lhs.encounter_global_angle, rhs.encounter_global_angle);
            const double time_error = std::abs(lhs.target_time_seconds - rhs.target_time_seconds);
            const double score = theta_error * theta_error * 1e12 +
                                 alpha_error * alpha_error * 1e12 +
                                 (time_error / kMatchTimeToleranceSeconds) *
                                     (time_error / kMatchTimeToleranceSeconds);
            if (score < best_score) {
                best = static_cast<int>(i);
                best_score = score;
            }
        }
        if (best < 0) {
            continue;
        }
        const auto& rhs = candidate[static_cast<std::size_t>(best)];
        const double theta_error = angle_distance(lhs.target_arrival_true_anomaly, rhs.target_arrival_true_anomaly);
        const double alpha_error = angle_distance(lhs.encounter_global_angle, rhs.encounter_global_angle);
        const double time_error = std::abs(lhs.target_time_seconds - rhs.target_time_seconds);
        if (theta_error <= kMatchAngleTolerance &&
            alpha_error <= kMatchAngleTolerance &&
            time_error <= kMatchTimeToleranceSeconds) {
            used[static_cast<std::size_t>(best)] = true;
            stats.matched_branch_count += 1;
            stats.theta_prime_error_sum += theta_error;
            stats.theta_prime_error_max = std::max(stats.theta_prime_error_max, theta_error);
            stats.alpha_error_sum += alpha_error;
            stats.alpha_error_max = std::max(stats.alpha_error_max, alpha_error);
            stats.arrival_time_error_sum += time_error;
            stats.arrival_time_error_max = std::max(stats.arrival_time_error_max, time_error);
        }
    }
    stats.missing_count = stats.direct_branch_count - stats.matched_branch_count;
    stats.new_count = stats.candidate_branch_count - stats.matched_branch_count;
    return stats;
}

inline double percentile(std::vector<double> values, double p) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const double index = p * static_cast<double>(values.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(index));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(index));
    if (lo == hi) return values[lo];
    const double t = index - static_cast<double>(lo);
    return values[lo] * (1.0 - t) + values[hi] * t;
}

inline double mean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    double sum = 0.0;
    for (double value : values) sum += value;
    return sum / static_cast<double>(values.size());
}

inline std::string pair_name(const PairCase& pair) {
    return std::string(spaceship_cpp::planet_params::planet_name(pair.from)) + "_to_" +
           spaceship_cpp::planet_params::planet_name(pair.to);
}

}  // namespace problem1_5deg_test
