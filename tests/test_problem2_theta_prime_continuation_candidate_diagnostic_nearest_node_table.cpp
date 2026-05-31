#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kThetaSampleCount = 32;
constexpr double kGeometryDenominatorEpsilon = 1e-12;
constexpr double kBackSubstitutionAbsoluteAmbiguousTolerance = 1e-3;
constexpr double kBackSubstitutionRelativeAmbiguousTolerance = 1e-12;

struct GeometryDiagnostic {
    bool valid_geometry = false;
    std::string invalid_reason;
    double r0 = std::numeric_limits<double>::quiet_NaN();
    double r1 = std::numeric_limits<double>::quiet_NaN();
    double denominator = std::numeric_limits<double>::quiet_NaN();
    double e_prime = std::numeric_limits<double>::quiet_NaN();
    double R_prime = std::numeric_limits<double>::quiet_NaN();
    double B_encounter = std::numeric_limits<double>::quiet_NaN();
    double B_target = std::numeric_limits<double>::quiet_NaN();
    double encounter_radius_error = std::numeric_limits<double>::quiet_NaN();
    double target_radius_error = std::numeric_limits<double>::quiet_NaN();
};

struct Problem2ThetaBranch {
    bool problem1_valid = false;
    bool slingshot_valid = false;
    bool slingshot_boundary_ambiguous = false;
    bool search_residual_usable = false;

    std::string problem1_invalid_reason;
    std::string slingshot_invalid_reason;
    std::string boundary_ambiguous_reason;
    std::string search_residual_source;

    double theta_prime = 0.0;
    double theta_A = 0.0;

    double alpha = 0.0;
    int transfer_revolution = 0;
    int target_revolution = 0;

    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double problem1_residual_seconds = 0.0;

    double slingshot_residual = std::numeric_limits<double>::quiet_NaN();
    double outgoing_eccentricity = std::numeric_limits<double>::quiet_NaN();
    double outgoing_semi_latus_rectum = std::numeric_limits<double>::quiet_NaN();

    bool derivatives_available = false;
    double d_alpha_d_theta_A = std::numeric_limits<double>::quiet_NaN();

    std::string continuation_invalid_reason;
    double seed_encounter_global_angle = std::numeric_limits<double>::quiet_NaN();
};

struct SampleSummary {
    double theta_prime = 0.0;
    std::vector<Problem2ThetaBranch> branches;
    std::map<int, int> raw_count_by_k;
    bool table_fallback_used = false;
};

struct SignChangeCandidate {
    int interval_left = 0;
    bool origin_was_topology_change = false;
    double theta_left = 0.0;
    double theta_right = 0.0;
    int k = 0;
    int rank_in_k = 0;
    Problem2ThetaBranch left;
    Problem2ThetaBranch right;
};

struct ContinuationDiagnostic {
    bool valid = false;
    std::string invalid_reason;
    double delta_theta_prime = std::numeric_limits<double>::quiet_NaN();
    double alpha_seed = std::numeric_limits<double>::quiet_NaN();
    double initial_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double final_residual_seconds = std::numeric_limits<double>::quiet_NaN();
};

struct Summary {
    int initial_interval_count = 0;
    int initial_stable_interval_count = 0;
    int initial_topology_change_interval_count = 0;
    int adaptive_stable_subinterval_count = 0;
    int adaptive_topology_split_count = 0;
    int topology_transition_core_skipped_count = 0;
    int max_topology_recursion_depth_reached = 0;
    int raw_stable_interval_count = 0;
    int raw_topology_change_interval_count = 0;
    int raw_pair_count = 0;
    int boundary_ambiguous_pair_count = 0;
    int boundary_ambiguous_endpoint_count = 0;
    int boundary_ambiguous_candidate_count = 0;
    int residual_unusable_pair_count = 0;
    int residual_invalid_pair_count = 0;
    int sign_change_candidate_count = 0;
    int continuation_stable_candidate_count = 0;
    int two_sided_continuation_failed_count = 0;
    int one_sided_continuation_only_count = 0;
    int two_sided_continuation_mismatch_count = 0;
    int endpoint_derivatives_unavailable_count = 0;
    int continuation_refine_failed_count = 0;
    int table_fallback_count = 0;
    int seed_abnormal_count = 0;
    double max_seed_consistency_diff = 0.0;
};

struct BisectionResult {
    bool valid = false;
    std::string invalid_reason;
    SignChangeCandidate candidate;
    Problem2ThetaBranch root_branch;
    double theta_root = std::numeric_limits<double>::quiet_NaN();
    double last_theta = std::numeric_limits<double>::quiet_NaN();
    double last_residual = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    double final_width = std::numeric_limits<double>::quiet_NaN();
};

struct BisectionSummary {
    int sign_change_candidate_count = 0;
    int continuation_stable_candidate_count = 0;
    int bisection_attempt_count = 0;
    int bisection_success_count = 0;
    int dedup_success_count = 0;
    int bisection_failure_count = 0;
    int midpoint_continuation_not_stable_count = 0;
    int midpoint_boundary_ambiguous_count = 0;
    int midpoint_slingshot_invalid_count = 0;
    int midpoint_residual_unusable_count = 0;
    int width_converged_residual_too_large_count = 0;
    int max_iterations_residual_too_large_count = 0;
    int root_boundary_ambiguous_count = 0;
    int strict_root_count = 0;
    int relaxed_boundary_root_count = 0;
    double max_abs_slingshot_residual_at_root = 0.0;
    double max_abs_problem1_residual_seconds_at_root = 0.0;
};

struct BoundaryUsableResidualResult {
    bool valid = false;
    std::string invalid_reason;
    std::string source;
    double slingshot_residual = std::numeric_limits<double>::quiet_NaN();
    double outgoing_eccentricity = std::numeric_limits<double>::quiet_NaN();
    double outgoing_semi_latus_rectum = std::numeric_limits<double>::quiet_NaN();
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

int theta_sample_count_from_env() {
    if (const char* raw = std::getenv("PROBLEM2_THETA_SAMPLE_COUNT")) {
        if (*raw != '\0') {
            const int value = std::atoi(raw);
            if (value >= 4) {
                return value;
            }
        }
    }
    return kThetaSampleCount;
}

bool topology_adaptive_enabled_from_env() {
    if (const char* raw = std::getenv("PROBLEM2_TOPOLOGY_ADAPTIVE")) {
        return std::string(raw) == "1";
    }
    return false;
}

int topology_max_depth_from_env() {
    if (const char* raw = std::getenv("PROBLEM2_TOPOLOGY_MAX_DEPTH")) {
        if (*raw != '\0') {
            const int value = std::atoi(raw);
            if (value >= 0) {
                return value;
            }
        }
    }
    return 10;
}

double topology_epsilon_from_env() {
    if (const char* raw = std::getenv("PROBLEM2_TOPOLOGY_EPSILON")) {
        if (*raw != '\0') {
            const double value = std::atof(raw);
            if (is_finite(value) && value > 0.0) {
                return value;
            }
        }
    }
    return 1e-5;
}

long long theta_cache_key(double theta_prime) {
    return static_cast<long long>(std::llround(theta_prime / 1e-14));
}

std::string format_count_by_k(const std::map<int, int>& count_by_k) {
    std::ostringstream os;
    os << "{";
    bool first = true;
    for (const auto& [k, count] : count_by_k) {
        if (!first) {
            os << ",";
        }
        first = false;
        os << k << ":" << count;
    }
    os << "}";
    return os.str();
}

double wrapped_angle_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

GeometryDiagnostic compute_geometry_diagnostic(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
) {
    namespace problem2 = spaceship_cpp::problem2;
    GeometryDiagnostic g{};
    g.r0 = problem2::problem2_orbit_radius(R_J, e_J, phi);
    if (!is_finite(g.r0) || !(g.r0 > 0.0)) {
        g.invalid_reason = "invalid_encounter_radius";
        return g;
    }
    g.r1 = problem2::problem2_orbit_radius(R_K, e_K, alpha);
    if (!is_finite(g.r1) || !(g.r1 > 0.0)) {
        g.invalid_reason = "invalid_target_radius";
        return g;
    }
    const double cos_encounter = std::cos(phi - theta_prime);
    const double cos_target = std::cos(alpha - theta_prime);
    g.denominator = g.r0 * cos_encounter - g.r1 * cos_target;
    if (!is_finite(g.denominator) || std::abs(g.denominator) <= kGeometryDenominatorEpsilon) {
        g.invalid_reason = "geometry_denominator_too_small";
        return g;
    }
    g.e_prime = (g.r1 - g.r0) / g.denominator;
    if (!is_finite(g.e_prime)) {
        g.invalid_reason = "non_finite_outgoing_eccentricity";
        return g;
    }
    g.B_encounter = 1.0 + g.e_prime * cos_encounter;
    g.R_prime = g.r0 * g.B_encounter;
    if (!is_finite(g.R_prime) || !(g.R_prime > 0.0)) {
        g.invalid_reason = "non_positive_outgoing_semi_latus_rectum";
        return g;
    }
    if (!is_finite(g.B_encounter) || !(g.B_encounter > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }
    g.B_target = 1.0 + g.e_prime * cos_target;
    if (!is_finite(g.B_target) || !(g.B_target > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }
    const double encounter_reconstructed_radius = g.R_prime / g.B_encounter;
    g.encounter_radius_error = encounter_reconstructed_radius - g.r0;
    if (!is_finite(encounter_reconstructed_radius) || std::abs(g.encounter_radius_error) > 1e-10) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }
    const double target_reconstructed_radius = g.R_prime / g.B_target;
    g.target_radius_error = target_reconstructed_radius - g.r1;
    if (!is_finite(target_reconstructed_radius) || std::abs(g.target_radius_error) > 1e-10) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }
    g.valid_geometry = true;
    return g;
}

bool is_boundary_ambiguous(const std::string& reason, const GeometryDiagnostic& g) {
    if (reason != "outgoing_orbit_does_not_pass_target" &&
        reason != "outgoing_orbit_does_not_pass_encounter") {
        return false;
    }
    if (!is_finite(g.denominator) || !(std::abs(g.denominator) > kGeometryDenominatorEpsilon) ||
        !is_finite(g.R_prime) || !(g.R_prime > 0.0) ||
        !is_finite(g.B_encounter) || !(g.B_encounter > 0.0) ||
        !is_finite(g.B_target) || !(g.B_target > 0.0)) {
        return false;
    }
    const double encounter_abs = std::abs(g.encounter_radius_error);
    const double target_abs = std::abs(g.target_radius_error);
    const bool encounter_roundoff =
        (is_finite(encounter_abs) && encounter_abs <= kBackSubstitutionAbsoluteAmbiguousTolerance) ||
        (is_finite(encounter_abs) && is_finite(g.r0) && g.r0 > 0.0 &&
         encounter_abs / g.r0 <= kBackSubstitutionRelativeAmbiguousTolerance);
    const bool target_roundoff =
        (is_finite(target_abs) && target_abs <= kBackSubstitutionAbsoluteAmbiguousTolerance) ||
        (is_finite(target_abs) && is_finite(g.r1) && g.r1 > 0.0 &&
         target_abs / g.r1 <= kBackSubstitutionRelativeAmbiguousTolerance);
    return encounter_roundoff && target_roundoff;
}

BoundaryUsableResidualResult evaluate_problem2_slingshot_residual_relaxed_for_roundoff_test_only(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    const std::string& strict_invalid_reason
) {
    namespace problem2 = spaceship_cpp::problem2;
    BoundaryUsableResidualResult result{};
    const GeometryDiagnostic geometry = compute_geometry_diagnostic(R_J, e_J, R_K, e_K, phi, alpha, theta_prime);
    if (!is_boundary_ambiguous(strict_invalid_reason, geometry)) {
        result.invalid_reason = "not_roundoff_boundary_ambiguous";
        return result;
    }
    const auto residual = problem2::evaluate_problem2_slingshot_residual(
        phi,
        e_J,
        incoming_e,
        incoming_theta,
        geometry.e_prime,
        theta_prime);
    if (!residual.valid) {
        result.invalid_reason = residual.invalid_reason.empty() ? "relaxed_invariant_invalid" : residual.invalid_reason;
        return result;
    }
    result.valid = true;
    result.source = "boundary_ambiguous_roundoff";
    result.slingshot_residual = residual.residual;
    result.outgoing_eccentricity = geometry.e_prime;
    result.outgoing_semi_latus_rectum = geometry.R_prime;
    return result;
}

Problem2ThetaBranch convert_branch(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double phi,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    double theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem2 = spaceship_cpp::problem2;
    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);

    Problem2ThetaBranch out{};
    out.theta_prime = theta_prime;
    out.theta_A = theta_A;
    out.alpha = branch.target_arrival_true_anomaly;
    out.transfer_revolution = branch.transfer_revolution;
    out.target_revolution = branch.target_revolution;
    out.time_of_flight_seconds = branch.time_of_flight_seconds;
    out.target_time_seconds = branch.target_time_seconds;
    out.problem1_residual_seconds = branch.residual_seconds;
    out.derivatives_available = branch.derivatives_available;
    out.d_alpha_d_theta_A = branch.d_encounter_global_angle_d_theta_A;
    if (!branch.valid) {
        out.problem1_invalid_reason =
            branch.invalid_reason.empty() ? "problem1_branch_invalid" : branch.invalid_reason;
        return out;
    }
    out.problem1_valid = true;
    const GeometryDiagnostic geometry = compute_geometry_diagnostic(
        departure_params.orbit.p,
        departure_params.orbit.e,
        target_params.orbit.p,
        target_params.orbit.e,
        phi,
        out.alpha,
        theta_prime);
    const auto residual = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
        departure_params.orbit.p,
        departure_params.orbit.e,
        target_params.orbit.p,
        target_params.orbit.e,
        phi,
        out.alpha,
        incoming_e,
        incoming_theta,
        theta_prime);
    if (!residual.valid) {
        out.slingshot_invalid_reason =
            residual.invalid_reason.empty() ? "slingshot_invalid" : residual.invalid_reason;
        if (is_boundary_ambiguous(out.slingshot_invalid_reason, geometry)) {
            out.slingshot_boundary_ambiguous = true;
            out.boundary_ambiguous_reason = "geometry_back_substitution_roundoff";
            const auto relaxed = evaluate_problem2_slingshot_residual_relaxed_for_roundoff_test_only(
                departure_params.orbit.p,
                departure_params.orbit.e,
                target_params.orbit.p,
                target_params.orbit.e,
                phi,
                out.alpha,
                incoming_e,
                incoming_theta,
                theta_prime,
                out.slingshot_invalid_reason);
            if (relaxed.valid) {
                out.search_residual_usable = true;
                out.search_residual_source = relaxed.source;
                out.slingshot_residual = relaxed.slingshot_residual;
                out.outgoing_eccentricity = relaxed.outgoing_eccentricity;
                out.outgoing_semi_latus_rectum = relaxed.outgoing_semi_latus_rectum;
            } else {
                out.outgoing_eccentricity = geometry.e_prime;
                out.outgoing_semi_latus_rectum = geometry.R_prime;
            }
        }
        return out;
    }
    out.slingshot_valid = true;
    out.search_residual_usable = true;
    out.search_residual_source = "strict";
    out.slingshot_residual = residual.slingshot_residual;
    out.outgoing_eccentricity = residual.outgoing_eccentricity;
    out.outgoing_semi_latus_rectum = residual.outgoing_semi_latus_rectum;
    return out;
}

SampleSummary build_table_sample(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto query = problem1::query_problem1_from_2deg_nearest_node(
        loader, departure_planet, target_planet, phi, beta, theta_A, 1, 1, options);
    SampleSummary sample{};
    sample.theta_prime = theta_prime;
    sample.table_fallback_used = query.used_direct_solve_fallback;
    sample.branches.reserve(query.branches.size());
    for (const auto& branch : query.branches) {
        sample.branches.push_back(
            convert_branch(departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime, theta_A,
                           branch));
    }
    std::sort(sample.branches.begin(), sample.branches.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.time_of_flight_seconds) <
               std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.time_of_flight_seconds);
    });
    for (const auto& branch : sample.branches) {
        if (branch.problem1_valid) {
            sample.raw_count_by_k[branch.transfer_revolution] += 1;
        }
    }
    return sample;
}

std::map<int, std::vector<Problem2ThetaBranch>> raw_branches_by_k(const SampleSummary& sample) {
    std::map<int, std::vector<Problem2ThetaBranch>> by_k;
    for (const auto& branch : sample.branches) {
        if (branch.problem1_valid) {
            by_k[branch.transfer_revolution].push_back(branch);
        }
    }
    for (auto& [k, branches] : by_k) {
        (void)k;
        std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    }
    return by_k;
}

Problem2ThetaBranch continue_problem2_theta_branch_to_theta_prime_test_only(
    const Problem2ThetaBranch& endpoint,
    double theta_prime_new,
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    ContinuationDiagnostic* diagnostic,
    Summary* summary
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    (void)loader;
    Problem2ThetaBranch out{};
    out.theta_prime = theta_prime_new;
    out.transfer_revolution = endpoint.transfer_revolution;
    out.target_revolution = endpoint.target_revolution;
    if (!endpoint.derivatives_available || !is_finite(endpoint.d_alpha_d_theta_A)) {
        out.continuation_invalid_reason = "endpoint_derivatives_unavailable";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = out.continuation_invalid_reason;
        }
        if (summary != nullptr) {
            summary->endpoint_derivatives_unavailable_count += 1;
        }
        return out;
    }

    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto& target_params = planet_params::get_planet_params(target_planet);
    const double theta_A_new = normalize_angle_0_2pi(departure_state.theta_global - theta_prime_new);
    const double delta_theta_prime = normalize_angle_minus_pi_pi(theta_prime_new - endpoint.theta_prime);
    const double endpoint_encounter_global_angle = normalize_angle_0_2pi(endpoint.alpha + target_params.orbit.theta_0);
    const double alpha_seed = normalize_angle_0_2pi(
        endpoint_encounter_global_angle - endpoint.d_alpha_d_theta_A * delta_theta_prime);
    if (diagnostic != nullptr) {
        diagnostic->delta_theta_prime = delta_theta_prime;
        diagnostic->alpha_seed = alpha_seed;
    }

    const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
        departure_planet,
        target_planet,
        phi,
        beta,
        theta_A_new,
        endpoint.transfer_revolution,
        endpoint.target_revolution,
        alpha_seed,
        80,
        1e-2,
        1e-12,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    if (diagnostic != nullptr) {
        diagnostic->initial_residual_seconds = refined.diagnostic.initial_residual_seconds;
        diagnostic->final_residual_seconds = refined.diagnostic.final_residual_seconds;
    }
    if (!refined.valid || !refined.branch.valid) {
        out.continuation_invalid_reason = "continuation_refine_failed";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = out.continuation_invalid_reason;
        }
        if (summary != nullptr) {
            summary->continuation_refine_failed_count += 1;
        }
        return out;
    }

    const auto differentiated = problem1::attach_problem1_root_derivatives_with_mode(
        departure_planet,
        target_planet,
        phi,
        beta,
        theta_A_new,
        refined.branch,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    out = convert_branch(
        departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime_new, theta_A_new,
        differentiated);
    if (!out.problem1_valid) {
        out.continuation_invalid_reason = "continuation_refine_failed";
        if (diagnostic != nullptr) {
            diagnostic->invalid_reason = out.continuation_invalid_reason;
        }
        if (summary != nullptr) {
            summary->continuation_refine_failed_count += 1;
        }
        return out;
    }
    out.seed_encounter_global_angle = alpha_seed;
    if (diagnostic != nullptr) {
        diagnostic->valid = true;
    }
    return out;
}

bool residual_sign_changed(double lhs, double rhs) {
    return (lhs < 0.0 && rhs > 0.0) || (lhs > 0.0 && rhs < 0.0);
}

void print_candidate(const SignChangeCandidate& candidate);

void append_stable_interval_candidates(
    int interval_left,
    bool origin_was_topology_change,
    const SampleSummary& left,
    const SampleSummary& right,
    Summary* summary,
    std::vector<SignChangeCandidate>* candidates,
    bool print_candidates
) {
    const auto left_by_k = raw_branches_by_k(left);
    const auto right_by_k = raw_branches_by_k(right);
    for (const auto& [k, left_group] : left_by_k) {
        const auto right_it = right_by_k.find(k);
        if (right_it == right_by_k.end()) {
            continue;
        }
        const auto& right_group = right_it->second;
        const int pair_count = std::min<int>(left_group.size(), right_group.size());
        for (int rank = 0; rank < pair_count; ++rank) {
            const auto& a = left_group[static_cast<std::size_t>(rank)];
            const auto& b = right_group[static_cast<std::size_t>(rank)];
            summary->raw_pair_count += 1;
            if (a.slingshot_boundary_ambiguous || b.slingshot_boundary_ambiguous) {
                summary->boundary_ambiguous_pair_count += 1;
                summary->boundary_ambiguous_endpoint_count += (a.slingshot_boundary_ambiguous ? 1 : 0);
                summary->boundary_ambiguous_endpoint_count += (b.slingshot_boundary_ambiguous ? 1 : 0);
            }
            if (!a.search_residual_usable || !b.search_residual_usable) {
                summary->residual_unusable_pair_count += 1;
                if (!a.slingshot_valid || !b.slingshot_valid) {
                    summary->residual_invalid_pair_count += 1;
                }
                continue;
            }
            if (!residual_sign_changed(a.slingshot_residual, b.slingshot_residual)) {
                continue;
            }
            SignChangeCandidate candidate{};
            candidate.interval_left = interval_left;
            candidate.origin_was_topology_change = origin_was_topology_change;
            candidate.theta_left = left.theta_prime;
            candidate.theta_right = right.theta_prime;
            candidate.k = k;
            candidate.rank_in_k = rank;
            candidate.left = a;
            candidate.right = b;
            candidates->push_back(candidate);
            if (a.slingshot_boundary_ambiguous || b.slingshot_boundary_ambiguous) {
                summary->boundary_ambiguous_candidate_count += 1;
            }
            if (print_candidates) {
                print_candidate(candidate);
            }
        }
    }
}

void print_candidate(const SignChangeCandidate& candidate) {
    std::cout << "Problem2ContinuationCandidate"
              << " interval=" << candidate.interval_left << "-" << (candidate.interval_left + 1)
              << " theta_left=" << candidate.theta_left
              << " theta_right=" << candidate.theta_right
              << " k=" << candidate.k
              << " left_q=" << candidate.left.target_revolution
              << " right_q=" << candidate.right.target_revolution
              << " rank_in_k=" << candidate.rank_in_k
              << " left_alpha=" << candidate.left.alpha
              << " right_alpha=" << candidate.right.alpha
              << " left_time=" << candidate.left.time_of_flight_seconds
              << " right_time=" << candidate.right.time_of_flight_seconds
              << " left_residual=" << candidate.left.slingshot_residual
              << " right_residual=" << candidate.right.slingshot_residual
              << " left_e_prime=" << candidate.left.outgoing_eccentricity
              << " right_e_prime=" << candidate.right.outgoing_eccentricity
              << " left_search_residual_source=" << candidate.left.search_residual_source
              << " right_search_residual_source=" << candidate.right.search_residual_source
              << '\n';
}

bool same_continued_branch(
    const Problem2ThetaBranch& lhs,
    const Problem2ThetaBranch& rhs,
    double* alpha_diff,
    double* time_diff,
    double* e_prime_diff,
    double* residual_diff
) {
    *alpha_diff = wrapped_angle_distance(lhs.alpha, rhs.alpha);
    *time_diff = std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds);
    *e_prime_diff = std::abs(lhs.outgoing_eccentricity - rhs.outgoing_eccentricity);
    *residual_diff = std::abs(lhs.slingshot_residual - rhs.slingshot_residual);
    return lhs.transfer_revolution == rhs.transfer_revolution &&
           *alpha_diff <= 1e-5 &&
           *time_diff <= 1e3 &&
           *e_prime_diff <= 1e-4 &&
           *residual_diff <= 1e-4 &&
           lhs.search_residual_usable &&
           rhs.search_residual_usable &&
           is_finite(lhs.slingshot_residual) &&
           is_finite(rhs.slingshot_residual) &&
           is_finite(lhs.outgoing_eccentricity) &&
           is_finite(rhs.outgoing_eccentricity);
}

std::vector<SignChangeCandidate> collect_sign_change_candidates(
    const std::vector<SampleSummary>& samples,
    Summary* summary,
    bool print_candidates
) {
    std::vector<SignChangeCandidate> candidates;
    for (int i = 0; i + 1 < kThetaSampleCount; ++i) {
        const auto& left = samples[static_cast<std::size_t>(i)];
        const auto& right = samples[static_cast<std::size_t>(i + 1)];
        const bool stable = left.raw_count_by_k == right.raw_count_by_k;
        if (stable) {
            summary->raw_stable_interval_count += 1;
            append_stable_interval_candidates(i, false, left, right, summary, &candidates, print_candidates);
        } else {
            summary->raw_topology_change_interval_count += 1;
        }
    }
    summary->sign_change_candidate_count = static_cast<int>(candidates.size());
    return candidates;
}

std::vector<SignChangeCandidate> collect_sign_change_candidates_for_sample_count(
    const std::vector<SampleSummary>& samples,
    int theta_sample_count,
    Summary* summary,
    bool print_candidates
) {
    std::vector<SignChangeCandidate> candidates;
    for (int i = 0; i + 1 < theta_sample_count; ++i) {
        const auto& left = samples[static_cast<std::size_t>(i)];
        const auto& right = samples[static_cast<std::size_t>(i + 1)];
        const bool stable = left.raw_count_by_k == right.raw_count_by_k;
        if (stable) {
            summary->raw_stable_interval_count += 1;
            append_stable_interval_candidates(i, false, left, right, summary, &candidates, print_candidates);
        } else {
            summary->raw_topology_change_interval_count += 1;
        }
    }
    summary->sign_change_candidate_count = static_cast<int>(candidates.size());
    return candidates;
}

SampleSummary evaluate_cached_or_build(
    double theta_prime,
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    std::map<long long, SampleSummary>* cache,
    Summary* summary
) {
    const long long key = theta_cache_key(theta_prime);
    const auto found = cache->find(key);
    if (found != cache->end()) {
        return found->second;
    }
    auto sample = build_table_sample(
        loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
        theta_prime);
    if (sample.table_fallback_used) {
        summary->table_fallback_count += 1;
    }
    (*cache)[key] = sample;
    return sample;
}

void collect_sign_change_candidates_topology_adaptive(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    int origin_interval_left,
    double theta_left,
    double theta_right,
    const SampleSummary& left_sample,
    const SampleSummary& right_sample,
    int depth,
    int max_depth,
    double topology_epsilon,
    std::map<long long, SampleSummary>* cache,
    Summary* summary,
    std::vector<SignChangeCandidate>* candidates
) {
    summary->max_topology_recursion_depth_reached =
        std::max(summary->max_topology_recursion_depth_reached, depth);
    if (left_sample.raw_count_by_k == right_sample.raw_count_by_k) {
        summary->adaptive_stable_subinterval_count += 1;
        append_stable_interval_candidates(
            origin_interval_left, true, left_sample, right_sample, summary, candidates, false);
        return;
    }
    if (depth >= max_depth || std::abs(theta_right - theta_left) <= topology_epsilon) {
        summary->topology_transition_core_skipped_count += 1;
        std::cout << "Problem2TopologyTransitionCoreSkipped"
                  << " theta_left=" << theta_left
                  << " theta_right=" << theta_right
                  << " width=" << std::abs(theta_right - theta_left)
                  << " depth=" << depth
                  << " left_raw_count_by_k=" << format_count_by_k(left_sample.raw_count_by_k)
                  << " right_raw_count_by_k=" << format_count_by_k(right_sample.raw_count_by_k)
                  << " reason=topology_core_too_small_or_depth_limit\n";
        return;
    }

    summary->adaptive_topology_split_count += 1;
    const double theta_mid = 0.5 * (theta_left + theta_right);
    const SampleSummary mid_sample = evaluate_cached_or_build(
        theta_mid, loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
        incoming_theta, cache, summary);
    collect_sign_change_candidates_topology_adaptive(
        loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
        origin_interval_left, theta_left, theta_mid, left_sample, mid_sample, depth + 1, max_depth, topology_epsilon,
        cache, summary, candidates);
    collect_sign_change_candidates_topology_adaptive(
        loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
        origin_interval_left, theta_mid, theta_right, mid_sample, right_sample, depth + 1, max_depth, topology_epsilon,
        cache, summary, candidates);
}

bool midpoint_continuation_stable(
    const SignChangeCandidate& candidate,
    double theta_mid,
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    Summary* summary,
    Problem2ThetaBranch* selected_mid,
    std::string* invalid_reason
) {
    ContinuationDiagnostic left_diag{};
    ContinuationDiagnostic right_diag{};
    const auto left_mid = continue_problem2_theta_branch_to_theta_prime_test_only(
        candidate.left,
        theta_mid,
        loader,
        departure_planet,
        target_planet,
        encounter_time,
        phi,
        beta,
        incoming_e,
        incoming_theta,
        &left_diag,
        summary);
    const auto right_mid = continue_problem2_theta_branch_to_theta_prime_test_only(
        candidate.right,
        theta_mid,
        loader,
        departure_planet,
        target_planet,
        encounter_time,
        phi,
        beta,
        incoming_e,
        incoming_theta,
        &right_diag,
        summary);
    if (!left_mid.problem1_valid || !right_mid.problem1_valid) {
        *invalid_reason = "midpoint_continuation_not_stable";
        return false;
    }
    if (!left_mid.search_residual_usable || !right_mid.search_residual_usable) {
        *invalid_reason = "midpoint_residual_unusable";
        return false;
    }
    double alpha_diff = 0.0;
    double time_diff = 0.0;
    double e_prime_diff = 0.0;
    double residual_diff = 0.0;
    if (!same_continued_branch(left_mid, right_mid, &alpha_diff, &time_diff, &e_prime_diff, &residual_diff)) {
        *invalid_reason = "midpoint_continuation_not_stable";
        return false;
    }
    *selected_mid = std::abs(left_mid.problem1_residual_seconds) <= std::abs(right_mid.problem1_residual_seconds)
        ? left_mid
        : right_mid;
    return true;
}

BisectionResult bisect_candidate(
    const SignChangeCandidate& candidate,
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    Summary* continuation_summary
) {
    constexpr double kResidualTolerance = 1e-8;
    constexpr double kThetaTolerance = 1e-10;
    constexpr int kMaxIterations = 80;

    BisectionResult result{};
    result.candidate = candidate;
    double left_theta = candidate.theta_left;
    double right_theta = candidate.theta_right;
    Problem2ThetaBranch left_branch = candidate.left;
    Problem2ThetaBranch right_branch = candidate.right;
    for (int iteration = 0; iteration < kMaxIterations; ++iteration) {
        const double theta_mid = 0.5 * (left_theta + right_theta);
        Problem2ThetaBranch mid_branch{};
        std::string invalid_reason;
        result.iterations = iteration + 1;
        result.last_theta = theta_mid;
        result.final_width = right_theta - left_theta;
        SignChangeCandidate current_candidate = candidate;
        current_candidate.theta_left = left_theta;
        current_candidate.theta_right = right_theta;
        current_candidate.left = left_branch;
        current_candidate.right = right_branch;
        if (!midpoint_continuation_stable(
                current_candidate,
                theta_mid,
                loader,
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                continuation_summary,
                &mid_branch,
                &invalid_reason)) {
            result.invalid_reason = invalid_reason;
            return result;
        }
        result.last_residual = mid_branch.slingshot_residual;
        if (std::abs(mid_branch.slingshot_residual) <= kResidualTolerance) {
            result.valid = true;
            result.theta_root = theta_mid;
            result.root_branch = mid_branch;
            return result;
        }
        if ((right_theta - left_theta) <= kThetaTolerance) {
            result.invalid_reason = "width_converged_residual_too_large";
            return result;
        }
        if (residual_sign_changed(left_branch.slingshot_residual, mid_branch.slingshot_residual)) {
            right_theta = theta_mid;
            right_branch = mid_branch;
        } else {
            left_theta = theta_mid;
            left_branch = mid_branch;
        }
    }
    result.invalid_reason = "max_iterations_residual_too_large";
    return result;
}

bool duplicate_root(const BisectionResult& lhs, const BisectionResult& rhs) {
    return lhs.root_branch.transfer_revolution == rhs.root_branch.transfer_revolution &&
           std::abs(lhs.theta_root - rhs.theta_root) <= 1e-6 &&
           wrapped_angle_distance(lhs.root_branch.alpha, rhs.root_branch.alpha) <= 1e-6 &&
           std::abs(lhs.root_branch.time_of_flight_seconds - rhs.root_branch.time_of_flight_seconds) <= 1e3;
}

}  // namespace

#if !defined(PROBLEM2_THETA_PRIME_BISECTION_DIAGNOSTIC) && \
    !defined(PROBLEM2_THETA_PRIME_ROOT_VALIDATION_DIAGNOSTIC)
int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_continuation_candidate_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    try {
        const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
        problem1::Problem1NearestNodeQueryOptions options{};
        options.residual_tolerance_seconds = 1e-2;
        options.max_newton_iterations = 80;
        options.fallback_direct_solve = true;

        const auto departure_planet = planet_params::PlanetId::Earth;
        const auto target_planet = planet_params::PlanetId::Mars;
        const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
        const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
        const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
        const double phi = departure_state.varphi;
        const double beta = target_state.varphi;
        const double incoming_e = 0.3;
        const double incoming_theta = 0.4;

        std::vector<SampleSummary> samples;
        samples.reserve(kThetaSampleCount);
        Summary summary{};
        for (int i = 0; i < kThetaSampleCount; ++i) {
            const double theta_prime = kTwoPi * static_cast<double>(i) / static_cast<double>(kThetaSampleCount);
            samples.push_back(build_table_sample(
                loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                incoming_theta, theta_prime));
            if (samples.back().table_fallback_used) {
                summary.table_fallback_count += 1;
            }
        }

        std::vector<SignChangeCandidate> candidates;
        for (int i = 0; i + 1 < kThetaSampleCount; ++i) {
            const auto& left = samples[static_cast<std::size_t>(i)];
            const auto& right = samples[static_cast<std::size_t>(i + 1)];
            const bool stable = left.raw_count_by_k == right.raw_count_by_k;
            if (stable) {
                summary.raw_stable_interval_count += 1;
            } else {
                summary.raw_topology_change_interval_count += 1;
                continue;
            }

            const auto left_by_k = raw_branches_by_k(left);
            const auto right_by_k = raw_branches_by_k(right);
            for (const auto& [k, left_group] : left_by_k) {
                const auto right_it = right_by_k.find(k);
                if (right_it == right_by_k.end()) {
                    continue;
                }
                const auto& right_group = right_it->second;
                const int pair_count = std::min<int>(left_group.size(), right_group.size());
                for (int rank = 0; rank < pair_count; ++rank) {
                    const auto& a = left_group[static_cast<std::size_t>(rank)];
                    const auto& b = right_group[static_cast<std::size_t>(rank)];
                    summary.raw_pair_count += 1;
                    if (a.slingshot_boundary_ambiguous || b.slingshot_boundary_ambiguous) {
                        summary.boundary_ambiguous_pair_count += 1;
                        continue;
                    }
                    if (!a.slingshot_valid || !b.slingshot_valid) {
                        summary.residual_invalid_pair_count += 1;
                        continue;
                    }
                    if (!residual_sign_changed(a.slingshot_residual, b.slingshot_residual)) {
                        continue;
                    }
                    SignChangeCandidate candidate{};
                    candidate.interval_left = i;
                    candidate.theta_left = left.theta_prime;
                    candidate.theta_right = right.theta_prime;
                    candidate.k = k;
                    candidate.rank_in_k = rank;
                    candidate.left = a;
                    candidate.right = b;
                    candidates.push_back(candidate);
                    print_candidate(candidate);
                }
            }
        }
        summary.sign_change_candidate_count = static_cast<int>(candidates.size());

        const auto& target_params = planet_params::get_planet_params(target_planet);
        for (const auto& candidate : candidates) {
            const double theta_mid = 0.5 * (candidate.theta_left + candidate.theta_right);
            ContinuationDiagnostic left_diag{};
            ContinuationDiagnostic right_diag{};
            const auto left_mid = continue_problem2_theta_branch_to_theta_prime_test_only(
                candidate.left,
                theta_mid,
                loader,
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                &left_diag,
                &summary);
            const auto right_mid = continue_problem2_theta_branch_to_theta_prime_test_only(
                candidate.right,
                theta_mid,
                loader,
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                &right_diag,
                &summary);
            const bool left_valid = left_mid.problem1_valid;
            const bool right_valid = right_mid.problem1_valid;
            bool same_branch = false;
            double alpha_diff = std::numeric_limits<double>::quiet_NaN();
            double time_diff = std::numeric_limits<double>::quiet_NaN();
            double e_prime_diff = std::numeric_limits<double>::quiet_NaN();
            double residual_diff = std::numeric_limits<double>::quiet_NaN();
            if (left_valid && right_valid) {
                alpha_diff = wrapped_angle_distance(left_mid.alpha, right_mid.alpha);
                time_diff = std::abs(left_mid.time_of_flight_seconds - right_mid.time_of_flight_seconds);
                e_prime_diff = std::abs(left_mid.outgoing_eccentricity - right_mid.outgoing_eccentricity);
                residual_diff = std::abs(left_mid.slingshot_residual - right_mid.slingshot_residual);
                same_branch =
                    left_mid.transfer_revolution == right_mid.transfer_revolution &&
                    alpha_diff <= 1e-5 &&
                    time_diff <= 1e3 &&
                    e_prime_diff <= 1e-4 &&
                    residual_diff <= 1e-4 &&
                    left_mid.slingshot_valid &&
                    right_mid.slingshot_valid;
                if (same_branch) {
                    summary.continuation_stable_candidate_count += 1;
                } else {
                    summary.two_sided_continuation_mismatch_count += 1;
                }
            } else if (left_valid || right_valid) {
                summary.one_sided_continuation_only_count += 1;
            } else {
                summary.two_sided_continuation_failed_count += 1;
            }
            const double left_seed_consistency = left_valid
                ? wrapped_angle_distance(left_mid.seed_encounter_global_angle,
                                         normalize_angle_0_2pi(left_mid.alpha + target_params.orbit.theta_0))
                : std::numeric_limits<double>::quiet_NaN();
            const double right_seed_consistency = right_valid
                ? wrapped_angle_distance(right_mid.seed_encounter_global_angle,
                                         normalize_angle_0_2pi(right_mid.alpha + target_params.orbit.theta_0))
                : std::numeric_limits<double>::quiet_NaN();
            if (is_finite(left_seed_consistency)) {
                summary.max_seed_consistency_diff = std::max(summary.max_seed_consistency_diff, left_seed_consistency);
                if (left_seed_consistency > 0.2) {
                    summary.seed_abnormal_count += 1;
                }
            }
            if (is_finite(right_seed_consistency)) {
                summary.max_seed_consistency_diff =
                    std::max(summary.max_seed_consistency_diff, right_seed_consistency);
                if (right_seed_consistency > 0.2) {
                    summary.seed_abnormal_count += 1;
                }
            }
            std::cout << "Problem2ContinuationMidpointDiagnostic"
                      << " interval=" << candidate.interval_left << "-" << (candidate.interval_left + 1)
                      << " k=" << candidate.k
                      << " rank_in_k=" << candidate.rank_in_k
                      << " theta_mid=" << theta_mid
                      << " left_mid_valid=" << (left_valid ? 1 : 0)
                      << " left_mid_invalid_reason=" << left_mid.continuation_invalid_reason
                      << " right_mid_valid=" << (right_valid ? 1 : 0)
                      << " right_mid_invalid_reason=" << right_mid.continuation_invalid_reason
                      << " same_branch=" << (same_branch ? 1 : 0)
                      << " alpha_diff=" << alpha_diff
                      << " time_diff=" << time_diff
                      << " e_prime_diff=" << e_prime_diff
                      << " residual_diff=" << residual_diff
                      << " left_mid_slingshot_valid=" << (left_mid.slingshot_valid ? 1 : 0)
                      << " right_mid_slingshot_valid=" << (right_mid.slingshot_valid ? 1 : 0)
                      << " left_mid_boundary_ambiguous=" << (left_mid.slingshot_boundary_ambiguous ? 1 : 0)
                      << " right_mid_boundary_ambiguous=" << (right_mid.slingshot_boundary_ambiguous ? 1 : 0)
                      << " left_mid_search_residual_usable=" << (left_mid.search_residual_usable ? 1 : 0)
                      << " right_mid_search_residual_usable=" << (right_mid.search_residual_usable ? 1 : 0)
                      << " left_mid_search_residual_source=" << left_mid.search_residual_source
                      << " right_mid_search_residual_source=" << right_mid.search_residual_source
                      << " left_initial_residual_seconds=" << left_diag.initial_residual_seconds
                      << " right_initial_residual_seconds=" << right_diag.initial_residual_seconds
                      << " left_alpha_seed=" << left_diag.alpha_seed
                      << " right_alpha_seed=" << right_diag.alpha_seed
                      << '\n';
        }

        std::cout << "Problem2ContinuationCandidateDiagnosticSummary\n";
        std::cout << "theta_sample_count=" << kThetaSampleCount << '\n';
        std::cout << "raw_stable_interval_count=" << summary.raw_stable_interval_count << '\n';
        std::cout << "raw_topology_change_interval_count=" << summary.raw_topology_change_interval_count << '\n';
        std::cout << "raw_pair_count=" << summary.raw_pair_count << '\n';
        std::cout << "boundary_ambiguous_pair_count=" << summary.boundary_ambiguous_pair_count << '\n';
        std::cout << "boundary_ambiguous_endpoint_count=" << summary.boundary_ambiguous_endpoint_count << '\n';
        std::cout << "boundary_ambiguous_candidate_count=" << summary.boundary_ambiguous_candidate_count << '\n';
        std::cout << "residual_unusable_pair_count=" << summary.residual_unusable_pair_count << '\n';
        std::cout << "residual_invalid_pair_count=" << summary.residual_invalid_pair_count << '\n';
        std::cout << "sign_change_candidate_count=" << summary.sign_change_candidate_count << '\n';
        std::cout << "continuation_stable_candidate_count=" << summary.continuation_stable_candidate_count << '\n';
        std::cout << "two_sided_continuation_failed_count=" << summary.two_sided_continuation_failed_count << '\n';
        std::cout << "one_sided_continuation_only_count=" << summary.one_sided_continuation_only_count << '\n';
        std::cout << "two_sided_continuation_mismatch_count=" << summary.two_sided_continuation_mismatch_count << '\n';
        std::cout << "endpoint_derivatives_unavailable_count=" << summary.endpoint_derivatives_unavailable_count
                  << '\n';
        std::cout << "continuation_refine_failed_count=" << summary.continuation_refine_failed_count << '\n';
        std::cout << "table_fallback_count=" << summary.table_fallback_count << '\n';
        std::cout << "seed_abnormal_count=" << summary.seed_abnormal_count << '\n';
        std::cout << "max_seed_consistency_diff=" << summary.max_seed_consistency_diff << '\n';
        std::cout << "diagnostic_ok=1\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cout << "problem2_theta_prime_continuation_candidate_diagnostic_error=" << ex.what() << '\n';
        return 1;
    }
}
#elif defined(PROBLEM2_THETA_PRIME_BISECTION_DIAGNOSTIC)
int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_bisection_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    try {
        const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
        problem1::Problem1NearestNodeQueryOptions options{};
        options.residual_tolerance_seconds = 1e-2;
        options.max_newton_iterations = 80;
        options.fallback_direct_solve = true;

        const auto departure_planet = planet_params::PlanetId::Earth;
        const auto target_planet = planet_params::PlanetId::Mars;
        const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
        const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
        const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
        const double phi = departure_state.varphi;
        const double beta = target_state.varphi;
        const double incoming_e = 0.3;
        const double incoming_theta = 0.4;

        std::vector<SampleSummary> samples;
        samples.reserve(kThetaSampleCount);
        Summary candidate_summary{};
        for (int i = 0; i < kThetaSampleCount; ++i) {
            const double theta_prime = kTwoPi * static_cast<double>(i) / static_cast<double>(kThetaSampleCount);
            samples.push_back(build_table_sample(
                loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                incoming_theta, theta_prime));
            if (samples.back().table_fallback_used) {
                candidate_summary.table_fallback_count += 1;
            }
        }

        auto candidates = collect_sign_change_candidates(samples, &candidate_summary, false);
        std::vector<SignChangeCandidate> stable_candidates;
        Summary midpoint_probe_summary{};
        for (const auto& candidate : candidates) {
            const double theta_mid = 0.5 * (candidate.theta_left + candidate.theta_right);
            Problem2ThetaBranch mid_branch{};
            std::string invalid_reason;
            if (midpoint_continuation_stable(
                    candidate,
                    theta_mid,
                    loader,
                    departure_planet,
                    target_planet,
                    encounter_time,
                    phi,
                    beta,
                    incoming_e,
                    incoming_theta,
                    &midpoint_probe_summary,
                    &mid_branch,
                    &invalid_reason)) {
                stable_candidates.push_back(candidate);
            }
        }

        BisectionSummary bisection_summary{};
        bisection_summary.sign_change_candidate_count = static_cast<int>(candidates.size());
        bisection_summary.continuation_stable_candidate_count = static_cast<int>(stable_candidates.size());

        std::vector<BisectionResult> successes;
        for (const auto& candidate : stable_candidates) {
            bisection_summary.bisection_attempt_count += 1;
            Summary bisection_continuation_summary{};
            const BisectionResult result = bisect_candidate(
                candidate,
                loader,
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                &bisection_continuation_summary);
            if (result.valid) {
                bisection_summary.bisection_success_count += 1;
                bisection_summary.max_abs_slingshot_residual_at_root = std::max(
                    bisection_summary.max_abs_slingshot_residual_at_root,
                    std::abs(result.root_branch.slingshot_residual));
                bisection_summary.max_abs_problem1_residual_seconds_at_root = std::max(
                    bisection_summary.max_abs_problem1_residual_seconds_at_root,
                    std::abs(result.root_branch.problem1_residual_seconds));
                if (result.root_branch.slingshot_boundary_ambiguous) {
                    bisection_summary.root_boundary_ambiguous_count += 1;
                }
                if (result.root_branch.search_residual_source == "strict") {
                    bisection_summary.strict_root_count += 1;
                } else if (result.root_branch.search_residual_source == "boundary_ambiguous_roundoff") {
                    bisection_summary.relaxed_boundary_root_count += 1;
                }
                successes.push_back(result);
                std::cout << "Problem2ThetaPrimeBisectionSolution\n";
                std::cout << "valid=1\n";
                std::cout << "interval=" << result.candidate.interval_left << "-"
                          << (result.candidate.interval_left + 1) << '\n';
                std::cout << "k=" << result.candidate.k << '\n';
                std::cout << "rank_in_k=" << result.candidate.rank_in_k << '\n';
                std::cout << "theta_root=" << result.theta_root << '\n';
                std::cout << "alpha_root=" << result.root_branch.alpha << '\n';
                std::cout << "target_revolution=" << result.root_branch.target_revolution << '\n';
                std::cout << "time_of_flight_seconds=" << result.root_branch.time_of_flight_seconds << '\n';
                std::cout << "outgoing_eccentricity=" << result.root_branch.outgoing_eccentricity << '\n';
                std::cout << "outgoing_semi_latus_rectum=" << result.root_branch.outgoing_semi_latus_rectum << '\n';
                std::cout << "slingshot_residual=" << result.root_branch.slingshot_residual << '\n';
                std::cout << "problem1_residual_seconds=" << result.root_branch.problem1_residual_seconds << '\n';
                std::cout << "search_residual_source=" << result.root_branch.search_residual_source << '\n';
                std::cout << "root_boundary_ambiguous=" << (result.root_branch.slingshot_boundary_ambiguous ? 1 : 0)
                          << '\n';
                std::cout << "iterations=" << result.iterations << '\n';
                std::cout << "final_width=" << result.final_width << '\n';
            } else {
                bisection_summary.bisection_failure_count += 1;
                if (result.invalid_reason == "midpoint_continuation_not_stable") {
                    bisection_summary.midpoint_continuation_not_stable_count += 1;
                } else if (result.invalid_reason == "midpoint_boundary_ambiguous") {
                    bisection_summary.midpoint_boundary_ambiguous_count += 1;
                } else if (result.invalid_reason == "midpoint_slingshot_invalid") {
                    bisection_summary.midpoint_slingshot_invalid_count += 1;
                } else if (result.invalid_reason == "midpoint_residual_unusable") {
                    bisection_summary.midpoint_residual_unusable_count += 1;
                } else if (result.invalid_reason == "width_converged_residual_too_large") {
                    bisection_summary.width_converged_residual_too_large_count += 1;
                } else if (result.invalid_reason == "max_iterations_residual_too_large") {
                    bisection_summary.max_iterations_residual_too_large_count += 1;
                }
                std::cout << "Problem2ThetaPrimeBisectionFailure\n";
                std::cout << "interval=" << result.candidate.interval_left << "-"
                          << (result.candidate.interval_left + 1) << '\n';
                std::cout << "k=" << result.candidate.k << '\n';
                std::cout << "rank_in_k=" << result.candidate.rank_in_k << '\n';
                std::cout << "invalid_reason=" << result.invalid_reason << '\n';
                std::cout << "last_theta=" << result.last_theta << '\n';
                std::cout << "last_residual=" << result.last_residual << '\n';
                std::cout << "iterations=" << result.iterations << '\n';
                std::cout << "final_width=" << result.final_width << '\n';
            }
        }

        std::vector<BisectionResult> deduped;
        for (const auto& root : successes) {
            bool duplicate = false;
            for (const auto& existing : deduped) {
                if (duplicate_root(root, existing)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                deduped.push_back(root);
            }
        }
        bisection_summary.dedup_success_count = static_cast<int>(deduped.size());

        std::cout << "Problem2ThetaPrimeBisectionNearestNodeTableSummary\n";
        std::cout << "theta_sample_count=" << kThetaSampleCount << '\n';
        std::cout << "sign_change_candidate_count=" << bisection_summary.sign_change_candidate_count << '\n';
        std::cout << "boundary_ambiguous_endpoint_count="
                  << candidate_summary.boundary_ambiguous_endpoint_count << '\n';
        std::cout << "boundary_ambiguous_candidate_count="
                  << candidate_summary.boundary_ambiguous_candidate_count << '\n';
        std::cout << "residual_unusable_pair_count=" << candidate_summary.residual_unusable_pair_count << '\n';
        std::cout << "continuation_stable_candidate_count="
                  << bisection_summary.continuation_stable_candidate_count << '\n';
        std::cout << "bisection_attempt_count=" << bisection_summary.bisection_attempt_count << '\n';
        std::cout << "bisection_success_count=" << bisection_summary.bisection_success_count << '\n';
        std::cout << "dedup_success_count=" << bisection_summary.dedup_success_count << '\n';
        std::cout << "bisection_failure_count=" << bisection_summary.bisection_failure_count << '\n';
        std::cout << "midpoint_continuation_not_stable_count="
                  << bisection_summary.midpoint_continuation_not_stable_count << '\n';
        std::cout << "midpoint_boundary_ambiguous_count="
                  << bisection_summary.midpoint_boundary_ambiguous_count << '\n';
        std::cout << "midpoint_slingshot_invalid_count="
                  << bisection_summary.midpoint_slingshot_invalid_count << '\n';
        std::cout << "midpoint_residual_unusable_count="
                  << bisection_summary.midpoint_residual_unusable_count << '\n';
        std::cout << "width_converged_residual_too_large_count="
                  << bisection_summary.width_converged_residual_too_large_count << '\n';
        std::cout << "max_iterations_residual_too_large_count="
                  << bisection_summary.max_iterations_residual_too_large_count << '\n';
        std::cout << "root_boundary_ambiguous_count=" << bisection_summary.root_boundary_ambiguous_count << '\n';
        std::cout << "strict_root_count=" << bisection_summary.strict_root_count << '\n';
        std::cout << "relaxed_boundary_root_count=" << bisection_summary.relaxed_boundary_root_count << '\n';
        std::cout << "max_abs_slingshot_residual_at_root="
                  << bisection_summary.max_abs_slingshot_residual_at_root << '\n';
        std::cout << "max_abs_problem1_residual_seconds_at_root="
                  << bisection_summary.max_abs_problem1_residual_seconds_at_root << '\n';
        std::cout << "table_fallback_count=" << candidate_summary.table_fallback_count << '\n';
        const bool ok = bisection_summary.sign_change_candidate_count > 0;
        std::cout << "diagnostic_ok=" << (ok ? 1 : 0) << '\n';
        return ok ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cout << "problem2_theta_prime_bisection_nearest_node_table_error=" << ex.what() << '\n';
        return 1;
    }
}
#else
int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_bisection_root_validation_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    try {
        const int theta_sample_count = theta_sample_count_from_env();
        const bool topology_adaptive_enabled = topology_adaptive_enabled_from_env();
        const int topology_max_depth = topology_max_depth_from_env();
        const double topology_epsilon = topology_epsilon_from_env();
        const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
        problem1::Problem1NearestNodeQueryOptions options{};
        options.residual_tolerance_seconds = 1e-2;
        options.max_newton_iterations = 80;
        options.fallback_direct_solve = true;

        const auto departure_planet = planet_params::PlanetId::Earth;
        const auto target_planet = planet_params::PlanetId::Mars;
        const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
        const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
        const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
        const double phi = departure_state.varphi;
        const double beta = target_state.varphi;
        const double incoming_e = 0.3;
        const double incoming_theta = 0.4;

        std::vector<SampleSummary> samples;
        samples.reserve(static_cast<std::size_t>(theta_sample_count));
        Summary candidate_summary{};
        std::map<long long, SampleSummary> sample_cache;
        for (int i = 0; i < theta_sample_count; ++i) {
            const double theta_prime =
                kTwoPi * static_cast<double>(i) / static_cast<double>(theta_sample_count);
            samples.push_back(build_table_sample(
                loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                incoming_theta, theta_prime));
            sample_cache[theta_cache_key(theta_prime)] = samples.back();
            if (samples.back().table_fallback_used) {
                candidate_summary.table_fallback_count += 1;
            }
        }

        candidate_summary.initial_interval_count = theta_sample_count - 1;
        for (int i = 0; i + 1 < theta_sample_count; ++i) {
            if (samples[static_cast<std::size_t>(i)].raw_count_by_k ==
                samples[static_cast<std::size_t>(i + 1)].raw_count_by_k) {
                candidate_summary.initial_stable_interval_count += 1;
            } else {
                candidate_summary.initial_topology_change_interval_count += 1;
            }
        }

        std::vector<SignChangeCandidate> candidates;
        if (topology_adaptive_enabled) {
            for (int i = 0; i + 1 < theta_sample_count; ++i) {
                const auto& left = samples[static_cast<std::size_t>(i)];
                const auto& right = samples[static_cast<std::size_t>(i + 1)];
                if (left.raw_count_by_k == right.raw_count_by_k) {
                    candidate_summary.raw_stable_interval_count += 1;
                    append_stable_interval_candidates(i, false, left, right, &candidate_summary, &candidates, false);
                } else {
                    candidate_summary.raw_topology_change_interval_count += 1;
                    collect_sign_change_candidates_topology_adaptive(
                        loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                        incoming_theta, i, left.theta_prime, right.theta_prime, left, right, 0, topology_max_depth,
                        topology_epsilon, &sample_cache, &candidate_summary, &candidates);
                }
            }
            candidate_summary.sign_change_candidate_count = static_cast<int>(candidates.size());
        } else {
            candidates = collect_sign_change_candidates_for_sample_count(
                samples, theta_sample_count, &candidate_summary, false);
        }
        std::vector<SignChangeCandidate> stable_candidates;
        Summary midpoint_probe_summary{};
        for (const auto& candidate : candidates) {
            const double theta_mid = 0.5 * (candidate.theta_left + candidate.theta_right);
            Problem2ThetaBranch mid_branch{};
            std::string invalid_reason;
            if (midpoint_continuation_stable(
                    candidate,
                    theta_mid,
                    loader,
                    departure_planet,
                    target_planet,
                    encounter_time,
                    phi,
                    beta,
                    incoming_e,
                    incoming_theta,
                    &midpoint_probe_summary,
                    &mid_branch,
                    &invalid_reason)) {
                stable_candidates.push_back(candidate);
            }
        }

        BisectionSummary bisection_summary{};
        bisection_summary.sign_change_candidate_count = static_cast<int>(candidates.size());
        bisection_summary.continuation_stable_candidate_count = static_cast<int>(stable_candidates.size());

        std::vector<BisectionResult> successes;
        for (const auto& candidate : stable_candidates) {
            bisection_summary.bisection_attempt_count += 1;
            Summary bisection_continuation_summary{};
            const BisectionResult result = bisect_candidate(
                candidate,
                loader,
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                &bisection_continuation_summary);
            if (!result.valid) {
                bisection_summary.bisection_failure_count += 1;
                continue;
            }
            bisection_summary.bisection_success_count += 1;
            bisection_summary.max_abs_slingshot_residual_at_root = std::max(
                bisection_summary.max_abs_slingshot_residual_at_root,
                std::abs(result.root_branch.slingshot_residual));
            bisection_summary.max_abs_problem1_residual_seconds_at_root = std::max(
                bisection_summary.max_abs_problem1_residual_seconds_at_root,
                std::abs(result.root_branch.problem1_residual_seconds));
            if (result.root_branch.search_residual_source == "strict") {
                bisection_summary.strict_root_count += 1;
            } else if (result.root_branch.search_residual_source == "boundary_ambiguous_roundoff") {
                bisection_summary.relaxed_boundary_root_count += 1;
            }
            successes.push_back(result);
        }

        std::vector<BisectionResult> deduped;
        for (const auto& root : successes) {
            bool duplicate = false;
            for (const auto& existing : deduped) {
                if (duplicate_root(root, existing)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                deduped.push_back(root);
            }
        }
        bisection_summary.dedup_success_count = static_cast<int>(deduped.size());

        int validation_success_count = 0;
        int validation_failure_count = 0;
        double max_alpha_diff_direct = 0.0;
        double max_time_diff_direct = 0.0;
        double max_abs_direct_residual = 0.0;
        int root_index = 0;
        for (const auto& root : deduped) {
            const auto& root_branch = root.root_branch;
            std::cout << "Problem2RootValidationInput"
                      << " root_index=" << root_index
                      << " theta_root=" << root.theta_root
                      << " alpha_root=" << root_branch.alpha
                      << " k=" << root_branch.transfer_revolution
                      << " q=" << root_branch.target_revolution
                      << " time_of_flight_seconds=" << root_branch.time_of_flight_seconds
                      << " outgoing_eccentricity=" << root_branch.outgoing_eccentricity
                      << " outgoing_semi_latus_rectum=" << root_branch.outgoing_semi_latus_rectum
                      << " slingshot_residual=" << root_branch.slingshot_residual
                      << " problem1_residual_seconds=" << root_branch.problem1_residual_seconds
                      << " search_residual_source=" << root_branch.search_residual_source
                      << " root_boundary_ambiguous=" << (root_branch.slingshot_boundary_ambiguous ? 1 : 0)
                      << '\n';

            const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - root.theta_root);
            const auto direct_problem1 = problem1::solve_problem1_from_departure_anomalies(
                departure_planet, target_planet, phi, beta, theta_A, 1, 1);
            bool matched = false;
            Problem2ThetaBranch best{};
            double best_score = std::numeric_limits<double>::infinity();
            for (const auto& branch : direct_problem1) {
                if (!branch.valid || branch.transfer_revolution != root_branch.transfer_revolution) {
                    continue;
                }
                const auto converted = convert_branch(
                    departure_planet, target_planet, phi, incoming_e, incoming_theta, root.theta_root, theta_A,
                    branch);
                const double alpha_diff = wrapped_angle_distance(converted.alpha, root_branch.alpha);
                const double time_diff =
                    std::abs(converted.time_of_flight_seconds - root_branch.time_of_flight_seconds);
                const double score = alpha_diff + time_diff / 1.0e8;
                if (score < best_score) {
                    best_score = score;
                    best = converted;
                    matched = true;
                }
            }
            const double alpha_diff = matched ? wrapped_angle_distance(best.alpha, root_branch.alpha)
                                              : std::numeric_limits<double>::quiet_NaN();
            const double time_diff = matched
                ? std::abs(best.time_of_flight_seconds - root_branch.time_of_flight_seconds)
                : std::numeric_limits<double>::quiet_NaN();
            const bool root_validated = matched &&
                                        alpha_diff <= 1e-7 &&
                                        time_diff <= 1e3 &&
                                        best.search_residual_usable &&
                                        std::abs(best.slingshot_residual) <= 1e-7;
            if (root_validated) {
                validation_success_count += 1;
                max_alpha_diff_direct = std::max(max_alpha_diff_direct, alpha_diff);
                max_time_diff_direct = std::max(max_time_diff_direct, time_diff);
                max_abs_direct_residual = std::max(max_abs_direct_residual, std::abs(best.slingshot_residual));
            } else {
                validation_failure_count += 1;
            }
            std::cout << "Problem2RootDirectValidation"
                      << " root_index=" << root_index
                      << " matched=" << (matched ? 1 : 0)
                      << " matched_k=" << best.transfer_revolution
                      << " matched_q=" << best.target_revolution
                      << " alpha_diff=" << alpha_diff
                      << " time_diff_seconds=" << time_diff
                      << " problem1_residual_seconds=" << best.problem1_residual_seconds
                      << " slingshot_residual_direct=" << best.slingshot_residual
                      << " slingshot_residual_source_direct=" << best.search_residual_source
                      << " root_validated=" << (root_validated ? 1 : 0)
                      << '\n';
            root_index += 1;
        }

        std::cout << "Problem2TopologyAdaptiveSearchSummary\n";
        std::cout << "theta_sample_count=" << theta_sample_count << '\n';
        std::cout << "topology_adaptive_enabled=" << (topology_adaptive_enabled ? 1 : 0) << '\n';
        std::cout << "initial_interval_count=" << candidate_summary.initial_interval_count << '\n';
        std::cout << "initial_stable_interval_count=" << candidate_summary.initial_stable_interval_count << '\n';
        std::cout << "initial_topology_change_interval_count="
                  << candidate_summary.initial_topology_change_interval_count << '\n';
        std::cout << "adaptive_stable_subinterval_count="
                  << candidate_summary.adaptive_stable_subinterval_count << '\n';
        std::cout << "adaptive_topology_split_count=" << candidate_summary.adaptive_topology_split_count << '\n';
        std::cout << "topology_transition_core_skipped_count="
                  << candidate_summary.topology_transition_core_skipped_count << '\n';
        std::cout << "max_topology_recursion_depth_reached="
                  << candidate_summary.max_topology_recursion_depth_reached << '\n';
        std::cout << "sign_change_candidate_count=" << bisection_summary.sign_change_candidate_count << '\n';
        std::cout << "continuation_stable_candidate_count="
                  << bisection_summary.continuation_stable_candidate_count << '\n';
        std::cout << "bisection_attempt_count=" << bisection_summary.bisection_attempt_count << '\n';
        std::cout << "bisection_success_count=" << bisection_summary.bisection_success_count << '\n';
        std::cout << "dedup_success_count=" << bisection_summary.dedup_success_count << '\n';
        std::cout << "strict_root_count=" << bisection_summary.strict_root_count << '\n';
        std::cout << "relaxed_boundary_root_count=" << bisection_summary.relaxed_boundary_root_count << '\n';
        std::cout << "direct_validation_success_count=" << validation_success_count << '\n';
        std::cout << "direct_validation_failure_count=" << validation_failure_count << '\n';
        std::cout << "table_fallback_count=" << candidate_summary.table_fallback_count << '\n';
        std::cout << "max_abs_direct_slingshot_residual=" << max_abs_direct_residual << '\n';

        std::cout << "Problem2ThetaSampleScalingSummary\n";
        std::cout << "theta_sample_count=" << theta_sample_count << '\n';
        std::cout << "sign_change_candidate_count=" << bisection_summary.sign_change_candidate_count << '\n';
        std::cout << "continuation_stable_candidate_count="
                  << bisection_summary.continuation_stable_candidate_count << '\n';
        std::cout << "bisection_attempt_count=" << bisection_summary.bisection_attempt_count << '\n';
        std::cout << "bisection_success_count=" << bisection_summary.bisection_success_count << '\n';
        std::cout << "dedup_success_count=" << bisection_summary.dedup_success_count << '\n';
        std::cout << "strict_root_count=" << bisection_summary.strict_root_count << '\n';
        std::cout << "relaxed_boundary_root_count=" << bisection_summary.relaxed_boundary_root_count << '\n';
        std::cout << "max_abs_slingshot_residual_at_root="
                  << bisection_summary.max_abs_slingshot_residual_at_root << '\n';
        std::cout << "table_fallback_count=" << candidate_summary.table_fallback_count << '\n';

        std::cout << "Problem2RootValidationSummary\n";
        std::cout << "theta_sample_count=" << theta_sample_count << '\n';
        std::cout << "dedup_root_count=" << bisection_summary.dedup_success_count << '\n';
        std::cout << "direct_validation_success_count=" << validation_success_count << '\n';
        std::cout << "direct_validation_failure_count=" << validation_failure_count << '\n';
        std::cout << "max_alpha_diff_direct=" << max_alpha_diff_direct << '\n';
        std::cout << "max_time_diff_direct_seconds=" << max_time_diff_direct << '\n';
        std::cout << "max_abs_direct_slingshot_residual=" << max_abs_direct_residual << '\n';
        std::cout << "strict_root_count=" << bisection_summary.strict_root_count << '\n';
        std::cout << "relaxed_boundary_root_count=" << bisection_summary.relaxed_boundary_root_count << '\n';
        std::cout << "table_fallback_count=" << candidate_summary.table_fallback_count << '\n';
        std::cout << "diagnostic_ok=1\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cout << "problem2_theta_prime_bisection_root_validation_error=" << ex.what() << '\n';
        return 1;
    }
}
#endif
