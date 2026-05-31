#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace spaceship_cpp::problem2 {

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kGeometryDenominatorEpsilon = 1e-12;

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

int chunk_index_for_linear_index(const problem1::Problem1RootTable2DegLoader& loader, long long linear_index) {
    const auto& chunks = loader.chunks();
    for (std::size_t i = 0; i < chunks.size(); ++i) {
        if (linear_index >= chunks[i].start_node && linear_index < chunks[i].end_node) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

struct GeometryCheck {
    bool valid = false;
    std::string invalid_reason;
    double r0 = std::numeric_limits<double>::quiet_NaN();
    double r1 = std::numeric_limits<double>::quiet_NaN();
    double denominator = std::numeric_limits<double>::quiet_NaN();
    double e_prime = std::numeric_limits<double>::quiet_NaN();
    double p_prime = std::numeric_limits<double>::quiet_NaN();
    double encounter_factor = std::numeric_limits<double>::quiet_NaN();
    double target_factor = std::numeric_limits<double>::quiet_NaN();
    double encounter_radius_error = std::numeric_limits<double>::quiet_NaN();
    double target_radius_error = std::numeric_limits<double>::quiet_NaN();
};

struct ThetaBranch {
    bool problem1_valid = false;
    bool slingshot_valid = false;
    bool boundary_ambiguous = false;
    bool residual_usable = false;

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
    Problem2ResidualSource residual_source = Problem2ResidualSource::Strict;
};

struct Sample {
    double theta_prime = 0.0;
    std::vector<ThetaBranch> branches;
    std::map<int, int> raw_count_by_k;
    bool table_fallback_used = false;
};

struct Candidate {
    int interval_left = 0;
    bool origin_was_topology_change = false;
    double theta_left = 0.0;
    double theta_right = 0.0;
    int k = 0;
    int rank_in_k = 0;
    ThetaBranch left;
    ThetaBranch right;
};

struct BisectionResult {
    bool valid = false;
    std::string invalid_reason;
    Candidate candidate;
    ThetaBranch root_branch;
    double theta_root = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    double final_width = std::numeric_limits<double>::quiet_NaN();
};

enum class Problem2SampleBuildPhase {
    Initial,
    Adaptive
};

void accumulate_nearest_node_query_profile(
    const problem1::Problem1NearestNodeQueryProfile& query_profile,
    Problem2SampleBuildPhase phase,
    Problem2GravityAssistSolverProfile* profile
);

ThetaBranch continue_branch(
    const ThetaBranch& endpoint,
    double theta_prime_new,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    Problem2GravityAssistSolverProfile* profile
);

const Sample& choose_midpoint_seed_sample(
    const Sample& left_sample,
    const Sample& right_sample,
    const Problem2GravityAssistSolverOptions& options
);

bool sample_passes_midpoint_routeA_validation(const Sample& sample);

Sample build_midpoint_sample_from_endpoint_routeA(
    const Sample& seed_sample,
    double theta_mid,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    Problem2GravityAssistSolverProfile* profile
);

long long theta_cache_key(double theta_prime) {
    return static_cast<long long>(std::llround(theta_prime / 1e-14));
}

double wrapped_angle_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

bool residual_sign_changed(double lhs, double rhs) {
    return (lhs < 0.0 && rhs > 0.0) || (lhs > 0.0 && rhs < 0.0);
}

GeometryCheck compute_geometry(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
) {
    GeometryCheck g{};
    g.r0 = problem2_orbit_radius(R_J, e_J, phi);
    if (!is_finite(g.r0) || !(g.r0 > 0.0)) {
        g.invalid_reason = "invalid_encounter_radius";
        return g;
    }
    g.r1 = problem2_orbit_radius(R_K, e_K, alpha);
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

    g.encounter_factor = 1.0 + g.e_prime * cos_encounter;
    g.p_prime = g.r0 * g.encounter_factor;
    if (!is_finite(g.p_prime) || !(g.p_prime > 0.0)) {
        g.invalid_reason = "non_positive_outgoing_semi_latus_rectum";
        return g;
    }
    if (!is_finite(g.encounter_factor) || !(g.encounter_factor > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }

    g.target_factor = 1.0 + g.e_prime * cos_target;
    if (!is_finite(g.target_factor) || !(g.target_factor > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }

    const double encounter_radius = g.p_prime / g.encounter_factor;
    g.encounter_radius_error = encounter_radius - g.r0;
    if (!is_finite(encounter_radius) || std::abs(g.encounter_radius_error) > 1e-10) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }
    const double target_radius = g.p_prime / g.target_factor;
    g.target_radius_error = target_radius - g.r1;
    if (!is_finite(target_radius) || std::abs(g.target_radius_error) > 1e-10) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }

    g.valid = true;
    return g;
}

bool is_boundary_ambiguous(
    const std::string& strict_reason,
    const GeometryCheck& g,
    const Problem2GravityAssistSolverOptions& options
) {
    if (!options.allow_roundoff_boundary_relaxed_residual) {
        return false;
    }
    if (strict_reason != "outgoing_orbit_does_not_pass_target" &&
        strict_reason != "outgoing_orbit_does_not_pass_encounter") {
        return false;
    }
    if (!is_finite(g.denominator) || !(std::abs(g.denominator) > kGeometryDenominatorEpsilon) ||
        !is_finite(g.p_prime) || !(g.p_prime > 0.0) ||
        !is_finite(g.encounter_factor) || !(g.encounter_factor > 0.0) ||
        !is_finite(g.target_factor) || !(g.target_factor > 0.0)) {
        return false;
    }

    const auto roundoff_ok = [&](double abs_error, double radius) {
        return (is_finite(abs_error) && abs_error <= options.boundary_roundoff_abs_tolerance_m) ||
               (is_finite(abs_error) && is_finite(radius) && radius > 0.0 &&
                abs_error / radius <= options.boundary_roundoff_rel_tolerance);
    };
    return roundoff_ok(std::abs(g.encounter_radius_error), g.r0) &&
           roundoff_ok(std::abs(g.target_radius_error), g.r1);
}

ThetaBranch convert_branch(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double phi,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    double theta_A,
    const problem1::Problem1SolutionBranch& branch,
    const Problem2GravityAssistSolverOptions& options
) {
    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);

    ThetaBranch out{};
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
        return out;
    }
    out.problem1_valid = true;

    const GeometryCheck geometry = compute_geometry(
        departure_params.orbit.p,
        departure_params.orbit.e,
        target_params.orbit.p,
        target_params.orbit.e,
        phi,
        out.alpha,
        theta_prime);
    const auto strict = evaluate_problem2_slingshot_residual_from_theta_alpha(
        departure_params.orbit.p,
        departure_params.orbit.e,
        target_params.orbit.p,
        target_params.orbit.e,
        phi,
        out.alpha,
        incoming_e,
        incoming_theta,
        theta_prime);

    if (strict.valid) {
        out.slingshot_valid = true;
        out.residual_usable = true;
        out.residual_source = Problem2ResidualSource::Strict;
        out.slingshot_residual = strict.slingshot_residual;
        out.outgoing_eccentricity = strict.outgoing_eccentricity;
        out.outgoing_semi_latus_rectum = strict.outgoing_semi_latus_rectum;
        return out;
    }

    const std::string strict_reason = strict.invalid_reason.empty() ? "slingshot_invalid" : strict.invalid_reason;
    if (!is_boundary_ambiguous(strict_reason, geometry, options)) {
        return out;
    }

    out.boundary_ambiguous = true;
    const auto relaxed = evaluate_problem2_slingshot_residual(
        phi,
        departure_params.orbit.e,
        incoming_e,
        incoming_theta,
        geometry.e_prime,
        theta_prime);
    if (!relaxed.valid) {
        out.outgoing_eccentricity = geometry.e_prime;
        out.outgoing_semi_latus_rectum = geometry.p_prime;
        return out;
    }

    out.residual_usable = true;
    out.residual_source = Problem2ResidualSource::BoundaryAmbiguousRoundoff;
    out.slingshot_residual = relaxed.residual;
    out.outgoing_eccentricity = geometry.e_prime;
    out.outgoing_semi_latus_rectum = geometry.p_prime;
    return out;
}

Sample build_sample(
    const problem1::Problem1RootTable2DegLoader& loader,
    const problem1::Problem1NearestNodeQueryOptions& query_options,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    Problem2SampleBuildPhase phase,
    const Problem2GravityAssistSolverOptions& options,
    Problem2GravityAssistSolverProfile* profile
) {
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    if (profile != nullptr) {
        profile->total_sample_build_count += 1;
    }

    std::vector<problem1::Problem1SolutionBranch> problem1_branches;
    bool table_fallback_used = false;
    problem1::Problem1NearestNodeQueryResult query{};
    double direct_solve_ms = 0.0;
    double direct_attach_ms = 0.0;

    if (options.use_problem1_solve_for_problem2_seed) {
        const auto solve_start = Clock::now();
        problem1_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            phi,
            beta,
            theta_A,
            options.max_transfer_revolution,
            options.max_target_revolution);
        direct_solve_ms = elapsed_ms(solve_start, Clock::now());

        const auto attach_start = Clock::now();
        for (auto& branch : problem1_branches) {
            if (!branch.valid) {
                continue;
            }
            const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
                departure_planet,
                target_planet,
                phi,
                beta,
                theta_A,
                branch,
                problem1::Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
            if (attached.valid) {
                branch = attached;
            }
        }
        direct_attach_ms = elapsed_ms(attach_start, Clock::now());
        if (profile != nullptr) {
            profile->direct_problem1_seed_solve_count += 1;
            profile->direct_problem1_seed_branch_count += static_cast<int>(problem1_branches.size());
            profile->direct_problem1_seed_solve_ms += direct_solve_ms;
            profile->direct_problem1_seed_derivative_attach_ms += direct_attach_ms;
            for (const auto& branch : problem1_branches) {
                if (branch.valid && branch.derivatives_available) {
                    profile->direct_problem1_seed_derivative_attach_count += 1;
                }
            }
        }
    } else {
        if (profile != nullptr) {
            profile->problem1_nearest_node_query_count += 1;
        }
        query = problem1::query_problem1_from_2deg_nearest_node(
            loader,
            departure_planet,
            target_planet,
            phi,
            beta,
            theta_A,
            options.max_transfer_revolution,
            options.max_target_revolution,
            query_options);
        accumulate_nearest_node_query_profile(query.profile, phase, profile);
        problem1_branches = query.branches;
        table_fallback_used = query.used_direct_solve_fallback;
    }

    if (profile != nullptr && options.adaptive_interval_trace_enabled) {
        Problem2GravityAssistSolverProfile::ThetaSampleTrace trace{};
        trace.theta_prime = theta_prime;
        trace.query_nu_A = phi;
        trace.query_nu_B = beta;
        trace.query_theta_A = theta_A;
        trace.adaptive = phase == Problem2SampleBuildPhase::Adaptive;
        trace.from_cache = false;
        if (options.use_problem1_solve_for_problem2_seed) {
            trace.nearest_linear_index = -1;
            trace.nearest_nu_A_index = -1;
            trace.nearest_nu_B_index = -1;
            trace.nearest_theta_A_index = -1;
            trace.chunk_index = -1;
            trace.total_ms = direct_solve_ms + direct_attach_ms;
            trace.fallback_direct_solve_ms = direct_solve_ms;
            trace.loaded_branch_count = static_cast<int>(problem1_branches.size());
            trace.seed_branch_count = static_cast<int>(problem1_branches.size());
            trace.derivative_attach_attempt_count = static_cast<int>(problem1_branches.size());
            trace.fallback_direct_solve_count = 1;
        } else {
            trace.nearest_linear_index = query.nearest_linear_index;
            trace.nearest_nu_A_index = query.nearest_nu_A_index;
            trace.nearest_nu_B_index = query.nearest_nu_B_index;
            trace.nearest_theta_A_index = query.nearest_theta_A_index;
            trace.chunk_index = chunk_index_for_linear_index(loader, query.nearest_linear_index);
            trace.total_ms = query.profile.total_ms;
            trace.table_offset_build_ms = query.profile.table_offset_build_ms;
            trace.table_node_read_ms = query.profile.table_node_read_ms;
            trace.seed_generation_ms = query.profile.seed_generation_ms;
            trace.refine_ms = query.profile.refine_ms;
            trace.derivative_attach_ms = query.profile.derivative_attach_ms;
            trace.dedup_ms = query.profile.dedup_ms;
            trace.fallback_direct_solve_ms = query.profile.fallback_direct_solve_ms;
            trace.loaded_branch_count = query.profile.table_loaded_branch_count;
            trace.seed_branch_count = query.profile.seed_branch_count;
            trace.refine_attempt_count = query.profile.refine_attempt_count;
            trace.refine_success_count = query.profile.refine_success_count;
            trace.derivative_attach_attempt_count = query.profile.derivative_attach_attempt_count;
            trace.fallback_direct_solve_count = query.profile.fallback_direct_solve_count;
        }
        profile->theta_sample_traces.push_back(trace);
    }

    Sample sample{};
    sample.theta_prime = theta_prime;
    sample.table_fallback_used = table_fallback_used;
    sample.branches.reserve(problem1_branches.size());
    for (const auto& branch : problem1_branches) {
        sample.branches.push_back(convert_branch(
            departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime, theta_A, branch, options));
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

int valid_branch_count(const Sample& sample) {
    int count = 0;
    for (const auto& branch : sample.branches) {
        count += branch.problem1_valid ? 1 : 0;
    }
    return count;
}

std::map<int, std::vector<ThetaBranch>> raw_branches_by_k(const Sample& sample) {
    std::map<int, std::vector<ThetaBranch>> by_k;
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

void append_stable_interval_candidates(
    int interval_left,
    bool origin_was_topology_change,
    const Sample& left,
    const Sample& right,
    Problem2GravityAssistSolverSummary* summary,
    std::vector<Candidate>* candidates
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
            if (!a.residual_usable || !b.residual_usable) {
                continue;
            }
            if (!residual_sign_changed(a.slingshot_residual, b.slingshot_residual)) {
                continue;
            }
            Candidate candidate{};
            candidate.interval_left = interval_left;
            candidate.origin_was_topology_change = origin_was_topology_change;
            candidate.theta_left = left.theta_prime;
            candidate.theta_right = right.theta_prime;
            candidate.k = k;
            candidate.rank_in_k = rank;
            candidate.left = a;
            candidate.right = b;
            candidates->push_back(candidate);
        }
    }
    summary->sign_change_candidate_count = static_cast<int>(candidates->size());
}

Sample evaluate_cached_or_build(
    double theta_prime,
    const problem1::Problem1RootTable2DegLoader& loader,
    const problem1::Problem1NearestNodeQueryOptions& query_options,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    std::map<long long, Sample>* cache,
    Problem2GravityAssistSolverSummary* summary,
    Problem2GravityAssistSolverProfile* profile,
    bool* from_cache = nullptr
) {
    const long long key = theta_cache_key(theta_prime);
    const auto found = cache->find(key);
    if (found != cache->end()) {
        if (profile != nullptr) {
            profile->cached_sample_hit_count += 1;
            profile->adaptive_cache_hit_count += 1;
        }
        if (from_cache != nullptr) {
            *from_cache = true;
        }
        return found->second;
    }
    if (profile != nullptr) {
        profile->adaptive_midpoint_sample_count += 1;
        profile->adaptive_cache_miss_count += 1;
    }
    if (from_cache != nullptr) {
        *from_cache = false;
    }
    Sample sample = build_sample(
        loader, query_options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
        incoming_theta, theta_prime, Problem2SampleBuildPhase::Adaptive, options, profile);
    if (sample.table_fallback_used) {
        summary->table_fallback_count += 1;
    }
    (*cache)[key] = sample;
    return sample;
}

void collect_topology_adaptive(
    const problem1::Problem1RootTable2DegLoader& loader,
    const problem1::Problem1NearestNodeQueryOptions& query_options,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    int origin_interval_left,
    double theta_left,
    double theta_right,
    const Sample& left_sample,
    const Sample& right_sample,
    int depth,
    const Problem2GravityAssistSolverOptions& options,
    std::map<long long, Sample>* cache,
    Problem2GravityAssistSolverSummary* summary,
    Problem2GravityAssistSolverProfile* profile,
    std::vector<Candidate>* candidates
) {
    const auto trace_start = Clock::now();
    Problem2GravityAssistSolverProfile::AdaptiveIntervalTrace trace{};
    const bool trace_enabled = profile != nullptr && options.adaptive_interval_trace_enabled;
    if (profile != nullptr) {
        profile->adaptive_interval_count += 1;
        profile->adaptive_left_endpoint_query_count += 1;
        profile->adaptive_right_endpoint_query_count += 1;
    }
    if (trace_enabled) {
        trace.interval_id = profile->adaptive_interval_count;
        trace.depth = depth;
        trace.theta_left = theta_left;
        trace.theta_right = theta_right;
        trace.theta_mid = 0.5 * (theta_left + theta_right);
        trace.left_sample_built = true;
        trace.right_sample_built = true;
        trace.left_branch_count = valid_branch_count(left_sample);
        trace.right_branch_count = valid_branch_count(right_sample);
    }
    summary->max_topology_recursion_depth_reached =
        std::max(summary->max_topology_recursion_depth_reached, depth);
    if (left_sample.raw_count_by_k == right_sample.raw_count_by_k) {
        summary->adaptive_stable_subinterval_count += 1;
        if (profile != nullptr) {
            profile->adaptive_stable_interval_count += 1;
        }
        if (trace_enabled) {
            trace.classified_stable = true;
            profile->adaptive_interval_traces.push_back(trace);
            profile->adaptive_interval_trace_overhead_ms += elapsed_ms(trace_start, Clock::now());
        }
        append_stable_interval_candidates(origin_interval_left, true, left_sample, right_sample, summary, candidates);
        return;
    }
    if (depth >= options.topology_max_depth || std::abs(theta_right - theta_left) <= options.topology_epsilon) {
        summary->topology_transition_core_skipped_count += 1;
        if (profile != nullptr) {
            profile->adaptive_topology_change_interval_count += 1;
        }
        if (trace_enabled) {
            trace.classified_topology_change = true;
            trace.classified_discontinuous = true;
            trace.subdivision_reason = depth >= options.topology_max_depth ? "max_depth_reached" : "epsilon_reached";
            profile->adaptive_interval_traces.push_back(trace);
            profile->adaptive_interval_trace_overhead_ms += elapsed_ms(trace_start, Clock::now());
        }
        return;
    }

    summary->adaptive_topology_split_count += 1;
    if (profile != nullptr) {
        profile->adaptive_topology_change_interval_count += 1;
        profile->adaptive_subdivided_interval_count += 1;
        profile->adaptive_midpoint_query_count += 1;
    }
    const double theta_mid = 0.5 * (theta_left + theta_right);
    const auto sample_start = Clock::now();
    const int query_count_before = profile != nullptr ? profile->problem1_nearest_node_query_count : 0;
    const double query_ms_before = profile != nullptr ? profile->nearest_query_total_ms : 0.0;
    bool mid_from_cache = false;
    Sample mid_sample{};
    bool routeA_midpoint_used = false;
    if (options.adaptive_midpoint_routeA_from_endpoint_enabled &&
        options.adaptive_midpoint_seed_policy != 0) {
        if (profile != nullptr) {
            profile->adaptive_midpoint_routeA_attempt_count += 1;
        }
        Problem2GravityAssistSolverProfile::RouteAMidpointTrace route_trace{};
        route_trace.theta_mid = theta_mid;
        route_trace.seed_policy = options.adaptive_midpoint_seed_policy;
        const auto routeA_start = Clock::now();
        const Sample& seed_sample = choose_midpoint_seed_sample(left_sample, right_sample, options);
        route_trace.seed_endpoint_left = &seed_sample == &left_sample;
        route_trace.seed_valid_branch_count = valid_branch_count(seed_sample);
        const int continuation_count_before = profile != nullptr ? profile->continuation_refine_count : 0;
        mid_sample = build_midpoint_sample_from_endpoint_routeA(
            seed_sample,
            theta_mid,
            departure_planet,
            target_planet,
            encounter_time,
            phi,
            beta,
            incoming_e,
            incoming_theta,
            options,
            profile);
        route_trace.routeA_attempt_count_for_midpoint =
            profile != nullptr ? profile->continuation_refine_count - continuation_count_before : 0;
        route_trace.routeA_success_count_for_midpoint = valid_branch_count(mid_sample);
        if (profile != nullptr) {
            route_trace.routeA_total_ms_for_midpoint = elapsed_ms(routeA_start, Clock::now());
            profile->adaptive_midpoint_routeA_ms += route_trace.routeA_total_ms_for_midpoint;
        }
        const auto validation_start = Clock::now();
        const bool validation_passed = sample_passes_midpoint_routeA_validation(mid_sample);
        route_trace.validation_ms_for_midpoint = elapsed_ms(validation_start, Clock::now());
        if (validation_passed) {
            routeA_midpoint_used = true;
            if (profile != nullptr) {
                profile->adaptive_midpoint_routeA_success_count += 1;
                profile->adaptive_midpoint_nearest_query_avoided_count += 1;
            }
            (*cache)[theta_cache_key(theta_mid)] = mid_sample;
        } else {
            if (profile != nullptr) {
                profile->adaptive_midpoint_routeA_validation_failure_count += 1;
                profile->adaptive_midpoint_routeA_fallback_count += 1;
            }
            route_trace.routeA_validation_failure_count_for_midpoint = 1;
            route_trace.fallback_used = true;
        }
        if (profile != nullptr && options.adaptive_interval_trace_enabled) {
            profile->routeA_midpoint_traces.push_back(route_trace);
        }
    }
    if (!routeA_midpoint_used) {
        const auto fallback_start = Clock::now();
        mid_sample = evaluate_cached_or_build(
            theta_mid, loader, query_options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
            incoming_theta, options, cache, summary, profile, &mid_from_cache);
        if (profile != nullptr &&
            options.adaptive_midpoint_routeA_from_endpoint_enabled &&
            options.adaptive_midpoint_seed_policy != 0) {
            const double fallback_ms = elapsed_ms(fallback_start, Clock::now());
            profile->adaptive_midpoint_fallback_nearest_query_ms += fallback_ms;
            if (options.adaptive_interval_trace_enabled && !profile->routeA_midpoint_traces.empty()) {
                auto& route_trace = profile->routeA_midpoint_traces.back();
                if (route_trace.fallback_used && std::abs(route_trace.theta_mid - theta_mid) <= 1e-14) {
                    route_trace.fallback_nearest_query_ms_for_midpoint = fallback_ms;
                }
            }
        }
    }
    if (profile != nullptr) {
        profile->topology_adaptive_sampling_ms += elapsed_ms(sample_start, Clock::now());
    }
    if (trace_enabled) {
        trace.midpoint_sample_built = true;
        trace.midpoint_from_cache = mid_from_cache;
        trace.midpoint_branch_count = valid_branch_count(mid_sample);
        trace.classified_topology_change = true;
        trace.subdivided = true;
        trace.subdivision_reason = "raw_count_by_k_mismatch";
        trace.nearest_query_count_for_interval =
            profile->problem1_nearest_node_query_count - query_count_before;
        trace.nearest_query_ms_for_interval = profile->nearest_query_total_ms - query_ms_before;
        profile->adaptive_interval_traces.push_back(trace);
        profile->adaptive_interval_trace_overhead_ms += elapsed_ms(trace_start, Clock::now());
    }
    collect_topology_adaptive(
        loader, query_options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
        incoming_theta, origin_interval_left, theta_left, theta_mid, left_sample, mid_sample, depth + 1, options,
        cache, summary, profile, candidates);
    collect_topology_adaptive(
        loader, query_options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
        incoming_theta, origin_interval_left, theta_mid, theta_right, mid_sample, right_sample, depth + 1, options,
        cache, summary, profile, candidates);
}

ThetaBranch continue_branch(
    const ThetaBranch& endpoint,
    double theta_prime_new,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    Problem2GravityAssistSolverProfile* profile
) {
    ThetaBranch out{};
    out.theta_prime = theta_prime_new;
    out.transfer_revolution = endpoint.transfer_revolution;
    out.target_revolution = endpoint.target_revolution;
    if (!endpoint.derivatives_available || !is_finite(endpoint.d_alpha_d_theta_A)) {
        return out;
    }

    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto& target_params = planet_params::get_planet_params(target_planet);
    const double theta_A_new = normalize_angle_0_2pi(departure_state.theta_global - theta_prime_new);
    const double delta_theta_prime = normalize_angle_minus_pi_pi(theta_prime_new - endpoint.theta_prime);
    const double endpoint_global_alpha = normalize_angle_0_2pi(endpoint.alpha + target_params.orbit.theta_0);
    const double alpha_seed = normalize_angle_0_2pi(endpoint_global_alpha - endpoint.d_alpha_d_theta_A * delta_theta_prime);

    if (profile != nullptr) {
        profile->continuation_refine_count += 1;
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
        options.problem1_max_newton_iterations,
        options.problem1_residual_tolerance_seconds,
        1e-12,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    if (!refined.valid || !refined.branch.valid) {
        return out;
    }

    if (profile != nullptr) {
        profile->derivative_attach_count += 1;
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
    return convert_branch(
        departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime_new, theta_A_new,
        differentiated, options);
}

bool same_continued_branch(const ThetaBranch& lhs, const ThetaBranch& rhs) {
    const double alpha_diff = wrapped_angle_distance(lhs.alpha, rhs.alpha);
    const double time_diff = std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds);
    const double e_prime_diff = std::abs(lhs.outgoing_eccentricity - rhs.outgoing_eccentricity);
    const double residual_diff = std::abs(lhs.slingshot_residual - rhs.slingshot_residual);
    return lhs.transfer_revolution == rhs.transfer_revolution &&
           alpha_diff <= 1e-5 &&
           time_diff <= 1e3 &&
           e_prime_diff <= 1e-4 &&
           residual_diff <= 1e-4 &&
           lhs.residual_usable &&
           rhs.residual_usable &&
           is_finite(lhs.slingshot_residual) &&
           is_finite(rhs.slingshot_residual) &&
           is_finite(lhs.outgoing_eccentricity) &&
           is_finite(rhs.outgoing_eccentricity);
}

bool midpoint_continuation_stable(
    const Candidate& candidate,
    double theta_mid,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    ThetaBranch* selected_mid,
    Problem2GravityAssistSolverProfile* profile
) {
    const auto left_mid = continue_branch(
        candidate.left, theta_mid, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
        incoming_theta, options, profile);
    const auto right_mid = continue_branch(
        candidate.right, theta_mid, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
        incoming_theta, options, profile);
    if (!left_mid.problem1_valid || !right_mid.problem1_valid) {
        return false;
    }
    if (!same_continued_branch(left_mid, right_mid)) {
        return false;
    }
    *selected_mid = std::abs(left_mid.problem1_residual_seconds) <= std::abs(right_mid.problem1_residual_seconds)
        ? left_mid
        : right_mid;
    return true;
}

double sample_quality_score(const Sample& sample) {
    double score = 0.0;
    for (const auto& branch : sample.branches) {
        if (!branch.problem1_valid || !branch.residual_usable || !is_finite(branch.slingshot_residual)) {
            continue;
        }
        score += std::abs(branch.problem1_residual_seconds) + std::abs(branch.slingshot_residual);
    }
    return score;
}

const Sample& choose_midpoint_seed_sample(
    const Sample& left_sample,
    const Sample& right_sample,
    const Problem2GravityAssistSolverOptions& options
) {
    if (options.adaptive_midpoint_seed_policy == 2) {
        return left_sample;
    }
    if (options.adaptive_midpoint_seed_policy == 3) {
        return right_sample;
    }
    if (options.adaptive_midpoint_seed_policy == 4) {
        return sample_quality_score(left_sample) <= sample_quality_score(right_sample) ? left_sample : right_sample;
    }
    const int left_count = valid_branch_count(left_sample);
    const int right_count = valid_branch_count(right_sample);
    if (left_count != right_count) {
        return left_count > right_count ? left_sample : right_sample;
    }
    return left_sample;
}

bool sample_passes_midpoint_routeA_validation(const Sample& sample) {
    if (sample.branches.empty()) {
        return false;
    }
    for (const auto& branch : sample.branches) {
        if (!branch.problem1_valid ||
            !branch.residual_usable ||
            !branch.derivatives_available ||
            !is_finite(branch.alpha) ||
            !is_finite(branch.time_of_flight_seconds) ||
            !is_finite(branch.outgoing_eccentricity) ||
            !is_finite(branch.outgoing_semi_latus_rectum) ||
            !is_finite(branch.slingshot_residual)) {
            return false;
        }
    }
    return true;
}

Sample build_midpoint_sample_from_endpoint_routeA(
    const Sample& seed_sample,
    double theta_mid,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    Problem2GravityAssistSolverProfile* profile
) {
    Sample sample{};
    sample.theta_prime = theta_mid;
    sample.branches.reserve(seed_sample.branches.size());
    for (const auto& branch : seed_sample.branches) {
        if (!branch.problem1_valid || !branch.derivatives_available) {
            continue;
        }
        auto continued = continue_branch(
            branch,
            theta_mid,
            departure_planet,
            target_planet,
            encounter_time,
            phi,
            beta,
            incoming_e,
            incoming_theta,
            options,
            profile);
        if (!continued.problem1_valid ||
            std::abs(continued.problem1_residual_seconds) > options.problem1_residual_tolerance_seconds) {
            continue;
        }
        sample.branches.push_back(continued);
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

BisectionResult bisect_candidate(
    const Candidate& candidate,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options,
    Problem2GravityAssistSolverProfile* profile
) {
    BisectionResult result{};
    result.candidate = candidate;
    double left_theta = candidate.theta_left;
    double right_theta = candidate.theta_right;
    ThetaBranch left_branch = candidate.left;
    ThetaBranch right_branch = candidate.right;
    for (int iteration = 0; iteration < options.bisection_max_iterations; ++iteration) {
        if (profile != nullptr) {
            profile->bisection_iteration_count += 1;
        }
        const double theta_mid = 0.5 * (left_theta + right_theta);
        result.iterations = iteration + 1;
        result.final_width = right_theta - left_theta;

        Candidate current = candidate;
        current.theta_left = left_theta;
        current.theta_right = right_theta;
        current.left = left_branch;
        current.right = right_branch;

        ThetaBranch mid_branch{};
        if (!midpoint_continuation_stable(
                current, theta_mid, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                incoming_theta, options, &mid_branch, profile)) {
            result.invalid_reason = "midpoint_continuation_not_stable";
            return result;
        }
        if (std::abs(mid_branch.slingshot_residual) <= options.bisection_residual_tolerance) {
            result.valid = true;
            result.theta_root = theta_mid;
            result.root_branch = mid_branch;
            return result;
        }
        if ((right_theta - left_theta) <= options.bisection_theta_tolerance) {
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

void accumulate_nearest_node_query_profile(
    const problem1::Problem1NearestNodeQueryProfile& query_profile,
    Problem2SampleBuildPhase phase,
    Problem2GravityAssistSolverProfile* profile
) {
    if (profile == nullptr) {
        return;
    }
    profile->nearest_query_table_loaded_branch_count += query_profile.table_loaded_branch_count;
    profile->nearest_query_seed_branch_count += query_profile.seed_branch_count;
    profile->nearest_query_refine_attempt_count += query_profile.refine_attempt_count;
    profile->nearest_query_refine_success_count += query_profile.refine_success_count;
    profile->nearest_query_refine_fail_count += query_profile.refine_fail_count;
    profile->nearest_query_derivative_attach_attempt_count += query_profile.derivative_attach_attempt_count;
    profile->nearest_query_derivative_attach_success_count += query_profile.derivative_attach_success_count;
    profile->nearest_query_derivative_attach_fail_count += query_profile.derivative_attach_fail_count;
    profile->nearest_query_fallback_direct_solve_count += query_profile.fallback_direct_solve_count;
    profile->nearest_query_fallback_direct_branch_count += query_profile.fallback_direct_branch_count;
    profile->cell_vertex_routeA_attempt_count += query_profile.cell_vertex_routeA_attempt_count;
    profile->cell_vertex_routeA_success_count += query_profile.cell_vertex_routeA_success_count;
    profile->cell_vertex_routeA_failure_count += query_profile.cell_vertex_routeA_failure_count;
    profile->cell_vertex_routeA_no_valid_vertex_count += query_profile.cell_vertex_routeA_no_valid_vertex_count;
    profile->cell_vertex_routeA_no_valid_branch_count += query_profile.cell_vertex_routeA_no_valid_branch_count;
    profile->cell_vertex_routeA_validation_failure_count += query_profile.cell_vertex_routeA_validation_failure_count;
    profile->cell_vertex_loaded_vertex_count += query_profile.cell_vertex_loaded_vertex_count;
    profile->cell_vertex_selected_vertex_count += query_profile.cell_vertex_selected_vertex_count;
    profile->cell_vertex_attempted_branch_count += query_profile.cell_vertex_attempted_branch_count;
    profile->direct_fallback_avoided_count += query_profile.direct_fallback_avoided_count;
    profile->direct_fallback_skipped_count += query_profile.direct_fallback_skipped_count;

    profile->nearest_query_table_load_ms += query_profile.table_load_ms;
    profile->nearest_query_table_offset_build_ms += query_profile.table_offset_build_ms;
    profile->nearest_query_node_read_ms += query_profile.table_node_read_ms;
    profile->nearest_query_seed_generation_ms += query_profile.seed_generation_ms;
    profile->nearest_query_refine_ms += query_profile.refine_ms;
    profile->nearest_query_derivative_attach_ms += query_profile.derivative_attach_ms;
    profile->nearest_query_dedup_ms += query_profile.dedup_ms;
    profile->nearest_query_fallback_direct_solve_ms += query_profile.fallback_direct_solve_ms;
    profile->cell_vertex_selection_ms += query_profile.cell_vertex_selection_ms;
    profile->cell_vertex_routeA_ms += query_profile.cell_vertex_routeA_ms;
    profile->cell_vertex_validation_ms += query_profile.cell_vertex_validation_ms;
    profile->nearest_query_total_ms += query_profile.total_ms;

    if (phase == Problem2SampleBuildPhase::Initial) {
        profile->initial_nearest_query_count += 1;
        profile->initial_nearest_query_refine_attempt_count += query_profile.refine_attempt_count;
        profile->initial_nearest_query_table_load_ms += query_profile.table_load_ms;
        profile->initial_nearest_query_refine_ms += query_profile.refine_ms;
        profile->initial_nearest_query_derivative_attach_ms += query_profile.derivative_attach_ms;
        profile->initial_nearest_query_total_ms += query_profile.total_ms;
    } else {
        profile->adaptive_nearest_query_count += 1;
        profile->adaptive_nearest_query_refine_attempt_count += query_profile.refine_attempt_count;
        profile->adaptive_nearest_query_table_load_ms += query_profile.table_load_ms;
        profile->adaptive_nearest_query_refine_ms += query_profile.refine_ms;
        profile->adaptive_nearest_query_derivative_attach_ms += query_profile.derivative_attach_ms;
        profile->adaptive_nearest_query_total_ms += query_profile.total_ms;
    }
}

Problem2GravityAssistSolution make_solution(const BisectionResult& root) {
    Problem2GravityAssistSolution solution{};
    solution.valid = true;
    solution.theta_prime = root.theta_root;
    solution.alpha = root.root_branch.alpha;
    solution.transfer_revolution = root.root_branch.transfer_revolution;
    solution.target_revolution = root.root_branch.target_revolution;
    solution.time_of_flight_seconds = root.root_branch.time_of_flight_seconds;
    solution.target_time_seconds = root.root_branch.target_time_seconds;
    solution.outgoing_eccentricity = root.root_branch.outgoing_eccentricity;
    solution.outgoing_semi_latus_rectum = root.root_branch.outgoing_semi_latus_rectum;
    solution.slingshot_residual = root.root_branch.slingshot_residual;
    solution.problem1_residual_seconds = root.root_branch.problem1_residual_seconds;
    solution.residual_source = root.root_branch.residual_source;
    solution.boundary_ambiguous = root.root_branch.boundary_ambiguous;
    solution.bisection_iterations = root.iterations;
    solution.final_theta_width = root.final_width;
    solution.origin_was_topology_change = root.candidate.origin_was_topology_change;
    return solution;
}

}  // namespace

Problem2GravityAssistSolverResult solve_problem2_gravity_assist_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options
) {
    const auto total_start = Clock::now();
    Problem2GravityAssistSolverResult result{};
    result.summary.theta_sample_count = options.theta_sample_count;

    try {
        if (options.theta_sample_count < 4) {
            result.error_message = "theta_sample_count_must_be_at_least_4";
            result.profile.total_ms = elapsed_ms(total_start, Clock::now());
            return result;
        }
        if (options.max_transfer_revolution < 0 || options.max_target_revolution < 0) {
            result.error_message = "max_revolution_must_be_non_negative";
            result.profile.total_ms = elapsed_ms(total_start, Clock::now());
            return result;
        }
        if (!(options.bisection_max_iterations > 0)) {
            result.error_message = "bisection_max_iterations_must_be_positive";
            result.profile.total_ms = elapsed_ms(total_start, Clock::now());
            return result;
        }

        problem1::Problem1NearestNodeQueryOptions query_options{};
        query_options.residual_tolerance_seconds = options.problem1_residual_tolerance_seconds;
        query_options.max_newton_iterations = options.problem1_max_newton_iterations;
        query_options.fallback_direct_solve = true;
        query_options.fallback_mode = options.problem1_fallback_mode;
        query_options.cell_vertex_routeA_enabled = options.problem1_cell_vertex_routeA_enabled;
        query_options.cell_vertex_seed_policy = options.problem1_cell_vertex_seed_policy;
        query_options.cell_vertex_top_k_vertices = options.problem1_cell_vertex_top_k_vertices;
        query_options.cell_vertex_max_branch_attempts = options.problem1_cell_vertex_max_branch_attempts;
        query_options.scout_no_direct_fallback = options.problem1_scout_no_direct_fallback;

        const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
        const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
        const double phi = departure_state.varphi;
        const double beta = target_state.varphi;

        std::vector<Sample> samples;
        samples.reserve(static_cast<std::size_t>(options.theta_sample_count));
        std::map<long long, Sample> sample_cache;
        const auto initial_sampling_start = Clock::now();
        for (int i = 0; i < options.theta_sample_count; ++i) {
            const double theta_prime = kTwoPi * static_cast<double>(i) / static_cast<double>(options.theta_sample_count);
            result.profile.initial_sample_count += 1;
            samples.push_back(build_sample(
                loader, query_options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                incoming_theta, theta_prime, Problem2SampleBuildPhase::Initial, options, &result.profile));
            sample_cache[theta_cache_key(theta_prime)] = samples.back();
            if (samples.back().table_fallback_used) {
                result.summary.table_fallback_count += 1;
            }
        }
        result.profile.initial_sampling_ms += elapsed_ms(initial_sampling_start, Clock::now());

        for (int i = 0; i + 1 < options.theta_sample_count; ++i) {
            if (samples[static_cast<std::size_t>(i)].raw_count_by_k ==
                samples[static_cast<std::size_t>(i + 1)].raw_count_by_k) {
                result.summary.initial_stable_interval_count += 1;
            } else {
                result.summary.initial_topology_change_interval_count += 1;
            }
        }

        std::vector<Candidate> candidates;
        const auto candidate_collection_start = Clock::now();
        if (options.topology_adaptive_enabled) {
            for (int i = 0; i + 1 < options.theta_sample_count; ++i) {
                const auto& left = samples[static_cast<std::size_t>(i)];
                const auto& right = samples[static_cast<std::size_t>(i + 1)];
                if (left.raw_count_by_k == right.raw_count_by_k) {
                    append_stable_interval_candidates(i, false, left, right, &result.summary, &candidates);
                } else {
                    collect_topology_adaptive(
                        loader, query_options, departure_planet, target_planet, encounter_time, phi, beta,
                        incoming_e, incoming_theta, i, left.theta_prime, right.theta_prime, left, right, 0, options,
                        &sample_cache, &result.summary, &result.profile, &candidates);
                }
            }
        } else {
            for (int i = 0; i + 1 < options.theta_sample_count; ++i) {
                const auto& left = samples[static_cast<std::size_t>(i)];
                const auto& right = samples[static_cast<std::size_t>(i + 1)];
                if (left.raw_count_by_k == right.raw_count_by_k) {
                    append_stable_interval_candidates(i, false, left, right, &result.summary, &candidates);
                }
            }
        }
        result.summary.sign_change_candidate_count = static_cast<int>(candidates.size());
        result.profile.candidate_collection_ms += elapsed_ms(candidate_collection_start, Clock::now());

        std::vector<Candidate> stable_candidates;
        const auto midpoint_continuation_start = Clock::now();
        for (const auto& candidate : candidates) {
            const double theta_mid = 0.5 * (candidate.theta_left + candidate.theta_right);
            ThetaBranch mid_branch{};
            if (midpoint_continuation_stable(
                    candidate, theta_mid, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                    incoming_theta, options, &mid_branch, &result.profile)) {
                stable_candidates.push_back(candidate);
            } else {
                result.summary.continuation_failure_count += 1;
            }
        }
        result.summary.continuation_stable_candidate_count = static_cast<int>(stable_candidates.size());
        result.profile.midpoint_continuation_ms += elapsed_ms(midpoint_continuation_start, Clock::now());

        std::vector<BisectionResult> successes;
        const auto bisection_start = Clock::now();
        for (const auto& candidate : stable_candidates) {
            result.summary.bisection_attempt_count += 1;
            const BisectionResult bisection = bisect_candidate(
                candidate, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
                options, &result.profile);
            if (!bisection.valid) {
                result.summary.bisection_failure_count += 1;
                continue;
            }
            result.summary.bisection_success_count += 1;
            successes.push_back(bisection);
        }
        result.profile.bisection_ms += elapsed_ms(bisection_start, Clock::now());

        std::vector<BisectionResult> deduped;
        const auto dedup_start = Clock::now();
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
        result.profile.dedup_ms += elapsed_ms(dedup_start, Clock::now());

        result.summary.dedup_success_count = static_cast<int>(deduped.size());
        result.solutions.reserve(deduped.size());
        for (const auto& root : deduped) {
            const auto solution = make_solution(root);
            result.summary.max_abs_slingshot_residual_at_root = std::max(
                result.summary.max_abs_slingshot_residual_at_root, std::abs(solution.slingshot_residual));
            result.summary.max_abs_problem1_residual_seconds_at_root = std::max(
                result.summary.max_abs_problem1_residual_seconds_at_root,
                std::abs(solution.problem1_residual_seconds));
            if (solution.residual_source == Problem2ResidualSource::Strict) {
                result.summary.strict_root_count += 1;
            } else {
                result.summary.relaxed_boundary_root_count += 1;
            }
            result.solutions.push_back(solution);
        }

        result.ok = true;
        result.profile.total_ms = elapsed_ms(total_start, Clock::now());
        return result;
    } catch (const std::exception& ex) {
        result.ok = false;
        result.error_message = ex.what();
        result.profile.total_ms = elapsed_ms(total_start, Clock::now());
        return result;
    }
}

}  // namespace spaceship_cpp::problem2
