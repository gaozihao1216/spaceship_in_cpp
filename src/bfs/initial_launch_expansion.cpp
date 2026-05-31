#include "spaceship_cpp/bfs/initial_launch_expansion.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

struct ResolvedPeriapsisAngle {
    double angle = 0.0;
    bool flipped = false;
    double p_error = 0.0;
};

double transfer_p_from_angle(double radius, double eccentricity, double global_true_anomaly, double periapsis_angle) {
    return radius * (1.0 + eccentricity * std::cos(global_true_anomaly - periapsis_angle));
}

ResolvedPeriapsisAngle resolve_transfer_periapsis_angle(
    double raw_periapsis_angle,
    const planet_params::PlanetState& departure_state,
    const problem1::Problem1SolutionBranch& branch
) {
    ResolvedPeriapsisAngle resolved{};
    resolved.angle = normalize_angle_0_2pi(raw_periapsis_angle);
    if (!is_finite(branch.transfer_p) || !is_finite(branch.transfer_e) || branch.transfer_e == 0.0) {
        resolved.p_error = 0.0;
        return resolved;
    }
    const double normalized_raw = normalize_angle_0_2pi(raw_periapsis_angle);
    const double flipped = normalize_angle_0_2pi(raw_periapsis_angle + kPi);
    const double raw_p = transfer_p_from_angle(
        departure_state.radius, branch.transfer_e, departure_state.theta_global, normalized_raw);
    const double flipped_p = transfer_p_from_angle(
        departure_state.radius, branch.transfer_e, departure_state.theta_global, flipped);
    const double raw_error = std::abs(raw_p - branch.transfer_p);
    const double flipped_error = std::abs(flipped_p - branch.transfer_p);
    if (flipped_error < raw_error) {
        resolved.angle = flipped;
        resolved.flipped = true;
        resolved.p_error = flipped_error;
    } else {
        resolved.angle = normalized_raw;
        resolved.flipped = false;
        resolved.p_error = raw_error;
    }
    return resolved;
}

void update_v_inf_bounds(InitialLaunchExpansionResult* result, double launch_v_inf) {
    if (result->accepted_candidate_count == 0) {
        result->min_launch_v_inf = launch_v_inf;
        result->max_launch_v_inf = launch_v_inf;
        return;
    }
    result->min_launch_v_inf = std::min(result->min_launch_v_inf, launch_v_inf);
    result->max_launch_v_inf = std::max(result->max_launch_v_inf, launch_v_inf);
}

}  // namespace

InitialLaunchExpansionResult expand_initial_launch_with_problem1_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    double initial_theta,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const InitialLaunchExpansionOptions& options
) {
    InitialLaunchExpansionResult result{};
    if (!is_finite(launch_time) || !is_finite(initial_theta)) {
        result.error_message = "invalid_launch_time_or_initial_theta";
        return result;
    }
    if (options.max_transfer_revolution < 0 || options.max_target_revolution < 0) {
        result.error_message = "max_revolution_must_be_non_negative";
        return result;
    }

    result.target_planet_count = static_cast<int>(allowed_first_targets.size());
    const auto departure_state = planet_params::planet_state_at_time(launch_planet, launch_time);
    const double phi = departure_state.varphi;
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - initial_theta);

    problem1::Problem1NearestNodeQueryOptions query_options{};
    query_options.residual_tolerance_seconds = options.problem1_residual_tolerance_seconds;
    query_options.max_newton_iterations = options.problem1_max_newton_iterations;
    query_options.fallback_direct_solve = true;

    for (const auto target_planet : allowed_first_targets) {
        if (!planet_params::is_valid_planet_id(target_planet)) {
            std::ostringstream os;
            os << "invalid_allowed_first_target_raw=" << planet_params::planet_id_raw_value(target_planet);
            result.target_error_messages.push_back(os.str());
            continue;
        }
        const auto target_state = planet_params::planet_state_at_time(target_planet, launch_time);
        const double beta = target_state.varphi;
        result.problem1_query_count += 1;
        problem1::Problem1NearestNodeQueryResult query{};
        try {
            query = problem1::query_problem1_from_2deg_nearest_node(
                loader,
                launch_planet,
                target_planet,
                phi,
                beta,
                theta_A,
                options.max_transfer_revolution,
                options.max_target_revolution,
                query_options);
        } catch (const std::exception& ex) {
            std::ostringstream os;
            os << "initial_problem1_query_exception"
               << " target_planet=" << planet_params::planet_name(target_planet)
               << " target_planet_raw=" << planet_params::planet_id_raw_value(target_planet)
               << " message=" << ex.what();
            result.target_error_messages.push_back(os.str());
            continue;
        }

        result.raw_branch_count += static_cast<int>(query.branches.size());
        for (const auto& branch : query.branches) {
            if (!branch.valid || !is_finite(branch.transfer_e) || !(branch.transfer_e >= 0.0) ||
                !is_finite(branch.transfer_p) || !(branch.transfer_p > 0.0) ||
                !is_finite(branch.time_of_flight_seconds) || !(branch.time_of_flight_seconds > 0.0)) {
                continue;
            }

            const auto resolved_periapsis =
                resolve_transfer_periapsis_angle(initial_theta, departure_state, branch);
            const double launch_v_inf = trajectory::relative_speed_to_planet(
                launch_planet, launch_time, branch.transfer_e, resolved_periapsis.angle);
            if (!is_finite(launch_v_inf) || launch_v_inf < 0.0) {
                continue;
            }
            if (launch_v_inf > options.max_launch_v_inf) {
                result.launch_v_inf_pruned_count += 1;
                continue;
            }

            InitialLaunchCandidate candidate{};
            candidate.valid = true;
            candidate.launch_planet = launch_planet;
            candidate.target_planet = target_planet;
            candidate.launch_time = launch_time;
            candidate.arrival_time = launch_time + branch.time_of_flight_seconds;
            candidate.transfer_time_seconds = branch.time_of_flight_seconds;
            candidate.initial_theta = normalize_angle_0_2pi(initial_theta);
            candidate.alpha = branch.target_arrival_true_anomaly;
            candidate.transfer_eccentricity = branch.transfer_e;
            candidate.transfer_periapsis_angle = resolved_periapsis.angle;
            candidate.transfer_revolution = branch.transfer_revolution;
            candidate.target_revolution = branch.target_revolution;
            candidate.launch_v_inf = launch_v_inf;
            candidate.problem1_residual_seconds = branch.residual_seconds;
            candidate.periapsis_flipped_from_input_theta = resolved_periapsis.flipped;
            candidate.transfer_p_resolution_error = resolved_periapsis.p_error;

            candidate.edge.valid = true;
            candidate.edge.from_planet = launch_planet;
            candidate.edge.to_planet = target_planet;
            candidate.edge.departure_time = launch_time;
            candidate.edge.arrival_time = candidate.arrival_time;
            candidate.edge.transfer_time_seconds = branch.time_of_flight_seconds;
            candidate.edge.outgoing_e = branch.transfer_e;
            candidate.edge.outgoing_p = branch.transfer_p;
            candidate.edge.outgoing_theta = resolved_periapsis.angle;
            candidate.edge.theta_prime = resolved_periapsis.angle;
            candidate.edge.alpha = branch.target_arrival_true_anomaly;
            candidate.edge.transfer_revolution = branch.transfer_revolution;
            candidate.edge.target_revolution = branch.target_revolution;
            candidate.edge.problem1_residual_seconds = branch.residual_seconds;

            candidate.next_state.valid = true;
            candidate.next_state.current_planet = target_planet;
            candidate.next_state.current_time = candidate.edge.arrival_time;
            candidate.next_state.incoming_e = branch.transfer_e;
            candidate.next_state.incoming_theta = resolved_periapsis.angle;
            candidate.next_state.depth = 1;
            candidate.next_state.accumulated_time_seconds = candidate.edge.transfer_time_seconds;
            candidate.next_state.launch_v_inf = launch_v_inf;
            candidate.next_state.accumulated_score =
                launch_v_inf +
                options.time_weight_m_per_s_per_day * (candidate.edge.transfer_time_seconds / 86400.0);

            update_v_inf_bounds(&result, launch_v_inf);
            result.accepted_candidate_count += 1;
            if (resolved_periapsis.flipped) {
                result.periapsis_flipped_selected_count += 1;
            } else {
                result.periapsis_raw_selected_count += 1;
            }
            result.max_transfer_p_resolution_error = std::max(
                result.max_transfer_p_resolution_error, resolved_periapsis.p_error);
            result.candidates.push_back(candidate);
        }
    }

    result.ok = true;
    return result;
}

double compute_arrival_v_inf_from_state(const TrajectorySearchState& state) {
    if (!state.valid) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return trajectory::compute_arrival_v_inf_at_planet(
        state.current_planet,
        state.current_time,
        state.incoming_e,
        state.incoming_theta);
}

}  // namespace spaceship_cpp::bfs
