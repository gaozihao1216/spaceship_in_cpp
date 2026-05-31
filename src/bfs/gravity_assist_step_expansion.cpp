#include "spaceship_cpp/bfs/gravity_assist_step_expansion.hpp"

#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_expansion.hpp"

namespace spaceship_cpp::bfs {

TrajectorySearchExpansionResult expand_trajectory_state_by_gravity_assist(
    const problem1::Problem1RootTable2DegLoader& loader,
    const TrajectorySearchState& state,
    planet_params::PlanetId target_planet,
    const problem2::Problem2GravityAssistSolverOptions& options
) {
    TrajectorySearchExpansionResult result{};
    if (!state.valid) {
        result.ok = false;
        result.error_message = "invalid_trajectory_search_state";
        return result;
    }
    if (!planet_params::is_valid_planet_id(state.current_planet)) {
        result.ok = false;
        result.error_message = "invalid_current_planet";
        return result;
    }
    if (!planet_params::is_valid_planet_id(target_planet)) {
        result.ok = false;
        result.error_message = "invalid_target_planet";
        return result;
    }

    problem2::Problem2GravityAssistExpansionInput input{};
    input.encounter_planet = state.current_planet;
    input.target_planet = target_planet;
    input.encounter_time = state.current_time;
    input.incoming_e = state.incoming_e;
    input.incoming_theta = global_periapsis_angle_to_problem2_local(
        state.current_planet,
        state.incoming_theta);

    const auto expansion = problem2::expand_one_gravity_assist_step_with_table_profiled(loader, input, options);
    result.problem2_solver_profile = expansion.solver_profile;
    if (!expansion.ok) {
        result.ok = false;
        result.error_message = "problem2_gravity_assist_expansion_failed";
        return result;
    }
    result.edges.reserve(expansion.states.size());
    result.next_states.reserve(expansion.states.size());

    for (const auto& expanded : expansion.states) {
        TrajectorySearchEdge edge{};
        edge.valid = expanded.valid;
        edge.from_planet = expanded.current_planet;
        edge.to_planet = expanded.next_planet;
        edge.departure_time = expanded.departure_time;
        edge.arrival_time = expanded.arrival_time;
        edge.transfer_time_seconds = expanded.transfer_time_seconds;
        edge.outgoing_e = expanded.outgoing_e;
        edge.outgoing_p = expanded.outgoing_p;
        edge.outgoing_theta = problem2_local_periapsis_angle_to_global(
            state.current_planet,
            expanded.theta_prime);
        edge.theta_prime = expanded.theta_prime;
        edge.alpha = expanded.alpha;
        edge.transfer_revolution = expanded.transfer_revolution;
        edge.target_revolution = expanded.target_revolution;
        edge.slingshot_residual = expanded.slingshot_residual;
        edge.problem1_residual_seconds = expanded.problem1_residual_seconds;
        edge.residual_source = expanded.residual_source;
        edge.boundary_ambiguous = expanded.boundary_ambiguous;
        edge.origin_was_topology_change = expanded.origin_was_topology_change;
        result.edges.push_back(edge);

        TrajectorySearchState next_state{};
        next_state.valid = edge.valid;
        next_state.current_planet = edge.to_planet;
        next_state.current_time = edge.arrival_time;
        next_state.incoming_e = edge.outgoing_e;
        next_state.incoming_theta = edge.outgoing_theta;
        next_state.depth = state.depth + 1;
        next_state.accumulated_time_seconds = state.accumulated_time_seconds + edge.transfer_time_seconds;
        next_state.launch_v_inf = state.launch_v_inf;
        next_state.accumulated_score = state.accumulated_score;
        result.next_states.push_back(next_state);
    }

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
