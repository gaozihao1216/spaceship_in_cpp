#include "spaceship_cpp/problem2/problem2_gravity_assist_expansion.hpp"

#include "spaceship_cpp/problem2/problem2_gravity_assist_transition.hpp"

#include <algorithm>

namespace spaceship_cpp::problem2 {

std::vector<Problem2GravityAssistExpandedState> expand_one_gravity_assist_step_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    const Problem2GravityAssistExpansionInput& input,
    const Problem2GravityAssistSolverOptions& options
) {
    const auto transition_result = generate_problem2_gravity_assist_transitions(
        loader,
        input.encounter_planet,
        input.target_planet,
        input.encounter_time,
        input.incoming_e,
        input.incoming_theta,
        options);
    if (!transition_result.ok) {
        return {};
    }

    std::vector<Problem2GravityAssistExpandedState> states;
    states.reserve(transition_result.candidates.size());
    for (const auto& candidate : transition_result.candidates) {
        Problem2GravityAssistExpandedState state{};
        state.valid = candidate.valid;
        state.current_planet = candidate.encounter_planet;
        state.next_planet = candidate.target_planet;
        state.departure_time = candidate.encounter_time;
        state.arrival_time = candidate.arrival_time;
        state.transfer_time_seconds = candidate.transfer_time_seconds;
        state.outgoing_e = candidate.outgoing_eccentricity;
        state.outgoing_p = candidate.outgoing_semi_latus_rectum;
        state.outgoing_theta = candidate.theta_prime;
        state.theta_prime = candidate.theta_prime;
        state.alpha = candidate.alpha;
        state.transfer_revolution = candidate.transfer_revolution;
        state.target_revolution = candidate.target_revolution;
        state.slingshot_residual = candidate.slingshot_residual;
        state.problem1_residual_seconds = candidate.problem1_residual_seconds;
        state.residual_source = candidate.residual_source;
        state.boundary_ambiguous = candidate.boundary_ambiguous;
        state.origin_was_topology_change = candidate.origin_was_topology_change;
        states.push_back(state);
    }

    std::sort(states.begin(), states.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.arrival_time < rhs.arrival_time;
    });
    return states;
}

Problem2GravityAssistExpansionResult expand_one_gravity_assist_step_with_table_profiled(
    const problem1::Problem1RootTable2DegLoader& loader,
    const Problem2GravityAssistExpansionInput& input,
    const Problem2GravityAssistSolverOptions& options
) {
    Problem2GravityAssistExpansionResult result{};
    const auto transition_result = generate_problem2_gravity_assist_transitions(
        loader,
        input.encounter_planet,
        input.target_planet,
        input.encounter_time,
        input.incoming_e,
        input.incoming_theta,
        options);
    result.ok = transition_result.ok;
    result.solver_summary = transition_result.solver_summary;
    result.solver_profile = transition_result.solver_profile;
    if (!transition_result.ok) {
        return result;
    }

    result.states.reserve(transition_result.candidates.size());
    for (const auto& candidate : transition_result.candidates) {
        Problem2GravityAssistExpandedState state{};
        state.valid = candidate.valid;
        state.current_planet = candidate.encounter_planet;
        state.next_planet = candidate.target_planet;
        state.departure_time = candidate.encounter_time;
        state.arrival_time = candidate.arrival_time;
        state.transfer_time_seconds = candidate.transfer_time_seconds;
        state.outgoing_e = candidate.outgoing_eccentricity;
        state.outgoing_p = candidate.outgoing_semi_latus_rectum;
        state.outgoing_theta = candidate.theta_prime;
        state.theta_prime = candidate.theta_prime;
        state.alpha = candidate.alpha;
        state.transfer_revolution = candidate.transfer_revolution;
        state.target_revolution = candidate.target_revolution;
        state.slingshot_residual = candidate.slingshot_residual;
        state.problem1_residual_seconds = candidate.problem1_residual_seconds;
        state.residual_source = candidate.residual_source;
        state.boundary_ambiguous = candidate.boundary_ambiguous;
        state.origin_was_topology_change = candidate.origin_was_topology_change;
        result.states.push_back(state);
    }

    std::sort(result.states.begin(), result.states.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.arrival_time < rhs.arrival_time;
    });
    return result;
}

}  // namespace spaceship_cpp::problem2
