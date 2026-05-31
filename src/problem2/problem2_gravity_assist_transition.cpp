#include "spaceship_cpp/problem2/problem2_gravity_assist_transition.hpp"

#include <algorithm>
#include <tuple>

namespace spaceship_cpp::problem2 {

Problem2GravityAssistTransitionResult generate_problem2_gravity_assist_transitions(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId encounter_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options
) {
    Problem2GravityAssistTransitionResult result{};

    const auto solver_result = solve_problem2_gravity_assist_with_table(
        loader,
        encounter_planet,
        target_planet,
        encounter_time,
        incoming_e,
        incoming_theta,
        options);
    result.solver_summary = solver_result.summary;
    result.solver_profile = solver_result.profile;
    if (!solver_result.ok) {
        result.ok = false;
        result.error_message = solver_result.error_message;
        return result;
    }

    result.candidates.reserve(solver_result.solutions.size());
    for (const auto& solution : solver_result.solutions) {
        Problem2GravityAssistTransitionCandidate candidate{};
        candidate.valid = solution.valid;
        candidate.encounter_planet = encounter_planet;
        candidate.target_planet = target_planet;
        candidate.encounter_time = encounter_time;
        candidate.transfer_time_seconds = solution.time_of_flight_seconds;
        candidate.arrival_time = encounter_time + solution.time_of_flight_seconds;
        candidate.theta_prime = solution.theta_prime;
        candidate.alpha = solution.alpha;
        candidate.transfer_revolution = solution.transfer_revolution;
        candidate.target_revolution = solution.target_revolution;
        candidate.outgoing_eccentricity = solution.outgoing_eccentricity;
        candidate.outgoing_semi_latus_rectum = solution.outgoing_semi_latus_rectum;
        candidate.slingshot_residual = solution.slingshot_residual;
        candidate.problem1_residual_seconds = solution.problem1_residual_seconds;
        candidate.residual_source = solution.residual_source;
        candidate.boundary_ambiguous = solution.boundary_ambiguous;
        candidate.origin_was_topology_change = solution.origin_was_topology_change;
        result.candidates.push_back(candidate);
    }

    std::sort(result.candidates.begin(), result.candidates.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.arrival_time, lhs.target_planet, lhs.transfer_revolution, lhs.target_revolution) <
               std::tie(rhs.arrival_time, rhs.target_planet, rhs.transfer_revolution, rhs.target_revolution);
    });

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::problem2
