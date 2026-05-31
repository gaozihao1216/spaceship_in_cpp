#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::problem2 {

struct Problem2GravityAssistTransitionCandidate {
    bool valid = false;

    planet_params::PlanetId encounter_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId target_planet = planet_params::PlanetId::Mercury;

    double encounter_time = 0.0;
    double transfer_time_seconds = 0.0;
    // Uses the transfer time of flight from the Problem 1 branch. target_time_seconds is
    // the matching target-orbit sheet time and should agree up to the Problem 1 residual.
    double arrival_time = 0.0;

    double theta_prime = 0.0;
    double alpha = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;

    double slingshot_residual = 0.0;
    double problem1_residual_seconds = 0.0;

    Problem2ResidualSource residual_source = Problem2ResidualSource::Strict;
    bool boundary_ambiguous = false;

    bool origin_was_topology_change = false;
};

struct Problem2GravityAssistTransitionResult {
    bool ok = false;
    std::string error_message;

    std::vector<Problem2GravityAssistTransitionCandidate> candidates;
    Problem2GravityAssistSolverSummary solver_summary;
    Problem2GravityAssistSolverProfile solver_profile;
};

Problem2GravityAssistTransitionResult generate_problem2_gravity_assist_transitions(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId encounter_planet,
    planet_params::PlanetId target_planet,
    double encounter_time,
    double incoming_e,
    double incoming_theta,
    const Problem2GravityAssistSolverOptions& options
);

}  // namespace spaceship_cpp::problem2
