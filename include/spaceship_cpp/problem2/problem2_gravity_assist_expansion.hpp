#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"

#include <vector>

namespace spaceship_cpp::problem2 {

struct Problem2GravityAssistExpansionInput {
    planet_params::PlanetId encounter_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId target_planet = planet_params::PlanetId::Mercury;

    double encounter_time = 0.0;

    double incoming_e = 0.0;
    double incoming_theta = 0.0;
};

struct Problem2GravityAssistExpandedState {
    bool valid = false;

    planet_params::PlanetId current_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId next_planet = planet_params::PlanetId::Mercury;

    double departure_time = 0.0;
    double arrival_time = 0.0;
    double transfer_time_seconds = 0.0;

    double outgoing_e = 0.0;
    double outgoing_p = 0.0;
    double outgoing_theta = 0.0;

    double theta_prime = 0.0;
    double alpha = 0.0;

    int transfer_revolution = 0;
    int target_revolution = 0;

    double slingshot_residual = 0.0;
    double problem1_residual_seconds = 0.0;

    Problem2ResidualSource residual_source = Problem2ResidualSource::Strict;
    bool boundary_ambiguous = false;
    bool origin_was_topology_change = false;
};

std::vector<Problem2GravityAssistExpandedState> expand_one_gravity_assist_step_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    const Problem2GravityAssistExpansionInput& input,
    const Problem2GravityAssistSolverOptions& options
);

struct Problem2GravityAssistExpansionResult {
    bool ok = false;
    std::vector<Problem2GravityAssistExpandedState> states;
    Problem2GravityAssistSolverSummary solver_summary;
    Problem2GravityAssistSolverProfile solver_profile;
};

Problem2GravityAssistExpansionResult expand_one_gravity_assist_step_with_table_profiled(
    const problem1::Problem1RootTable2DegLoader& loader,
    const Problem2GravityAssistExpansionInput& input,
    const Problem2GravityAssistSolverOptions& options
);

}  // namespace spaceship_cpp::problem2
