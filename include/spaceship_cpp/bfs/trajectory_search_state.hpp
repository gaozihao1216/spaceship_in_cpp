#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct TrajectorySearchState {
    bool valid = false;

    planet_params::PlanetId current_planet = planet_params::PlanetId::Mercury;

    double current_time = 0.0;

    double incoming_e = 0.0;
    double incoming_theta = 0.0;

    int depth = 0;

    double accumulated_time_seconds = 0.0;
    double launch_v_inf = 0.0;
    double accumulated_score = 0.0;
};

struct TrajectorySearchEdge {
    bool valid = false;

    planet_params::PlanetId from_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId to_planet = planet_params::PlanetId::Mercury;

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

    bool boundary_ambiguous = false;
    bool origin_was_topology_change = false;
};

struct TrajectorySearchExpansionResult {
    bool ok = false;
    std::string error_message;

    std::vector<TrajectorySearchEdge> edges;
    std::vector<TrajectorySearchState> next_states;
};

}  // namespace spaceship_cpp::bfs
