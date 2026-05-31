#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

namespace spaceship_cpp::bfs {

double global_periapsis_angle_to_problem2_local(
    planet_params::PlanetId flyby_planet,
    double global_periapsis_angle
);

double problem2_local_periapsis_angle_to_global(
    planet_params::PlanetId flyby_planet,
    double local_periapsis_angle
);

}  // namespace spaceship_cpp::bfs
