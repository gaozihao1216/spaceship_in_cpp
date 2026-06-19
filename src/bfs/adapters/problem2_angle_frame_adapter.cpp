#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"

#include "spaceship_cpp/common/common.hpp"

namespace spaceship_cpp::bfs {

double global_periapsis_angle_to_problem2_local(
    planet_params::PlanetId flyby_planet,
    double global_periapsis_angle
) {
    const auto& planet = planet_params::get_planet_params(flyby_planet);
    return common::normalize_angle_0_2pi(global_periapsis_angle - planet.orbit.theta_0);
}

double problem2_local_periapsis_angle_to_global(
    planet_params::PlanetId flyby_planet,
    double local_periapsis_angle
) {
    const auto& planet = planet_params::get_planet_params(flyby_planet);
    return common::normalize_angle_0_2pi(local_periapsis_angle + planet.orbit.theta_0);
}

}  // namespace spaceship_cpp::bfs
