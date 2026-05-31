#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace trajectory = spaceship_cpp::trajectory;

    std::cout << std::setprecision(12) << std::scientific;

    const auto earth = planet_params::PlanetId::Earth;
    const auto& earth_params = planet_params::get_planet_params(earth);
    const double time = 0.17 * planet_params::planet_orbital_period(earth);
    const double earth_self_relative_speed = trajectory::relative_speed_to_planet(
        earth, time, earth_params.orbit.e, earth_params.orbit.theta_0);
    assert(std::isfinite(earth_self_relative_speed));
    assert(earth_self_relative_speed <= 1e-9);

    const double artificial_relative_speed = trajectory::relative_speed_to_planet(earth, time, 0.3, 0.4);
    assert(std::isfinite(artificial_relative_speed));
    assert(artificial_relative_speed >= 0.0);

    std::cout << "OrbitVelocityHelperSummary\n";
    std::cout << "earth_self_relative_speed=" << earth_self_relative_speed << '\n';
    std::cout << "artificial_relative_speed=" << artificial_relative_speed << '\n';
    std::cout << "velocity_helper_ok=1\n";
    return 0;
}
