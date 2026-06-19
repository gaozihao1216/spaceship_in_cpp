/*
 * 文件作用：测试轨道速度辅助函数。
 * 主要工作：验证行星速度、转移轨道速度和相对速度计算的基本物理一致性。
 */
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

    // 中文说明：验证地球相对自身（同轨道参数）的相对速度应接近零，排除自洽性错误。
    const auto earth = planet_params::PlanetId::Earth;
    const auto& earth_params = planet_params::get_planet_params(earth);
    const double time = 0.17 * planet_params::planet_orbital_period(earth);
    const double earth_self_relative_speed = trajectory::relative_speed_to_planet(
        earth, time, earth_params.orbit.e, earth_params.orbit.theta_0);
    assert(std::isfinite(earth_self_relative_speed));
    assert(earth_self_relative_speed <= 1e-9);

    // 中文说明：验证人造转移轨道（e/theta_0 与地球不同）相对地球速度为有限非负值。
    const double artificial_relative_speed = trajectory::relative_speed_to_planet(earth, time, 0.3, 0.4);
    assert(std::isfinite(artificial_relative_speed));
    assert(artificial_relative_speed >= 0.0);

    std::cout << "OrbitVelocityHelperSummary\n";
    std::cout << "earth_self_relative_speed=" << earth_self_relative_speed << '\n';
    std::cout << "artificial_relative_speed=" << artificial_relative_speed << '\n';
    std::cout << "velocity_helper_ok=1\n";
    return 0;
}
