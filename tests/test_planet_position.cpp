/*
 * 文件作用：测试行星随时间的位置计算。
 * 主要工作：验证 J2000 附近的半径、真近点角和全局角归一化行为。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <cassert>
#include <cmath>

namespace {

bool approx_equal(double a, double b, double eps = 1e-9) {
    return std::abs(a - b) <= eps;
}

bool angle_close(double a, double b, double eps = 1e-9) {
    return std::abs(spaceship_cpp::common::normalize_angle_minus_pi_pi(a - b)) <= eps;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;

    for (const auto& planet : planet_params::all_planet_params()) {
        const double varphi_at_epoch = planet_params::planet_true_anomaly_at_time(planet.id, 0.0);
        assert(angle_close(varphi_at_epoch, planet.orbit.varphi_0));

        const planet_params::PlanetState state = planet_params::planet_state_at_time(planet.id, 0.0);
        assert(angle_close(state.theta_global, common::normalize_angle_0_2pi(planet.orbit.theta_0 + planet.orbit.varphi_0)));

        const double expected_radius = planet.orbit.p / (1.0 + planet.orbit.e * std::cos(state.varphi));
        assert(approx_equal(state.radius, expected_radius, 1e-6));
        assert(approx_equal(state.x, state.radius * std::cos(state.theta_global), 1e-6));
        assert(approx_equal(state.y, state.radius * std::sin(state.theta_global), 1e-6));

        const double period = planet_params::planet_orbital_period(planet.id);
        assert(angle_close(planet_params::planet_true_anomaly_at_time(planet.id, period), planet.orbit.varphi_0, 1e-8));
        assert(angle_close(planet_params::planet_true_anomaly_at_time(planet.id, -period), planet.orbit.varphi_0, 1e-8));

        const double large_positive = planet_params::planet_true_anomaly_at_time(planet.id, 1000.0 * period);
        const double large_negative = planet_params::planet_true_anomaly_at_time(planet.id, -1000.0 * period);
        assert(angle_close(large_positive, planet.orbit.varphi_0, 1e-9));
        assert(angle_close(large_negative, planet.orbit.varphi_0, 1e-9));

        const double quarter_period = 0.25 * period;
        const double quarter_reference = planet_params::planet_true_anomaly_at_time(planet.id, quarter_period);
        const double quarter_large = planet_params::planet_true_anomaly_at_time(planet.id, (1000.0 + 0.25) * period);
        assert(angle_close(quarter_reference, quarter_large, 1e-9));

        const double negative_quarter_reference = planet_params::planet_true_anomaly_at_time(planet.id, -quarter_period);
        const double negative_quarter_large = planet_params::planet_true_anomaly_at_time(planet.id, -(1000.0 + 0.25) * period);
        assert(angle_close(negative_quarter_reference, negative_quarter_large, 1e-9));

        for (const double time : {0.0, 0.1 * period, 0.5 * period, -0.3 * period}) {
            const double dt_reduced = std::remainder(time, period);
            const double orbit_scale = std::sqrt(
                planet_params::get_solar_system_physical_params().GM_sun /
                (planet.orbit.p * planet.orbit.p * planet.orbit.p));
            const double F_expected =
                common::orbit_F(planet.orbit.e, planet.orbit.varphi_0) + dt_reduced * orbit_scale;
            const double f = planet_params::planet_true_anomaly_at_time(planet.id, time);
            assert(f >= 0.0);
            assert(f < common::kTwoPi);

            const double phase_error = common::orbit_F(planet.orbit.e, f) - F_expected;
            const double F_period = common::kTwoPi / std::pow(1.0 - planet.orbit.e * planet.orbit.e, 1.5);
            const double phase_error_reduced = std::remainder(phase_error, F_period);
            assert(std::abs(phase_error_reduced) < 1e-8);
        }
    }

    const double earth_mean_motion = planet_params::planet_mean_motion(planet_params::PlanetId::Earth);
    assert(earth_mean_motion > 0.0);
    const double earth_period = planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
    assert(earth_period > 3.0e7);
    assert(earth_period < 3.3e7);

    return 0;
}
