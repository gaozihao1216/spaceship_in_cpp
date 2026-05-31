#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <cassert>
#include <cmath>
#include <cstring>

namespace {

bool approx_equal(double a, double b, double abs_tol, double rel_tol) {
    const double diff = std::abs(a - b);
    if (diff <= abs_tol) {
        return true;
    }
    return diff <= rel_tol * std::max(std::abs(a), std::abs(b));
}

}  // namespace

int main() {
    constexpr double kTwoPi = 6.28318530717958647692;
    namespace planet_params = spaceship_cpp::planet_params;

    const auto& planets = planet_params::all_planet_params();
    assert(planets.size() == 8);

    for (const auto& planet : planets) {
        assert(planet.name != nullptr);
        assert(planet.name[0] != '\0');
        assert(planet.orbit.p > 0.0);
        assert(planet.orbit.e >= 0.0);
        assert(std::isfinite(planet.orbit.theta_0));
        assert(std::isfinite(planet.orbit.varphi_0));
        assert(planet.orbit.theta_0 >= 0.0);
        assert(planet.orbit.theta_0 < kTwoPi);
        assert(planet.orbit.varphi_0 >= 0.0);
        assert(planet.orbit.varphi_0 < kTwoPi);
        assert(planet.physical.radius > 0.0);
        assert(planet.physical.GM > 0.0);
    }

    assert(planet_params::get_solar_system_physical_params().GM_sun > 0.0);

    assert(std::strcmp(planet_params::get_planet_params(planet_params::PlanetId::Earth).name, "Earth") == 0);
    assert(std::strcmp(planet_params::planet_name(planet_params::PlanetId::Mercury), "Mercury") == 0);
    assert(std::strcmp(planet_params::kPlanetParamsEpochName, "J2000") == 0);
    assert(planet_params::kPlanetParamsEpochJulianDate == 2451545.0);
    assert(std::strcmp(planet_params::kPlanetParamsEpochIso, "2000-01-01T12:00:00Z") == 0);
    assert(approx_equal(
        planet_params::get_solar_system_physical_params().GM_sun,
        1.32712440041279419e20,
        1.0,
        1e-15));

    const auto& earth = planet_params::get_planet_params(planet_params::PlanetId::Earth);
    const auto& mercury = planet_params::get_planet_params(planet_params::PlanetId::Mercury);
    const auto& jupiter = planet_params::get_planet_params(planet_params::PlanetId::Jupiter);
    const auto& neptune = planet_params::get_planet_params(planet_params::PlanetId::Neptune);

    assert(approx_equal(earth.orbit.theta_0, 1.796601474049, 1e-12, 1e-12));
    assert(approx_equal(earth.orbit.varphi_0, 6.238548481796, 1e-12, 1e-12));
    assert(approx_equal(mercury.orbit.p, 5.546046912930e10, 1e-3, 1e-12));
    assert(approx_equal(neptune.orbit.varphi_0, 4.519493198198, 1e-12, 1e-12));
    assert(approx_equal(earth.physical.GM, 3.98600435507e14, 1e-3, 1e-12));
    assert(approx_equal(jupiter.physical.GM, 1.26712764100e17, 1e-3, 1e-12));

    return 0;
}
