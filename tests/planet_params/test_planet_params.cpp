/*
 * 文件作用：测试行星参数表。
 * 主要工作：检查行星编号、名称、轨道参数和太阳系物理常量的有效性。
 */
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <cassert>
#include <cmath>
#include <cstring>

namespace {

// 中文说明：提供绝对/相对容差浮点比较，用于核对 J2000 参考轨道与 GM 常量。
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

    // 中文说明：遍历八大行星，验证名称非空、轨道要素（p,e,theta_0,varphi_0）在合法范围内且物理 GM/半径为正。
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

    // 中文说明：验证太阳引力参数 GM_sun 为有限正数。
    assert(planet_params::get_solar_system_physical_params().GM_sun > 0.0);

    // 中文说明：验证行星名称查询、历元常量（J2000/JD/ISO）与 GM_sun 参考值正确。
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

    // 中文说明：抽查 Earth/Mercury/Jupiter/Neptune 的 J2000 轨道角、p 与 GM 与已知参考数据一致。
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
