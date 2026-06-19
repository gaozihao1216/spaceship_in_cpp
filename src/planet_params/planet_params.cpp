/*
 * 文件作用：实现行星参数表和行星状态计算。
 * 主要工作：保存太阳系与行星轨道常量，并按时间求解行星半径、真近点角和全局角。
 */
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"

#include <array>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace spaceship_cpp::planet_params {

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::orbit_F;
using spaceship_cpp::common::orbit_F_xi_derivative;

// Radius values are approximate mean radii retained from the previous dataset.
constexpr PlanetParams kMercury{
    PlanetId::Mercury,
    "Mercury",
    {
        // JPL Table 1 J2000:
        // a = 0.38709927 au, e = 0.20563593,
        // L = 252.25032350 deg, long.peri. = 77.45779628 deg
        5.546046912930e10,
        0.20563593,
        1.351893576425,
        3.080381544911,
    },
    {
        2439700.0,
        2.2031868551e13,
    },
};

constexpr PlanetParams kVenus{
    PlanetId::Venus,
    "Venus",
    {
        // JPL Table 1 J2000:
        // a = 0.72333566 au, e = 0.00677672,
        // L = 181.97909950 deg, long.peri. = 131.60246718 deg
        1.082045051326e11,
        0.00677672,
        2.296896356039,
        0.889734178645,
    },
    {
        6051800.0,
        3.24858592000e14,
    },
};

constexpr PlanetParams kEarth{
    PlanetId::Earth,
    "Earth",
    {
        // JPL Table 1 J2000 (Earth uses EM Bary as a temporary approximation):
        // a = 1.00000261 au, e = 0.01671123,
        // L = 100.46457166 deg, long.peri. = 102.93768193 deg
        1.495564835609e11,
        0.01671123,
        1.796601474049,
        6.238548481796,
    },
    {
        6371000.0,
        3.98600435507e14,
    },
};

constexpr PlanetParams kMars{
    PlanetId::Mars,
    "Mars",
    {
        // JPL Table 1 J2000:
        // a = 1.52371034 au, e = 0.09339410,
        // L = -4.55343205 deg, long.peri. = -23.94362959 deg
        2.259555920295e11,
        0.09339410,
        5.865290135956,
        0.408009787494,
    },
    {
        3389500.0,
        4.2828375816e13,
    },
};

constexpr PlanetParams kJupiter{
    PlanetId::Jupiter,
    "Jupiter",
    {
        // JPL Table 1 J2000:
        // a = 5.20288700 au, e = 0.04838624,
        // L = 34.39644051 deg, long.peri. = 14.72847983 deg
        7.765185432069e11,
        0.04838624,
        0.257060466847,
        0.377796289891,
    },
    {
        69911000.0,
        1.26712764100e17,
    },
};

constexpr PlanetParams kSaturn{
    PlanetId::Saturn,
    "Saturn",
    {
        // JPL Table 1 J2000:
        // a = 9.53667594 au, e = 0.05386179,
        // L = 49.95424423 deg, long.peri. = 92.59887831 deg
        1.422527523057e12,
        0.05386179,
        1.616155310163,
        5.462200522686,
    },
    {
        58232000.0,
        3.79405848418e16,
    },
};

constexpr PlanetParams kUranus{
    PlanetId::Uranus,
    "Uranus",
    {
        // JPL Table 1 J2000:
        // a = 19.18916464 au, e = 0.04725744,
        // L = 313.23810451 deg, long.peri. = 170.95427630 deg
        2.864247228412e12,
        0.04725744,
        2.983714991799,
        2.538527378245,
    },
    {
        25362000.0,
        5.79455640000e15,
    },
};

constexpr PlanetParams kNeptune{
    PlanetId::Neptune,
    "Neptune",
    {
        // JPL Table 1 J2000:
        // a = 30.06992276 au, e = 0.00859048,
        // L = -55.12002969 deg, long.peri. = 44.96476227 deg
        4.498064451788e12,
        0.00859048,
        0.784783148988,
        4.519493198198,
    },
    {
        24622000.0,
        6.83652710058e15,
    },
};

constexpr std::array<PlanetParams, 8> kAllPlanetParams{
    kMercury,
    kVenus,
    kEarth,
    kMars,
    kJupiter,
    kSaturn,
    kUranus,
    kNeptune,
};

constexpr SolarSystemPhysicalParams kSolarSystemPhysicalParams{
    1.32712440041279419e20,
};

// 由半通径 p 和偏心率 e 计算半长轴 a = p/(1-e²)；
// 解决平均角速度和轨道周期公式需要半长轴的问题。
double semi_major_axis_from_orbit(const PlanetOrbitParams& orbit) {
    return orbit.p / (1.0 - orbit.e * orbit.e);
}

// 将绝对时间折叠到一个轨道周期内；解决按时间推进真近点角时
// 大时间偏移需等价到 [0, T) 的问题。
double reduce_time_by_orbital_period(double time_seconds, double period_seconds) {
    if (!is_finite(time_seconds) || !is_finite(period_seconds) || !(period_seconds > 0.0)) {
        throw std::invalid_argument("invalid time or orbital period");
    }
    return std::remainder(time_seconds, period_seconds);
}

// 已知目标轨道积分值 F，用牛顿法反解椭圆轨道真近点角 xi；
// 解决"给定时刻行星应在哪个相位"的核心逆问题。
double solve_true_anomaly_from_orbit_F(double e, double target_F) {
    if (!is_finite(e) || !is_finite(target_F)) {
        throw std::invalid_argument("solve_true_anomaly_from_orbit_F requires finite inputs");
    }
    if (!(e >= 0.0 && e < 1.0)) {
        throw std::invalid_argument("solve_true_anomaly_from_orbit_F supports only elliptic orbits");
    }

    const double delta = 1.0 - e * e;
    double xi = target_F * delta * std::sqrt(delta);

    for (int iteration = 0; iteration < 30; ++iteration) {
        const double residual = orbit_F(e, xi) - target_F;
        const double derivative = orbit_F_xi_derivative(e, xi);
        double step = residual / derivative;
        if (step > spaceship_cpp::common::kPi) {
            step = spaceship_cpp::common::kPi;
        } else if (step < -spaceship_cpp::common::kPi) {
            step = -spaceship_cpp::common::kPi;
        }
        xi -= step;
        if (!is_finite(xi)) {
            throw std::runtime_error("solve_true_anomaly_from_orbit_F produced non-finite iterate");
        }
        if (std::abs(step) < 1e-13) {
            return xi;
        }
    }

    if (!is_finite(xi)) {
        throw std::runtime_error("solve_true_anomaly_from_orbit_F failed");
    }
    return xi;
}

// 启动时校验单颗行星静态参数的物理合法性；防止错误常数进入运行时计算。
void validate_planet_params(const PlanetParams& planet) {
    if (planet.name == nullptr || planet.name[0] == '\0') {
        throw std::runtime_error("planet name must be non-empty");
    }
    if (!(planet.orbit.p > 0.0)) {
        throw std::runtime_error("planet semi-latus rectum must be positive");
    }
    if (!(planet.orbit.e >= 0.0)) {
        throw std::runtime_error("planet eccentricity must be non-negative");
    }
    if (!is_finite(planet.orbit.theta_0)) {
        throw std::runtime_error("planet theta_0 must be finite");
    }
    if (!is_finite(planet.orbit.varphi_0)) {
        throw std::runtime_error("planet varphi_0 must be finite");
    }
    if (!(planet.orbit.theta_0 >= 0.0 && planet.orbit.theta_0 < kTwoPi)) {
        throw std::runtime_error("planet theta_0 must be in [0, 2pi)");
    }
    if (!(planet.orbit.varphi_0 >= 0.0 && planet.orbit.varphi_0 < kTwoPi)) {
        throw std::runtime_error("planet varphi_0 must be in [0, 2pi)");
    }
    if (!(planet.physical.radius > 0.0)) {
        throw std::runtime_error("planet radius must be positive");
    }
    if (!(planet.physical.GM > 0.0)) {
        throw std::runtime_error("planet GM must be positive");
    }
}

// 启动时校验太阳引力参数 GM_sun 为正。
void validate_solar_system_physical_params(const SolarSystemPhysicalParams& params) {
    if (!(params.GM_sun > 0.0)) {
        throw std::runtime_error("GM_sun must be positive");
    }
}

// PlanetId 枚举 → 数组下标；解决八大行星参数表的 O(1) 索引访问。
std::size_t planet_index(PlanetId id) {
    switch (id) {
        case PlanetId::Mercury:
            return 0;
        case PlanetId::Venus:
            return 1;
        case PlanetId::Earth:
            return 2;
        case PlanetId::Mars:
            return 3;
        case PlanetId::Jupiter:
            return 4;
        case PlanetId::Saturn:
            return 5;
        case PlanetId::Uranus:
            return 6;
        case PlanetId::Neptune:
            return 7;
    }
    throw std::invalid_argument("unknown PlanetId");
}

// 程序启动时自动触发行星参数校验（静态初始化）。
struct PlanetParamsValidationRunner {
    PlanetParamsValidationRunner() {
        for (const PlanetParams& planet : kAllPlanetParams) {
            validate_planet_params(planet);
        }
        validate_solar_system_physical_params(kSolarSystemPhysicalParams);
    }
};

const PlanetParamsValidationRunner kValidationRunner{};

}  // namespace

// 按 PlanetId 返回对应行星的静态轨道/物理参数。
const PlanetParams& get_planet_params(PlanetId id) {
    (void)kValidationRunner;
    return kAllPlanetParams.at(planet_index(id));
}

// 返回太阳标准引力参数 GM_sun。
const SolarSystemPhysicalParams& get_solar_system_physical_params() {
    (void)kValidationRunner;
    return kSolarSystemPhysicalParams;
}

// 返回全部 8 颗行星的参数数组，供遍历和诊断使用。
const std::array<PlanetParams, 8>& all_planet_params() {
    (void)kValidationRunner;
    return kAllPlanetParams;
}

// 返回行星英文名称字符串。
const char* planet_name(PlanetId id) {
    return get_planet_params(id).name;
}

// 判断 PlanetId 是否为已知的八大行星之一。
bool is_valid_planet_id(PlanetId id) {
    switch (id) {
        case PlanetId::Mercury:
        case PlanetId::Venus:
        case PlanetId::Earth:
        case PlanetId::Mars:
        case PlanetId::Jupiter:
        case PlanetId::Saturn:
        case PlanetId::Uranus:
        case PlanetId::Neptune:
            return true;
    }
    return false;
}

// 返回 PlanetId 的底层整数值（0=Mercury, …, 7=Neptune）。
int planet_id_raw_value(PlanetId id) {
    return static_cast<int>(id);
}

// 计算行星平均角速度 n = √(GM/a³)；解决按时间线性推进轨道相位的问题。
double planet_mean_motion(PlanetId id) {
    const PlanetParams& planet = get_planet_params(id);
    const double semi_major_axis = semi_major_axis_from_orbit(planet.orbit);
    return std::sqrt(get_solar_system_physical_params().GM_sun /
        (semi_major_axis * semi_major_axis * semi_major_axis));
}

// 计算行星轨道周期 T = 2π/n。
double planet_orbital_period(PlanetId id) {
    return kTwoPi / planet_mean_motion(id);
}

// 计算指定 J2000 秒偏移时刻的行星真近点角（归一化到 [0,2π)）；
// 解决 Problem 1/2 需要知道行星在任意时刻相位的问题。
double planet_true_anomaly_at_time(PlanetId id, double time_seconds_since_j2000) {
    if (!is_finite(time_seconds_since_j2000)) {
        throw std::invalid_argument("time_seconds_since_j2000 must be finite");
    }

    const PlanetParams& planet = get_planet_params(id);
    const double period = planet_orbital_period(id);
    const double reduced_time = reduce_time_by_orbital_period(time_seconds_since_j2000, period);
    const double orbit_scale = std::sqrt(
        get_solar_system_physical_params().GM_sun /
        (planet.orbit.p * planet.orbit.p * planet.orbit.p));
    const double target_F = orbit_F(planet.orbit.e, planet.orbit.varphi_0) + reduced_time * orbit_scale;
    const double varphi_unwrapped = solve_true_anomaly_from_orbit_F(planet.orbit.e, target_F);
    return normalize_angle_0_2pi(varphi_unwrapped);
}

// 由真近点角计算日心半径 r = p/(1+e·cos(φ))。
double planet_radius_at_true_anomaly(PlanetId id, double varphi) {
    if (!is_finite(varphi)) {
        throw std::invalid_argument("varphi must be finite");
    }

    const PlanetParams& planet = get_planet_params(id);
    return planet.orbit.p / (1.0 + planet.orbit.e * std::cos(varphi));
}

// 一站式计算指定时刻的行星状态（真近点角、全局角、半径、平面坐标）；
// 解决调用方需同时获取多个几何量的问题。
PlanetState planet_state_at_time(PlanetId id, double time_seconds_since_j2000) {
    const PlanetParams& planet = get_planet_params(id);
    const double varphi = planet_true_anomaly_at_time(id, time_seconds_since_j2000);
    const double theta_global = normalize_angle_0_2pi(planet.orbit.theta_0 + varphi);
    const double radius = planet_radius_at_true_anomaly(id, varphi);
    return PlanetState{
        id,
        time_seconds_since_j2000,
        varphi,
        theta_global,
        radius,
        radius * std::cos(theta_global),
        radius * std::sin(theta_global),
    };
}

}  // namespace spaceship_cpp::planet_params
