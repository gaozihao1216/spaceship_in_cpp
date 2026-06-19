/*
 * 文件作用：声明行星参数、行星编号和状态查询接口。
 * 主要工作：保存太阳系物理常量、行星轨道参数，并计算指定时刻的行星位置。
 */
#pragma once

#include <array>

namespace spaceship_cpp::planet_params {

enum class PlanetId {
    // 八大行星按太阳系从内到外排序，底层整数值也按此顺序使用。
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
};

struct PlanetOrbitParams {
    // 半通径 p 和偏心率 e 描述轨道形状。
    double p;
    double e;
    // theta_0 是轨道近日点全局方向，varphi_0 是 J2000 时刻真近点角。
    double theta_0;
    double varphi_0;
};

struct PlanetPhysicalParams {
    // 行星半径和标准引力参数，用于飞掠可行性计算。
    double radius;
    double GM;
};

struct PlanetParams {
    // 单颗行星的静态参数集合。
    PlanetId id;
    const char* name;
    PlanetOrbitParams orbit;
    PlanetPhysicalParams physical;
};

struct PlanetState {
    // 指定时刻下行星在太阳中心平面模型中的状态。
    PlanetId id;
    double time_seconds_since_j2000;
    double varphi;
    double theta_global;
    double radius;
    double x;
    double y;
};

struct SolarSystemPhysicalParams {
    // 太阳标准引力参数，用于所有日心轨道速度和时间尺度计算。
    double GM_sun;
};

constexpr const char* kPlanetParamsEpochName = "J2000";
// 工程标签：用于说明时间零点是 J2000，不试图区分精密天文 TT/TDB 时间尺度。
constexpr const char* kPlanetParamsEpochIso = "2000-01-01T12:00:00Z";
constexpr double kPlanetParamsEpochJulianDate = 2451545.0;

// 按 PlanetId 取静态行星参数，非法 id 会抛出异常。
const PlanetParams& get_planet_params(PlanetId id);

// 取太阳系全局物理参数。
const SolarSystemPhysicalParams& get_solar_system_physical_params();

// 返回完整行星参数表，便于遍历测试和诊断。
const std::array<PlanetParams, 8>& all_planet_params();

// 行星 id 到名称/整数值的辅助转换。
const char* planet_name(PlanetId id);

bool is_valid_planet_id(PlanetId id);

int planet_id_raw_value(PlanetId id);

// 平均角速度和轨道周期，解决按时间推进真近点角的问题。
double planet_mean_motion(PlanetId id);

double planet_orbital_period(PlanetId id);

// 计算指定时刻的行星真近点角。
double planet_true_anomaly_at_time(PlanetId id, double time_seconds_since_j2000);

// 根据真近点角计算行星日心半径。
double planet_radius_at_true_anomaly(PlanetId id, double varphi);

// 一站式计算指定时刻的行星半径、局部角、全局角和平面坐标。
PlanetState planet_state_at_time(PlanetId id, double time_seconds_since_j2000);

}  // namespace spaceship_cpp::planet_params
