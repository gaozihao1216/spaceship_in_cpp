#pragma once

#include <array>

namespace spaceship_cpp::planet_params {

enum class PlanetId {
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
    double p;
    double e;
    double theta_0;
    double varphi_0;
};

struct PlanetPhysicalParams {
    double radius;
    double GM;
};

struct PlanetParams {
    PlanetId id;
    const char* name;
    PlanetOrbitParams orbit;
    PlanetPhysicalParams physical;
};

struct PlanetState {
    PlanetId id;
    double time_seconds_since_j2000;
    double varphi;
    double theta_global;
    double radius;
    double x;
    double y;
};

struct SolarSystemPhysicalParams {
    double GM_sun;
};

constexpr const char* kPlanetParamsEpochName = "J2000";
// Engineering label for the J2000 epoch. This is not intended to model the
// exact TT/TDB astronomical timescale distinction.
constexpr const char* kPlanetParamsEpochIso = "2000-01-01T12:00:00Z";
constexpr double kPlanetParamsEpochJulianDate = 2451545.0;

const PlanetParams& get_planet_params(PlanetId id);

const SolarSystemPhysicalParams& get_solar_system_physical_params();

const std::array<PlanetParams, 8>& all_planet_params();

const char* planet_name(PlanetId id);

bool is_valid_planet_id(PlanetId id);

int planet_id_raw_value(PlanetId id);

double planet_mean_motion(PlanetId id);

double planet_orbital_period(PlanetId id);

double planet_true_anomaly_at_time(PlanetId id, double time_seconds_since_j2000);

double planet_radius_at_true_anomaly(PlanetId id, double varphi);

PlanetState planet_state_at_time(PlanetId id, double time_seconds_since_j2000);

}  // namespace spaceship_cpp::planet_params
