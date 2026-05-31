#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <limits>

namespace spaceship_cpp::trajectory {

enum class FlybyPhysicalFilterMode {
    Disabled,
    ObserveOnly,
    Enforce
};

struct CanonicalOrbitETheta {
    bool valid = false;
    double eccentricity = 0.0;
    double periapsis_angle = 0.0;
};

struct FlybyPhysicalFeasibilityOptions {
    bool enabled = true;
    FlybyPhysicalFilterMode mode = FlybyPhysicalFilterMode::Enforce;

    double min_flyby_altitude_m = 300000.0;

    double turn_angle_tolerance_rad = 1e-10;
    double relative_speed_tolerance = 1e-6;
};

struct FlybyPhysicalFeasibilityResult {
    bool valid = false;

    bool feasible = false;

    planet_params::PlanetId planet;

    double v_infinity_in = 0.0;
    double v_infinity_out = 0.0;
    double v_infinity_mismatch = 0.0;

    double turn_angle_rad = 0.0;
    double max_turn_angle_rad = 0.0;

    double required_periapsis_radius_m = std::numeric_limits<double>::infinity();
    double min_allowed_periapsis_radius_m = std::numeric_limits<double>::infinity();

    bool rejected_by_vinf_mismatch = false;
    bool rejected_by_turn_angle = false;
    bool rejected_by_periapsis_radius = false;
};

CanonicalOrbitETheta canonicalize_orbit_e_theta(double eccentricity, double periapsis_angle);

double compute_max_flyby_turn_angle_rad(
    double planet_mu,
    double min_periapsis_radius,
    double v_inf
);

double compute_required_flyby_periapsis_radius(
    double planet_mu,
    double v_inf,
    double turn_angle_rad
);

FlybyPhysicalFeasibilityResult evaluate_flyby_physical_feasibility(
    planet_params::PlanetId flyby_planet,
    double flyby_time,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta,
    const FlybyPhysicalFeasibilityOptions& options
);

}  // namespace spaceship_cpp::trajectory
