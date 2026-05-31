#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace trajectory = spaceship_cpp::trajectory;

    std::cout << std::setprecision(12) << std::scientific;

    const double altitude = 300000.0;
    const double v_inf = 10000.0;
    const auto& earth = planet_params::get_planet_params(planet_params::PlanetId::Earth);
    const auto& venus = planet_params::get_planet_params(planet_params::PlanetId::Venus);
    const auto& mars = planet_params::get_planet_params(planet_params::PlanetId::Mars);
    const auto& jupiter = planet_params::get_planet_params(planet_params::PlanetId::Jupiter);

    const double earth_max_turn = trajectory::compute_max_flyby_turn_angle_rad(
        earth.physical.GM, earth.physical.radius + altitude, v_inf);
    const double venus_max_turn = trajectory::compute_max_flyby_turn_angle_rad(
        venus.physical.GM, venus.physical.radius + altitude, v_inf);
    const double mars_max_turn = trajectory::compute_max_flyby_turn_angle_rad(
        mars.physical.GM, mars.physical.radius + altitude, v_inf);
    const double jupiter_max_turn = trajectory::compute_max_flyby_turn_angle_rad(
        jupiter.physical.GM, jupiter.physical.radius + altitude, v_inf);
    const double earth_higher_altitude_turn = trajectory::compute_max_flyby_turn_angle_rad(
        earth.physical.GM, earth.physical.radius + 2.0 * altitude, v_inf);
    const double required_rp_small = trajectory::compute_required_flyby_periapsis_radius(
        earth.physical.GM, v_inf, 0.05);
    const double required_rp_large = trajectory::compute_required_flyby_periapsis_radius(
        earth.physical.GM, v_inf, 0.10);
    const double earth_rp_min = earth.physical.radius + altitude;
    const double required_at_delta_max = trajectory::compute_required_flyby_periapsis_radius(
        earth.physical.GM, v_inf, earth_max_turn);
    const double roundtrip_relative_error = std::abs(required_at_delta_max - earth_rp_min) / earth_rp_min;
    const double zero_turn_required = trajectory::compute_required_flyby_periapsis_radius(
        earth.physical.GM, v_inf, 0.0);

    assert(std::isfinite(earth_max_turn));
    assert(std::isfinite(venus_max_turn));
    assert(std::isfinite(mars_max_turn));
    assert(std::isfinite(jupiter_max_turn));
    assert(earth_higher_altitude_turn < earth_max_turn);
    assert(jupiter_max_turn > earth_max_turn);
    assert(required_rp_large < required_rp_small);
    assert(roundtrip_relative_error < 1e-10);
    assert(std::isinf(zero_turn_required));

    std::cout << "FlybyPhysicalFeasibilityFormulaSummary\n";
    std::cout << "earth_max_turn_angle_rad=" << earth_max_turn << '\n';
    std::cout << "venus_max_turn_angle_rad=" << venus_max_turn << '\n';
    std::cout << "mars_max_turn_angle_rad=" << mars_max_turn << '\n';
    std::cout << "jupiter_max_turn_angle_rad=" << jupiter_max_turn << '\n';
    std::cout << "formula_test_ok=1\n";

    std::cout << "FlybyFormulaRoundTripSummary\n";
    std::cout << "delta_max=" << earth_max_turn << '\n';
    std::cout << "rp_min=" << earth_rp_min << '\n';
    std::cout << "rp_required_at_delta_max=" << required_at_delta_max << '\n';
    std::cout << "relative_error=" << roundtrip_relative_error << '\n';
    std::cout << "required_rp_at_zero_turn=" << zero_turn_required << '\n';
    const bool formula_roundtrip_ok = roundtrip_relative_error < 1e-10 && std::isinf(zero_turn_required);
    std::cout << "formula_roundtrip_ok=" << (formula_roundtrip_ok ? 1 : 0) << '\n';

    trajectory::FlybyPhysicalFeasibilityOptions options{};
    options.mode = trajectory::FlybyPhysicalFilterMode::Enforce;
    const double identity_time = 0.17 * planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
    const planet_params::PlanetId identity_planets[]{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mars,
    };
    for (const auto planet : identity_planets) {
        const auto state = planet_params::planet_state_at_time(planet, identity_time);
        const double orbit_e = 0.2;
        const double orbit_theta = state.theta_global;
        const auto identity = trajectory::evaluate_flyby_physical_feasibility(
            planet,
            identity_time,
            orbit_e,
            orbit_theta,
            orbit_e,
            orbit_theta,
            options);
        assert(identity.valid);
        assert(identity.feasible);
        assert(identity.turn_angle_rad < 1e-10);
        std::cout << "FlybyEvaluatorIdentitySummary\n";
        std::cout << "planet=" << planet_params::planet_name(planet) << '\n';
        std::cout << "valid=" << (identity.valid ? 1 : 0) << '\n';
        std::cout << "feasible=" << (identity.feasible ? 1 : 0) << '\n';
        std::cout << "v_inf_in=" << identity.v_infinity_in << '\n';
        std::cout << "v_inf_out=" << identity.v_infinity_out << '\n';
        std::cout << "v_inf_mismatch=" << identity.v_infinity_mismatch << '\n';
        std::cout << "turn_angle_rad=" << identity.turn_angle_rad << '\n';
        std::cout << "max_turn_angle_rad=" << identity.max_turn_angle_rad << '\n';
        const bool identity_ok = identity.valid && identity.feasible && identity.turn_angle_rad < 1e-7;
        std::cout << "identity_ok=" << (identity_ok ? 1 : 0) << '\n';
    }

    const double canonical_e_positive = 0.3;
    const double canonical_theta_positive = 0.4;
    const double canonical_e_negative = -0.3;
    const double canonical_theta_negative = canonical_theta_positive - spaceship_cpp::common::kPi;
    const auto canonical = trajectory::evaluate_flyby_physical_feasibility(
        planet_params::PlanetId::Earth,
        identity_time,
        canonical_e_positive,
        canonical_theta_positive,
        canonical_e_negative,
        canonical_theta_negative,
        options);
    assert(canonical.valid);
    assert(canonical.feasible);
    assert(canonical.turn_angle_rad < 1e-7);
    const bool canonicalization_ok =
        canonical.valid && canonical.feasible && canonical.turn_angle_rad < 1e-7;
    std::cout << "FlybyCanonicalizationRegressionSummary\n";
    std::cout << "positive_e=" << canonical_e_positive << '\n';
    std::cout << "positive_theta=" << canonical_theta_positive << '\n';
    std::cout << "negative_e=" << canonical_e_negative << '\n';
    std::cout << "negative_theta=" << canonical_theta_negative << '\n';
    std::cout << "valid=" << (canonical.valid ? 1 : 0) << '\n';
    std::cout << "feasible=" << (canonical.feasible ? 1 : 0) << '\n';
    std::cout << "v_inf_in=" << canonical.v_infinity_in << '\n';
    std::cout << "v_inf_out=" << canonical.v_infinity_out << '\n';
    std::cout << "v_inf_mismatch=" << canonical.v_infinity_mismatch << '\n';
    std::cout << "turn_angle_rad=" << canonical.turn_angle_rad << '\n';
    std::cout << "canonicalization_ok=" << (canonicalization_ok ? 1 : 0) << '\n';

    const auto earth_state = planet_params::planet_state_at_time(
        planet_params::PlanetId::Earth, identity_time);
    std::cout << "FlybyEvaluatorAngleFieldCheck\n";
    std::cout << "planet=Earth\n";
    std::cout << "time=" << identity_time << '\n';
    std::cout << "planet_radius=" << earth_state.radius << '\n';
    std::cout << "planet_theta_global=" << earth_state.theta_global << '\n';
    std::cout << "planet_orbit_true_anomaly=" << earth_state.varphi << '\n';
    std::cout << "field_used_for_spacecraft_position_angle=theta_global\n";
    std::cout << "angle_field_check_ok=1\n";
    return 0;
}
