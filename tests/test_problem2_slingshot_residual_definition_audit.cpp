#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct DirectSlingshotVelocityResidual {
    bool valid = false;
    double v_inf_in = std::numeric_limits<double>::quiet_NaN();
    double v_inf_out = std::numeric_limits<double>::quiet_NaN();
    double v_inf_mismatch = std::numeric_limits<double>::quiet_NaN();
    double squared_residual = std::numeric_limits<double>::quiet_NaN();
    double normalized_squared_residual = std::numeric_limits<double>::quiet_NaN();
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

double norm(spaceship_cpp::trajectory::CartesianVelocity velocity) {
    return std::hypot(velocity.vx, velocity.vy);
}

spaceship_cpp::trajectory::CartesianVelocity subtract(
    spaceship_cpp::trajectory::CartesianVelocity lhs,
    spaceship_cpp::trajectory::CartesianVelocity rhs
) {
    return spaceship_cpp::trajectory::CartesianVelocity{lhs.vx - rhs.vx, lhs.vy - rhs.vy};
}

DirectSlingshotVelocityResidual compute_direct_slingshot_velocity_residual(
    spaceship_cpp::planet_params::PlanetId flyby_planet,
    double flyby_time,
    double incoming_e,
    double incoming_theta,
    double outgoing_e,
    double outgoing_theta
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace trajectory = spaceship_cpp::trajectory;
    DirectSlingshotVelocityResidual result{};
    const auto incoming = trajectory::canonicalize_orbit_e_theta(incoming_e, incoming_theta);
    const auto outgoing = trajectory::canonicalize_orbit_e_theta(outgoing_e, outgoing_theta);
    if (!incoming.valid || !outgoing.valid || !std::isfinite(flyby_time)) {
        return result;
    }
    const auto& solar = planet_params::get_solar_system_physical_params();
    const auto state = planet_params::planet_state_at_time(flyby_planet, flyby_time);
    const auto planet_velocity = trajectory::compute_planet_heliocentric_velocity(flyby_planet, flyby_time);
    const auto v_in = trajectory::compute_heliocentric_velocity_on_orbit(
        solar.GM_sun, state.radius, state.theta_global, incoming.eccentricity, incoming.periapsis_angle);
    const auto v_out = trajectory::compute_heliocentric_velocity_on_orbit(
        solar.GM_sun, state.radius, state.theta_global, outgoing.eccentricity, outgoing.periapsis_angle);
    const auto u_in = subtract(v_in, planet_velocity);
    const auto u_out = subtract(v_out, planet_velocity);
    result.v_inf_in = norm(u_in);
    result.v_inf_out = norm(u_out);
    if (!std::isfinite(result.v_inf_in) || !std::isfinite(result.v_inf_out)) {
        return result;
    }
    const double in2 = result.v_inf_in * result.v_inf_in;
    const double out2 = result.v_inf_out * result.v_inf_out;
    result.v_inf_mismatch = std::abs(result.v_inf_out - result.v_inf_in);
    result.squared_residual = out2 - in2;
    result.normalized_squared_residual = result.squared_residual / std::max({1.0, in2, out2});
    result.valid = std::isfinite(result.normalized_squared_residual);
    return result;
}

double orbit_radius_from_elements(double p, double e, double phi, double theta) {
    const double denominator = 1.0 + e * std::cos(phi - theta);
    if (!std::isfinite(denominator) || denominator == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return p / denominator;
}

std::string best_interpretation_name(int index) {
    switch (index) {
    case 0:
        return "outgoing_theta";
    case 1:
        return "outgoing_theta_plus_pi";
    case 2:
        return "theta_prime";
    case 3:
        return "theta_prime_plus_pi";
    case 4:
        return "canonical_outgoing_theta";
    case 5:
        return "negative_e_equivalent";
    case 6:
        return "theta_prime_plus_planet_theta0";
    default:
        return "unknown";
    }
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;
    namespace trajectory = spaceship_cpp::trajectory;

    std::cout << std::setprecision(12) << std::scientific;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_slingshot_residual_definition_audit_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const std::vector<planet_params::PlanetId> sequence{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };
    const double launch_time =
        0.17 * planet_params::planet_orbital_period(planet_params::PlanetId::Earth) - 315583.195659;
    const double theta = 2.572077523548;

    bfs::InitialLaunchExpansionOptions initial_options{};
    initial_options.max_transfer_revolution = 1;
    initial_options.max_target_revolution = 1;
    initial_options.max_launch_v_inf = 7000.0;

    problem2::Problem2GravityAssistSolverOptions problem2_options{};
    problem2_options.theta_sample_count = 64;
    problem2_options.topology_adaptive_enabled = true;
    problem2_options.max_transfer_revolution = 1;
    problem2_options.max_target_revolution = 1;

    bfs::FixedSequenceThetaSearchOptions sequence_options{};
    sequence_options.beam_width = 10;
    sequence_options.max_launch_v_inf = 7000.0;
    sequence_options.flyby_physical_options.enabled = true;
    sequence_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::ObserveOnly;

    const auto search = bfs::run_fixed_sequence_theta_search_with_table(
        loader, launch_time, theta, sequence, initial_options, problem2_options, sequence_options);
    assert(search.ok);

    int edge_count = 0;
    int small_slingshot_large_direct_count = 0;
    int fixed_by_theta_plus_pi_count = 0;
    int geometry_radius_bad_count = 0;
    double max_edge_slingshot = 0.0;
    double max_direct_mismatch = 0.0;
    double max_direct_normalized_squared = 0.0;
    double min_best_interpretation_mismatch = std::numeric_limits<double>::infinity();
    double max_best_interpretation_mismatch = 0.0;
    double max_geometry_relative_error = 0.0;
    double max_analytic_direct_difference = 0.0;

    for (std::size_t node_index = 1; node_index < search.nodes.size(); ++node_index) {
        const auto& node = search.nodes[node_index];
        if (!node.valid || node.parent_index <= 0 || !node.incoming_edge.valid) {
            continue;
        }
        const auto& parent = search.nodes[static_cast<std::size_t>(node.parent_index)];
        const auto& edge = node.incoming_edge;
        const auto planet_state = planet_params::planet_state_at_time(
            parent.state.current_planet, parent.state.current_time);
        const auto& planet = planet_params::get_planet_params(parent.state.current_planet);

        const auto direct = compute_direct_slingshot_velocity_residual(
            parent.state.current_planet,
            parent.state.current_time,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            edge.outgoing_e,
            edge.outgoing_theta);
        const auto canonical = trajectory::canonicalize_orbit_e_theta(edge.outgoing_e, edge.outgoing_theta);
        const auto analytic = problem2::evaluate_problem2_slingshot_residual(
            planet_state.varphi,
            planet.orbit.e,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            edge.outgoing_e,
            edge.outgoing_theta);

        const double alternatives_e[7] = {
            edge.outgoing_e,
            edge.outgoing_e,
            edge.outgoing_e,
            edge.outgoing_e,
            canonical.eccentricity,
            -canonical.eccentricity,
            edge.outgoing_e,
        };
        const double alternatives_theta[7] = {
            edge.outgoing_theta,
            edge.outgoing_theta + common::kPi,
            edge.theta_prime,
            edge.theta_prime + common::kPi,
            canonical.periapsis_angle,
            canonical.periapsis_angle + common::kPi,
            edge.theta_prime + planet.orbit.theta_0,
        };
        double mismatches[7] = {};
        int best_index = 0;
        double best_mismatch = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 7; ++i) {
            const auto trial = compute_direct_slingshot_velocity_residual(
                parent.state.current_planet,
                parent.state.current_time,
                parent.state.incoming_e,
                parent.state.incoming_theta,
                alternatives_e[i],
                alternatives_theta[i]);
            mismatches[i] = trial.valid ? trial.v_inf_mismatch : std::numeric_limits<double>::infinity();
            if (mismatches[i] < best_mismatch) {
                best_mismatch = mismatches[i];
                best_index = i;
            }
        }

        const double r_reconstructed = orbit_radius_from_elements(
            edge.outgoing_p, edge.outgoing_e, planet_state.theta_global, edge.outgoing_theta);
        const double radius_error = r_reconstructed - planet_state.radius;
        const double radius_relative_error =
            std::abs(radius_error) / std::max(1.0, planet_state.radius);
        const double analytic_direct_difference = analytic.valid
            ? analytic.residual - direct.normalized_squared_residual
            : std::numeric_limits<double>::quiet_NaN();

        edge_count += 1;
        max_edge_slingshot = std::max(max_edge_slingshot, std::abs(edge.slingshot_residual));
        if (direct.valid) {
            max_direct_mismatch = std::max(max_direct_mismatch, direct.v_inf_mismatch);
            max_direct_normalized_squared =
                std::max(max_direct_normalized_squared, std::abs(direct.normalized_squared_residual));
            if (std::abs(edge.slingshot_residual) <= 1e-8 && direct.v_inf_mismatch > 1.0) {
                small_slingshot_large_direct_count += 1;
            }
        }
        min_best_interpretation_mismatch = std::min(min_best_interpretation_mismatch, best_mismatch);
        max_best_interpretation_mismatch = std::max(max_best_interpretation_mismatch, best_mismatch);
        if (best_index == 1 || best_index == 3) {
            fixed_by_theta_plus_pi_count += 1;
        }
        if (!std::isfinite(radius_relative_error) || radius_relative_error > 1e-8) {
            geometry_radius_bad_count += 1;
        }
        if (std::isfinite(analytic_direct_difference)) {
            max_analytic_direct_difference =
                std::max(max_analytic_direct_difference, std::abs(analytic_direct_difference));
        }
        max_geometry_relative_error = std::max(max_geometry_relative_error, radius_relative_error);

        std::cout << "Problem2SlingshotResidualDefinitionAuditEdge\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "depth=" << node.state.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(parent.state.current_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
        std::cout << "parent_time=" << parent.state.current_time << '\n';
        std::cout << "incoming_e=" << parent.state.incoming_e << '\n';
        std::cout << "incoming_theta=" << parent.state.incoming_theta << '\n';
        std::cout << "outgoing_e=" << edge.outgoing_e << '\n';
        std::cout << "outgoing_theta=" << edge.outgoing_theta << '\n';
        std::cout << "theta_prime=" << edge.theta_prime << '\n';
        std::cout << "alpha=" << edge.alpha << '\n';
        std::cout << "k=" << edge.transfer_revolution << '\n';
        std::cout << "q=" << edge.target_revolution << '\n';
        std::cout << "edge_slingshot_residual=" << edge.slingshot_residual << '\n';
        std::cout << "direct_v_inf_in=" << direct.v_inf_in << '\n';
        std::cout << "direct_v_inf_out=" << direct.v_inf_out << '\n';
        std::cout << "direct_v_inf_mismatch=" << direct.v_inf_mismatch << '\n';
        std::cout << "direct_squared_residual=" << direct.squared_residual << '\n';
        std::cout << "direct_normalized_squared_residual=" << direct.normalized_squared_residual << '\n';
        std::cout << "problem1_residual_seconds=" << edge.problem1_residual_seconds << '\n';
        std::cout << "outgoing_p=" << edge.outgoing_p << '\n';
        std::cout << "outgoing_semi_latus_rectum=" << edge.outgoing_p << '\n';
        std::cout << "edge_outgoing_theta_equals_theta_prime="
                  << (std::abs(common::normalize_angle_minus_pi_pi(edge.outgoing_theta - edge.theta_prime)) < 1e-12
                          ? 1 : 0) << '\n';
        std::cout << "edge_outgoing_theta_minus_theta_prime="
                  << common::normalize_angle_minus_pi_pi(edge.outgoing_theta - edge.theta_prime) << '\n';

        std::cout << "Problem2ThetaPrimeInterpretationAudit\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "mismatch_outgoing_theta=" << mismatches[0] << '\n';
        std::cout << "mismatch_outgoing_theta_plus_pi=" << mismatches[1] << '\n';
        std::cout << "mismatch_theta_prime=" << mismatches[2] << '\n';
        std::cout << "mismatch_theta_prime_plus_pi=" << mismatches[3] << '\n';
        std::cout << "mismatch_canonical_outgoing_theta=" << mismatches[4] << '\n';
        std::cout << "mismatch_negative_e_equivalent=" << mismatches[5] << '\n';
        std::cout << "mismatch_theta_prime_plus_planet_theta0=" << mismatches[6] << '\n';
        std::cout << "best_interpretation=" << best_interpretation_name(best_index) << '\n';
        std::cout << "best_mismatch=" << best_mismatch << '\n';

        std::cout << "Problem2AnalyticResidualAudit\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "analytic_residual_from_problem2=" << (analytic.valid ? analytic.residual : std::numeric_limits<double>::quiet_NaN()) << '\n';
        std::cout << "edge_slingshot_residual=" << edge.slingshot_residual << '\n';
        std::cout << "direct_normalized_squared_residual=" << direct.normalized_squared_residual << '\n';
        std::cout << "analytic_direct_difference=" << analytic_direct_difference << '\n';

        std::cout << "Problem2OutgoingOrbitGeometryAudit\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "r_flyby=" << planet_state.radius << '\n';
        std::cout << "phi_flyby=" << planet_state.theta_global << '\n';
        std::cout << "outgoing_p=" << edge.outgoing_p << '\n';
        std::cout << "outgoing_e=" << edge.outgoing_e << '\n';
        std::cout << "outgoing_theta=" << edge.outgoing_theta << '\n';
        std::cout << "reconstructed_r_from_outgoing_orbit=" << r_reconstructed << '\n';
        std::cout << "geometry_radius_error=" << radius_error << '\n';
        std::cout << "geometry_radius_relative_error=" << radius_relative_error << '\n';
    }

    std::cout << "Problem2SlingshotResidualDefinitionAuditSummary\n";
    std::cout << "edge_count=" << edge_count << '\n';
    std::cout << "max_edge_slingshot_residual=" << max_edge_slingshot << '\n';
    std::cout << "max_direct_v_inf_mismatch=" << max_direct_mismatch << '\n';
    std::cout << "max_direct_normalized_squared_residual=" << max_direct_normalized_squared << '\n';
    std::cout << "min_best_theta_interpretation_mismatch=" << min_best_interpretation_mismatch << '\n';
    std::cout << "max_best_theta_interpretation_mismatch=" << max_best_interpretation_mismatch << '\n';
    std::cout << "edge_count_fixed_by_theta_plus_pi=" << fixed_by_theta_plus_pi_count << '\n';
    std::cout << "edge_count_geometry_radius_bad=" << geometry_radius_bad_count << '\n';
    std::cout << "max_geometry_radius_relative_error=" << max_geometry_relative_error << '\n';
    std::cout << "max_analytic_direct_difference=" << max_analytic_direct_difference << '\n';
    std::cout << "problem2_residual_matches_direct_velocity="
              << (small_slingshot_large_direct_count == 0 ? 1 : 0) << '\n';
    std::cout << "theta_interpretation_suspect="
              << (max_best_interpretation_mismatch < 1.0 && max_direct_mismatch > 1.0 ? 1 : 0) << '\n';
    std::cout << "outgoing_geometry_suspect=" << (geometry_radius_bad_count > 0 ? 1 : 0) << '\n';
    std::cout << "small_slingshot_large_direct_count=" << small_slingshot_large_direct_count << '\n';
    std::cout << "audit_ok=1\n";

    return 0;
}
