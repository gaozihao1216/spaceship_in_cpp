#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"
#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_expansion.hpp"
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
};

struct GeometryResidual {
    bool valid = false;
    double reconstructed_radius = std::numeric_limits<double>::quiet_NaN();
    double radius_error = std::numeric_limits<double>::quiet_NaN();
    double radius_relative_error = std::numeric_limits<double>::quiet_NaN();
};

struct RequeryBest {
    bool valid = false;
    spaceship_cpp::problem2::Problem2GravityAssistExpandedState state;
    double theta_prime_global = std::numeric_limits<double>::quiet_NaN();
    DirectSlingshotVelocityResidual direct;
    GeometryResidual geometry;
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
    result.v_inf_mismatch = std::abs(result.v_inf_out - result.v_inf_in);
    result.valid = std::isfinite(result.v_inf_mismatch);
    return result;
}

GeometryResidual compute_geometry_residual(
    spaceship_cpp::planet_params::PlanetId flyby_planet,
    double flyby_time,
    double outgoing_p,
    double outgoing_e,
    double outgoing_theta_global
) {
    GeometryResidual result{};
    const auto planet_state = spaceship_cpp::planet_params::planet_state_at_time(flyby_planet, flyby_time);
    const double denominator = 1.0 + outgoing_e * std::cos(planet_state.theta_global - outgoing_theta_global);
    if (!std::isfinite(denominator) || denominator == 0.0 || !std::isfinite(outgoing_p)) {
        return result;
    }
    result.reconstructed_radius = outgoing_p / denominator;
    result.radius_error = result.reconstructed_radius - planet_state.radius;
    result.radius_relative_error = std::abs(result.radius_error) / std::max(1.0, planet_state.radius);
    result.valid = std::isfinite(result.radius_relative_error);
    return result;
}

RequeryBest find_best_requery_candidate(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::bfs::FixedLaunchThetaSearchNode& parent,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverOptions& problem2_options
) {
    namespace bfs = spaceship_cpp::bfs;
    namespace problem2 = spaceship_cpp::problem2;
    RequeryBest best{};

    problem2::Problem2GravityAssistExpansionInput input{};
    input.encounter_planet = parent.state.current_planet;
    input.target_planet = target_planet;
    input.encounter_time = parent.state.current_time;
    input.incoming_e = parent.state.incoming_e;
    input.incoming_theta = bfs::global_periapsis_angle_to_problem2_local(
        parent.state.current_planet, parent.state.incoming_theta);

    const auto states = problem2::expand_one_gravity_assist_step_with_table(loader, input, problem2_options);
    for (const auto& state : states) {
        const double theta_global = bfs::problem2_local_periapsis_angle_to_global(
            parent.state.current_planet, state.theta_prime);
        const auto direct = compute_direct_slingshot_velocity_residual(
            parent.state.current_planet,
            parent.state.current_time,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            state.outgoing_e,
            theta_global);
        const auto geometry = compute_geometry_residual(
            parent.state.current_planet,
            parent.state.current_time,
            state.outgoing_p,
            state.outgoing_e,
            theta_global);
        if (!direct.valid || !geometry.valid) {
            continue;
        }
        if (!best.valid || direct.v_inf_mismatch < best.direct.v_inf_mismatch) {
            best.valid = true;
            best.state = state;
            best.theta_prime_global = theta_global;
            best.direct = direct;
            best.geometry = geometry;
        }
    }
    return best;
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
        std::cout << "problem2_angle_frame_adapter_audit_skipped_missing_table\n";
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
    int output_converted_geometry_good_count = 0;
    int requery_geometry_good_count = 0;
    int requery_vinf_good_count = 0;
    double max_geometry_error_current = 0.0;
    double max_geometry_error_output_converted = 0.0;
    double max_geometry_error_requery = 0.0;
    double max_vinf_mismatch_current = 0.0;
    double max_vinf_mismatch_output_converted = 0.0;
    double max_vinf_mismatch_requery = 0.0;
    double min_vinf_mismatch_requery = std::numeric_limits<double>::infinity();

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

        const double incoming_theta_local = bfs::global_periapsis_angle_to_problem2_local(
            parent.state.current_planet, parent.state.incoming_theta);
        const double output_converted_theta_global = bfs::problem2_local_periapsis_angle_to_global(
            parent.state.current_planet, edge.theta_prime);

        const auto direct_current = compute_direct_slingshot_velocity_residual(
            parent.state.current_planet,
            parent.state.current_time,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            edge.outgoing_e,
            edge.outgoing_theta);
        const auto direct_output_converted = compute_direct_slingshot_velocity_residual(
            parent.state.current_planet,
            parent.state.current_time,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            edge.outgoing_e,
            output_converted_theta_global);
        const auto geometry_current = compute_geometry_residual(
            parent.state.current_planet,
            parent.state.current_time,
            edge.outgoing_p,
            edge.outgoing_e,
            edge.outgoing_theta);
        const auto geometry_output_converted = compute_geometry_residual(
            parent.state.current_planet,
            parent.state.current_time,
            edge.outgoing_p,
            edge.outgoing_e,
            output_converted_theta_global);
        const auto requery = find_best_requery_candidate(loader, parent, edge.to_planet, problem2_options);

        edge_count += 1;
        if (geometry_current.valid) {
            max_geometry_error_current =
                std::max(max_geometry_error_current, geometry_current.radius_relative_error);
        }
        if (geometry_output_converted.valid) {
            max_geometry_error_output_converted =
                std::max(max_geometry_error_output_converted, geometry_output_converted.radius_relative_error);
            if (geometry_output_converted.radius_relative_error <= 1e-8) {
                output_converted_geometry_good_count += 1;
            }
        }
        if (requery.valid) {
            max_geometry_error_requery =
                std::max(max_geometry_error_requery, requery.geometry.radius_relative_error);
            max_vinf_mismatch_requery = std::max(max_vinf_mismatch_requery, requery.direct.v_inf_mismatch);
            min_vinf_mismatch_requery = std::min(min_vinf_mismatch_requery, requery.direct.v_inf_mismatch);
            if (requery.geometry.radius_relative_error <= 1e-8) {
                requery_geometry_good_count += 1;
            }
            if (requery.direct.v_inf_mismatch <= 1.0) {
                requery_vinf_good_count += 1;
            }
        }
        if (direct_current.valid) {
            max_vinf_mismatch_current = std::max(max_vinf_mismatch_current, direct_current.v_inf_mismatch);
        }
        if (direct_output_converted.valid) {
            max_vinf_mismatch_output_converted =
                std::max(max_vinf_mismatch_output_converted, direct_output_converted.v_inf_mismatch);
        }

        std::cout << "Problem2AngleFrameAdapterAuditEdge\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "depth=" << node.state.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(parent.state.current_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
        std::cout << "planet_theta0=" << planet.orbit.theta_0 << '\n';
        std::cout << "planet_varphi=" << planet_state.varphi << '\n';
        std::cout << "planet_theta_global=" << planet_state.theta_global << '\n';
        std::cout << "theta_global_minus_varphi="
                  << common::normalize_angle_0_2pi(planet_state.theta_global - planet_state.varphi) << '\n';
        std::cout << "incoming_theta_global=" << parent.state.incoming_theta << '\n';
        std::cout << "incoming_theta_local_candidate=" << incoming_theta_local << '\n';
        std::cout << "edge_theta_prime_raw=" << edge.theta_prime << '\n';
        std::cout << "outgoing_theta_current_global=" << edge.outgoing_theta << '\n';
        std::cout << "outgoing_theta_output_converted_global=" << output_converted_theta_global << '\n';
        std::cout << "geometry_error_current=" << geometry_current.radius_error << '\n';
        std::cout << "geometry_error_output_converted=" << geometry_output_converted.radius_error << '\n';
        std::cout << "vinf_mismatch_current=" << direct_current.v_inf_mismatch << '\n';
        std::cout << "vinf_mismatch_output_converted=" << direct_output_converted.v_inf_mismatch << '\n';

        std::cout << "Problem2AngleFrameAdapterRequeryAuditEdge\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "depth=" << node.state.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(parent.state.current_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
        std::cout << "incoming_theta_local_used=" << incoming_theta_local << '\n';
        std::cout << "requery_valid=" << (requery.valid ? 1 : 0) << '\n';
        std::cout << "theta_prime_local=" << (requery.valid ? requery.state.theta_prime : std::numeric_limits<double>::quiet_NaN()) << '\n';
        std::cout << "theta_prime_global_converted=" << requery.theta_prime_global << '\n';
        std::cout << "geometry_error_requery_converted=" << requery.geometry.radius_error << '\n';
        std::cout << "geometry_relative_error_requery_converted=" << requery.geometry.radius_relative_error << '\n';
        std::cout << "vinf_mismatch_requery_converted=" << requery.direct.v_inf_mismatch << '\n';
        std::cout << "slingshot_residual_requery="
                  << (requery.valid ? requery.state.slingshot_residual : std::numeric_limits<double>::quiet_NaN()) << '\n';
        std::cout << "problem1_residual_requery="
                  << (requery.valid ? requery.state.problem1_residual_seconds : std::numeric_limits<double>::quiet_NaN()) << '\n';
    }

    std::cout << "Problem2AngleFrameAdapterAuditSummary\n";
    std::cout << "edge_count=" << edge_count << '\n';
    std::cout << "output_converted_geometry_good_count=" << output_converted_geometry_good_count << '\n';
    std::cout << "requery_geometry_good_count=" << requery_geometry_good_count << '\n';
    std::cout << "requery_vinf_good_count=" << requery_vinf_good_count << '\n';
    std::cout << "max_geometry_error_current=" << max_geometry_error_current << '\n';
    std::cout << "max_geometry_error_output_converted=" << max_geometry_error_output_converted << '\n';
    std::cout << "max_geometry_error_requery=" << max_geometry_error_requery << '\n';
    std::cout << "max_vinf_mismatch_current=" << max_vinf_mismatch_current << '\n';
    std::cout << "max_vinf_mismatch_output_converted=" << max_vinf_mismatch_output_converted << '\n';
    std::cout << "min_vinf_mismatch_requery=" << min_vinf_mismatch_requery << '\n';
    std::cout << "max_vinf_mismatch_requery=" << max_vinf_mismatch_requery << '\n';
    std::cout << "requery_conversion_fixes_geometry="
              << (requery_geometry_good_count == edge_count && edge_count > 0 ? 1 : 0) << '\n';
    std::cout << "requery_conversion_fixes_vinf="
              << (requery_vinf_good_count == edge_count && edge_count > 0 ? 1 : 0) << '\n';
    std::cout << "audit_ok=1\n";

    return 0;
}
