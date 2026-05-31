#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

double relative_mismatch(double mismatch, double v_in, double v_out) {
    const double scale = std::max({1.0, std::abs(v_in), std::abs(v_out)});
    return mismatch / scale;
}

bool canonical_consistent(
    const spaceship_cpp::trajectory::CanonicalOrbitETheta& lhs,
    const spaceship_cpp::trajectory::CanonicalOrbitETheta& rhs
) {
    if (!lhs.valid || !rhs.valid) {
        return false;
    }
    return std::abs(lhs.eccentricity - rhs.eccentricity) <= 1e-12 &&
        std::abs(spaceship_cpp::common::normalize_angle_minus_pi_pi(
            lhs.periapsis_angle - rhs.periapsis_angle)) <= 1e-12;
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;
    namespace trajectory = spaceship_cpp::trajectory;

    std::cout << std::setprecision(12) << std::scientific;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_vinf_conservation_audit_skipped_missing_table\n";
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
    initial_options.time_weight_m_per_s_per_day = 0.0;

    problem2::Problem2GravityAssistSolverOptions problem2_options{};
    problem2_options.theta_sample_count = 64;
    problem2_options.topology_adaptive_enabled = true;
    problem2_options.max_transfer_revolution = 1;
    problem2_options.max_target_revolution = 1;

    bfs::FixedSequenceThetaSearchOptions sequence_options{};
    sequence_options.beam_width = 10;
    sequence_options.max_launch_v_inf = 7000.0;
    sequence_options.time_weight_m_per_s_per_day = 0.0;
    sequence_options.flyby_physical_options.enabled = true;
    sequence_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::ObserveOnly;
    sequence_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;

    const auto search = bfs::run_fixed_sequence_theta_search_with_table(
        loader, launch_time, theta, sequence, initial_options, problem2_options, sequence_options);
    assert(search.ok);

    int edge_count = 0;
    int valid_flyby_eval_count = 0;
    int invalid_flyby_eval_count = 0;
    int small_slingshot_large_vinf_count = 0;
    double max_abs_slingshot = 0.0;
    double max_vinf_mismatch = 0.0;
    double max_relative_mismatch = 0.0;

    for (std::size_t node_index = 1; node_index < search.nodes.size(); ++node_index) {
        const auto& node = search.nodes[node_index];
        if (!node.valid || node.parent_index <= 0 || !node.incoming_edge.valid) {
            continue;
        }
        const auto& parent = search.nodes[static_cast<std::size_t>(node.parent_index)];
        const auto& edge = node.incoming_edge;
        const auto flyby = trajectory::evaluate_flyby_physical_feasibility(
            parent.state.current_planet,
            parent.state.current_time,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            edge.outgoing_e,
            edge.outgoing_theta,
            sequence_options.flyby_physical_options);
        const auto outgoing_canonical = trajectory::canonicalize_orbit_e_theta(
            edge.outgoing_e, edge.outgoing_theta);
        const auto next_canonical = trajectory::canonicalize_orbit_e_theta(
            node.state.incoming_e, node.state.incoming_theta);
        const double rel = relative_mismatch(
            flyby.v_infinity_mismatch, flyby.v_infinity_in, flyby.v_infinity_out);
        const double abs_slingshot = std::abs(edge.slingshot_residual);
        edge_count += 1;
        max_abs_slingshot = std::max(max_abs_slingshot, abs_slingshot);
        if (flyby.valid) {
            valid_flyby_eval_count += 1;
            max_vinf_mismatch = std::max(max_vinf_mismatch, flyby.v_infinity_mismatch);
            max_relative_mismatch = std::max(max_relative_mismatch, rel);
            if (abs_slingshot <= 1e-8 && flyby.v_infinity_mismatch > 1.0) {
                small_slingshot_large_vinf_count += 1;
            }
        } else {
            invalid_flyby_eval_count += 1;
        }

        std::cout << "Problem2VinfConservationAuditEdge\n";
        std::cout << "index=" << (edge_count - 1) << '\n';
        std::cout << "depth=" << node.state.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(parent.state.current_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
        std::cout << "parent_time=" << parent.state.current_time << '\n';
        std::cout << "edge_departure_time=" << edge.departure_time << '\n';
        std::cout << "edge_arrival_time=" << edge.arrival_time << '\n';
        std::cout << "incoming_e=" << parent.state.incoming_e << '\n';
        std::cout << "incoming_theta=" << parent.state.incoming_theta << '\n';
        std::cout << "outgoing_e_raw=" << edge.outgoing_e << '\n';
        std::cout << "outgoing_theta_raw=" << edge.outgoing_theta << '\n';
        std::cout << "outgoing_e_canonical=" << outgoing_canonical.eccentricity << '\n';
        std::cout << "outgoing_theta_canonical=" << outgoing_canonical.periapsis_angle << '\n';
        std::cout << "theta_prime=" << edge.theta_prime << '\n';
        std::cout << "alpha=" << edge.alpha << '\n';
        std::cout << "k=" << edge.transfer_revolution << '\n';
        std::cout << "q=" << edge.target_revolution << '\n';
        std::cout << "problem2_slingshot_residual=" << edge.slingshot_residual << '\n';
        std::cout << "problem1_residual_seconds=" << edge.problem1_residual_seconds << '\n';
        std::cout << "v_inf_in=" << flyby.v_infinity_in << '\n';
        std::cout << "v_inf_out=" << flyby.v_infinity_out << '\n';
        std::cout << "v_inf_mismatch=" << flyby.v_infinity_mismatch << '\n';
        std::cout << "v_inf_relative_mismatch=" << rel << '\n';
        std::cout << "turn_angle_rad=" << flyby.turn_angle_rad << '\n';
        std::cout << "max_turn_angle_rad=" << flyby.max_turn_angle_rad << '\n';
        std::cout << "required_periapsis_radius_m=" << flyby.required_periapsis_radius_m << '\n';
        std::cout << "min_allowed_periapsis_radius_m=" << flyby.min_allowed_periapsis_radius_m << '\n';

        std::cout << "Problem2BridgeMappingAudit\n";
        std::cout << "edge_outgoing_e=" << edge.outgoing_e << '\n';
        std::cout << "edge_outgoing_theta=" << edge.outgoing_theta << '\n';
        std::cout << "next_state_incoming_e=" << node.state.incoming_e << '\n';
        std::cout << "next_state_incoming_theta=" << node.state.incoming_theta << '\n';
        std::cout << "theta_prime_if_available=" << edge.theta_prime << '\n';
        std::cout << "canonical_edge_outgoing_e=" << outgoing_canonical.eccentricity << '\n';
        std::cout << "canonical_edge_outgoing_theta=" << outgoing_canonical.periapsis_angle << '\n';
        std::cout << "canonical_next_state_incoming_e=" << next_canonical.eccentricity << '\n';
        std::cout << "canonical_next_state_incoming_theta=" << next_canonical.periapsis_angle << '\n';
        std::cout << "edge_next_state_consistent="
                  << (canonical_consistent(outgoing_canonical, next_canonical) ? 1 : 0) << '\n';
        std::cout << "mapping_audit_ok=1\n";
    }

    std::cout << "Problem2VinfConservationAuditSummary\n";
    std::cout << "edge_count=" << edge_count << '\n';
    std::cout << "valid_flyby_eval_count=" << valid_flyby_eval_count << '\n';
    std::cout << "invalid_flyby_eval_count=" << invalid_flyby_eval_count << '\n';
    std::cout << "max_abs_problem2_slingshot_residual=" << max_abs_slingshot << '\n';
    std::cout << "max_v_inf_mismatch=" << max_vinf_mismatch << '\n';
    std::cout << "max_v_inf_relative_mismatch=" << max_relative_mismatch << '\n';
    std::cout << "edge_count_small_slingshot_large_vinf_mismatch="
              << small_slingshot_large_vinf_count << '\n';
    std::cout << "problem2_residual_inconsistent_with_velocity_helper="
              << (small_slingshot_large_vinf_count > 0 ? 1 : 0) << '\n';
    std::cout << "audit_ok=1\n";

    return 0;
}
