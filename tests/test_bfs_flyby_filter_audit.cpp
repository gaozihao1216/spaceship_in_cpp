#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"
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
#include <sstream>
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

long long sum_reject_invalid(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.reject_flyby_evaluator_invalid_count;
    }
    return total;
}

long long sum_reject_physical(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.reject_flyby_physical_infeasible_count;
    }
    return total;
}

long long sum_would_reject_invalid(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.would_reject_flyby_evaluator_invalid_count;
    }
    return total;
}

long long sum_would_reject_physical(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.would_reject_flyby_physical_infeasible_count;
    }
    return total;
}

long long sum_rejected_by_vinf_mismatch(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.rejected_by_vinf_mismatch_count;
    }
    return total;
}

long long sum_rejected_by_turn_angle(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.rejected_by_turn_angle_count;
    }
    return total;
}

long long sum_rejected_by_periapsis(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.rejected_by_periapsis_radius_count;
    }
    return total;
}

std::string sequence_string(const std::vector<spaceship_cpp::planet_params::PlanetId>& sequence) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::ostringstream os;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        if (i > 0) {
            os << "->";
        }
        os << planet_params::planet_name(sequence[i]);
    }
    return os.str();
}

std::vector<spaceship_cpp::planet_params::PlanetId> sequence_from_path(
    spaceship_cpp::planet_params::PlanetId launch_planet,
    const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& path
) {
    std::vector<spaceship_cpp::planet_params::PlanetId> sequence;
    sequence.push_back(launch_planet);
    for (const auto& edge : path) {
        sequence.push_back(edge.to_planet);
    }
    return sequence;
}

double mismatch_tolerance(
    const spaceship_cpp::trajectory::FlybyPhysicalFeasibilityResult& flyby,
    const spaceship_cpp::trajectory::FlybyPhysicalFeasibilityOptions& options
) {
    const double v_inf = 0.5 * (flyby.v_infinity_in + flyby.v_infinity_out);
    return options.relative_speed_tolerance * std::max(1.0, v_inf);
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
        std::cout << "bfs_flyby_filter_audit_skipped_missing_table\n";
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

    bfs::FixedSequenceThetaSearchOptions observe_options{};
    observe_options.beam_width = 10;
    observe_options.max_launch_v_inf = 7000.0;
    observe_options.time_weight_m_per_s_per_day = 0.0;
    observe_options.flyby_physical_options.enabled = true;
    observe_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::ObserveOnly;
    observe_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;

    auto enforce_options = observe_options;
    enforce_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::Enforce;

    const auto observe = bfs::run_fixed_sequence_theta_search_with_table(
        loader, launch_time, theta, sequence, initial_options, problem2_options, observe_options);
    const auto enforce = bfs::run_fixed_sequence_theta_search_with_table(
        loader, launch_time, theta, sequence, initial_options, problem2_options, enforce_options);
    assert(observe.ok);
    assert(enforce.ok);

    int printed = 0;
    for (std::size_t node_index = 1; node_index < observe.nodes.size() && printed < 50; ++node_index) {
        const auto& node = observe.nodes[node_index];
        if (!node.valid || node.parent_index <= 0 || !node.incoming_edge.valid) {
            continue;
        }
        const auto& parent = observe.nodes[static_cast<std::size_t>(node.parent_index)];
        const auto& edge = node.incoming_edge;
        const auto flyby = trajectory::evaluate_flyby_physical_feasibility(
            parent.state.current_planet,
            parent.state.current_time,
            parent.state.incoming_e,
            parent.state.incoming_theta,
            edge.outgoing_e,
            edge.outgoing_theta,
            observe_options.flyby_physical_options);
        const double tolerance = mismatch_tolerance(flyby, observe_options.flyby_physical_options);
        std::cout << "BfsFlybyAuditEdge\n";
        std::cout << "index=" << printed << '\n';
        std::cout << "depth=" << node.state.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(parent.state.current_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
        std::cout << "parent_time=" << parent.state.current_time << '\n';
        std::cout << "incoming_e=" << parent.state.incoming_e << '\n';
        std::cout << "incoming_theta=" << parent.state.incoming_theta << '\n';
        std::cout << "outgoing_e=" << edge.outgoing_e << '\n';
        std::cout << "outgoing_theta=" << edge.outgoing_theta << '\n';
        std::cout << "flyby_valid=" << (flyby.valid ? 1 : 0) << '\n';
        std::cout << "flyby_feasible=" << (flyby.feasible ? 1 : 0) << '\n';
        std::cout << "v_inf_in=" << flyby.v_infinity_in << '\n';
        std::cout << "v_inf_out=" << flyby.v_infinity_out << '\n';
        std::cout << "v_inf_mismatch=" << flyby.v_infinity_mismatch << '\n';
        std::cout << "v_inf_mismatch_tolerance=" << tolerance << '\n';
        std::cout << "turn_angle_rad=" << flyby.turn_angle_rad << '\n';
        std::cout << "max_turn_angle_rad=" << flyby.max_turn_angle_rad << '\n';
        std::cout << "turn_angle_margin_rad=" << (flyby.max_turn_angle_rad - flyby.turn_angle_rad) << '\n';
        std::cout << "required_periapsis_radius_m=" << flyby.required_periapsis_radius_m << '\n';
        std::cout << "min_allowed_periapsis_radius_m=" << flyby.min_allowed_periapsis_radius_m << '\n';
        std::cout << "periapsis_margin_m="
                  << (flyby.required_periapsis_radius_m - flyby.min_allowed_periapsis_radius_m) << '\n';
        std::cout << "rejected_by_vinf_mismatch=" << (flyby.rejected_by_vinf_mismatch ? 1 : 0) << '\n';
        std::cout << "rejected_by_turn_angle=" << (flyby.rejected_by_turn_angle ? 1 : 0) << '\n';
        std::cout << "rejected_by_periapsis_radius=" << (flyby.rejected_by_periapsis_radius ? 1 : 0) << '\n';
        if (!flyby.valid) {
            std::cout << "BfsFlybyAuditInvalidEdge\n";
            std::cout << "depth=" << node.state.depth << '\n';
            std::cout << "from_planet=" << planet_params::planet_name(parent.state.current_planet) << '\n';
            std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
            std::cout << "incoming_e=" << parent.state.incoming_e << '\n';
            std::cout << "incoming_theta=" << parent.state.incoming_theta << '\n';
            std::cout << "outgoing_e=" << edge.outgoing_e << '\n';
            std::cout << "outgoing_theta=" << edge.outgoing_theta << '\n';
            std::cout << "reason_hint=nonfinite_or_nonpositive_velocity_after_canonicalization\n";
        }
        printed += 1;
    }

    for (const auto& stats : observe.problem2_solve_branch_stats) {
        std::cout << "BfsFlybyAuditSolveSummary\n";
        std::cout << "depth=" << stats.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(stats.from_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(stats.to_planet) << '\n';
        std::cout << "raw_edge_count=" << stats.raw_edge_count << '\n';
        std::cout << "accepted_edge_count=" << stats.accepted_edge_count << '\n';
        std::cout << "flyby_evaluated_count=" << stats.flyby_evaluated_count << '\n';
        std::cout << "flyby_valid_count=" << stats.flyby_valid_count << '\n';
        std::cout << "flyby_feasible_count=" << stats.flyby_feasible_count << '\n';
        std::cout << "flyby_infeasible_count=" << stats.flyby_infeasible_count << '\n';
        std::cout << "flyby_invalid_count=" << stats.flyby_invalid_count << '\n';
        std::cout << "would_reject_flyby_evaluator_invalid_count="
                  << stats.would_reject_flyby_evaluator_invalid_count << '\n';
        std::cout << "would_reject_flyby_physical_infeasible_count="
                  << stats.would_reject_flyby_physical_infeasible_count << '\n';
        std::cout << "rejected_by_vinf_mismatch_count=" << stats.rejected_by_vinf_mismatch_count << '\n';
        std::cout << "rejected_by_turn_angle_count=" << stats.rejected_by_turn_angle_count << '\n';
        std::cout << "rejected_by_periapsis_radius_count=" << stats.rejected_by_periapsis_radius_count << '\n';
        std::cout << "min_v_inf_in=" << stats.min_v_inf_in << '\n';
        std::cout << "max_v_inf_in=" << stats.max_v_inf_in << '\n';
        std::cout << "min_v_inf_out=" << stats.min_v_inf_out << '\n';
        std::cout << "max_v_inf_out=" << stats.max_v_inf_out << '\n';
        std::cout << "min_v_inf_mismatch=" << stats.min_v_inf_mismatch << '\n';
        std::cout << "max_v_inf_mismatch=" << stats.max_v_inf_mismatch << '\n';
        std::cout << "min_turn_angle_rad=" << stats.min_flyby_turn_angle_rad << '\n';
        std::cout << "max_turn_angle_rad=" << stats.max_flyby_turn_angle_rad << '\n';
        std::cout << "min_max_turn_angle_rad=" << stats.min_flyby_max_turn_angle_rad << '\n';
        std::cout << "max_max_turn_angle_rad=" << stats.max_flyby_max_turn_angle_rad << '\n';
        std::cout << "min_required_periapsis_radius_m=" << stats.min_flyby_required_periapsis_radius_m << '\n';
        std::cout << "max_required_periapsis_radius_m=" << stats.max_flyby_required_periapsis_radius_m << '\n';
        std::cout << "min_allowed_periapsis_radius_m=" << stats.min_flyby_allowed_periapsis_radius_m << '\n';
        std::cout << "max_allowed_periapsis_radius_m=" << stats.max_flyby_allowed_periapsis_radius_m << '\n';
    }

    std::vector<planet_params::PlanetId> observe_best_sequence;
    if (observe.best_terminal_solution.valid) {
        observe_best_sequence = sequence_from_path(
            sequence.front(),
            bfs::reconstruct_fixed_launch_theta_path(observe, observe.best_terminal_solution.node_index));
    }

    std::cout << "BfsFlybyFilterModeComparison\n";
    std::cout << "mode_observe_best_valid=" << (observe.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "mode_observe_best_total_delta_v="
              << (observe.best_terminal_solution.valid
                      ? observe.best_terminal_solution.total_delta_v
                      : std::numeric_limits<double>::infinity()) << '\n';
    std::cout << "mode_observe_best_sequence=" << sequence_string(observe_best_sequence) << '\n';
    std::cout << "mode_observe_terminal_solution_count=" << observe.terminal_solutions.size() << '\n';
    std::cout << "mode_observe_would_reject_invalid_count=" << sum_would_reject_invalid(observe) << '\n';
    std::cout << "mode_observe_would_reject_physical_count=" << sum_would_reject_physical(observe) << '\n';
    std::cout << "mode_observe_rejected_by_vinf_mismatch_count="
              << sum_rejected_by_vinf_mismatch(observe) << '\n';
    std::cout << "mode_observe_rejected_by_turn_angle_count="
              << sum_rejected_by_turn_angle(observe) << '\n';
    std::cout << "mode_observe_rejected_by_periapsis_radius_count="
              << sum_rejected_by_periapsis(observe) << '\n';
    std::cout << "mode_enforce_best_valid=" << (enforce.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "mode_enforce_best_total_delta_v="
              << (enforce.best_terminal_solution.valid
                      ? enforce.best_terminal_solution.total_delta_v
                      : std::numeric_limits<double>::infinity()) << '\n';
    std::cout << "mode_enforce_terminal_solution_count=" << enforce.terminal_solutions.size() << '\n';
    std::cout << "mode_enforce_reject_invalid_count=" << sum_reject_invalid(enforce) << '\n';
    std::cout << "mode_enforce_reject_physical_count=" << sum_reject_physical(enforce) << '\n';
    std::cout << "comparison_ok=1\n";

    return 0;
}
