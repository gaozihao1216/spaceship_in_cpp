#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

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

long long total_raw_edges(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.raw_edge_count;
    }
    return total;
}

long long total_accepted_edges(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_solve_branch_stats) {
        total += stats.accepted_edge_count;
    }
    return total;
}

long long total_beam_pruned(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.problem2_layer_branch_stats) {
        total += stats.beam_pruned_count;
    }
    return total;
}

void print_solve_stats(
    int index,
    const spaceship_cpp::bfs::Problem2SolveBranchStats& stats
) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::cout << "BfsFixedSequenceProblem2SolveBranchStats\n";
    std::cout << "index=" << index << '\n';
    std::cout << "depth=" << stats.depth << '\n';
    std::cout << "parent_node_index=" << stats.parent_node_index << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(stats.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(stats.to_planet) << '\n';
    std::cout << "parent_time=" << stats.parent_time << '\n';
    std::cout << "parent_incoming_e=" << stats.parent_incoming_e << '\n';
    std::cout << "parent_incoming_theta=" << stats.parent_incoming_theta << '\n';
    std::cout << "raw_edge_count=" << stats.raw_edge_count << '\n';
    std::cout << "accepted_edge_count=" << stats.accepted_edge_count << '\n';
    std::cout << "reject_invalid_edge_count=" << stats.reject_invalid_edge_count << '\n';
    std::cout << "reject_invalid_next_state_count=" << stats.reject_invalid_next_state_count << '\n';
    std::cout << "reject_nonfinite_arrival_time_count=" << stats.reject_nonfinite_arrival_time_count << '\n';
    std::cout << "reject_nonfinite_incoming_orbit_count=" << stats.reject_nonfinite_incoming_orbit_count << '\n';
    std::cout << "reject_non_increasing_time_count=" << stats.reject_non_increasing_time_count << '\n';
    std::cout << "reject_non_positive_transfer_time_count=" << stats.reject_non_positive_transfer_time_count << '\n';
    std::cout << "reject_non_positive_outgoing_p_count=" << stats.reject_non_positive_outgoing_p_count << '\n';
    std::cout << "reject_slingshot_residual_count=" << stats.reject_slingshot_residual_count << '\n';
    std::cout << "reject_problem1_residual_count=" << stats.reject_problem1_residual_count << '\n';
    std::cout << "reject_arrival_time_count=" << stats.reject_arrival_time_count << '\n';
    std::cout << "reject_flyby_physical_infeasible_count="
              << stats.reject_flyby_physical_infeasible_count << '\n';
    std::cout << "min_transfer_time_seconds=" << stats.min_transfer_time_seconds << '\n';
    std::cout << "max_transfer_time_seconds=" << stats.max_transfer_time_seconds << '\n';
    std::cout << "min_abs_slingshot_residual=" << stats.min_abs_slingshot_residual << '\n';
    std::cout << "max_abs_slingshot_residual=" << stats.max_abs_slingshot_residual << '\n';
    std::cout << "min_abs_problem1_residual_seconds=" << stats.min_abs_problem1_residual_seconds << '\n';
    std::cout << "max_abs_problem1_residual_seconds=" << stats.max_abs_problem1_residual_seconds << '\n';
    std::cout << "min_flyby_required_periapsis_radius_m="
              << stats.min_flyby_required_periapsis_radius_m << '\n';
    std::cout << "max_flyby_required_periapsis_radius_m="
              << stats.max_flyby_required_periapsis_radius_m << '\n';
    std::cout << "min_flyby_turn_angle_rad=" << stats.min_flyby_turn_angle_rad << '\n';
    std::cout << "max_flyby_turn_angle_rad=" << stats.max_flyby_turn_angle_rad << '\n';
    std::cout << "min_flyby_max_turn_angle_rad=" << stats.min_flyby_max_turn_angle_rad << '\n';
    std::cout << "max_flyby_max_turn_angle_rad=" << stats.max_flyby_max_turn_angle_rad << '\n';
    std::cout << "solve_wall_time_ms=" << stats.solve_wall_time_ms << '\n';
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "bfs_fixed_sequence_problem2_branch_profile_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time_offset_seconds = -315583.195659;
    const double launch_time =
        0.17 * planet_params::planet_orbital_period(launch_planet) + launch_time_offset_seconds;
    const double theta = 2.572077523548;
    const std::vector<planet_params::PlanetId> sequence{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };

    bfs::InitialLaunchExpansionOptions initial_options{};
    initial_options.max_transfer_revolution = 1;
    initial_options.max_target_revolution = 1;
    initial_options.max_launch_v_inf = 7000.0;
    initial_options.time_weight_m_per_s_per_day = 0.0;

    problem2::Problem2GravityAssistSolverOptions problem2_options{};
    problem2_options.theta_sample_count = 64;
    problem2_options.topology_adaptive_enabled = true;

    auto run_with_beam = [&](int beam_width, bool flyby_filter_enabled) {
        bfs::FixedSequenceThetaSearchOptions sequence_options{};
        sequence_options.beam_width = beam_width;
        sequence_options.max_launch_v_inf = 7000.0;
        sequence_options.time_weight_m_per_s_per_day = 0.0;
        sequence_options.max_abs_slingshot_residual = 1e-8;
        sequence_options.max_abs_problem1_residual_seconds = 1e-2;
        sequence_options.flyby_physical_options.enabled = flyby_filter_enabled;
        return bfs::run_fixed_sequence_theta_search_with_table(
            loader,
            launch_time,
            theta,
            sequence,
            initial_options,
            problem2_options,
            sequence_options);
    };

    const auto result = run_with_beam(10, false);
    const auto flyby_enabled = run_with_beam(10, true);
    assert(result.ok);
    assert(flyby_enabled.ok);
    assert(!result.problem2_solve_branch_stats.empty());
    assert(!result.problem2_layer_branch_stats.empty());

    std::cout << "BfsFixedSequenceProblem2BranchProfileSummary\n";
    std::cout << "sequence=" << sequence_string(sequence) << '\n';
    std::cout << "launch_time_offset_seconds=" << launch_time_offset_seconds << '\n';
    std::cout << "theta=" << theta << '\n';
    std::cout << "initial_candidate_count=" << result.initial_candidate_count << '\n';
    std::cout << "terminal_solution_count=" << result.terminal_solutions.size() << '\n';
    std::cout << "best_valid=" << (result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "best_total_delta_v=" << result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_terminal_solution.launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_terminal_solution.arrival_v_inf << '\n';
    std::cout << "problem2_solve_count=" << result.problem2_solve_branch_stats.size() << '\n';
    std::cout << "problem2_layer_count=" << result.problem2_layer_branch_stats.size() << '\n';
    std::cout << "total_raw_edge_count=" << total_raw_edges(result) << '\n';
    std::cout << "total_accepted_edge_count=" << total_accepted_edges(result) << '\n';
    std::cout << "total_beam_pruned_count=" << total_beam_pruned(result) << '\n';
    std::cout << "profile_ok=1\n";

    for (const auto& stats : result.problem2_layer_branch_stats) {
        std::cout << "BfsFixedSequenceProblem2LayerBranchStats\n";
        std::cout << "depth=" << stats.depth << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(stats.from_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(stats.to_planet) << '\n';
        std::cout << "input_state_count=" << stats.input_state_count << '\n';
        std::cout << "problem2_solve_count=" << stats.problem2_solve_count << '\n';
        std::cout << "raw_edge_count_total=" << stats.raw_edge_count_total << '\n';
        std::cout << "accepted_edge_count_total=" << stats.accepted_edge_count_total << '\n';
        std::cout << "output_state_count_after_beam=" << stats.output_state_count_after_beam << '\n';
        std::cout << "beam_pruned_count=" << stats.beam_pruned_count << '\n';
        std::cout << "raw_edge_count_min=" << stats.raw_edge_count_min << '\n';
        std::cout << "raw_edge_count_mean=" << stats.raw_edge_count_mean << '\n';
        std::cout << "raw_edge_count_max=" << stats.raw_edge_count_max << '\n';
        std::cout << "accepted_edge_count_min=" << stats.accepted_edge_count_min << '\n';
        std::cout << "accepted_edge_count_mean=" << stats.accepted_edge_count_mean << '\n';
        std::cout << "accepted_edge_count_max=" << stats.accepted_edge_count_max << '\n';
        std::cout << "layer_wall_time_ms=" << stats.layer_wall_time_ms << '\n';
    }

    const int solve_print_count = std::min<int>(50, static_cast<int>(result.problem2_solve_branch_stats.size()));
    for (int i = 0; i < solve_print_count; ++i) {
        print_solve_stats(i, result.problem2_solve_branch_stats[static_cast<std::size_t>(i)]);
    }
    const int flyby_solve_print_count =
        std::min<int>(50, static_cast<int>(flyby_enabled.problem2_solve_branch_stats.size()));
    for (int i = 0; i < flyby_solve_print_count; ++i) {
        print_solve_stats(i, flyby_enabled.problem2_solve_branch_stats[static_cast<std::size_t>(i)]);
    }

    long long flyby_reject_count = 0;
    for (const auto& stats : flyby_enabled.problem2_solve_branch_stats) {
        flyby_reject_count += stats.reject_flyby_physical_infeasible_count;
    }

    std::cout << "BfsFixedSequenceProblem2FlybyFilterComparison\n";
    std::cout << "filter_disabled_best_valid=" << (result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "filter_disabled_best_total_delta_v=" << result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "filter_disabled_solve_count=" << result.problem2_solve_branch_stats.size() << '\n';
    std::cout << "filter_disabled_raw_edges=" << total_raw_edges(result) << '\n';
    std::cout << "filter_disabled_accepted_edges=" << total_accepted_edges(result) << '\n';
    std::cout << "filter_enabled_best_valid=" << (flyby_enabled.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "filter_enabled_best_total_delta_v=" << flyby_enabled.best_terminal_solution.total_delta_v << '\n';
    std::cout << "filter_enabled_solve_count=" << flyby_enabled.problem2_solve_branch_stats.size() << '\n';
    std::cout << "filter_enabled_raw_edges=" << total_raw_edges(flyby_enabled) << '\n';
    std::cout << "filter_enabled_accepted_edges=" << total_accepted_edges(flyby_enabled) << '\n';
    std::cout << "filter_enabled_reject_flyby_physical_count=" << flyby_reject_count << '\n';
    std::cout << "comparison_ok=1\n";

    const auto beam5 = run_with_beam(5, false);
    const auto beam3 = run_with_beam(3, false);
    assert(beam5.ok);
    assert(beam3.ok);
    std::cout << "BfsFixedSequenceProblem2BeamWidthComparison\n";
    std::cout << "beam_width_a=10\n";
    std::cout << "solve_count_a=" << result.problem2_solve_branch_stats.size() << '\n';
    std::cout << "best_total_delta_v_a=" << result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "beam_width_b=5\n";
    std::cout << "solve_count_b=" << beam5.problem2_solve_branch_stats.size() << '\n';
    std::cout << "best_total_delta_v_b=" << beam5.best_terminal_solution.total_delta_v << '\n';
    std::cout << "beam_width_c=3\n";
    std::cout << "solve_count_c=" << beam3.problem2_solve_branch_stats.size() << '\n';
    std::cout << "best_total_delta_v_c=" << beam3.best_terminal_solution.total_delta_v << '\n';
    std::cout << "comparison_ok=1\n";

    return 0;
}
