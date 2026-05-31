#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

int open_depth_from_env() {
    if (const char* raw = std::getenv("BFS_OPEN_DEPTH")) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return 3;
}

int open_beam_width_from_env() {
    if (const char* raw = std::getenv("BFS_OPEN_BEAM_WIDTH")) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return 10;
}

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

long long total_raw_edges(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.raw_edge_count_total;
    }
    return total;
}

long long total_accepted_edges(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.accepted_edge_count_total;
    }
    return total;
}

long long total_flyby_rejects(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.reject_flyby_physical_infeasible_count;
    }
    return total;
}

long long total_beam_pruned(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.beam_pruned_count;
    }
    return total;
}

std::vector<spaceship_cpp::planet_params::PlanetId> sequence_from_path(
    spaceship_cpp::planet_params::PlanetId launch_planet,
    const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& path
) {
    std::vector<spaceship_cpp::planet_params::PlanetId> sequence;
    sequence.push_back(launch_planet);
    for (const auto& edge : path) {
        if (edge.valid) {
            sequence.push_back(edge.to_planet);
        }
    }
    return sequence;
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

void print_summary(
    spaceship_cpp::planet_params::PlanetId launch_planet,
    bool beam_enabled,
    int max_depth,
    const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result,
    double wall_time_ms
) {
    std::vector<spaceship_cpp::planet_params::PlanetId> best_sequence;
    if (result.best_terminal_solution.valid) {
        best_sequence = sequence_from_path(
            launch_planet,
            spaceship_cpp::bfs::reconstruct_fixed_launch_theta_path(
                result,
                result.best_terminal_solution.node_index));
    }
    std::cout << "BfsOpenSearchDepthWidthFlybySummary\n";
    std::cout << "beam_enabled=" << (beam_enabled ? 1 : 0) << '\n';
    std::cout << "max_depth=" << max_depth << '\n';
    std::cout << "initial_candidate_count=" << result.initial_candidate_count << '\n';
    std::cout << "terminal_solution_count=" << result.terminal_solutions.size() << '\n';
    std::cout << "best_valid=" << (result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "best_total_delta_v=" << result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_terminal_solution.launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_terminal_solution.arrival_v_inf << '\n';
    std::cout << "best_sequence=" << sequence_string(best_sequence) << '\n';
    std::cout << "node_count=" << result.nodes.size() << '\n';
    std::cout << "total_raw_edge_count=" << total_raw_edges(result) << '\n';
    std::cout << "total_accepted_edge_count=" << total_accepted_edges(result) << '\n';
    std::cout << "total_flyby_reject_count=" << total_flyby_rejects(result) << '\n';
    std::cout << "total_beam_pruned_count=" << total_beam_pruned(result) << '\n';
    std::cout << "wall_time_ms=" << wall_time_ms << '\n';
    std::cout << "profile_ok=1\n";
}

void print_layers(
    bool beam_enabled,
    const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result
) {
    for (const auto& stats : result.depth_width_stats) {
        std::cout << "BfsOpenSearchDepthWidthLayerStats\n";
        std::cout << "beam_enabled=" << (beam_enabled ? 1 : 0) << '\n';
        std::cout << "depth=" << stats.depth << '\n';
        std::cout << "input_state_count=" << stats.input_state_count << '\n';
        std::cout << "attempted_expansion_count=" << stats.attempted_expansion_count << '\n';
        std::cout << "raw_edge_count_total=" << stats.raw_edge_count_total << '\n';
        std::cout << "accepted_edge_count_total=" << stats.accepted_edge_count_total << '\n';
        std::cout << "reject_flyby_physical_infeasible_count="
                  << stats.reject_flyby_physical_infeasible_count << '\n';
        std::cout << "output_state_count_before_beam=" << stats.output_state_count_before_beam << '\n';
        std::cout << "output_state_count_after_beam=" << stats.output_state_count_after_beam << '\n';
        std::cout << "beam_pruned_count=" << stats.beam_pruned_count << '\n';
        std::cout << "layer_wall_time_ms=" << stats.layer_wall_time_ms << '\n';
    }
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
        std::cout << "bfs_open_search_depth_width_with_flyby_filter_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double launch_time_offset_seconds = -315583.195659;
    const double launch_time =
        0.17 * planet_params::planet_orbital_period(launch_planet) + launch_time_offset_seconds;
    const double initial_theta = 2.572077523548;
    const int max_depth = open_depth_from_env();
    const int beam_width = open_beam_width_from_env();

    const std::vector<planet_params::PlanetId> allowed_first_targets{
        planet_params::PlanetId::Mercury,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mars,
    };
    const std::vector<planet_params::PlanetId> allowed_transfer_planets{
        planet_params::PlanetId::Mercury,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
    };

    bfs::InitialLaunchExpansionOptions initial_options{};
    initial_options.max_transfer_revolution = 1;
    initial_options.max_target_revolution = 1;
    initial_options.max_launch_v_inf = 7000.0;
    initial_options.time_weight_m_per_s_per_day = 0.0;

    problem2::Problem2GravityAssistSolverOptions problem2_options{};
    problem2_options.theta_sample_count = 64;
    problem2_options.topology_adaptive_enabled = true;

    auto run_search = [&](bool enable_beam_pruning) {
        bfs::FixedLaunchThetaSearchOptions search_options{};
        search_options.max_depth = max_depth;
        search_options.beam_width = beam_width;
        search_options.enable_beam_pruning = enable_beam_pruning;
        search_options.flyby_physical_options.enabled = true;
        search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
        search_options.max_launch_v_inf = 7000.0;
        search_options.time_weight_m_per_s_per_day = 0.0;
        search_options.continue_after_reaching_terminal = true;
        const auto start = Clock::now();
        auto result = bfs::run_fixed_launch_theta_beam_search_with_table(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            initial_theta,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options);
        const double wall_time_ms = elapsed_ms(start, Clock::now());
        return std::pair<bfs::FixedLaunchThetaSearchResult, double>{std::move(result), wall_time_ms};
    };

    auto [beam_result, beam_wall_time_ms] = run_search(true);
    auto [no_beam_result, no_beam_wall_time_ms] = run_search(false);
    assert(beam_result.ok);
    assert(no_beam_result.ok);

    print_summary(launch_planet, true, max_depth, beam_result, beam_wall_time_ms);
    print_layers(true, beam_result);
    print_summary(launch_planet, false, max_depth, no_beam_result, no_beam_wall_time_ms);
    print_layers(false, no_beam_result);

    std::cout << "BfsOpenSearchBeamVsNoBeamComparison\n";
    std::cout << "max_depth=" << max_depth << '\n';
    std::cout << "beam_enabled_node_count=" << beam_result.nodes.size() << '\n';
    std::cout << "beam_enabled_terminal_solution_count=" << beam_result.terminal_solutions.size() << '\n';
    std::cout << "beam_enabled_best_valid=" << (beam_result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "beam_enabled_best_total_delta_v=" << beam_result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "beam_enabled_total_raw_edges=" << total_raw_edges(beam_result) << '\n';
    std::cout << "beam_enabled_total_accepted_edges=" << total_accepted_edges(beam_result) << '\n';
    std::cout << "beam_enabled_total_flyby_reject_count=" << total_flyby_rejects(beam_result) << '\n';
    std::cout << "beam_enabled_total_beam_pruned_count=" << total_beam_pruned(beam_result) << '\n';
    std::cout << "beam_enabled_wall_time_ms=" << beam_wall_time_ms << '\n';
    std::cout << "beam_disabled_node_count=" << no_beam_result.nodes.size() << '\n';
    std::cout << "beam_disabled_terminal_solution_count=" << no_beam_result.terminal_solutions.size() << '\n';
    std::cout << "beam_disabled_best_valid=" << (no_beam_result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "beam_disabled_best_total_delta_v=" << no_beam_result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "beam_disabled_total_raw_edges=" << total_raw_edges(no_beam_result) << '\n';
    std::cout << "beam_disabled_total_accepted_edges=" << total_accepted_edges(no_beam_result) << '\n';
    std::cout << "beam_disabled_total_flyby_reject_count=" << total_flyby_rejects(no_beam_result) << '\n';
    std::cout << "beam_disabled_wall_time_ms=" << no_beam_wall_time_ms << '\n';
    std::cout << "comparison_ok=1\n";

    return 0;
}
