#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/bfs/theta_candidate_selector.hpp"
#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"
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
#include <vector>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

void print_layer_summaries(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    for (const auto& summary : result.layer_summaries) {
        std::cout << "BfsFixedLaunchThetaSearchLayerSummary\n";
        std::cout << "depth=" << summary.depth << '\n';
        std::cout << "input_state_count=" << summary.input_state_count << '\n';
        std::cout << "attempted_expansion_count=" << summary.attempted_expansion_count << '\n';
        std::cout << "generated_edge_count=" << summary.generated_edge_count << '\n';
        std::cout << "accepted_edge_count=" << summary.accepted_edge_count << '\n';
        std::cout << "terminal_solution_count_after_layer="
                  << summary.terminal_solution_count_after_layer << '\n';
        std::cout << "output_state_count=" << summary.output_state_count << '\n';
        std::cout << "layer_wall_time_ms=" << summary.layer_wall_time_ms << '\n';
    }
}

void print_path(
    const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& path
) {
    namespace planet_params = spaceship_cpp::planet_params;
    for (std::size_t i = 0; i < path.size(); ++i) {
        const auto& edge = path[i];
        std::cout << "BfsFixedLaunchThetaBestPathEdge\n";
        std::cout << "index=" << i << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(edge.from_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(edge.to_planet) << '\n';
        std::cout << "departure_time=" << edge.departure_time << '\n';
        std::cout << "arrival_time=" << edge.arrival_time << '\n';
        std::cout << "transfer_time_seconds=" << edge.transfer_time_seconds << '\n';
        std::cout << "outgoing_e=" << edge.outgoing_e << '\n';
        std::cout << "outgoing_theta=" << edge.outgoing_theta << '\n';
        std::cout << "transfer_revolution=" << edge.transfer_revolution << '\n';
        std::cout << "target_revolution=" << edge.target_revolution << '\n';
        std::cout << "slingshot_residual=" << edge.slingshot_residual << '\n';
        std::cout << "problem1_residual_seconds=" << edge.problem1_residual_seconds << '\n';
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
        std::cout << "bfs_fixed_launch_theta_search_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet);
    const double initial_theta = 2.783800156931;
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

    bfs::FixedLaunchThetaSearchOptions search_options{};
    search_options.max_depth = 3;
    search_options.beam_width = 10;
    search_options.max_launch_v_inf = 7000.0;
    search_options.time_weight_m_per_s_per_day = 0.0;
    search_options.continue_after_reaching_terminal = true;

    const auto result = bfs::run_fixed_launch_theta_beam_search_with_table(
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
    assert(result.ok);
    assert(result.initial_candidate_count > 0);
    assert(!result.nodes.empty());
    assert(!result.layer_summaries.empty());
    assert(result.layer_summaries.front().output_state_count <= search_options.beam_width);

    std::vector<bfs::TrajectorySearchEdge> best_path;
    if (!result.terminal_solutions.empty()) {
        assert(result.best_terminal_solution.valid);
        assert(std::isfinite(result.best_terminal_solution.total_delta_v));
        assert(result.best_terminal_solution.launch_v_inf <= 7000.0);
        best_path = bfs::reconstruct_fixed_launch_theta_path(
            result, result.best_terminal_solution.node_index);
        assert(static_cast<int>(best_path.size()) <= search_options.max_depth);
        assert(!best_path.empty());
        assert(best_path.front().from_planet == launch_planet);
        assert(best_path.back().to_planet == terminal_planet);
    }

    std::cout << "BfsFixedLaunchThetaSearchSummary\n";
    std::cout << "initial_theta=" << result.initial_theta << '\n';
    std::cout << "initial_candidate_count=" << result.initial_candidate_count << '\n';
    std::cout << "launch_v_inf_pruned_count=" << result.launch_v_inf_pruned_count << '\n';
    std::cout << "node_count=" << result.nodes.size() << '\n';
    std::cout << "terminal_solution_count=" << result.terminal_solutions.size() << '\n';
    std::cout << "best_solution_valid=" << (result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "best_score=" << result.best_terminal_solution.score << '\n';
    std::cout << "best_total_delta_v=" << result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_terminal_solution.launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_terminal_solution.arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds="
              << result.best_terminal_solution.total_flight_time_seconds << '\n';
    std::cout << "fixed_theta_search_ok=1\n";
    print_layer_summaries(result);
    if (!best_path.empty()) {
        print_path(best_path);
    }

    bfs::ThetaLaunchFeasibilityScoutOptions scout_options{};
    scout_options.theta_scout_count = 720;
    scout_options.max_launch_v_inf = 7000.0;
    scout_options.near_v_inf_buffer = 1000.0;
    scout_options.refine_interval_boundaries = true;
    const auto scout = bfs::scout_theta_launch_feasibility_with_problem1_table(
        loader, launch_planet, launch_time, allowed_first_targets, scout_options);
    assert(scout.ok);

    bfs::ThetaCandidateSelectorOptions selector_options{};
    selector_options.max_theta_candidates = 50;
    selector_options.max_samples_per_valid_interval = 30;
    const auto selection = bfs::select_theta_candidates_from_launch_scout(scout, selector_options);
    assert(selection.ok);
    const int theta_candidate_tested_count = std::min(3, selection.candidate_count);
    int terminal_theta_count = 0;
    double best_theta = 0.0;
    double best_score = std::numeric_limits<double>::infinity();
    for (int i = 0; i < theta_candidate_tested_count; ++i) {
        const double theta = selection.candidates[static_cast<std::size_t>(i)].theta;
        const auto theta_result = bfs::run_fixed_launch_theta_beam_search_with_table(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            theta,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options);
        assert(theta_result.ok);
        if (theta_result.best_terminal_solution.valid) {
            terminal_theta_count += 1;
            if (theta_result.best_terminal_solution.score < best_score) {
                best_score = theta_result.best_terminal_solution.score;
                best_theta = theta;
            }
        }
    }

    std::cout << "BfsFixedLaunchThetaSelectorIntegrationSummary\n";
    std::cout << "theta_candidate_tested_count=" << theta_candidate_tested_count << '\n';
    std::cout << "terminal_theta_count=" << terminal_theta_count << '\n';
    std::cout << "best_theta=" << best_theta << '\n';
    std::cout << "best_score=" << best_score << '\n';
    std::cout << "selector_integration_ok=1\n";
    return 0;
}
