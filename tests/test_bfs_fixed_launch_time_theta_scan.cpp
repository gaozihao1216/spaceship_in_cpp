#include "spaceship_cpp/bfs/fixed_launch_time_theta_scan.hpp"
#include "spaceship_cpp/bfs/theta_candidate_selector.hpp"
#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <cassert>
#include <algorithm>
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
    return "/tmp/problem1_root_table_2deg_full";
}

int scan_max_from_env() {
    if (const char* raw = std::getenv("BFS_THETA_SCAN_MAX")) {
        if (*raw != '\0') {
            return std::max(1, std::atoi(raw));
        }
    }
    return 20;
}

bool comparison_enabled_from_env() {
    if (const char* raw = std::getenv("BFS_THETA_SCAN_COMPARE")) {
        return std::atoi(raw) != 0;
    }
    return true;
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

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "bfs_fixed_launch_time_theta_scan_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet);
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
    selector_options.include_near_valid_intervals = false;
    const auto selection = bfs::select_theta_candidates_from_launch_scout(scout, selector_options);
    assert(selection.ok);

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

    bfs::FixedLaunchTimeThetaScanOptions scan_options{};
    scan_options.max_theta_to_scan = scan_max_from_env();
    scan_options.theta_selection_mode = bfs::ThetaScanSelectionMode::Hybrid;
    scan_options.hybrid_lowest_vinf_count = 6;
    scan_options.hybrid_stratified_count = 14;
    scan_options.top_solution_count = 10;
    scan_options.top_sequence_count = 10;
    const auto result = bfs::scan_fixed_launch_time_theta_candidates_with_table(
        loader,
        launch_planet,
        terminal_planet,
        launch_time,
        selection.candidates,
        allowed_first_targets,
        allowed_transfer_planets,
        initial_options,
        problem2_options,
        search_options,
        scan_options);
    assert(result.ok);
    assert(result.theta_scanned_count > 0);
    assert(result.theta_scanned_count <= scan_options.max_theta_to_scan);
    if (result.terminal_theta_count > 0) {
        assert(result.best_solution.valid);
        assert(!result.best_solution.path.empty());
        assert(result.best_solution.sequence.front() == launch_planet);
        assert(result.best_solution.sequence.back() == terminal_planet);
        assert(result.best_solution.terminal_solution.launch_v_inf <= 7000.0);
        assert(std::isfinite(result.best_solution.terminal_solution.total_delta_v));
    }

    double min_scanned_theta = std::numeric_limits<double>::infinity();
    double max_scanned_theta = -std::numeric_limits<double>::infinity();
    for (const auto& summary : result.theta_summaries) {
        min_scanned_theta = std::min(min_scanned_theta, summary.theta);
        max_scanned_theta = std::max(max_scanned_theta, summary.theta);
    }
    const double scanned_theta_width = result.theta_summaries.empty()
        ? 0.0
        : max_scanned_theta - min_scanned_theta;

    std::cout << "BfsFixedLaunchTimeThetaScanSelectionSummary\n";
    std::cout << "selection_mode=" << bfs::theta_scan_selection_mode_name(scan_options.theta_selection_mode) << '\n';
    std::cout << "max_theta_to_scan=" << scan_options.max_theta_to_scan << '\n';
    std::cout << "hybrid_lowest_vinf_count=" << scan_options.hybrid_lowest_vinf_count << '\n';
    std::cout << "hybrid_stratified_count=" << scan_options.hybrid_stratified_count << '\n';
    std::cout << "theta_scanned_count=" << result.theta_scanned_count << '\n';
    std::cout << "min_scanned_theta=" << min_scanned_theta << '\n';
    std::cout << "max_scanned_theta=" << max_scanned_theta << '\n';
    std::cout << "scanned_theta_width=" << scanned_theta_width << '\n';

    std::cout << "BfsFixedLaunchTimeThetaScanSummary\n";
    std::cout << "theta_candidate_count=" << result.theta_candidate_count << '\n';
    std::cout << "theta_scanned_count=" << result.theta_scanned_count << '\n';
    std::cout << "terminal_theta_count=" << result.terminal_theta_count << '\n';
    std::cout << "total_terminal_solution_count=" << result.total_terminal_solution_count << '\n';
    std::cout << "top_solution_count=" << result.top_solutions.size() << '\n';
    std::cout << "top_sequence_count=" << result.top_sequences.size() << '\n';
    std::cout << "best_solution_valid=" << (result.best_solution.valid ? 1 : 0) << '\n';
    std::cout << "best_theta=" << result.best_solution.theta << '\n';
    std::cout << "best_score=" << result.best_solution.terminal_solution.score << '\n';
    std::cout << "best_total_delta_v=" << result.best_solution.terminal_solution.total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_solution.terminal_solution.launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_solution.terminal_solution.arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds="
              << result.best_solution.terminal_solution.total_flight_time_seconds << '\n';
    std::cout << "scan_ok=1\n";

    for (std::size_t i = 0; i < result.theta_summaries.size(); ++i) {
        const auto& summary = result.theta_summaries[i];
        std::cout << "BfsFixedLaunchTimeThetaSummary\n";
        std::cout << "index=" << i << '\n';
        std::cout << "theta=" << summary.theta << '\n';
        std::cout << "estimated_min_launch_v_inf=" << summary.estimated_min_launch_v_inf << '\n';
        std::cout << "initial_candidate_count=" << summary.initial_candidate_count << '\n';
        std::cout << "terminal_solution_count=" << summary.terminal_solution_count << '\n';
        std::cout << "node_count=" << summary.node_count << '\n';
        std::cout << "has_terminal=" << (summary.has_terminal ? 1 : 0) << '\n';
        std::cout << "best_score=" << summary.best_score << '\n';
        std::cout << "best_total_delta_v=" << summary.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << summary.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << summary.best_arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds=" << summary.best_total_flight_time_seconds << '\n';
        std::cout << "sequence=" << sequence_string(summary.best_sequence) << '\n';
    }

    for (std::size_t i = 0; i < result.top_sequences.size(); ++i) {
        const auto& sequence = result.top_sequences[i];
        std::cout << "BfsFixedLaunchTimeTopSequence\n";
        std::cout << "rank=" << i << '\n';
        std::cout << "sequence=" << sequence_string(sequence.sequence) << '\n';
        std::cout << "occurrence_count=" << sequence.occurrence_count << '\n';
        std::cout << "best_theta=" << sequence.best_theta << '\n';
        std::cout << "best_score=" << sequence.best_score << '\n';
        std::cout << "best_total_delta_v=" << sequence.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << sequence.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << sequence.best_arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds=" << sequence.best_total_flight_time_seconds << '\n';
    }

    for (std::size_t i = 0; i < result.top_solutions.size(); ++i) {
        const auto& solution = result.top_solutions[i];
        std::cout << "BfsFixedLaunchTimeTopSolution\n";
        std::cout << "rank=" << i << '\n';
        std::cout << "theta=" << solution.theta << '\n';
        std::cout << "sequence=" << sequence_string(solution.sequence) << '\n';
        std::cout << "score=" << solution.terminal_solution.score << '\n';
        std::cout << "total_delta_v=" << solution.terminal_solution.total_delta_v << '\n';
        std::cout << "launch_v_inf=" << solution.terminal_solution.launch_v_inf << '\n';
        std::cout << "arrival_v_inf=" << solution.terminal_solution.arrival_v_inf << '\n';
        std::cout << "total_flight_time_seconds=" << solution.terminal_solution.total_flight_time_seconds << '\n';
        std::cout << "path_length=" << solution.path.size() << '\n';
    }

    if (comparison_enabled_from_env()) {
        auto lowest_options = scan_options;
        lowest_options.max_theta_to_scan = std::min(scan_options.max_theta_to_scan, 10);
        lowest_options.theta_selection_mode = bfs::ThetaScanSelectionMode::LowestEstimatedVInf;
        const auto lowest_result = bfs::scan_fixed_launch_time_theta_candidates_with_table(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            selection.candidates,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options,
            lowest_options);
        assert(lowest_result.ok);

        std::cout << "BfsFixedLaunchTimeThetaScanModeComparison\n";
        std::cout << "lowest_mode_best_score=" << lowest_result.best_solution.terminal_solution.score << '\n';
        std::cout << "lowest_mode_best_theta=" << lowest_result.best_solution.theta << '\n';
        std::cout << "lowest_mode_best_sequence=" << sequence_string(lowest_result.best_solution.sequence) << '\n';
        std::cout << "hybrid_mode_best_score=" << result.best_solution.terminal_solution.score << '\n';
        std::cout << "hybrid_mode_best_theta=" << result.best_solution.theta << '\n';
        std::cout << "hybrid_mode_best_sequence=" << sequence_string(result.best_solution.sequence) << '\n';
        std::cout << "mode_comparison_ok=1\n";
    }
    return 0;
}
