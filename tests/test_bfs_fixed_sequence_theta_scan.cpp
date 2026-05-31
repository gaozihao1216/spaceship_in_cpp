#include "spaceship_cpp/bfs/fixed_sequence_theta_scan.hpp"
#include "spaceship_cpp/bfs/theta_candidate_selector.hpp"
#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
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

void print_sequence_scan_summary(
    const std::vector<spaceship_cpp::planet_params::PlanetId>& sequence,
    const spaceship_cpp::bfs::FixedLaunchTimeThetaScanResult& result,
    double wall_time_ms
) {
    std::cout << "BfsFixedSequenceThetaScanSummary\n";
    std::cout << "sequence=" << sequence_string(sequence) << '\n';
    std::cout << "theta_scanned_count=" << result.theta_scanned_count << '\n';
    std::cout << "terminal_theta_count=" << result.terminal_theta_count << '\n';
    std::cout << "best_solution_valid=" << (result.best_solution.valid ? 1 : 0) << '\n';
    std::cout << "best_theta=" << result.best_solution.theta << '\n';
    std::cout << "best_score=" << result.best_solution.terminal_solution.score << '\n';
    std::cout << "best_total_delta_v=" << result.best_solution.terminal_solution.total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_solution.terminal_solution.launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_solution.terminal_solution.arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds="
              << result.best_solution.terminal_solution.total_flight_time_seconds << '\n';
    std::cout << "wall_time_ms=" << wall_time_ms << '\n';
    std::cout << "sequence_scan_ok=" << (result.ok && result.best_solution.valid ? 1 : 0) << '\n';
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
        std::cout << "bfs_fixed_sequence_theta_scan_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet);
    const std::vector<planet_params::PlanetId> allowed_first_targets{
        planet_params::PlanetId::Mercury,
        planet_params::PlanetId::Venus,
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

    bfs::FixedSequenceThetaSearchOptions sequence_options{};
    sequence_options.beam_width = 10;
    sequence_options.max_launch_v_inf = 7000.0;
    sequence_options.time_weight_m_per_s_per_day = 0.0;

    bfs::FixedLaunchTimeThetaScanOptions scan_options{};
    scan_options.max_theta_to_scan = 20;
    scan_options.theta_selection_mode = bfs::ThetaScanSelectionMode::Hybrid;
    scan_options.hybrid_lowest_vinf_count = 6;
    scan_options.hybrid_stratified_count = 14;
    scan_options.top_solution_count = 10;
    scan_options.top_sequence_count = 2;

    const std::vector<planet_params::PlanetId> sequence_a{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };
    const std::vector<planet_params::PlanetId> sequence_b{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mercury,
    };

    const auto start_a = Clock::now();
    const auto result_a = bfs::scan_fixed_sequence_over_theta_candidates_with_table(
        loader,
        launch_time,
        selection.candidates,
        sequence_a,
        initial_options,
        problem2_options,
        sequence_options,
        scan_options);
    const double wall_a_ms = elapsed_ms(start_a, Clock::now());
    assert(result_a.ok);

    const auto start_b = Clock::now();
    const auto result_b = bfs::scan_fixed_sequence_over_theta_candidates_with_table(
        loader,
        launch_time,
        selection.candidates,
        sequence_b,
        initial_options,
        problem2_options,
        sequence_options,
        scan_options);
    const double wall_b_ms = elapsed_ms(start_b, Clock::now());
    assert(result_b.ok);

    print_sequence_scan_summary(sequence_a, result_a, wall_a_ms);
    print_sequence_scan_summary(sequence_b, result_b, wall_b_ms);

    const bool a_valid = result_a.best_solution.valid;
    const bool b_valid = result_b.best_solution.valid;
    std::string best_sequence = "none";
    if (a_valid && (!b_valid ||
        result_a.best_solution.terminal_solution.total_delta_v <=
            result_b.best_solution.terminal_solution.total_delta_v)) {
        best_sequence = sequence_string(sequence_a);
    } else if (b_valid) {
        best_sequence = sequence_string(sequence_b);
    }

    std::cout << "BfsFixedSequenceThetaScanComparison\n";
    std::cout << "sequence_a=" << sequence_string(sequence_a) << '\n';
    std::cout << "sequence_a_best_total_delta_v=" << result_a.best_solution.terminal_solution.total_delta_v << '\n';
    std::cout << "sequence_b=" << sequence_string(sequence_b) << '\n';
    std::cout << "sequence_b_best_total_delta_v=" << result_b.best_solution.terminal_solution.total_delta_v << '\n';
    std::cout << "best_sequence=" << best_sequence << '\n';
    std::cout << "comparison_ok=" << ((a_valid || b_valid) ? 1 : 0) << '\n';

    return 0;
}
