#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
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
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct CaseSummary {
    int max_depth = 0;
    int max_kq = 0;
    int theta_scanned_count = 0;
    int terminal_theta_count = 0;
    int total_terminal_solution_count = 0;

    bool best_valid = false;
    double best_theta = 0.0;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();
    std::vector<spaceship_cpp::planet_params::PlanetId> best_sequence;
    spaceship_cpp::bfs::FixedLaunchThetaSearchResult best_search;

    long long total_raw_edge_count = 0;
    long long total_accepted_edge_count = 0;
    long long total_flyby_reject_count = 0;
    long long total_beam_pruned_count = 0;

    int max_layer_output_width_seen = 0;
    double total_wall_time_ms = 0.0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

int theta_max_from_env() {
    if (const char* raw = std::getenv("BFS_SENSITIVITY_THETA_MAX")) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return 20;
}

bool include_kq3_from_env() {
    if (const char* raw = std::getenv("BFS_SENSITIVITY_INCLUDE_KQ3")) {
        return std::atoi(raw) != 0;
    }
    return false;
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

int max_output_width(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    int width = 0;
    for (const auto& stats : result.depth_width_stats) {
        width = std::max(width, stats.output_state_count_after_beam);
    }
    return width;
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

void print_case(const CaseSummary& summary) {
    std::cout << "BfsPhysicalDepthKqSensitivityCase\n";
    std::cout << "max_depth=" << summary.max_depth << '\n';
    std::cout << "max_kq=" << summary.max_kq << '\n';
    std::cout << "theta_scanned_count=" << summary.theta_scanned_count << '\n';
    std::cout << "terminal_theta_count=" << summary.terminal_theta_count << '\n';
    std::cout << "total_terminal_solution_count=" << summary.total_terminal_solution_count << '\n';
    std::cout << "best_valid=" << (summary.best_valid ? 1 : 0) << '\n';
    std::cout << "best_theta=" << summary.best_theta << '\n';
    std::cout << "best_sequence=" << sequence_string(summary.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << summary.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << summary.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << summary.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << summary.best_total_flight_time_seconds << '\n';
    std::cout << "total_raw_edge_count=" << summary.total_raw_edge_count << '\n';
    std::cout << "total_accepted_edge_count=" << summary.total_accepted_edge_count << '\n';
    std::cout << "total_flyby_reject_count=" << summary.total_flyby_reject_count << '\n';
    std::cout << "total_beam_pruned_count=" << summary.total_beam_pruned_count << '\n';
    std::cout << "max_layer_output_width_seen=" << summary.max_layer_output_width_seen << '\n';
    std::cout << "total_wall_time_ms=" << summary.total_wall_time_ms << '\n';
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
        std::cout << "bfs_physical_search_depth_kq_sensitivity_skipped_missing_table\n";
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

    const int theta_scanned_count =
        std::min(theta_max_from_env(), static_cast<int>(selection.candidates.size()));
    const std::vector<int> depth_values{4, 5, 6};
    std::vector<int> max_kq_values{1, 2};
    if (include_kq3_from_env()) {
        max_kq_values.push_back(3);
    }

    std::vector<CaseSummary> cases;
    CaseSummary global_best{};

    for (const int max_depth : depth_values) {
        for (const int max_kq : max_kq_values) {
            bfs::InitialLaunchExpansionOptions initial_options{};
            initial_options.max_transfer_revolution = max_kq;
            initial_options.max_target_revolution = max_kq;
            initial_options.max_launch_v_inf = 7000.0;
            initial_options.time_weight_m_per_s_per_day = 0.0;

            problem2::Problem2GravityAssistSolverOptions problem2_options{};
            problem2_options.theta_sample_count = 64;
            problem2_options.topology_adaptive_enabled = true;
            problem2_options.max_transfer_revolution = max_kq;
            problem2_options.max_target_revolution = max_kq;

            bfs::FixedLaunchThetaSearchOptions search_options{};
            search_options.max_depth = max_depth;
            search_options.beam_width = 3;
            search_options.enable_beam_pruning = true;
            search_options.flyby_physical_options.enabled = true;
            search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
            search_options.max_launch_v_inf = 7000.0;
            search_options.time_weight_m_per_s_per_day = 0.0;
            search_options.continue_after_reaching_terminal = true;

            CaseSummary summary{};
            summary.max_depth = max_depth;
            summary.max_kq = max_kq;
            summary.theta_scanned_count = theta_scanned_count;

            for (int theta_index = 0; theta_index < theta_scanned_count; ++theta_index) {
                const auto& candidate = selection.candidates[static_cast<std::size_t>(theta_index)];
                const auto start = Clock::now();
                auto search = bfs::run_fixed_launch_theta_beam_search_with_table(
                    loader,
                    launch_planet,
                    terminal_planet,
                    launch_time,
                    candidate.theta,
                    allowed_first_targets,
                    allowed_transfer_planets,
                    initial_options,
                    problem2_options,
                    search_options);
                const double wall_time_ms = elapsed_ms(start, Clock::now());
                assert(search.ok);

                summary.total_wall_time_ms += wall_time_ms;
                summary.total_raw_edge_count += total_raw_edges(search);
                summary.total_accepted_edge_count += total_accepted_edges(search);
                summary.total_flyby_reject_count += total_flyby_rejects(search);
                summary.total_beam_pruned_count += total_beam_pruned(search);
                summary.max_layer_output_width_seen =
                    std::max(summary.max_layer_output_width_seen, max_output_width(search));
                summary.total_terminal_solution_count += static_cast<int>(search.terminal_solutions.size());

                if (search.best_terminal_solution.valid) {
                    summary.terminal_theta_count += 1;
                    if (!summary.best_valid ||
                        search.best_terminal_solution.total_delta_v < summary.best_total_delta_v) {
                        const auto path = bfs::reconstruct_fixed_launch_theta_path(
                            search, search.best_terminal_solution.node_index);
                        summary.best_valid = true;
                        summary.best_theta = candidate.theta;
                        summary.best_total_delta_v = search.best_terminal_solution.total_delta_v;
                        summary.best_launch_v_inf = search.best_terminal_solution.launch_v_inf;
                        summary.best_arrival_v_inf = search.best_terminal_solution.arrival_v_inf;
                        summary.best_total_flight_time_seconds =
                            search.best_terminal_solution.total_flight_time_seconds;
                        summary.best_sequence = sequence_from_path(launch_planet, path);
                        summary.best_search = search;
                    }
                }
            }

            print_case(summary);
            if (summary.best_valid &&
                (!global_best.best_valid || summary.best_total_delta_v < global_best.best_total_delta_v)) {
                global_best = summary;
            }
            cases.push_back(std::move(summary));
        }
    }

    std::cout << "BfsPhysicalDepthKqSensitivitySummary\n";
    std::cout << "case_count=" << cases.size() << '\n';
    std::cout << "theta_candidate_count=" << selection.candidates.size() << '\n';
    std::cout << "theta_scanned_count=" << theta_scanned_count << '\n';
    std::cout << "any_terminal_found=" << (global_best.best_valid ? 1 : 0) << '\n';
    std::cout << "best_max_depth=" << global_best.max_depth << '\n';
    std::cout << "best_max_kq=" << global_best.max_kq << '\n';
    std::cout << "best_theta=" << global_best.best_theta << '\n';
    std::cout << "best_sequence=" << sequence_string(global_best.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << global_best.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << global_best.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << global_best.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << global_best.best_total_flight_time_seconds << '\n';
    std::cout << "best_case_max_layer_output_width=" << global_best.max_layer_output_width_seen << '\n';
    std::cout << "sensitivity_ok=1\n";

    if (global_best.best_valid) {
        bfs::InitialLaunchExpansionOptions initial_options{};
        initial_options.max_transfer_revolution = global_best.max_kq;
        initial_options.max_target_revolution = global_best.max_kq;
        initial_options.max_launch_v_inf = 7000.0;
        initial_options.time_weight_m_per_s_per_day = 0.0;

        problem2::Problem2GravityAssistSolverOptions problem2_options{};
        problem2_options.theta_sample_count = 64;
        problem2_options.topology_adaptive_enabled = true;
        problem2_options.max_transfer_revolution = global_best.max_kq;
        problem2_options.max_target_revolution = global_best.max_kq;

        bfs::FixedLaunchThetaSearchOptions search_options{};
        search_options.max_depth = global_best.max_depth;
        search_options.beam_width = 3;
        search_options.enable_beam_pruning = false;
        search_options.flyby_physical_options.enabled = true;
        search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
        search_options.max_launch_v_inf = 7000.0;
        search_options.time_weight_m_per_s_per_day = 0.0;
        search_options.continue_after_reaching_terminal = true;

        const auto start = Clock::now();
        const auto no_beam = bfs::run_fixed_launch_theta_beam_search_with_table(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            global_best.best_theta,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options);
        const double wall_time_ms = elapsed_ms(start, Clock::now());
        assert(no_beam.ok);

        std::cout << "BfsPhysicalDepthKqNoBeamValidation\n";
        std::cout << "max_depth=" << global_best.max_depth << '\n';
        std::cout << "max_kq=" << global_best.max_kq << '\n';
        std::cout << "theta=" << global_best.best_theta << '\n';
        std::cout << "beam_best_total_delta_v=" << global_best.best_total_delta_v << '\n';
        std::cout << "no_beam_best_valid=" << (no_beam.best_terminal_solution.valid ? 1 : 0) << '\n';
        std::cout << "no_beam_best_total_delta_v="
                  << (no_beam.best_terminal_solution.valid
                          ? no_beam.best_terminal_solution.total_delta_v
                          : std::numeric_limits<double>::infinity()) << '\n';
        std::cout << "beam_node_count=" << global_best.best_search.nodes.size() << '\n';
        std::cout << "no_beam_node_count=" << no_beam.nodes.size() << '\n';
        std::cout << "no_beam_max_layer_output_width=" << max_output_width(no_beam) << '\n';
        std::cout << "no_beam_wall_time_ms=" << wall_time_ms << '\n';
        std::cout << "no_beam_validation_ok=1\n";
    }

    assert(theta_scanned_count > 0);
    return 0;
}
