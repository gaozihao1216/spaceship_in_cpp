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
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct SearchRun {
    spaceship_cpp::bfs::FixedLaunchThetaSearchResult result;
    double wall_time_ms = 0.0;
};

struct ThetaRunSummary {
    int index = 0;
    spaceship_cpp::bfs::ThetaCandidate candidate;
    SearchRun run;
    std::vector<spaceship_cpp::planet_params::PlanetId> best_sequence;
    int max_output_width = 0;
};

struct SequenceAggregate {
    std::vector<spaceship_cpp::planet_params::PlanetId> sequence;
    int occurrence_count = 0;
    double best_theta = 0.0;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

int scan_max_from_env() {
    if (const char* raw = std::getenv("BFS_PHYSICAL_THETA_SCAN_MAX")) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return 20;
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

SearchRun run_search(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId launch_planet,
    spaceship_cpp::planet_params::PlanetId terminal_planet,
    double launch_time,
    double theta,
    const std::vector<spaceship_cpp::planet_params::PlanetId>& allowed_first_targets,
    const std::vector<spaceship_cpp::planet_params::PlanetId>& allowed_transfer_planets,
    const spaceship_cpp::bfs::InitialLaunchExpansionOptions& initial_options,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverOptions& problem2_options,
    spaceship_cpp::bfs::FixedLaunchThetaSearchOptions search_options
) {
    const auto start = Clock::now();
    spaceship_cpp::bfs::FixedLaunchThetaSearchResult result;
    try {
        result = spaceship_cpp::bfs::run_fixed_launch_theta_beam_search_with_table(
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
    } catch (const std::exception& ex) {
        result.ok = false;
        result.error_message = ex.what();
    }
    const double wall_time_ms = elapsed_ms(start, Clock::now());
    return SearchRun{std::move(result), wall_time_ms};
}

double best_total_delta_v_or_inf(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    if (!result.best_terminal_solution.valid) {
        return std::numeric_limits<double>::infinity();
    }
    return result.best_terminal_solution.total_delta_v;
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
        std::cout << "bfs_open_search_physical_theta_scan_skipped_missing_table\n";
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

    bfs::FixedLaunchThetaSearchOptions base_search_options{};
    base_search_options.max_depth = 4;
    base_search_options.beam_width = 10;
    base_search_options.enable_beam_pruning = true;
    base_search_options.flyby_physical_options.enabled = true;
    base_search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
    base_search_options.max_launch_v_inf = 7000.0;
    base_search_options.time_weight_m_per_s_per_day = 0.0;
    base_search_options.continue_after_reaching_terminal = true;

    const int theta_scanned_count =
        std::min(scan_max_from_env(), static_cast<int>(selection.candidates.size()));
    std::vector<ThetaRunSummary> theta_runs;
    theta_runs.reserve(static_cast<std::size_t>(theta_scanned_count));

    int terminal_theta_count = 0;
    int best_index = -1;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    long long scan_total_raw_edges = 0;
    long long scan_total_accepted_edges = 0;
    long long scan_total_flyby_rejects = 0;
    long long scan_total_beam_pruned = 0;
    int scan_max_layer_output_width = 0;
    double scan_total_wall_time_ms = 0.0;

    for (int i = 0; i < theta_scanned_count; ++i) {
        const auto& candidate = selection.candidates[static_cast<std::size_t>(i)];
        auto run = run_search(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            candidate.theta,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            base_search_options);
        if (!run.result.ok) {
            theta_runs.push_back(ThetaRunSummary{i, candidate, std::move(run), {}, 0});
            const auto& stored = theta_runs.back();
            std::cout << "PhysicalThetaScanThetaSummary\n";
            std::cout << "index=" << stored.index << '\n';
            std::cout << "theta=" << stored.candidate.theta << '\n';
            std::cout << "estimated_min_launch_v_inf=" << stored.candidate.estimated_min_launch_v_inf << '\n';
            std::cout << "initial_candidate_count=0\n";
            std::cout << "terminal_solution_count=0\n";
            std::cout << "best_valid=0\n";
            std::cout << "best_total_delta_v=" << std::numeric_limits<double>::infinity() << '\n';
            std::cout << "best_launch_v_inf=" << std::numeric_limits<double>::infinity() << '\n';
            std::cout << "best_arrival_v_inf=" << std::numeric_limits<double>::infinity() << '\n';
            std::cout << "best_total_flight_time_seconds=" << std::numeric_limits<double>::infinity() << '\n';
            std::cout << "node_count=0\n";
            std::cout << "total_raw_edge_count=0\n";
            std::cout << "total_accepted_edge_count=0\n";
            std::cout << "total_flyby_reject_count=0\n";
            std::cout << "total_beam_pruned_count=0\n";
            std::cout << "max_output_width=0\n";
            std::cout << "wall_time_ms=" << stored.run.wall_time_ms << '\n';
            std::cout << "best_sequence=\n";
            std::cout << "error_message=" << stored.run.result.error_message << '\n';
            continue;
        }

        std::vector<planet_params::PlanetId> best_sequence;
        if (run.result.best_terminal_solution.valid) {
            const auto path = bfs::reconstruct_fixed_launch_theta_path(
                run.result, run.result.best_terminal_solution.node_index);
            best_sequence = sequence_from_path(launch_planet, path);
        }

        const int width = max_output_width(run.result);
        scan_total_raw_edges += total_raw_edges(run.result);
        scan_total_accepted_edges += total_accepted_edges(run.result);
        scan_total_flyby_rejects += total_flyby_rejects(run.result);
        scan_total_beam_pruned += total_beam_pruned(run.result);
        scan_max_layer_output_width = std::max(scan_max_layer_output_width, width);
        scan_total_wall_time_ms += run.wall_time_ms;

        if (run.result.best_terminal_solution.valid) {
            terminal_theta_count += 1;
            if (run.result.best_terminal_solution.total_delta_v < best_total_delta_v) {
                best_total_delta_v = run.result.best_terminal_solution.total_delta_v;
                best_index = i;
            }
        }

        theta_runs.push_back(ThetaRunSummary{i, candidate, std::move(run), std::move(best_sequence), width});
        const auto& stored = theta_runs.back();
        std::cout << "PhysicalThetaScanThetaSummary\n";
        std::cout << "index=" << stored.index << '\n';
        std::cout << "theta=" << stored.candidate.theta << '\n';
        std::cout << "estimated_min_launch_v_inf=" << stored.candidate.estimated_min_launch_v_inf << '\n';
        std::cout << "initial_candidate_count=" << stored.run.result.initial_candidate_count << '\n';
        std::cout << "terminal_solution_count=" << stored.run.result.terminal_solutions.size() << '\n';
        std::cout << "best_valid=" << (stored.run.result.best_terminal_solution.valid ? 1 : 0) << '\n';
        std::cout << "best_total_delta_v=" << best_total_delta_v_or_inf(stored.run.result) << '\n';
        std::cout << "best_launch_v_inf=" << stored.run.result.best_terminal_solution.launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << stored.run.result.best_terminal_solution.arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds="
                  << stored.run.result.best_terminal_solution.total_flight_time_seconds << '\n';
        std::cout << "node_count=" << stored.run.result.nodes.size() << '\n';
        std::cout << "total_raw_edge_count=" << total_raw_edges(stored.run.result) << '\n';
        std::cout << "total_accepted_edge_count=" << total_accepted_edges(stored.run.result) << '\n';
        std::cout << "total_flyby_reject_count=" << total_flyby_rejects(stored.run.result) << '\n';
        std::cout << "total_beam_pruned_count=" << total_beam_pruned(stored.run.result) << '\n';
        std::cout << "max_output_width=" << stored.max_output_width << '\n';
        std::cout << "wall_time_ms=" << stored.run.wall_time_ms << '\n';
        std::cout << "best_sequence=" << sequence_string(stored.best_sequence) << '\n';
        for (const auto& message : stored.run.result.initial_target_error_messages) {
            std::cout << "initial_target_error=" << message << '\n';
        }
    }

    std::map<std::string, SequenceAggregate> sequence_aggregates;
    for (const auto& theta_run : theta_runs) {
        for (const auto& terminal : theta_run.run.result.terminal_solutions) {
            if (!terminal.valid) {
                continue;
            }
            const auto path = bfs::reconstruct_fixed_launch_theta_path(theta_run.run.result, terminal.node_index);
            const auto sequence = sequence_from_path(launch_planet, path);
            const auto key = sequence_string(sequence);
            auto& aggregate = sequence_aggregates[key];
            if (aggregate.occurrence_count == 0) {
                aggregate.sequence = sequence;
            }
            aggregate.occurrence_count += 1;
            if (terminal.total_delta_v < aggregate.best_total_delta_v) {
                aggregate.best_theta = theta_run.candidate.theta;
                aggregate.best_total_delta_v = terminal.total_delta_v;
                aggregate.best_launch_v_inf = terminal.launch_v_inf;
                aggregate.best_arrival_v_inf = terminal.arrival_v_inf;
                aggregate.best_total_flight_time_seconds = terminal.total_flight_time_seconds;
            }
        }
    }
    std::vector<SequenceAggregate> top_sequences;
    for (const auto& [unused, aggregate] : sequence_aggregates) {
        (void)unused;
        top_sequences.push_back(aggregate);
    }
    std::sort(top_sequences.begin(), top_sequences.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.best_total_delta_v != rhs.best_total_delta_v) {
            return lhs.best_total_delta_v < rhs.best_total_delta_v;
        }
        return lhs.occurrence_count > rhs.occurrence_count;
    });

    const bool best_valid = best_index >= 0;
    const auto* best_run = best_valid ? &theta_runs[static_cast<std::size_t>(best_index)] : nullptr;
    std::cout << "BfsOpenSearchPhysicalThetaScanSummary\n";
    std::cout << "theta_candidate_count=" << selection.candidates.size() << '\n';
    std::cout << "theta_scanned_count=" << theta_scanned_count << '\n';
    std::cout << "terminal_theta_count=" << terminal_theta_count << '\n';
    std::cout << "best_valid=" << (best_valid ? 1 : 0) << '\n';
    std::cout << "best_theta=" << (best_valid ? best_run->candidate.theta : 0.0) << '\n';
    std::cout << "best_sequence=" << (best_valid ? sequence_string(best_run->best_sequence) : std::string{}) << '\n';
    std::cout << "best_total_delta_v=" << (best_valid ? best_run->run.result.best_terminal_solution.total_delta_v
                                                       : std::numeric_limits<double>::infinity()) << '\n';
    std::cout << "best_launch_v_inf=" << (best_valid ? best_run->run.result.best_terminal_solution.launch_v_inf
                                                     : std::numeric_limits<double>::infinity()) << '\n';
    std::cout << "best_arrival_v_inf=" << (best_valid ? best_run->run.result.best_terminal_solution.arrival_v_inf
                                                      : std::numeric_limits<double>::infinity()) << '\n';
    std::cout << "best_total_flight_time_seconds="
              << (best_valid ? best_run->run.result.best_terminal_solution.total_flight_time_seconds
                             : std::numeric_limits<double>::infinity()) << '\n';
    std::cout << "total_raw_edge_count=" << scan_total_raw_edges << '\n';
    std::cout << "total_accepted_edge_count=" << scan_total_accepted_edges << '\n';
    std::cout << "total_flyby_reject_count=" << scan_total_flyby_rejects << '\n';
    std::cout << "total_beam_pruned_count=" << scan_total_beam_pruned << '\n';
    std::cout << "max_layer_output_width_seen=" << scan_max_layer_output_width << '\n';
    std::cout << "total_wall_time_ms=" << scan_total_wall_time_ms << '\n';
    std::cout << "physical_theta_scan_ok=1\n";

    const std::size_t top_sequence_count = std::min<std::size_t>(10, top_sequences.size());
    for (std::size_t i = 0; i < top_sequence_count; ++i) {
        const auto& sequence = top_sequences[i];
        std::cout << "BfsOpenSearchPhysicalTopSequence\n";
        std::cout << "rank=" << i << '\n';
        std::cout << "sequence=" << sequence_string(sequence.sequence) << '\n';
        std::cout << "occurrence_count=" << sequence.occurrence_count << '\n';
        std::cout << "best_theta=" << sequence.best_theta << '\n';
        std::cout << "best_total_delta_v=" << sequence.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << sequence.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << sequence.best_arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds=" << sequence.best_total_flight_time_seconds << '\n';
    }

    std::vector<int> no_beam_indices;
    if (best_valid) {
        no_beam_indices.push_back(best_index);
        for (const auto& theta_run : theta_runs) {
            if (static_cast<int>(no_beam_indices.size()) >= 4) {
                break;
            }
            if (theta_run.index == best_index || !theta_run.run.result.best_terminal_solution.valid) {
                continue;
            }
            no_beam_indices.push_back(theta_run.index);
        }
    } else {
        for (int i = 0; i < std::min(3, theta_scanned_count); ++i) {
            no_beam_indices.push_back(i);
        }
    }
    std::sort(no_beam_indices.begin(), no_beam_indices.end());
    no_beam_indices.erase(std::unique(no_beam_indices.begin(), no_beam_indices.end()), no_beam_indices.end());

    for (const int index : no_beam_indices) {
        auto search_options = base_search_options;
        search_options.enable_beam_pruning = false;
        const auto& beam_theta = theta_runs[static_cast<std::size_t>(index)];
        auto no_beam = run_search(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            beam_theta.candidate.theta,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options);
        if (!no_beam.result.ok) {
            std::cout << "BfsOpenSearchPhysicalBeamComparison\n";
            std::cout << "theta=" << beam_theta.candidate.theta << '\n';
            std::cout << "beam_enabled_best_valid="
                      << (beam_theta.run.result.best_terminal_solution.valid ? 1 : 0) << '\n';
            std::cout << "beam_enabled_best_total_delta_v="
                      << best_total_delta_v_or_inf(beam_theta.run.result) << '\n';
            std::cout << "beam_enabled_node_count=" << beam_theta.run.result.nodes.size() << '\n';
            std::cout << "beam_enabled_max_output_width=" << beam_theta.max_output_width << '\n';
            std::cout << "beam_enabled_total_beam_pruned_count=" << total_beam_pruned(beam_theta.run.result) << '\n';
            std::cout << "beam_disabled_best_valid=0\n";
            std::cout << "beam_disabled_best_total_delta_v=" << std::numeric_limits<double>::infinity() << '\n';
            std::cout << "beam_disabled_node_count=0\n";
            std::cout << "beam_disabled_max_output_width=0\n";
            std::cout << "beam_disabled_wall_time_ms=" << no_beam.wall_time_ms << '\n';
            std::cout << "beam_disabled_exploded=0\n";
            std::cout << "comparison_ok=0\n";
            continue;
        }

        std::cout << "BfsOpenSearchPhysicalBeamComparison\n";
        std::cout << "theta=" << beam_theta.candidate.theta << '\n';
        std::cout << "beam_enabled_best_valid="
                  << (beam_theta.run.result.best_terminal_solution.valid ? 1 : 0) << '\n';
        std::cout << "beam_enabled_best_total_delta_v="
                  << best_total_delta_v_or_inf(beam_theta.run.result) << '\n';
        std::cout << "beam_enabled_node_count=" << beam_theta.run.result.nodes.size() << '\n';
        std::cout << "beam_enabled_max_output_width=" << beam_theta.max_output_width << '\n';
        std::cout << "beam_enabled_total_beam_pruned_count=" << total_beam_pruned(beam_theta.run.result) << '\n';
        std::cout << "beam_disabled_best_valid="
                  << (no_beam.result.best_terminal_solution.valid ? 1 : 0) << '\n';
        std::cout << "beam_disabled_best_total_delta_v=" << best_total_delta_v_or_inf(no_beam.result) << '\n';
        std::cout << "beam_disabled_node_count=" << no_beam.result.nodes.size() << '\n';
        std::cout << "beam_disabled_max_output_width=" << max_output_width(no_beam.result) << '\n';
        std::cout << "beam_disabled_wall_time_ms=" << no_beam.wall_time_ms << '\n';
        std::cout << "beam_disabled_exploded=0\n";
        std::cout << "comparison_ok=1\n";
    }

    assert(theta_scanned_count > 0);
    return 0;
}
