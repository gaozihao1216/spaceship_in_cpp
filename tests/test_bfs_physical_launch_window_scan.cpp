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
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct TimeSummary {
    int index = 0;
    double launch_time = 0.0;
    double launch_time_offset_seconds = 0.0;
    int theta_candidate_count = 0;
    int theta_scanned_count = 0;
    int valid_interval_count = 0;
    int near_valid_interval_count = 0;
    int terminal_theta_count = 0;
    int total_terminal_solution_count = 0;
    bool best_valid = false;
    double best_theta = 0.0;
    std::vector<spaceship_cpp::planet_params::PlanetId> best_sequence;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();
    long long total_raw_edge_count = 0;
    long long total_accepted_edge_count = 0;
    long long total_flyby_reject_count = 0;
    long long total_beam_pruned_count = 0;
    int max_layer_output_width_seen = 0;
    double wall_time_ms = 0.0;
    spaceship_cpp::bfs::FixedLaunchThetaSearchResult best_search;
};

struct SequenceAggregate {
    std::vector<spaceship_cpp::planet_params::PlanetId> sequence;
    int occurrence_count = 0;
    double best_launch_time = 0.0;
    double best_launch_time_offset_seconds = 0.0;
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

int env_int(const char* name, int default_value) {
    if (const char* raw = std::getenv(name)) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return default_value;
}

double env_double(const char* name, double default_value) {
    if (const char* raw = std::getenv(name)) {
        char* end = nullptr;
        const double value = std::strtod(raw, &end);
        if (end != raw && std::isfinite(value) && value > 0.0) {
            return value;
        }
    }
    return default_value;
}

bool env_bool(const char* name) {
    if (const char* raw = std::getenv(name)) {
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

std::vector<spaceship_cpp::bfs::ThetaCandidate> select_stratified_theta_candidates(
    std::vector<spaceship_cpp::bfs::ThetaCandidate> candidates,
    int max_count,
    bool include_near_valid
) {
    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(), [include_near_valid](const auto& candidate) {
            const bool usable_source = candidate.from_valid_interval ||
                (include_near_valid && candidate.from_near_valid_interval);
            return !candidate.valid || !usable_source || !std::isfinite(candidate.theta);
        }),
        candidates.end());
    std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.theta < rhs.theta;
    });
    if (static_cast<int>(candidates.size()) <= max_count) {
        return candidates;
    }

    std::vector<spaceship_cpp::bfs::ThetaCandidate> selected;
    selected.reserve(static_cast<std::size_t>(max_count));
    for (int i = 0; i < max_count; ++i) {
        const double position =
            static_cast<double>(i) * static_cast<double>(candidates.size() - 1) / static_cast<double>(max_count - 1);
        const std::size_t index = static_cast<std::size_t>(std::llround(position));
        selected.push_back(candidates[index]);
    }
    selected.erase(
        std::unique(selected.begin(), selected.end(), [](const auto& lhs, const auto& rhs) {
            return std::abs(lhs.theta - rhs.theta) < 1e-10;
        }),
        selected.end());
    return selected;
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

void print_time_summary(const TimeSummary& summary) {
    std::cout << "BfsPhysicalLaunchWindowTimeSummary\n";
    std::cout << "index=" << summary.index << '\n';
    std::cout << "launch_time=" << summary.launch_time << '\n';
    std::cout << "launch_time_offset_seconds=" << summary.launch_time_offset_seconds << '\n';
    std::cout << "theta_candidate_count=" << summary.theta_candidate_count << '\n';
    std::cout << "theta_scanned_count=" << summary.theta_scanned_count << '\n';
    std::cout << "valid_interval_count=" << summary.valid_interval_count << '\n';
    std::cout << "near_valid_interval_count=" << summary.near_valid_interval_count << '\n';
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
    std::cout << "wall_time_ms=" << summary.wall_time_ms << '\n';
    std::cout.flush();
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
        std::cout << "bfs_physical_launch_window_scan_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double earth_period = planet_params::planet_orbital_period(launch_planet);
    const double launch_time_center = 0.17 * earth_period;
    const int launch_time_sample_count = env_int("BFS_LAUNCH_WINDOW_SAMPLE_COUNT", 21);
    const double half_width_period = env_double("BFS_LAUNCH_WINDOW_HALF_WIDTH_PERIOD", 0.10);
    const double launch_time_half_width = half_width_period * earth_period;
    const int theta_scan_max_per_time = env_int("BFS_LAUNCH_WINDOW_THETA_MAX", 20);
    const bool include_near_valid = env_bool("BFS_LAUNCH_WINDOW_INCLUDE_NEAR_VALID");

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
    problem2_options.max_transfer_revolution = 1;
    problem2_options.max_target_revolution = 1;

    bfs::FixedLaunchThetaSearchOptions search_options{};
    search_options.max_depth = 5;
    search_options.beam_width = 3;
    search_options.enable_beam_pruning = true;
    search_options.flyby_physical_options.enabled = true;
    search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
    search_options.max_launch_v_inf = 7000.0;
    search_options.time_weight_m_per_s_per_day = 0.0;
    search_options.continue_after_reaching_terminal = true;

    std::vector<TimeSummary> time_summaries;
    time_summaries.reserve(static_cast<std::size_t>(launch_time_sample_count));
    std::map<std::string, SequenceAggregate> sequence_aggregates;
    TimeSummary global_best{};

    int time_with_valid_theta_count = 0;
    int time_with_terminal_count = 0;
    int total_theta_scanned_count = 0;
    int total_terminal_solution_count = 0;
    long long global_total_raw_edges = 0;
    long long global_total_accepted_edges = 0;
    long long global_total_flyby_rejects = 0;
    long long global_total_beam_pruned = 0;
    int global_max_layer_output_width = 0;
    double total_wall_time_ms = 0.0;

    for (int time_index = 0; time_index < launch_time_sample_count; ++time_index) {
        const double fraction = launch_time_sample_count == 1
            ? 0.5
            : static_cast<double>(time_index) / static_cast<double>(launch_time_sample_count - 1);
        const double launch_time = launch_time_center - launch_time_half_width + 2.0 * launch_time_half_width * fraction;
        const auto time_start = Clock::now();

        bfs::ThetaLaunchFeasibilityScoutOptions scout_options{};
        scout_options.theta_scout_count = 720;
        scout_options.max_launch_v_inf = 7000.0;
        scout_options.near_v_inf_buffer = 1000.0;
        scout_options.refine_interval_boundaries = true;
        scout_options.max_transfer_revolution = 1;
        scout_options.max_target_revolution = 1;
        const auto scout = bfs::scout_theta_launch_feasibility_with_problem1_table(
            loader, launch_planet, launch_time, allowed_first_targets, scout_options);
        assert(scout.ok);

        bfs::ThetaCandidateSelectorOptions selector_options{};
        selector_options.max_theta_candidates = 30;
        selector_options.min_samples_per_valid_interval = 5;
        selector_options.max_samples_per_valid_interval = 20;
        selector_options.include_interval_boundaries = true;
        selector_options.include_interval_midpoint = true;
        selector_options.include_interval_minimum = true;
        selector_options.include_near_valid_intervals = include_near_valid;
        const auto selection = bfs::select_theta_candidates_from_launch_scout(scout, selector_options);
        assert(selection.ok);

        auto theta_candidates = select_stratified_theta_candidates(
            selection.candidates, theta_scan_max_per_time, include_near_valid);

        TimeSummary summary{};
        summary.index = time_index;
        summary.launch_time = launch_time;
        summary.launch_time_offset_seconds = launch_time - launch_time_center;
        summary.theta_candidate_count = static_cast<int>(selection.candidates.size());
        summary.theta_scanned_count = static_cast<int>(theta_candidates.size());
        summary.valid_interval_count = static_cast<int>(scout.valid_intervals.size());
        summary.near_valid_interval_count = static_cast<int>(scout.near_valid_intervals.size());
        if (!theta_candidates.empty()) {
            time_with_valid_theta_count += 1;
        }

        for (const auto& candidate : theta_candidates) {
            const auto search_start = Clock::now();
            bfs::FixedLaunchThetaSearchResult search;
            try {
                search = bfs::run_fixed_launch_theta_beam_search_with_table(
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
            } catch (const std::exception& ex) {
                search.ok = false;
                search.error_message = ex.what();
            }
            const double search_wall_time_ms = elapsed_ms(search_start, Clock::now());
            if (!search.ok) {
                (void)search_wall_time_ms;
                continue;
            }
            (void)search_wall_time_ms;

            (void)search_wall_time_ms;
            summary.total_raw_edge_count += total_raw_edges(search);
            summary.total_accepted_edge_count += total_accepted_edges(search);
            summary.total_flyby_reject_count += total_flyby_rejects(search);
            summary.total_beam_pruned_count += total_beam_pruned(search);
            summary.max_layer_output_width_seen =
                std::max(summary.max_layer_output_width_seen, max_output_width(search));
            summary.total_terminal_solution_count += static_cast<int>(search.terminal_solutions.size());

            if (search.best_terminal_solution.valid) {
                summary.terminal_theta_count += 1;
                for (const auto& terminal : search.terminal_solutions) {
                    if (!terminal.valid) {
                        continue;
                    }
                    const auto path = bfs::reconstruct_fixed_launch_theta_path(search, terminal.node_index);
                    const auto sequence = sequence_from_path(launch_planet, path);
                    auto& aggregate = sequence_aggregates[sequence_string(sequence)];
                    if (aggregate.occurrence_count == 0) {
                        aggregate.sequence = sequence;
                    }
                    aggregate.occurrence_count += 1;
                    if (terminal.total_delta_v < aggregate.best_total_delta_v) {
                        aggregate.best_launch_time = launch_time;
                        aggregate.best_launch_time_offset_seconds = summary.launch_time_offset_seconds;
                        aggregate.best_theta = candidate.theta;
                        aggregate.best_total_delta_v = terminal.total_delta_v;
                        aggregate.best_launch_v_inf = terminal.launch_v_inf;
                        aggregate.best_arrival_v_inf = terminal.arrival_v_inf;
                        aggregate.best_total_flight_time_seconds = terminal.total_flight_time_seconds;
                    }
                }

                if (!summary.best_valid ||
                    search.best_terminal_solution.total_delta_v < summary.best_total_delta_v) {
                    const auto path = bfs::reconstruct_fixed_launch_theta_path(
                        search, search.best_terminal_solution.node_index);
                    summary.best_valid = true;
                    summary.best_theta = candidate.theta;
                    summary.best_sequence = sequence_from_path(launch_planet, path);
                    summary.best_total_delta_v = search.best_terminal_solution.total_delta_v;
                    summary.best_launch_v_inf = search.best_terminal_solution.launch_v_inf;
                    summary.best_arrival_v_inf = search.best_terminal_solution.arrival_v_inf;
                    summary.best_total_flight_time_seconds =
                        search.best_terminal_solution.total_flight_time_seconds;
                    summary.best_search = search;
                }
            }
        }

        summary.wall_time_ms = elapsed_ms(time_start, Clock::now());
        print_time_summary(summary);

        if (summary.best_valid) {
            time_with_terminal_count += 1;
            if (!global_best.best_valid || summary.best_total_delta_v < global_best.best_total_delta_v) {
                global_best = summary;
            }
        }
        total_theta_scanned_count += summary.theta_scanned_count;
        total_terminal_solution_count += summary.total_terminal_solution_count;
        global_total_raw_edges += summary.total_raw_edge_count;
        global_total_accepted_edges += summary.total_accepted_edge_count;
        global_total_flyby_rejects += summary.total_flyby_reject_count;
        global_total_beam_pruned += summary.total_beam_pruned_count;
        global_max_layer_output_width =
            std::max(global_max_layer_output_width, summary.max_layer_output_width_seen);
        total_wall_time_ms += summary.wall_time_ms;
        time_summaries.push_back(std::move(summary));
    }

    std::cout << "BfsPhysicalLaunchWindowScanSummary\n";
    std::cout << "launch_time_sample_count=" << launch_time_sample_count << '\n';
    std::cout << "launch_time_center=" << launch_time_center << '\n';
    std::cout << "launch_time_half_width_seconds=" << launch_time_half_width << '\n';
    std::cout << "theta_scan_max_per_time=" << theta_scan_max_per_time << '\n';
    std::cout << "time_with_valid_theta_count=" << time_with_valid_theta_count << '\n';
    std::cout << "time_with_terminal_count=" << time_with_terminal_count << '\n';
    std::cout << "total_theta_scanned_count=" << total_theta_scanned_count << '\n';
    std::cout << "total_terminal_solution_count=" << total_terminal_solution_count << '\n';
    std::cout << "best_valid=" << (global_best.best_valid ? 1 : 0) << '\n';
    std::cout << "best_launch_time=" << global_best.launch_time << '\n';
    std::cout << "best_launch_time_offset_seconds=" << global_best.launch_time_offset_seconds << '\n';
    std::cout << "best_theta=" << global_best.best_theta << '\n';
    std::cout << "best_sequence=" << sequence_string(global_best.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << global_best.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << global_best.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << global_best.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << global_best.best_total_flight_time_seconds << '\n';
    std::cout << "global_max_layer_output_width_seen=" << global_max_layer_output_width << '\n';
    std::cout << "total_raw_edge_count=" << global_total_raw_edges << '\n';
    std::cout << "total_accepted_edge_count=" << global_total_accepted_edges << '\n';
    std::cout << "total_flyby_reject_count=" << global_total_flyby_rejects << '\n';
    std::cout << "total_beam_pruned_count=" << global_total_beam_pruned << '\n';
    std::cout << "total_wall_time_ms=" << total_wall_time_ms << '\n';
    std::cout << "physical_launch_window_scan_ok=1\n";

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
    const std::size_t top_count = std::min<std::size_t>(10, top_sequences.size());
    for (std::size_t i = 0; i < top_count; ++i) {
        const auto& sequence = top_sequences[i];
        std::cout << "BfsPhysicalLaunchWindowTopSequence\n";
        std::cout << "rank=" << i << '\n';
        std::cout << "sequence=" << sequence_string(sequence.sequence) << '\n';
        std::cout << "occurrence_count=" << sequence.occurrence_count << '\n';
        std::cout << "best_launch_time=" << sequence.best_launch_time << '\n';
        std::cout << "best_launch_time_offset_seconds=" << sequence.best_launch_time_offset_seconds << '\n';
        std::cout << "best_theta=" << sequence.best_theta << '\n';
        std::cout << "best_total_delta_v=" << sequence.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << sequence.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << sequence.best_arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds=" << sequence.best_total_flight_time_seconds << '\n';
    }

    if (global_best.best_valid) {
        auto no_beam_options = search_options;
        no_beam_options.enable_beam_pruning = false;
        const auto start = Clock::now();
        bfs::FixedLaunchThetaSearchResult no_beam;
        try {
            no_beam = bfs::run_fixed_launch_theta_beam_search_with_table(
                loader,
                launch_planet,
                terminal_planet,
                global_best.launch_time,
                global_best.best_theta,
                allowed_first_targets,
                allowed_transfer_planets,
                initial_options,
                problem2_options,
                no_beam_options);
        } catch (const std::exception& ex) {
            no_beam.ok = false;
            no_beam.error_message = ex.what();
        }
        const double wall_time_ms = elapsed_ms(start, Clock::now());
        std::cout << "BfsPhysicalLaunchWindowNoBeamValidation\n";
        std::cout << "launch_time=" << global_best.launch_time << '\n';
        std::cout << "launch_time_offset_seconds=" << global_best.launch_time_offset_seconds << '\n';
        std::cout << "theta=" << global_best.best_theta << '\n';
        std::cout << "beam_best_total_delta_v=" << global_best.best_total_delta_v << '\n';
        std::cout << "no_beam_best_valid=" << (no_beam.ok && no_beam.best_terminal_solution.valid ? 1 : 0) << '\n';
        std::cout << "no_beam_best_total_delta_v="
                  << (no_beam.ok && no_beam.best_terminal_solution.valid
                          ? no_beam.best_terminal_solution.total_delta_v
                          : std::numeric_limits<double>::infinity()) << '\n';
        std::cout << "beam_node_count=" << global_best.best_search.nodes.size() << '\n';
        std::cout << "no_beam_node_count=" << no_beam.nodes.size() << '\n';
        std::cout << "no_beam_max_layer_output_width=" << max_output_width(no_beam) << '\n';
        std::cout << "no_beam_wall_time_ms=" << wall_time_ms << '\n';
        std::cout << "no_beam_validation_ok=" << (no_beam.ok ? 1 : 0) << '\n';
    }

    assert(launch_time_sample_count > 0);
    return 0;
}
