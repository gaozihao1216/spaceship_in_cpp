#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct RunProfile {
    int max_depth = 0;
    spaceship_cpp::bfs::FixedLaunchThetaSearchResult result;
    double wall_time_ms = 0.0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

bool env_flag_enabled(const char* name) {
    if (const char* raw = std::getenv(name)) {
        return std::string(raw) == "1";
    }
    return false;
}

int int_from_env(const char* name, int default_value) {
    if (const char* raw = std::getenv(name)) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return default_value;
}

double double_from_env(const char* name, double default_value) {
    if (const char* raw = std::getenv(name)) {
        const double value = std::atof(raw);
        if (value > 0.0) {
            return value;
        }
    }
    return default_value;
}

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

long long total_attempted_expansions(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.attempted_expansion_count;
    }
    return total;
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

long long total_zero_raw_edge_solves(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.zero_raw_edge_solve_count;
    }
    return total;
}

long long total_zero_accepted_edge_solves(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    long long total = 0;
    for (const auto& stats : result.depth_width_stats) {
        total += stats.zero_accepted_edge_solve_count;
    }
    return total;
}

int max_layer_input_width(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    int width = 0;
    for (const auto& stats : result.depth_width_stats) {
        width = std::max(width, stats.input_state_count);
    }
    return width;
}

int max_layer_output_width(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    int width = 0;
    for (const auto& stats : result.depth_width_stats) {
        width = std::max(width, stats.output_state_count_after_beam);
    }
    return width;
}

double mean_solve_wall_time_ms(const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    double layer_time_ms = 0.0;
    for (const auto& stats : result.depth_width_stats) {
        layer_time_ms += stats.layer_wall_time_ms;
    }
    const long long attempted = total_attempted_expansions(result);
    if (attempted <= 0) {
        return std::numeric_limits<double>::infinity();
    }
    return layer_time_ms / static_cast<double>(attempted);
}

double estimated_depth7_wall_time_ms_if_branching_constant(
    const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result,
    int max_depth,
    double wall_time_ms
) {
    if (max_depth >= 7) {
        return wall_time_ms;
    }
    if (result.depth_width_stats.size() < 2) {
        return std::numeric_limits<double>::infinity();
    }
    const auto& last = result.depth_width_stats.back();
    const double branching = last.input_state_count > 0
        ? static_cast<double>(last.output_state_count_after_beam) / static_cast<double>(last.input_state_count)
        : 0.0;
    if (!std::isfinite(branching) || branching <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return wall_time_ms * std::pow(branching, 7 - max_depth);
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

std::string best_sequence_string(
    spaceship_cpp::planet_params::PlanetId launch_planet,
    const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result
) {
    if (!result.best_terminal_solution.valid) {
        return "";
    }
    return sequence_string(sequence_from_path(
        launch_planet,
        spaceship_cpp::bfs::reconstruct_fixed_launch_theta_path(
            result,
            result.best_terminal_solution.node_index)));
}

double ratio_or_inf(double numerator, double denominator) {
    if (denominator == 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return numerator / denominator;
}

const RunProfile* find_profile(const std::vector<RunProfile>& profiles, int depth) {
    for (const auto& profile : profiles) {
        if (profile.max_depth == depth) {
            return &profile;
        }
    }
    return nullptr;
}

void print_layers(int max_depth, const spaceship_cpp::bfs::FixedLaunchThetaSearchResult& result) {
    for (const auto& stats : result.depth_width_stats) {
        std::cout << "BfsPhysicalNoBeamDepthLayerProfile\n";
        std::cout << "max_depth=" << max_depth << '\n';
        std::cout << "depth=" << stats.depth << '\n';
        std::cout << "input_state_count=" << stats.input_state_count << '\n';
        std::cout << "attempted_expansion_count=" << stats.attempted_expansion_count << '\n';
        std::cout << "raw_edge_count_total=" << stats.raw_edge_count_total << '\n';
        std::cout << "accepted_edge_count_total=" << stats.accepted_edge_count_total << '\n';
        std::cout << "reject_flyby_physical_infeasible_count="
                  << stats.reject_flyby_physical_infeasible_count << '\n';
        std::cout << "output_state_count_before_beam=" << stats.output_state_count_before_beam << '\n';
        std::cout << "output_state_count_after_beam=" << stats.output_state_count_after_beam << '\n';
        std::cout << "zero_raw_edge_solve_count=" << stats.zero_raw_edge_solve_count << '\n';
        std::cout << "zero_accepted_edge_solve_count=" << stats.zero_accepted_edge_solve_count << '\n';
        std::cout << "terminal_solution_count_at_depth=" << stats.terminal_solution_count_at_depth << '\n';
        std::cout << "layer_wall_time_ms=" << stats.layer_wall_time_ms << '\n';
        if (stats.output_state_count_before_beam != stats.output_state_count_after_beam) {
            std::cout << "BfsPhysicalNoBeamDepthLayerWarning\n";
            std::cout << "max_depth=" << max_depth << '\n';
            std::cout << "depth=" << stats.depth << '\n';
            std::cout << "message=before_after_beam_width_mismatch_in_no_beam_profile\n";
        }
    }
}

void print_summary(
    spaceship_cpp::planet_params::PlanetId launch_planet,
    const RunProfile& profile
) {
    const auto& result = profile.result;
    std::cout << "BfsPhysicalNoBeamDepthProfileSummary\n";
    std::cout << "max_depth=" << profile.max_depth << '\n';
    std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "stopped_by_node_limit=" << (result.stopped_by_node_limit ? 1 : 0) << '\n';
    std::cout << "stopped_by_wall_time_limit=" << (result.stopped_by_wall_time_limit ? 1 : 0) << '\n';
    std::cout << "node_count=" << result.nodes.size() << '\n';
    std::cout << "frontier_final_count=" << result.frontier_node_indices.size() << '\n';
    std::cout << "terminal_solution_count=" << result.terminal_solutions.size() << '\n';
    std::cout << "best_valid=" << (result.best_terminal_solution.valid ? 1 : 0) << '\n';
    std::cout << "best_sequence=" << best_sequence_string(launch_planet, result) << '\n';
    std::cout << "best_total_delta_v=" << result.best_terminal_solution.total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_terminal_solution.launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_terminal_solution.arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds="
              << result.best_terminal_solution.total_flight_time_seconds << '\n';
    std::cout << "total_problem2_solve_count=" << total_attempted_expansions(result) << '\n';
    std::cout << "total_raw_edge_count=" << total_raw_edges(result) << '\n';
    std::cout << "total_accepted_edge_count=" << total_accepted_edges(result) << '\n';
    std::cout << "total_flyby_reject_count=" << total_flyby_rejects(result) << '\n';
    std::cout << "total_zero_raw_edge_solve_count=" << total_zero_raw_edge_solves(result) << '\n';
    std::cout << "total_zero_accepted_edge_solve_count=" << total_zero_accepted_edge_solves(result) << '\n';
    std::cout << "max_layer_input_width=" << max_layer_input_width(result) << '\n';
    std::cout << "max_layer_output_width=" << max_layer_output_width(result) << '\n';
    std::cout << "total_wall_time_ms=" << profile.wall_time_ms << '\n';
    std::cout << "mean_expansion_attempt_wall_time_ms=" << mean_solve_wall_time_ms(result) << '\n';
    std::cout << "mean_problem2_solve_wall_time_ms=" << mean_solve_wall_time_ms(result) << '\n';
    std::cout << "estimated_depth7_wall_time_ms_if_branching_constant="
              << estimated_depth7_wall_time_ms_if_branching_constant(
                     result, profile.max_depth, profile.wall_time_ms)
              << '\n';
    if (!result.ok) {
        std::cout << "error_message=" << result.error_message << '\n';
    }
}

void print_growth_comparison(const std::vector<RunProfile>& profiles) {
    const auto value_for_depth = [&](int depth, auto getter) -> double {
        const RunProfile* profile = find_profile(profiles, depth);
        if (!profile) {
            return -1.0;
        }
        return getter(*profile);
    };
    const double d4_nodes = value_for_depth(4, [](const RunProfile& p) {
        return static_cast<double>(p.result.nodes.size());
    });
    const double d5_nodes = value_for_depth(5, [](const RunProfile& p) {
        return static_cast<double>(p.result.nodes.size());
    });
    const double d6_nodes = value_for_depth(6, [](const RunProfile& p) {
        return static_cast<double>(p.result.nodes.size());
    });
    const double d7_nodes = value_for_depth(7, [](const RunProfile& p) {
        return static_cast<double>(p.result.nodes.size());
    });
    const double d4_wall = value_for_depth(4, [](const RunProfile& p) { return p.wall_time_ms; });
    const double d5_wall = value_for_depth(5, [](const RunProfile& p) { return p.wall_time_ms; });
    const double d6_wall = value_for_depth(6, [](const RunProfile& p) { return p.wall_time_ms; });
    const double d7_wall = value_for_depth(7, [](const RunProfile& p) { return p.wall_time_ms; });

    std::cout << "BfsPhysicalNoBeamDepthGrowthComparison\n";
    std::cout << "depth4_node_count=" << d4_nodes << '\n';
    std::cout << "depth5_node_count=" << d5_nodes << '\n';
    std::cout << "depth6_node_count=" << d6_nodes << '\n';
    std::cout << "depth7_node_count=" << d7_nodes << '\n';
    std::cout << "depth4_wall_time_ms=" << d4_wall << '\n';
    std::cout << "depth5_wall_time_ms=" << d5_wall << '\n';
    std::cout << "depth6_wall_time_ms=" << d6_wall << '\n';
    std::cout << "depth7_wall_time_ms=" << d7_wall << '\n';
    std::cout << "node_growth_4_to_5=" << ratio_or_inf(d5_nodes, d4_nodes) << '\n';
    std::cout << "node_growth_5_to_6=" << ratio_or_inf(d6_nodes, d5_nodes) << '\n';
    std::cout << "node_growth_6_to_7=" << (d7_nodes < 0.0 ? -1.0 : ratio_or_inf(d7_nodes, d6_nodes)) << '\n';
    std::cout << "wall_time_growth_4_to_5=" << ratio_or_inf(d5_wall, d4_wall) << '\n';
    std::cout << "wall_time_growth_5_to_6=" << ratio_or_inf(d6_wall, d5_wall) << '\n';
    std::cout << "wall_time_growth_6_to_7=" << (d7_wall < 0.0 ? -1.0 : ratio_or_inf(d7_wall, d6_wall)) << '\n';
    std::cout << "depth_profile_ok=1\n";
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
        std::cout << "bfs_physical_no_beam_depth_profile_skipped_missing_table\n";
        std::cout << "ROOT_TABLE_2DEG_DIR=" << table_dir.string() << '\n';
        return 0;
    }

    const int max_nodes = int_from_env("BFS_NO_BEAM_MAX_NODES", 50000);
    const double max_wall_ms = double_from_env("BFS_NO_BEAM_MAX_WALL_MS", 300000.0);
    std::vector<int> max_depth_values{4, 5, 6};
    if (env_flag_enabled("BFS_NO_BEAM_INCLUDE_DEPTH7")) {
        max_depth_values.push_back(7);
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double launch_time_offset_seconds = 1262332.782636;
    const double launch_time =
        0.17 * planet_params::planet_orbital_period(launch_planet) + launch_time_offset_seconds;
    const double initial_theta = 6.169738905800;

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

    std::vector<RunProfile> profiles;
    for (const int max_depth : max_depth_values) {
        bfs::FixedLaunchThetaSearchOptions search_options{};
        search_options.max_depth = max_depth;
        search_options.enable_beam_pruning = false;
        search_options.beam_width = 3;
        search_options.flyby_physical_options.enabled = true;
        search_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::Enforce;
        search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
        search_options.max_launch_v_inf = 7000.0;
        search_options.time_weight_m_per_s_per_day = 0.0;
        search_options.continue_after_reaching_terminal = true;
        search_options.profile_max_node_count = max_nodes;
        search_options.profile_max_wall_time_ms = max_wall_ms;

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

        profiles.push_back(RunProfile{max_depth, std::move(result), wall_time_ms});
        print_summary(launch_planet, profiles.back());
        print_layers(max_depth, profiles.back().result);
    }

    print_growth_comparison(profiles);

    return 0;
}
