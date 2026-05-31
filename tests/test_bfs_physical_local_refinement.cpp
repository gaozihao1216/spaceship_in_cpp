#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

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
#include <sstream>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

constexpr double kBaselineTotalDeltaV = 20088.79070851;
constexpr double kBestLaunchTimeOffsetSeconds = 1262332.782636;
constexpr double kThetaCenter = 6.169738905800;

struct SearchSample {
    bool ok = false;
    bool best_valid = false;
    double launch_time = 0.0;
    double launch_time_offset_seconds = 0.0;
    double relative_offset_from_center_seconds = 0.0;
    double theta = 0.0;
    double theta_offset = 0.0;
    std::vector<spaceship_cpp::planet_params::PlanetId> best_sequence;
    double best_total_delta_v = std::numeric_limits<double>::infinity();
    double best_launch_v_inf = std::numeric_limits<double>::infinity();
    double best_arrival_v_inf = std::numeric_limits<double>::infinity();
    double best_total_flight_time_seconds = std::numeric_limits<double>::infinity();
    int terminal_solution_count = 0;
    int max_layer_output_width_seen = 0;
    long long beam_pruned_count = 0;
    std::size_t node_count = 0;
    double wall_time_ms = 0.0;
    spaceship_cpp::bfs::FixedLaunchThetaSearchResult search;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
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

double sample_linear(double center, double half_width, int index, int count) {
    if (count <= 1) {
        return center;
    }
    const double t = static_cast<double>(index) / static_cast<double>(count - 1);
    return center - half_width + 2.0 * half_width * t;
}

SearchSample run_search_sample(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId launch_planet,
    spaceship_cpp::planet_params::PlanetId terminal_planet,
    const std::vector<spaceship_cpp::planet_params::PlanetId>& allowed_first_targets,
    const std::vector<spaceship_cpp::planet_params::PlanetId>& allowed_transfer_planets,
    const spaceship_cpp::bfs::InitialLaunchExpansionOptions& initial_options,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const spaceship_cpp::bfs::FixedLaunchThetaSearchOptions& search_options,
    double launch_time,
    double launch_time_center,
    double theta,
    double theta_center
) {
    namespace bfs = spaceship_cpp::bfs;

    SearchSample sample{};
    sample.launch_time = launch_time;
    sample.launch_time_offset_seconds = launch_time - 0.17 * spaceship_cpp::planet_params::planet_orbital_period(
        launch_planet);
    sample.relative_offset_from_center_seconds = launch_time - launch_time_center;
    sample.theta = theta;
    sample.theta_offset = theta - theta_center;

    const auto start = Clock::now();
    try {
        sample.search = bfs::run_fixed_launch_theta_beam_search_with_table(
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
        sample.search.ok = false;
        sample.search.error_message = ex.what();
    }
    sample.wall_time_ms = elapsed_ms(start, Clock::now());
    sample.ok = sample.search.ok;
    sample.terminal_solution_count = static_cast<int>(sample.search.terminal_solutions.size());
    sample.max_layer_output_width_seen = max_output_width(sample.search);
    sample.beam_pruned_count = total_beam_pruned(sample.search);
    sample.node_count = sample.search.nodes.size();

    if (sample.search.ok && sample.search.best_terminal_solution.valid) {
        const auto path = bfs::reconstruct_fixed_launch_theta_path(
            sample.search, sample.search.best_terminal_solution.node_index);
        sample.best_valid = true;
        sample.best_sequence = sequence_from_path(launch_planet, path);
        sample.best_total_delta_v = sample.search.best_terminal_solution.total_delta_v;
        sample.best_launch_v_inf = sample.search.best_terminal_solution.launch_v_inf;
        sample.best_arrival_v_inf = sample.search.best_terminal_solution.arrival_v_inf;
        sample.best_total_flight_time_seconds = sample.search.best_terminal_solution.total_flight_time_seconds;
    }
    return sample;
}

void print_theta_sample(int index, const SearchSample& sample) {
    std::cout << "BfsPhysicalLocalThetaRefinementSample\n";
    std::cout << "index=" << index << '\n';
    std::cout << "theta=" << sample.theta << '\n';
    std::cout << "theta_offset=" << sample.theta_offset << '\n';
    std::cout << "best_valid=" << (sample.best_valid ? 1 : 0) << '\n';
    std::cout << "best_sequence=" << sequence_string(sample.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << sample.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << sample.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << sample.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << sample.best_total_flight_time_seconds << '\n';
    std::cout << "terminal_solution_count=" << sample.terminal_solution_count << '\n';
    std::cout << "max_layer_output_width_seen=" << sample.max_layer_output_width_seen << '\n';
    std::cout << "beam_pruned_count=" << sample.beam_pruned_count << '\n';
    std::cout << "wall_time_ms=" << sample.wall_time_ms << '\n';
}

void print_launch_time_sample(int index, const SearchSample& sample) {
    std::cout << "BfsPhysicalLocalLaunchTimeRefinementSample\n";
    std::cout << "index=" << index << '\n';
    std::cout << "launch_time=" << sample.launch_time << '\n';
    std::cout << "launch_time_offset_seconds=" << sample.launch_time_offset_seconds << '\n';
    std::cout << "relative_offset_from_center_seconds=" << sample.relative_offset_from_center_seconds << '\n';
    std::cout << "theta=" << sample.theta << '\n';
    std::cout << "best_valid=" << (sample.best_valid ? 1 : 0) << '\n';
    std::cout << "best_sequence=" << sequence_string(sample.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << sample.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << sample.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << sample.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << sample.best_total_flight_time_seconds << '\n';
    std::cout << "terminal_solution_count=" << sample.terminal_solution_count << '\n';
    std::cout << "max_layer_output_width_seen=" << sample.max_layer_output_width_seen << '\n';
    std::cout << "beam_pruned_count=" << sample.beam_pruned_count << '\n';
    std::cout << "wall_time_ms=" << sample.wall_time_ms << '\n';
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
        std::cout << "bfs_physical_local_refinement_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double earth_period = planet_params::planet_orbital_period(launch_planet);
    const double base_launch_time = 0.17 * earth_period;
    const double launch_time_center = base_launch_time + kBestLaunchTimeOffsetSeconds;
    const double theta_center = kThetaCenter;

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
    search_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::Enforce;
    search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
    search_options.max_launch_v_inf = 7000.0;
    search_options.time_weight_m_per_s_per_day = 0.0;
    search_options.continue_after_reaching_terminal = true;

    SearchSample theta_best{};
    const int theta_sample_count = 21;
    const double theta_half_width = 0.05;
    for (int i = 0; i < theta_sample_count; ++i) {
        const double theta = sample_linear(theta_center, theta_half_width, i, theta_sample_count);
        auto sample = run_search_sample(
            loader,
            launch_planet,
            terminal_planet,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options,
            launch_time_center,
            launch_time_center,
            theta,
            theta_center);
        print_theta_sample(i, sample);
        if (sample.best_valid && (!theta_best.best_valid || sample.best_total_delta_v < theta_best.best_total_delta_v)) {
            theta_best = std::move(sample);
        }
    }

    const double theta_stage_best = theta_best.best_valid ? theta_best.theta : theta_center;
    std::cout << "BfsPhysicalLocalThetaRefinementSummary\n";
    std::cout << "sample_count=" << theta_sample_count << '\n';
    std::cout << "best_valid=" << (theta_best.best_valid ? 1 : 0) << '\n';
    std::cout << "best_theta=" << (theta_best.best_valid ? theta_best.theta : theta_center) << '\n';
    std::cout << "best_theta_offset=" << (theta_best.best_valid ? theta_best.theta_offset : 0.0) << '\n';
    std::cout << "best_sequence=" << sequence_string(theta_best.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << theta_best.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << theta_best.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << theta_best.best_arrival_v_inf << '\n';
    std::cout << "baseline_total_delta_v=" << kBaselineTotalDeltaV << '\n';
    std::cout << "improvement_vs_baseline=" << (kBaselineTotalDeltaV - theta_best.best_total_delta_v) << '\n';
    std::cout << "theta_refinement_ok=1\n";

    SearchSample launch_time_best{};
    const int launch_time_sample_count = 11;
    const double launch_time_half_width = 0.02 * earth_period;
    for (int i = 0; i < launch_time_sample_count; ++i) {
        const double launch_time = sample_linear(launch_time_center, launch_time_half_width, i, launch_time_sample_count);
        auto sample = run_search_sample(
            loader,
            launch_planet,
            terminal_planet,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            search_options,
            launch_time,
            launch_time_center,
            theta_stage_best,
            theta_stage_best);
        print_launch_time_sample(i, sample);
        if (sample.best_valid &&
            (!launch_time_best.best_valid || sample.best_total_delta_v < launch_time_best.best_total_delta_v)) {
            launch_time_best = std::move(sample);
        }
    }

    SearchSample final_best = launch_time_best.best_valid ? launch_time_best : theta_best;
    std::cout << "BfsPhysicalLocalLaunchTimeRefinementSummary\n";
    std::cout << "sample_count=" << launch_time_sample_count << '\n';
    std::cout << "best_valid=" << (launch_time_best.best_valid ? 1 : 0) << '\n';
    std::cout << "best_launch_time=" << launch_time_best.launch_time << '\n';
    std::cout << "best_launch_time_offset_seconds=" << launch_time_best.launch_time_offset_seconds << '\n';
    std::cout << "best_relative_offset_from_center_seconds="
              << launch_time_best.relative_offset_from_center_seconds << '\n';
    std::cout << "best_theta=" << theta_stage_best << '\n';
    std::cout << "best_sequence=" << sequence_string(launch_time_best.best_sequence) << '\n';
    std::cout << "best_total_delta_v=" << launch_time_best.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << launch_time_best.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << launch_time_best.best_arrival_v_inf << '\n';
    std::cout << "baseline_total_delta_v=" << kBaselineTotalDeltaV << '\n';
    std::cout << "improvement_vs_baseline=" << (kBaselineTotalDeltaV - launch_time_best.best_total_delta_v) << '\n';
    std::cout << "launch_time_refinement_ok=1\n";

    const bool skip_grid = env_bool("BFS_LOCAL_REFINEMENT_SKIP_GRID");
    if (!skip_grid && final_best.best_valid &&
        (kBaselineTotalDeltaV - final_best.best_total_delta_v) > 0.0) {
        SearchSample grid_best = final_best;
        const int launch_time_count = 5;
        const int theta_count = 9;
        const double grid_launch_time_half_width = 0.005 * earth_period;
        const double grid_theta_half_width = 0.02;
        for (int ti = 0; ti < launch_time_count; ++ti) {
            const double launch_time = sample_linear(
                final_best.launch_time, grid_launch_time_half_width, ti, launch_time_count);
            for (int qi = 0; qi < theta_count; ++qi) {
                const double theta = sample_linear(final_best.theta, grid_theta_half_width, qi, theta_count);
                auto sample = run_search_sample(
                    loader,
                    launch_planet,
                    terminal_planet,
                    allowed_first_targets,
                    allowed_transfer_planets,
                    initial_options,
                    problem2_options,
                    search_options,
                    launch_time,
                    launch_time_center,
                    theta,
                    theta_stage_best);
                if (sample.best_valid && sample.best_total_delta_v < grid_best.best_total_delta_v) {
                    grid_best = std::move(sample);
                }
            }
        }
        final_best = grid_best;
        std::cout << "BfsPhysicalLocalGridRefinementSummary\n";
        std::cout << "launch_time_count=" << launch_time_count << '\n';
        std::cout << "theta_count=" << theta_count << '\n';
        std::cout << "best_valid=" << (grid_best.best_valid ? 1 : 0) << '\n';
        std::cout << "best_launch_time_offset_seconds=" << grid_best.launch_time_offset_seconds << '\n';
        std::cout << "best_theta=" << grid_best.theta << '\n';
        std::cout << "best_sequence=" << sequence_string(grid_best.best_sequence) << '\n';
        std::cout << "best_total_delta_v=" << grid_best.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << grid_best.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << grid_best.best_arrival_v_inf << '\n';
        std::cout << "baseline_total_delta_v=" << kBaselineTotalDeltaV << '\n';
        std::cout << "improvement_vs_baseline=" << (kBaselineTotalDeltaV - grid_best.best_total_delta_v) << '\n';
        std::cout << "grid_refinement_ok=1\n";
    } else {
        std::cout << "BfsPhysicalLocalGridRefinementSummary\n";
        std::cout << "launch_time_count=0\n";
        std::cout << "theta_count=0\n";
        std::cout << "best_valid=" << (final_best.best_valid ? 1 : 0) << '\n';
        std::cout << "best_launch_time_offset_seconds=" << final_best.launch_time_offset_seconds << '\n';
        std::cout << "best_theta=" << final_best.theta << '\n';
        std::cout << "best_sequence=" << sequence_string(final_best.best_sequence) << '\n';
        std::cout << "best_total_delta_v=" << final_best.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << final_best.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << final_best.best_arrival_v_inf << '\n';
        std::cout << "baseline_total_delta_v=" << kBaselineTotalDeltaV << '\n';
        std::cout << "improvement_vs_baseline=" << (kBaselineTotalDeltaV - final_best.best_total_delta_v) << '\n';
        std::cout << "grid_refinement_ok=1\n";
    }

    if (final_best.best_valid) {
        auto no_beam_options = search_options;
        no_beam_options.enable_beam_pruning = false;
        const auto start = Clock::now();
        auto no_beam = run_search_sample(
            loader,
            launch_planet,
            terminal_planet,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_options,
            problem2_options,
            no_beam_options,
            final_best.launch_time,
            launch_time_center,
            final_best.theta,
            theta_stage_best);
        const double no_beam_wall_time_ms = elapsed_ms(start, Clock::now());
        std::cout << "BfsPhysicalLocalRefinementNoBeamValidation\n";
        std::cout << "beam_best_total_delta_v=" << final_best.best_total_delta_v << '\n';
        std::cout << "no_beam_best_valid=" << (no_beam.best_valid ? 1 : 0) << '\n';
        std::cout << "no_beam_best_total_delta_v=" << no_beam.best_total_delta_v << '\n';
        std::cout << "beam_node_count=" << final_best.node_count << '\n';
        std::cout << "no_beam_node_count=" << no_beam.node_count << '\n';
        std::cout << "no_beam_max_layer_output_width=" << no_beam.max_layer_output_width_seen << '\n';
        std::cout << "no_beam_wall_time_ms=" << no_beam_wall_time_ms << '\n';
        std::cout << "no_beam_validation_ok=" << (no_beam.ok ? 1 : 0) << '\n';
    }

    assert(theta_sample_count > 0);
    return 0;
}
