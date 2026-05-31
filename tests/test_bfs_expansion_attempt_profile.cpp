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
#include <map>
#include <sstream>
#include <string>
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

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double mean(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    double total = 0.0;
    for (const double value : values) {
        total += value;
    }
    return total / static_cast<double>(values.size());
}

double percentile(std::vector<double> values, double q) {
    if (values.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    std::sort(values.begin(), values.end());
    const double pos = q * static_cast<double>(values.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    if (lo == hi) {
        return values[lo];
    }
    const double frac = pos - static_cast<double>(lo);
    return values[lo] * (1.0 - frac) + values[hi] * frac;
}

double max_value(const std::vector<double>& values) {
    if (values.empty()) {
        return std::numeric_limits<double>::infinity();
    }
    return *std::max_element(values.begin(), values.end());
}

struct Aggregate {
    int attempt_count = 0;
    int zero_raw_count = 0;
    int zero_accepted_count = 0;
    std::vector<double> outer_wall_ms;
    std::vector<double> problem2_total_ms;
    std::vector<double> nearest_query_count;
    std::vector<double> nearest_query_total_ms;
    std::vector<double> topology_adaptive_ms;
};

void add_attempt(Aggregate* aggregate, const spaceship_cpp::bfs::BfsExpansionAttemptProfile& attempt) {
    aggregate->attempt_count += 1;
    aggregate->zero_raw_count += attempt.zero_raw ? 1 : 0;
    aggregate->zero_accepted_count += attempt.zero_accepted ? 1 : 0;
    aggregate->outer_wall_ms.push_back(attempt.outer_expansion_wall_ms);
    aggregate->problem2_total_ms.push_back(attempt.problem2_solver_profile.total_ms);
    aggregate->nearest_query_count.push_back(
        static_cast<double>(attempt.problem2_solver_profile.problem1_nearest_node_query_count));
    aggregate->nearest_query_total_ms.push_back(attempt.problem2_solver_profile.nearest_query_total_ms);
    aggregate->topology_adaptive_ms.push_back(attempt.problem2_solver_profile.topology_adaptive_sampling_ms);
}

std::string pair_key(
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to
) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::ostringstream os;
    os << planet_params::planet_name(from) << "->" << planet_params::planet_name(to);
    return os.str();
}

void print_attempt(const spaceship_cpp::bfs::BfsExpansionAttemptProfile& attempt) {
    namespace planet_params = spaceship_cpp::planet_params;
    const auto& profile = attempt.problem2_solver_profile;
    std::cout << "BfsExpansionAttemptProfile\n";
    std::cout << "depth=" << attempt.depth << '\n';
    std::cout << "parent_index=" << attempt.parent_index << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(attempt.to_planet) << '\n';
    std::cout << "parent_time=" << attempt.parent_time << '\n';
    std::cout << "parent_incoming_e=" << attempt.parent_incoming_e << '\n';
    std::cout << "parent_incoming_theta=" << attempt.parent_incoming_theta << '\n';
    std::cout << "outer_expansion_wall_ms=" << attempt.outer_expansion_wall_ms << '\n';
    std::cout << "raw_edge_count=" << attempt.raw_edge_count << '\n';
    std::cout << "accepted_edge_count=" << attempt.accepted_edge_count << '\n';
    std::cout << "flyby_reject_count=" << attempt.flyby_reject_count << '\n';
    std::cout << "zero_raw=" << (attempt.zero_raw ? 1 : 0) << '\n';
    std::cout << "zero_accepted=" << (attempt.zero_accepted ? 1 : 0) << '\n';
    std::cout << "problem2_solver_total_ms=" << profile.total_ms << '\n';
    std::cout << "problem2_initial_sampling_ms=" << profile.initial_sampling_ms << '\n';
    std::cout << "problem2_topology_adaptive_sampling_ms=" << profile.topology_adaptive_sampling_ms << '\n';
    std::cout << "problem2_candidate_collection_ms=" << profile.candidate_collection_ms << '\n';
    std::cout << "problem2_midpoint_continuation_ms=" << profile.midpoint_continuation_ms << '\n';
    std::cout << "problem2_bisection_ms=" << profile.bisection_ms << '\n';
    std::cout << "problem2_dedup_ms=" << profile.dedup_ms << '\n';
    std::cout << "problem1_nearest_query_count=" << profile.problem1_nearest_node_query_count << '\n';
    std::cout << "problem1_nearest_query_table_load_ms=" << profile.nearest_query_table_load_ms << '\n';
    std::cout << "problem1_nearest_query_table_offset_build_ms="
              << profile.nearest_query_table_offset_build_ms << '\n';
    std::cout << "problem1_nearest_query_node_read_ms=" << profile.nearest_query_node_read_ms << '\n';
    std::cout << "problem1_nearest_query_refine_ms=" << profile.nearest_query_refine_ms << '\n';
    std::cout << "problem1_nearest_query_derivative_attach_ms="
              << profile.nearest_query_derivative_attach_ms << '\n';
    std::cout << "problem1_nearest_query_total_ms=" << profile.nearest_query_total_ms << '\n';
}

void print_aggregate_fields(const Aggregate& aggregate) {
    std::cout << "attempt_count=" << aggregate.attempt_count << '\n';
    std::cout << "outer_wall_ms_sum=";
    double sum = 0.0;
    for (const double value : aggregate.outer_wall_ms) {
        sum += value;
    }
    std::cout << sum << '\n';
    std::cout << "outer_wall_ms_mean=" << mean(aggregate.outer_wall_ms) << '\n';
    std::cout << "outer_wall_ms_p50=" << percentile(aggregate.outer_wall_ms, 0.50) << '\n';
    std::cout << "outer_wall_ms_p90=" << percentile(aggregate.outer_wall_ms, 0.90) << '\n';
    std::cout << "outer_wall_ms_max=" << max_value(aggregate.outer_wall_ms) << '\n';
    std::cout << "problem2_total_ms_mean=" << mean(aggregate.problem2_total_ms) << '\n';
    std::cout << "nearest_query_count_mean=" << mean(aggregate.nearest_query_count) << '\n';
    std::cout << "nearest_query_total_ms_mean=" << mean(aggregate.nearest_query_total_ms) << '\n';
    std::cout << "topology_adaptive_ms_mean=" << mean(aggregate.topology_adaptive_ms) << '\n';
    std::cout << "zero_raw_count=" << aggregate.zero_raw_count << '\n';
    std::cout << "zero_accepted_count=" << aggregate.zero_accepted_count << '\n';
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
        std::cout << "bfs_expansion_attempt_profile_skipped_missing_table\n";
        std::cout << "ROOT_TABLE_2DEG_DIR=" << table_dir.string() << '\n';
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const auto terminal_planet = planet_params::PlanetId::Mercury;
    const double launch_time =
        0.17 * planet_params::planet_orbital_period(launch_planet) + 1262332.782636;
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

    bfs::FixedLaunchThetaSearchOptions search_options{};
    search_options.max_depth = 5;
    search_options.enable_beam_pruning = false;
    search_options.beam_width = 3;
    search_options.flyby_physical_options.enabled = true;
    search_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::Enforce;
    search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
    search_options.max_launch_v_inf = 7000.0;
    search_options.time_weight_m_per_s_per_day = 0.0;
    search_options.continue_after_reaching_terminal = true;

    const auto start = Clock::now();
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
    const double total_wall_ms = elapsed_ms(start, Clock::now());

    Aggregate global;
    std::map<std::string, Aggregate> by_pair;
    std::map<std::string, std::pair<planet_params::PlanetId, planet_params::PlanetId>> pair_planets;
    for (const auto& attempt : result.expansion_attempt_profiles) {
        print_attempt(attempt);
        add_attempt(&global, attempt);
        const auto key = pair_key(attempt.from_planet, attempt.to_planet);
        add_attempt(&by_pair[key], attempt);
        pair_planets[key] = {attempt.from_planet, attempt.to_planet};
    }

    for (const auto& [key, aggregate] : by_pair) {
        const auto planets = pair_planets[key];
        std::cout << "BfsExpansionAttemptProfileByPair\n";
        std::cout << "from_planet=" << planet_params::planet_name(planets.first) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(planets.second) << '\n';
        print_aggregate_fields(aggregate);
    }

    std::vector<double> overhead;
    for (const auto& attempt : result.expansion_attempt_profiles) {
        overhead.push_back(attempt.outer_expansion_wall_ms - attempt.problem2_solver_profile.total_ms);
    }

    std::cout << "BfsExpansionAttemptProfileSummary\n";
    std::cout << "attempt_count=" << global.attempt_count << '\n';
    std::cout << "outer_wall_ms_mean=" << mean(global.outer_wall_ms) << '\n';
    std::cout << "outer_wall_ms_p50=" << percentile(global.outer_wall_ms, 0.50) << '\n';
    std::cout << "outer_wall_ms_p90=" << percentile(global.outer_wall_ms, 0.90) << '\n';
    std::cout << "outer_wall_ms_max=" << max_value(global.outer_wall_ms) << '\n';
    std::cout << "problem2_total_ms_mean=" << mean(global.problem2_total_ms) << '\n';
    std::cout << "problem2_total_ms_p50=" << percentile(global.problem2_total_ms, 0.50) << '\n';
    std::cout << "problem2_total_ms_p90=" << percentile(global.problem2_total_ms, 0.90) << '\n';
    std::cout << "problem2_total_ms_max=" << max_value(global.problem2_total_ms) << '\n';
    std::cout << "nearest_query_total_ms_mean=" << mean(global.nearest_query_total_ms) << '\n';
    std::cout << "topology_adaptive_ms_mean=" << mean(global.topology_adaptive_ms) << '\n';
    std::cout << "bfs_overhead_estimate_ms_mean=" << mean(overhead) << '\n';
    std::cout << "total_search_wall_ms=" << total_wall_ms << '\n';
    std::cout << "search_ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "profile_ok=1\n";
    return result.ok ? 0 : 1;
}
