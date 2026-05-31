#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct QuerySample {
    int source_attempt_index = -1;
    spaceship_cpp::planet_params::PlanetId from_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    double parent_time = 0.0;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
    bool adaptive = false;
};

struct RecallStats {
    int query_count = 0;
    int baseline_positive_query_count = 0;
    int query_level_miss_count = 0;
    int baseline_branch_count_total = 0;
    int missing_branch_count_total = 0;
    double baseline_ms_sum = 0.0;
    double scout_ms_sum = 0.0;
    int baseline_direct_fallback_query_count = 0;
    int scout_direct_fallback_skipped_count = 0;
    int scout_cell_vertex_success_count = 0;
    int scout_cell_vertex_failure_count = 0;
};

struct QueryComparison {
    int matched = 0;
    int missing = 0;
    int added = 0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

int max_queries_from_env() {
    if (const char* raw = std::getenv("SCOUT_QUERY_RECALL_MAX_QUERIES")) {
        if (*raw != '\0') {
            return std::max(1, std::atoi(raw));
        }
    }
    return std::numeric_limits<int>::max();
}

double angle_diff(double a, double b) {
    return std::abs(spaceship_cpp::common::normalize_angle_minus_pi_pi(a - b));
}

bool branch_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& a,
    const spaceship_cpp::problem1::Problem1SolutionBranch& b
) {
    return a.valid && b.valid &&
        a.transfer_revolution == b.transfer_revolution &&
        a.target_revolution == b.target_revolution &&
        angle_diff(a.encounter_global_angle, b.encounter_global_angle) <= 1e-8 &&
        std::abs(a.time_of_flight_seconds - b.time_of_flight_seconds) <= 1e-3;
}

QueryComparison compare_branches(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& baseline,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& scout
) {
    QueryComparison c{};
    std::vector<bool> used(scout.size(), false);
    for (const auto& b : baseline) {
        int best = -1;
        for (std::size_t i = 0; i < scout.size(); ++i) {
            if (!used[i] && branch_match(b, scout[i])) {
                best = static_cast<int>(i);
                break;
            }
        }
        if (best >= 0) {
            used[static_cast<std::size_t>(best)] = true;
            c.matched += 1;
        } else {
            c.missing += 1;
        }
    }
    int used_count = 0;
    for (const bool u : used) {
        used_count += u ? 1 : 0;
    }
    c.added = static_cast<int>(scout.size()) - used_count;
    return c;
}

std::string phase_name(bool adaptive) {
    return adaptive ? "adaptive" : "initial";
}

std::string pair_key(const QuerySample& sample) {
    return std::string(spaceship_cpp::planet_params::planet_name(sample.from_planet)) + "->" +
        spaceship_cpp::planet_params::planet_name(sample.to_planet);
}

std::string group_key(const QuerySample& sample) {
    return pair_key(sample) + "|" + phase_name(sample.adaptive);
}

std::vector<QuerySample> stratified_sample_queries(const std::vector<QuerySample>& queries, int max_queries) {
    if (static_cast<int>(queries.size()) <= max_queries) {
        return queries;
    }
    std::map<std::string, std::vector<int>> groups;
    for (std::size_t i = 0; i < queries.size(); ++i) {
        groups[group_key(queries[i])].push_back(static_cast<int>(i));
    }

    std::vector<int> selected;
    selected.reserve(static_cast<std::size_t>(max_queries));
    std::map<std::string, int> used_by_group;
    for (const auto& [key, indices] : groups) {
        const int quota = std::max(
            1,
            static_cast<int>(std::floor(
                static_cast<double>(max_queries) * static_cast<double>(indices.size()) /
                static_cast<double>(queries.size()))));
        const int take = std::min<int>(quota, indices.size());
        const double stride = static_cast<double>(indices.size()) / static_cast<double>(take);
        for (int j = 0; j < take && static_cast<int>(selected.size()) < max_queries; ++j) {
            const int pos = std::min<int>(indices.size() - 1, static_cast<int>(std::floor(stride * j)));
            selected.push_back(indices[static_cast<std::size_t>(pos)]);
            used_by_group[key] += 1;
        }
    }
    while (static_cast<int>(selected.size()) < max_queries) {
        bool added = false;
        for (const auto& [key, indices] : groups) {
            int& used = used_by_group[key];
            if (used < static_cast<int>(indices.size())) {
                selected.push_back(indices[static_cast<std::size_t>(used)]);
                used += 1;
                added = true;
                if (static_cast<int>(selected.size()) >= max_queries) {
                    break;
                }
            }
        }
        if (!added) {
            break;
        }
    }
    std::sort(selected.begin(), selected.end());
    std::vector<QuerySample> out;
    out.reserve(selected.size());
    for (const int index : selected) {
        out.push_back(queries[static_cast<std::size_t>(index)]);
    }
    return out;
}

void accumulate(
    RecallStats* stats,
    int baseline_count,
    int scout_count,
    const QueryComparison& comparison,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryResult& baseline,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryResult& scout
) {
    stats->query_count += 1;
    stats->baseline_branch_count_total += baseline_count;
    stats->missing_branch_count_total += comparison.missing;
    stats->baseline_ms_sum += baseline.profile.total_ms;
    stats->scout_ms_sum += scout.profile.total_ms;
    if (baseline_count > 0) {
        stats->baseline_positive_query_count += 1;
    }
    if (baseline_count > 0 && scout_count == 0) {
        stats->query_level_miss_count += 1;
    }
    if (baseline.profile.fallback_direct_solve_count > 0) {
        stats->baseline_direct_fallback_query_count += 1;
    }
    stats->scout_direct_fallback_skipped_count += scout.profile.direct_fallback_skipped_count;
    stats->scout_cell_vertex_success_count += scout.profile.cell_vertex_routeA_success_count;
    stats->scout_cell_vertex_failure_count += scout.profile.cell_vertex_routeA_failure_count;
}

void print_summary_values(const RecallStats& stats) {
    const double query_miss_rate = stats.baseline_positive_query_count > 0 ?
        static_cast<double>(stats.query_level_miss_count) /
            static_cast<double>(stats.baseline_positive_query_count) : 0.0;
    const double branch_miss_rate = stats.baseline_branch_count_total > 0 ?
        static_cast<double>(stats.missing_branch_count_total) /
            static_cast<double>(stats.baseline_branch_count_total) : 0.0;
    const double mean_baseline_ms = stats.query_count > 0 ?
        stats.baseline_ms_sum / static_cast<double>(stats.query_count) : 0.0;
    const double mean_scout_ms = stats.query_count > 0 ?
        stats.scout_ms_sum / static_cast<double>(stats.query_count) : 0.0;
    std::cout << "query_count=" << stats.query_count << '\n';
    std::cout << "baseline_positive_query_count=" << stats.baseline_positive_query_count << '\n';
    std::cout << "query_level_miss_count=" << stats.query_level_miss_count << '\n';
    std::cout << "query_level_miss_rate=" << query_miss_rate << '\n';
    std::cout << "baseline_branch_count_total=" << stats.baseline_branch_count_total << '\n';
    std::cout << "missing_branch_count_total=" << stats.missing_branch_count_total << '\n';
    std::cout << "branch_level_miss_rate=" << branch_miss_rate << '\n';
    std::cout << "mean_baseline_ms=" << mean_baseline_ms << '\n';
    std::cout << "mean_scout_ms=" << mean_scout_ms << '\n';
    std::cout << "speedup=" << (mean_scout_ms > 0.0 ? mean_baseline_ms / mean_scout_ms : 0.0) << '\n';
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
        std::cout << "problem1_scout_query_recall_rate_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);

    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet) + 1262332.782636;
    const double initial_theta = 6.169738905800;
    const std::vector<planet_params::PlanetId> allowed_first_targets{
        planet_params::PlanetId::Mercury, planet_params::PlanetId::Venus, planet_params::PlanetId::Mars};
    const std::vector<planet_params::PlanetId> allowed_transfer_planets{
        planet_params::PlanetId::Mercury, planet_params::PlanetId::Venus,
        planet_params::PlanetId::Earth, planet_params::PlanetId::Mars};

    bfs::InitialLaunchExpansionOptions initial_options{};
    initial_options.max_transfer_revolution = 1;
    initial_options.max_target_revolution = 1;
    initial_options.max_launch_v_inf = 7000.0;

    problem2::Problem2GravityAssistSolverOptions problem2_options{};
    problem2_options.theta_sample_count = 64;
    problem2_options.topology_adaptive_enabled = true;
    problem2_options.topology_max_depth = 10;
    problem2_options.max_transfer_revolution = 1;
    problem2_options.max_target_revolution = 1;
    problem2_options.adaptive_interval_trace_enabled = true;

    bfs::FixedLaunchThetaSearchOptions search_options{};
    search_options.max_depth = 5;
    search_options.enable_beam_pruning = false;
    search_options.beam_width = 3;
    search_options.flyby_physical_options.enabled = true;
    search_options.flyby_physical_options.mode = trajectory::FlybyPhysicalFilterMode::Enforce;
    search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
    search_options.max_launch_v_inf = 7000.0;
    search_options.continue_after_reaching_terminal = true;

    const auto search = bfs::run_fixed_launch_theta_beam_search_with_table(
        loader, planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury,
        launch_time, initial_theta, allowed_first_targets, allowed_transfer_planets,
        initial_options, problem2_options, search_options);
    if (!search.ok) {
        std::cout << "problem1_scout_query_recall_rate_search_failed\n";
        std::cout << "error_message=" << search.error_message << '\n';
        return 1;
    }

    std::vector<QuerySample> all_queries;
    for (std::size_t attempt_index = 0; attempt_index < search.expansion_attempt_profiles.size(); ++attempt_index) {
        const auto& attempt = search.expansion_attempt_profiles[attempt_index];
        for (const auto& trace : attempt.problem2_solver_profile.theta_sample_traces) {
            QuerySample sample{};
            sample.source_attempt_index = static_cast<int>(attempt_index);
            sample.from_planet = attempt.from_planet;
            sample.to_planet = attempt.to_planet;
            sample.parent_time = attempt.parent_time;
            sample.query_nu_A = trace.query_nu_A;
            sample.query_nu_B = trace.query_nu_B;
            sample.query_theta_A = trace.query_theta_A;
            sample.adaptive = trace.adaptive;
            all_queries.push_back(sample);
        }
    }

    const int max_queries = max_queries_from_env();
    const auto queries = stratified_sample_queries(all_queries, max_queries);
    std::cout << "Problem1ScoutQueryRecallCollection\n";
    std::cout << "collected_query_count=" << all_queries.size() << '\n';
    std::cout << "sampled_query_count=" << queries.size() << '\n';
    std::cout << "max_queries=" << max_queries << '\n';

    problem1::Problem1NearestNodeQueryOptions baseline_options{};
    baseline_options.fallback_mode = problem1::Problem1FallbackMode::DirectSolve;
    baseline_options.fallback_direct_solve = true;

    problem1::Problem1NearestNodeQueryOptions scout_options{};
    scout_options.fallback_mode = problem1::Problem1FallbackMode::ScoutNoDirectFallback;
    scout_options.fallback_direct_solve = true;
    scout_options.cell_vertex_routeA_enabled = true;
    scout_options.cell_vertex_seed_policy = 1;
    scout_options.cell_vertex_top_k_vertices = 2;
    scout_options.cell_vertex_max_branch_attempts = 8;
    scout_options.scout_no_direct_fallback = true;

    RecallStats global_stats{};
    std::map<std::string, RecallStats> pair_stats;
    std::map<std::string, RecallStats> phase_stats;

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const auto& sample = queries[i];
        const auto baseline = problem1::query_problem1_from_2deg_nearest_node(
            loader,
            sample.from_planet,
            sample.to_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            1,
            1,
            baseline_options);
        const auto scout = problem1::query_problem1_from_2deg_nearest_node(
            loader,
            sample.from_planet,
            sample.to_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            1,
            1,
            scout_options);
        const auto comparison = compare_branches(baseline.branches, scout.branches);
        const int baseline_count = static_cast<int>(baseline.branches.size());
        const int scout_count = static_cast<int>(scout.branches.size());
        const int query_level_miss = (baseline_count > 0 && scout_count == 0) ? 1 : 0;
        const double branch_miss_rate = baseline_count > 0 ?
            static_cast<double>(comparison.missing) / static_cast<double>(baseline_count) : 0.0;

        std::cout << "Problem1ScoutQueryRecallSample\n";
        std::cout << "index=" << i << '\n';
        std::cout << "source_attempt_index=" << sample.source_attempt_index << '\n';
        std::cout << "from_planet=" << planet_params::planet_name(sample.from_planet) << '\n';
        std::cout << "to_planet=" << planet_params::planet_name(sample.to_planet) << '\n';
        std::cout << "phase=" << phase_name(sample.adaptive) << '\n';
        std::cout << "parent_time=" << sample.parent_time << '\n';
        std::cout << "query_nu_A=" << sample.query_nu_A << '\n';
        std::cout << "query_nu_B=" << sample.query_nu_B << '\n';
        std::cout << "query_theta_A=" << sample.query_theta_A << '\n';
        std::cout << "baseline_branch_count=" << baseline_count << '\n';
        std::cout << "scout_branch_count=" << scout_count << '\n';
        std::cout << "matched_branch_count=" << comparison.matched << '\n';
        std::cout << "missing_branch_count=" << comparison.missing << '\n';
        std::cout << "new_branch_count=" << comparison.added << '\n';
        std::cout << "baseline_positive=" << (baseline_count > 0 ? 1 : 0) << '\n';
        std::cout << "scout_positive=" << (scout_count > 0 ? 1 : 0) << '\n';
        std::cout << "query_level_miss=" << query_level_miss << '\n';
        std::cout << "branch_level_miss_rate=" << branch_miss_rate << '\n';
        std::cout << "baseline_used_direct_fallback=" << (baseline.used_direct_solve_fallback ? 1 : 0) << '\n';
        std::cout << "scout_direct_fallback_skipped=" << scout.profile.direct_fallback_skipped_count << '\n';
        std::cout << "scout_cell_vertex_success=" << scout.profile.cell_vertex_routeA_success_count << '\n';
        std::cout << "scout_cell_vertex_failure=" << scout.profile.cell_vertex_routeA_failure_count << '\n';
        std::cout << "baseline_ms=" << baseline.profile.total_ms << '\n';
        std::cout << "scout_ms=" << scout.profile.total_ms << '\n';

        accumulate(&global_stats, baseline_count, scout_count, comparison, baseline, scout);
        accumulate(&pair_stats[pair_key(sample)], baseline_count, scout_count, comparison, baseline, scout);
        accumulate(&phase_stats[phase_name(sample.adaptive)], baseline_count, scout_count, comparison, baseline, scout);
    }

    double max_pair_query_miss_rate = 0.0;
    double max_pair_branch_miss_rate = 0.0;
    for (const auto& [key, stats] : pair_stats) {
        const auto arrow = key.find("->");
        std::cout << "Problem1ScoutQueryRecallPairSummary\n";
        std::cout << "from_planet=" << key.substr(0, arrow) << '\n';
        std::cout << "to_planet=" << key.substr(arrow + 2) << '\n';
        print_summary_values(stats);
        std::cout << "baseline_direct_fallback_query_count=" << stats.baseline_direct_fallback_query_count << '\n';
        std::cout << "scout_direct_fallback_skipped_count=" << stats.scout_direct_fallback_skipped_count << '\n';
        std::cout << "scout_cell_vertex_success_count=" << stats.scout_cell_vertex_success_count << '\n';
        std::cout << "scout_cell_vertex_failure_count=" << stats.scout_cell_vertex_failure_count << '\n';
        const double q = stats.baseline_positive_query_count > 0 ?
            static_cast<double>(stats.query_level_miss_count) /
                static_cast<double>(stats.baseline_positive_query_count) : 0.0;
        const double b = stats.baseline_branch_count_total > 0 ?
            static_cast<double>(stats.missing_branch_count_total) /
                static_cast<double>(stats.baseline_branch_count_total) : 0.0;
        max_pair_query_miss_rate = std::max(max_pair_query_miss_rate, q);
        max_pair_branch_miss_rate = std::max(max_pair_branch_miss_rate, b);
    }

    for (const auto& [phase, stats] : phase_stats) {
        std::cout << "Problem1ScoutQueryRecallPhaseSummary\n";
        std::cout << "phase=" << phase << '\n';
        print_summary_values(stats);
    }

    std::cout << "Problem1ScoutQueryRecallGlobalSummary\n";
    print_summary_values(global_stats);
    std::cout << "max_pair_query_level_miss_rate=" << max_pair_query_miss_rate << '\n';
    std::cout << "max_pair_branch_level_miss_rate=" << max_pair_branch_miss_rate << '\n';
    if (global_stats.query_level_miss_count == 0 && global_stats.baseline_positive_query_count > 0) {
        std::cout << "upper_95_query_miss_rate="
                  << (3.0 / static_cast<double>(global_stats.baseline_positive_query_count)) << '\n';
    } else {
        const double observed = global_stats.baseline_positive_query_count > 0 ?
            static_cast<double>(global_stats.query_level_miss_count) /
                static_cast<double>(global_stats.baseline_positive_query_count) : 0.0;
        std::cout << "upper_95_query_miss_rate=" << observed << '\n';
    }
    std::cout << "scout_query_recall_ok=1\n";
    return 0;
}
