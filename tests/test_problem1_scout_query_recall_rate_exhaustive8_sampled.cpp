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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct QueryRecord {
    spaceship_cpp::planet_params::PlanetId from_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    std::string phase = "initial";
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
    int baseline_branch_count_hint = 0;
    int scout_branch_count_hint = 0;
    bool baseline_positive = false;
    bool baseline_used_direct_fallback = false;
    bool query_level_miss_hint = false;
    double branch_level_miss_rate_hint = 0.0;
};

struct CompareResult {
    int matched = 0;
    int missing = 0;
    int added = 0;
};

struct Summary {
    int query_count = 0;
    int baseline_positive_query_count = 0;
    int query_level_miss_count = 0;
    int baseline_branch_count_total = 0;
    int missing_branch_count_total = 0;
    double baseline_ms_sum = 0.0;
    double config_ms_sum = 0.0;
    long long cell_vertex_attempted_branch_count = 0;
    int cell_vertex_routeA_success_count = 0;
    int cell_vertex_routeA_failure_count = 0;
    int direct_fallback_skipped_count = 0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

int env_int(const char* name, int fallback) {
    if (const char* raw = std::getenv(name)) {
        if (*raw != '\0') {
            return std::max(0, std::atoi(raw));
        }
    }
    return fallback;
}

spaceship_cpp::planet_params::PlanetId parse_planet(const std::string& name) {
    namespace planet_params = spaceship_cpp::planet_params;
    if (name == "Mercury") return planet_params::PlanetId::Mercury;
    if (name == "Venus") return planet_params::PlanetId::Venus;
    if (name == "Earth") return planet_params::PlanetId::Earth;
    if (name == "Mars") return planet_params::PlanetId::Mars;
    if (name == "Jupiter") return planet_params::PlanetId::Jupiter;
    if (name == "Saturn") return planet_params::PlanetId::Saturn;
    if (name == "Uranus") return planet_params::PlanetId::Uranus;
    return planet_params::PlanetId::Neptune;
}

std::string planet_name_string(spaceship_cpp::planet_params::PlanetId id) {
    return spaceship_cpp::planet_params::planet_name(id);
}

std::string pair_name(const QueryRecord& q) {
    return planet_name_string(q.from_planet) + "->" + planet_name_string(q.to_planet);
}

bool is_high_risk_pair(const QueryRecord& q) {
    const auto key = pair_name(q);
    return key == "Venus->Mars" ||
        key == "Mars->Earth" ||
        key == "Earth->Mercury" ||
        key == "Mars->Venus" ||
        key == "Venus->Earth" ||
        key == "Mars->Mercury" ||
        key == "Venus->Mercury";
}

int group_weight(const QueryRecord& q) {
    int weight = is_high_risk_pair(q) ? 8 : 1;
    if (pair_name(q) == "Earth->Mars") {
        weight = 3;
    }
    if (q.baseline_positive) {
        weight *= 3;
    }
    if (q.baseline_used_direct_fallback) {
        weight *= 4;
    }
    return weight;
}

std::string group_key(const QueryRecord& q) {
    return pair_name(q) + "|" + q.phase + "|" +
        (q.baseline_positive ? "positive" : "zero") + "|" +
        (q.baseline_used_direct_fallback ? "direct" : "table");
}

std::vector<int> deterministic_take(const std::vector<int>& indices, int take) {
    std::vector<int> out;
    if (take <= 0 || indices.empty()) {
        return out;
    }
    take = std::min<int>(take, indices.size());
    out.reserve(static_cast<std::size_t>(take));
    const double stride = static_cast<double>(indices.size()) / static_cast<double>(take);
    for (int i = 0; i < take; ++i) {
        const int pos = std::min<int>(indices.size() - 1, static_cast<int>(std::floor(stride * i)));
        out.push_back(indices[static_cast<std::size_t>(pos)]);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    for (int i = 0; static_cast<int>(out.size()) < take && i < static_cast<int>(indices.size()); ++i) {
        if (!std::binary_search(out.begin(), out.end(), indices[static_cast<std::size_t>(i)])) {
            out.push_back(indices[static_cast<std::size_t>(i)]);
            std::sort(out.begin(), out.end());
        }
    }
    return out;
}

std::vector<QueryRecord> sample_queries(const std::vector<QueryRecord>& records, int max_queries) {
    if (static_cast<int>(records.size()) <= max_queries) {
        return records;
    }
    std::map<std::string, std::vector<int>> groups;
    std::map<std::string, int> weights;
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto key = group_key(records[i]);
        groups[key].push_back(static_cast<int>(i));
        weights[key] = std::max(weights[key], group_weight(records[i]));
    }

    std::set<int> selected;
    for (const auto& [key, indices] : groups) {
        const int base_take = std::min<int>(20, indices.size());
        for (const int idx : deterministic_take(indices, base_take)) {
            if (static_cast<int>(selected.size()) < max_queries) {
                selected.insert(idx);
            }
        }
    }

    std::map<std::string, int> selected_by_group;
    for (const int idx : selected) {
        selected_by_group[group_key(records[static_cast<std::size_t>(idx)])] += 1;
    }
    while (static_cast<int>(selected.size()) < max_queries) {
        double best_score = -1.0;
        std::string best_key;
        for (const auto& [key, indices] : groups) {
            const int used = selected_by_group[key];
            if (used >= static_cast<int>(indices.size())) {
                continue;
            }
            const double score = static_cast<double>(weights[key]) *
                static_cast<double>(indices.size()) / static_cast<double>(used + 1);
            if (score > best_score) {
                best_score = score;
                best_key = key;
            }
        }
        if (best_key.empty()) {
            break;
        }
        const auto& indices = groups[best_key];
        const int used = selected_by_group[best_key];
        const double stride = static_cast<double>(indices.size()) / static_cast<double>(used + 1);
        int pos = std::min<int>(indices.size() - 1, static_cast<int>(std::floor(stride * used)));
        int idx = indices[static_cast<std::size_t>(pos)];
        while (selected.count(idx) != 0 && pos + 1 < static_cast<int>(indices.size())) {
            pos += 1;
            idx = indices[static_cast<std::size_t>(pos)];
        }
        if (selected.count(idx) == 0) {
            selected.insert(idx);
            selected_by_group[best_key] += 1;
        } else {
            selected_by_group[best_key] = static_cast<int>(indices.size());
        }
    }

    std::vector<QueryRecord> out;
    out.reserve(selected.size());
    for (const int idx : selected) {
        out.push_back(records[static_cast<std::size_t>(idx)]);
    }
    return out;
}

bool parse_bool_int(const std::string& raw) {
    return raw == "1" || raw == "true";
}

bool parse_existing_log(const std::filesystem::path& path, std::vector<QueryRecord>* records) {
    std::ifstream in(path);
    if (!in) {
        return false;
    }
    std::string line;
    QueryRecord current{};
    bool in_sample = false;
    bool have_from = false;
    bool have_to = false;
    bool have_phase = false;
    bool have_coords = false;
    auto flush = [&]() {
        if (in_sample && have_from && have_to && have_phase && have_coords) {
            records->push_back(current);
        }
        current = QueryRecord{};
        in_sample = false;
        have_from = have_to = have_phase = have_coords = false;
    };
    while (std::getline(in, line)) {
        if (line == "Problem1ScoutQueryRecallSample") {
            flush();
            in_sample = true;
            continue;
        }
        if (line.rfind("Problem1", 0) == 0) {
            flush();
            continue;
        }
        if (!in_sample) {
            continue;
        }
        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        const auto key = line.substr(0, eq);
        const auto value = line.substr(eq + 1);
        if (key == "from_planet") {
            current.from_planet = parse_planet(value);
            have_from = true;
        } else if (key == "to_planet") {
            current.to_planet = parse_planet(value);
            have_to = true;
        } else if (key == "phase") {
            current.phase = value;
            have_phase = true;
        } else if (key == "query_nu_A") {
            current.query_nu_A = std::stod(value);
        } else if (key == "query_nu_B") {
            current.query_nu_B = std::stod(value);
        } else if (key == "query_theta_A") {
            current.query_theta_A = std::stod(value);
            have_coords = true;
        } else if (key == "baseline_branch_count") {
            current.baseline_branch_count_hint = std::stoi(value);
        } else if (key == "scout_branch_count") {
            current.scout_branch_count_hint = std::stoi(value);
        } else if (key == "baseline_positive") {
            current.baseline_positive = parse_bool_int(value);
        } else if (key == "baseline_used_direct_fallback") {
            current.baseline_used_direct_fallback = parse_bool_int(value);
        } else if (key == "query_level_miss") {
            current.query_level_miss_hint = parse_bool_int(value);
        } else if (key == "branch_level_miss_rate") {
            current.branch_level_miss_rate_hint = std::stod(value);
        }
    }
    flush();
    return !records->empty();
}

std::vector<QueryRecord> regenerate_queries_from_bfs(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader
) {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem2 = spaceship_cpp::problem2;
    namespace trajectory = spaceship_cpp::trajectory;

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
    std::vector<QueryRecord> out;
    if (!search.ok) {
        return out;
    }
    for (const auto& attempt : search.expansion_attempt_profiles) {
        for (const auto& trace : attempt.problem2_solver_profile.theta_sample_traces) {
            QueryRecord q{};
            q.from_planet = attempt.from_planet;
            q.to_planet = attempt.to_planet;
            q.phase = trace.adaptive ? "adaptive" : "initial";
            q.query_nu_A = trace.query_nu_A;
            q.query_nu_B = trace.query_nu_B;
            q.query_theta_A = trace.query_theta_A;
            out.push_back(q);
        }
    }
    return out;
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

CompareResult compare_branches(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& baseline,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& config
) {
    CompareResult c{};
    std::vector<bool> used(config.size(), false);
    for (const auto& branch : baseline) {
        int best = -1;
        for (std::size_t i = 0; i < config.size(); ++i) {
            if (!used[i] && branch_match(branch, config[i])) {
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
    c.added = static_cast<int>(config.size()) - used_count;
    return c;
}

spaceship_cpp::problem1::Problem1NearestNodeQueryOptions make_options(
    const std::string& config_name
) {
    namespace problem1 = spaceship_cpp::problem1;
    problem1::Problem1NearestNodeQueryOptions options{};
    options.fallback_direct_solve = true;
    if (config_name == "Baseline") {
        options.fallback_mode = problem1::Problem1FallbackMode::DirectSolve;
        return options;
    }
    options.fallback_mode = problem1::Problem1FallbackMode::ScoutNoDirectFallback;
    options.cell_vertex_routeA_enabled = true;
    options.cell_vertex_seed_policy = 1;
    options.scout_no_direct_fallback = true;
    if (config_name == "Exhaustive8Scout") {
        options.cell_vertex_top_k_vertices = 8;
        options.cell_vertex_max_branch_attempts = 1000000;
    } else {
        options.cell_vertex_top_k_vertices = 2;
        options.cell_vertex_max_branch_attempts = 8;
    }
    return options;
}

void accumulate_summary(
    Summary* summary,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryResult& baseline,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryResult& config,
    const CompareResult& comparison
) {
    const int baseline_count = static_cast<int>(baseline.branches.size());
    const int config_count = static_cast<int>(config.branches.size());
    summary->query_count += 1;
    summary->baseline_branch_count_total += baseline_count;
    summary->missing_branch_count_total += comparison.missing;
    summary->baseline_ms_sum += baseline.profile.total_ms;
    summary->config_ms_sum += config.profile.total_ms;
    if (baseline_count > 0) {
        summary->baseline_positive_query_count += 1;
    }
    if (baseline_count > 0 && config_count == 0) {
        summary->query_level_miss_count += 1;
    }
    summary->cell_vertex_attempted_branch_count += config.profile.cell_vertex_attempted_branch_count;
    summary->cell_vertex_routeA_success_count += config.profile.cell_vertex_routeA_success_count;
    summary->cell_vertex_routeA_failure_count += config.profile.cell_vertex_routeA_failure_count;
    summary->direct_fallback_skipped_count += config.profile.direct_fallback_skipped_count;
}

double query_miss_rate(const Summary& s) {
    return s.baseline_positive_query_count > 0 ?
        static_cast<double>(s.query_level_miss_count) / static_cast<double>(s.baseline_positive_query_count) : 0.0;
}

double branch_miss_rate(const Summary& s) {
    return s.baseline_branch_count_total > 0 ?
        static_cast<double>(s.missing_branch_count_total) / static_cast<double>(s.baseline_branch_count_total) : 0.0;
}

double mean_baseline_ms(const Summary& s) {
    return s.query_count > 0 ? s.baseline_ms_sum / static_cast<double>(s.query_count) : 0.0;
}

double mean_config_ms(const Summary& s) {
    return s.query_count > 0 ? s.config_ms_sum / static_cast<double>(s.query_count) : 0.0;
}

void print_summary_fields(const Summary& s) {
    const double baseline_ms = mean_baseline_ms(s);
    const double config_ms = mean_config_ms(s);
    std::cout << "query_count=" << s.query_count << '\n';
    std::cout << "baseline_positive_query_count=" << s.baseline_positive_query_count << '\n';
    std::cout << "query_level_miss_count=" << s.query_level_miss_count << '\n';
    std::cout << "query_level_miss_rate=" << query_miss_rate(s) << '\n';
    std::cout << "baseline_branch_count_total=" << s.baseline_branch_count_total << '\n';
    std::cout << "missing_branch_count_total=" << s.missing_branch_count_total << '\n';
    std::cout << "branch_level_miss_rate=" << branch_miss_rate(s) << '\n';
    std::cout << "mean_baseline_ms=" << baseline_ms << '\n';
    std::cout << "mean_config_ms=" << config_ms << '\n';
    std::cout << "speedup_vs_baseline=" << (config_ms > 0.0 ? baseline_ms / config_ms : 0.0) << '\n';
}

void print_global_summary(const std::string& config_name, const Summary& s) {
    std::cout << "Problem1ScoutExhaustive8GlobalSummary\n";
    std::cout << "config_name=" << config_name << '\n';
    print_summary_fields(s);
    std::cout << "cell_vertex_attempted_branch_count=" << s.cell_vertex_attempted_branch_count << '\n';
    std::cout << "cell_vertex_routeA_success_count=" << s.cell_vertex_routeA_success_count << '\n';
    std::cout << "cell_vertex_routeA_failure_count=" << s.cell_vertex_routeA_failure_count << '\n';
    std::cout << "direct_fallback_skipped_count=" << s.direct_fallback_skipped_count << '\n';
    std::cout << "mean_attempted_branch_per_query="
              << (s.query_count > 0 ? static_cast<double>(s.cell_vertex_attempted_branch_count) /
                      static_cast<double>(s.query_count) : 0.0) << '\n';
    std::cout << "mean_attempted_branch_per_positive_query="
              << (s.baseline_positive_query_count > 0 ? static_cast<double>(s.cell_vertex_attempted_branch_count) /
                      static_cast<double>(s.baseline_positive_query_count) : 0.0) << '\n';
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem1_scout_query_recall_rate_exhaustive8_sampled_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);

    std::vector<QueryRecord> available;
    std::string source = "existing_log";
    if (!parse_existing_log("/tmp/problem1_scout_query_recall_rate.log", &available)) {
        source = "regenerated_bfs";
        available = regenerate_queries_from_bfs(loader);
    }
    if (available.empty()) {
        std::cout << "problem1_scout_query_recall_rate_exhaustive8_sampled_no_queries\n";
        return 1;
    }

    const int sample_max = env_int("SCOUT_EXHAUSTIVE8_SAMPLE_MAX_QUERIES", 5000);
    const int max_miss_examples = env_int("SCOUT_EXHAUSTIVE8_MAX_MISS_EXAMPLES", 50);
    const auto sampled = sample_queries(available, sample_max);
    std::set<std::string> groups;
    int positive_count = 0;
    int direct_count = 0;
    for (const auto& q : sampled) {
        groups.insert(group_key(q));
        positive_count += q.baseline_positive ? 1 : 0;
        direct_count += q.baseline_used_direct_fallback ? 1 : 0;
    }
    std::cout << "Problem1ScoutExhaustive8SampleCollection\n";
    std::cout << "source=" << source << '\n';
    std::cout << "available_query_count=" << available.size() << '\n';
    std::cout << "sampled_query_count=" << sampled.size() << '\n';
    std::cout << "sample_max_queries=" << sample_max << '\n';
    std::cout << "group_count=" << groups.size() << '\n';
    std::cout << "positive_query_count_in_sample=" << positive_count << '\n';
    std::cout << "direct_fallback_query_count_in_sample=" << direct_count << '\n';

    const auto baseline_options = make_options("Baseline");
    const auto current_options = make_options("CurrentScout");
    const auto exhaustive_options = make_options("Exhaustive8Scout");

    std::map<std::string, Summary> global_by_config;
    std::map<std::string, std::map<std::string, Summary>> pair_by_config;
    std::map<std::string, std::map<std::string, Summary>> phase_by_config;
    int miss_examples_printed = 0;

    for (std::size_t i = 0; i < sampled.size(); ++i) {
        const auto& q = sampled[i];
        const auto baseline = problem1::query_problem1_from_2deg_nearest_node(
            loader, q.from_planet, q.to_planet, q.query_nu_A, q.query_nu_B, q.query_theta_A, 1, 1, baseline_options);
        for (const auto& config : std::vector<std::pair<std::string, problem1::Problem1NearestNodeQueryOptions>>{
                 {"CurrentScout", current_options},
                 {"Exhaustive8Scout", exhaustive_options}}) {
            const auto result = problem1::query_problem1_from_2deg_nearest_node(
                loader, q.from_planet, q.to_planet, q.query_nu_A, q.query_nu_B, q.query_theta_A, 1, 1,
                config.second);
            const auto comparison = compare_branches(baseline.branches, result.branches);
            accumulate_summary(&global_by_config[config.first], baseline, result, comparison);
            accumulate_summary(&pair_by_config[config.first][pair_name(q)], baseline, result, comparison);
            accumulate_summary(&phase_by_config[config.first][q.phase], baseline, result, comparison);
            if (comparison.missing > 0 && miss_examples_printed < max_miss_examples) {
                std::cout << "Problem1ScoutExhaustive8MissExample\n";
                std::cout << "config_name=" << config.first << '\n';
                std::cout << "index=" << i << '\n';
                std::cout << "from_planet=" << planet_name_string(q.from_planet) << '\n';
                std::cout << "to_planet=" << planet_name_string(q.to_planet) << '\n';
                std::cout << "phase=" << q.phase << '\n';
                std::cout << "baseline_branch_count=" << baseline.branches.size() << '\n';
                std::cout << "config_branch_count=" << result.branches.size() << '\n';
                std::cout << "missing_branch_count=" << comparison.missing << '\n';
                std::cout << "baseline_used_direct_fallback=" << (baseline.used_direct_solve_fallback ? 1 : 0)
                          << '\n';
                std::cout << "config_direct_fallback_skipped=" << result.profile.direct_fallback_skipped_count
                          << '\n';
                std::cout << "config_cell_vertex_attempted_branch_count="
                          << result.profile.cell_vertex_attempted_branch_count << '\n';
                std::cout << "query_nu_A=" << q.query_nu_A << '\n';
                std::cout << "query_nu_B=" << q.query_nu_B << '\n';
                std::cout << "query_theta_A=" << q.query_theta_A << '\n';
                miss_examples_printed += 1;
            }
        }
    }

    for (const auto& config_name : {"CurrentScout", "Exhaustive8Scout"}) {
        const auto& global = global_by_config[config_name];
        print_global_summary(config_name, global);
        for (const auto& [pair, summary] : pair_by_config[config_name]) {
            const auto arrow = pair.find("->");
            std::cout << "Problem1ScoutExhaustive8PairSummary\n";
            std::cout << "config_name=" << config_name << '\n';
            std::cout << "from_planet=" << pair.substr(0, arrow) << '\n';
            std::cout << "to_planet=" << pair.substr(arrow + 2) << '\n';
            print_summary_fields(summary);
        }
        for (const auto& [phase, summary] : phase_by_config[config_name]) {
            std::cout << "Problem1ScoutExhaustive8PhaseSummary\n";
            std::cout << "config_name=" << config_name << '\n';
            std::cout << "phase=" << phase << '\n';
            print_summary_fields(summary);
        }
    }

    const auto& current = global_by_config["CurrentScout"];
    const auto& exhaustive = global_by_config["Exhaustive8Scout"];
    std::cout << "Problem1ScoutExhaustive8ComparisonSummary\n";
    std::cout << "sampled_query_count=" << sampled.size() << '\n';
    std::cout << "current_query_miss_rate=" << query_miss_rate(current) << '\n';
    std::cout << "exhaustive8_query_miss_rate=" << query_miss_rate(exhaustive) << '\n';
    std::cout << "current_branch_miss_rate=" << branch_miss_rate(current) << '\n';
    std::cout << "exhaustive8_branch_miss_rate=" << branch_miss_rate(exhaustive) << '\n';
    std::cout << "current_mean_ms=" << mean_config_ms(current) << '\n';
    std::cout << "exhaustive8_mean_ms=" << mean_config_ms(exhaustive) << '\n';
    std::cout << "baseline_mean_ms=" << mean_baseline_ms(exhaustive) << '\n';
    const bool improves = query_miss_rate(exhaustive) < query_miss_rate(current) &&
        branch_miss_rate(exhaustive) < branch_miss_rate(current);
    const bool faster = mean_config_ms(exhaustive) < mean_baseline_ms(exhaustive);
    const bool worth = improves && faster && query_miss_rate(exhaustive) < 0.05;
    std::cout << "exhaustive8_improves_recall=" << (improves ? 1 : 0) << '\n';
    std::cout << "exhaustive8_still_faster_than_baseline=" << (faster ? 1 : 0) << '\n';
    std::cout << "exhaustive8_worth_full_run=" << (worth ? 1 : 0) << '\n';
    std::cout << "problem1_scout_query_recall_rate_exhaustive8_sampled_ok=1\n";
    return 0;
}
