#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
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
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

struct SelectedCase {
    std::string name;
    spaceship_cpp::bfs::BfsExpansionAttemptProfile attempt;
};

struct ConfigCase {
    std::string name;
    spaceship_cpp::problem2::Problem2GravityAssistSolverOptions options;
};

struct RunRecord {
    std::string case_name;
    std::string config_name;
    spaceship_cpp::planet_params::PlanetId from_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::problem2::Problem2GravityAssistSolverResult result;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
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

bool same_pair(
    const spaceship_cpp::bfs::BfsExpansionAttemptProfile& attempt,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to
) {
    return attempt.from_planet == from && attempt.to_planet == to;
}

const spaceship_cpp::bfs::BfsExpansionAttemptProfile* find_first(
    const std::vector<spaceship_cpp::bfs::BfsExpansionAttemptProfile>& attempts,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    bool require_zero_raw,
    bool require_accepted
) {
    for (const auto& attempt : attempts) {
        if (!same_pair(attempt, from, to)) {
            continue;
        }
        if (require_zero_raw && !attempt.zero_raw) {
            continue;
        }
        if (require_accepted && attempt.accepted_edge_count <= 0) {
            continue;
        }
        return &attempt;
    }
    return nullptr;
}

void add_case_if_present(
    std::vector<SelectedCase>* cases,
    std::set<std::string>* names,
    const std::string& name,
    const spaceship_cpp::bfs::BfsExpansionAttemptProfile* attempt
) {
    if (attempt == nullptr || names->count(name) != 0) {
        return;
    }
    cases->push_back(SelectedCase{name, *attempt});
    names->insert(name);
}

std::vector<ConfigCase> make_configs() {
    namespace problem2 = spaceship_cpp::problem2;
    std::vector<ConfigCase> configs;

    problem2::Problem2GravityAssistSolverOptions baseline{};
    baseline.theta_sample_count = 64;
    baseline.topology_adaptive_enabled = true;
    baseline.topology_max_depth = 10;
    baseline.max_transfer_revolution = 1;
    baseline.max_target_revolution = 1;
    configs.push_back({"baseline_full", baseline});

    auto no_adaptive = baseline;
    no_adaptive.topology_adaptive_enabled = false;
    configs.push_back({"fast_no_adaptive", no_adaptive});

    for (const int depth : {2, 4, 6}) {
        auto limited = baseline;
        limited.topology_max_depth = depth;
        configs.push_back({"adaptive_depth" + std::to_string(depth), limited});
    }

    auto sample32 = baseline;
    sample32.theta_sample_count = 32;
    configs.push_back({"sample32_full_adaptive", sample32});

    auto cheap = baseline;
    cheap.theta_sample_count = 16;
    cheap.topology_adaptive_enabled = false;
    configs.push_back({"cheap_scout_16_no_adaptive", cheap});

    return configs;
}

void print_case_result(
    const SelectedCase& selected,
    const ConfigCase& config,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& result
) {
    namespace planet_params = spaceship_cpp::planet_params;
    const auto& profile = result.profile;
    const auto& summary = result.summary;
    std::cout << "Problem2EfficiencyAblationCase\n";
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(selected.attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(selected.attempt.to_planet) << '\n';
    std::cout << "config_name=" << config.name << '\n';
    std::cout << "theta_sample_count=" << config.options.theta_sample_count << '\n';
    std::cout << "topology_adaptive_enabled=" << (config.options.topology_adaptive_enabled ? 1 : 0) << '\n';
    std::cout << "topology_max_depth=" << config.options.topology_max_depth << '\n';
    std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "solution_count=" << result.solutions.size() << '\n';
    std::cout << "dedup_success_count=" << summary.dedup_success_count << '\n';
    std::cout << "max_abs_slingshot_residual=" << summary.max_abs_slingshot_residual_at_root << '\n';
    std::cout << "total_ms=" << profile.total_ms << '\n';
    std::cout << "initial_sampling_ms=" << profile.initial_sampling_ms << '\n';
    std::cout << "topology_adaptive_sampling_ms=" << profile.topology_adaptive_sampling_ms << '\n';
    std::cout << "candidate_collection_ms=" << profile.candidate_collection_ms << '\n';
    std::cout << "midpoint_continuation_ms=" << profile.midpoint_continuation_ms << '\n';
    std::cout << "bisection_ms=" << profile.bisection_ms << '\n';
    std::cout << "dedup_ms=" << profile.dedup_ms << '\n';
    std::cout << "nearest_query_count=" << profile.problem1_nearest_node_query_count << '\n';
    std::cout << "nearest_query_total_ms=" << profile.nearest_query_total_ms << '\n';
    std::cout << "nearest_query_table_load_ms=" << profile.nearest_query_table_load_ms << '\n';
    std::cout << "nearest_query_table_offset_build_ms=" << profile.nearest_query_table_offset_build_ms << '\n';
    std::cout << "nearest_query_node_read_ms=" << profile.nearest_query_node_read_ms << '\n';
    std::cout << "nearest_query_refine_ms=" << profile.nearest_query_refine_ms << '\n';
    std::cout << "nearest_query_derivative_attach_ms=" << profile.nearest_query_derivative_attach_ms << '\n';
    std::cout << "cached_sample_hit_count=" << profile.cached_sample_hit_count << '\n';
    std::cout << "total_sample_build_count=" << profile.total_sample_build_count << '\n';
    std::cout << "adaptive_midpoint_sample_count=" << profile.adaptive_midpoint_sample_count << '\n';
}

struct Comparison {
    int matched_count = 0;
    int missing_from_config_count = 0;
    int new_in_config_count = 0;
    double max_abs_arrival_time_diff = 0.0;
    double max_abs_alpha_diff = 0.0;
    double max_abs_theta_prime_diff = 0.0;
    bool matches_baseline = false;
};

Comparison compare_to_baseline(
    double encounter_time,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& baseline,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& candidate
) {
    constexpr double kArrivalToleranceSeconds = 1e-3;
    constexpr double kAngleTolerance = 1e-8;

    Comparison comparison{};
    std::vector<bool> used(candidate.size(), false);
    for (const auto& base : baseline) {
        int best_index = -1;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < candidate.size(); ++i) {
            if (used[i]) {
                continue;
            }
            const auto& other = candidate[i];
            const double arrival_diff =
                std::abs((encounter_time + base.time_of_flight_seconds) -
                         (encounter_time + other.time_of_flight_seconds));
            const double alpha_diff = std::abs(base.alpha - other.alpha);
            const double theta_diff = std::abs(base.theta_prime - other.theta_prime);
            const double score = arrival_diff + alpha_diff * 1e6 + theta_diff * 1e6;
            if (score < best_score) {
                best_score = score;
                best_index = static_cast<int>(i);
            }
        }
        if (best_index < 0) {
            comparison.missing_from_config_count += 1;
            continue;
        }
        const auto& matched = candidate[static_cast<std::size_t>(best_index)];
        const double arrival_diff =
            std::abs((encounter_time + base.time_of_flight_seconds) -
                     (encounter_time + matched.time_of_flight_seconds));
        const double alpha_diff = std::abs(base.alpha - matched.alpha);
        const double theta_diff = std::abs(base.theta_prime - matched.theta_prime);
        if (arrival_diff <= kArrivalToleranceSeconds &&
            alpha_diff <= kAngleTolerance &&
            theta_diff <= kAngleTolerance) {
            used[static_cast<std::size_t>(best_index)] = true;
            comparison.matched_count += 1;
            comparison.max_abs_arrival_time_diff =
                std::max(comparison.max_abs_arrival_time_diff, arrival_diff);
            comparison.max_abs_alpha_diff = std::max(comparison.max_abs_alpha_diff, alpha_diff);
            comparison.max_abs_theta_prime_diff =
                std::max(comparison.max_abs_theta_prime_diff, theta_diff);
        } else {
            comparison.missing_from_config_count += 1;
        }
    }
    int used_count = 0;
    for (const bool is_used : used) {
        used_count += is_used ? 1 : 0;
    }
    comparison.new_in_config_count = static_cast<int>(candidate.size()) - used_count;
    comparison.matches_baseline =
        comparison.missing_from_config_count == 0 &&
        comparison.new_in_config_count == 0 &&
        baseline.size() == candidate.size();
    return comparison;
}

void print_comparison(
    const SelectedCase& selected,
    const std::string& config_name,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& baseline,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& candidate
) {
    const auto comparison = compare_to_baseline(
        selected.attempt.parent_time,
        baseline.solutions,
        candidate.solutions);
    std::cout << "Problem2EfficiencyAblationComparison\n";
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "config_name=" << config_name << '\n';
    std::cout << "baseline_solution_count=" << baseline.solutions.size() << '\n';
    std::cout << "config_solution_count=" << candidate.solutions.size() << '\n';
    std::cout << "matched_count=" << comparison.matched_count << '\n';
    std::cout << "missing_from_config_count=" << comparison.missing_from_config_count << '\n';
    std::cout << "new_in_config_count=" << comparison.new_in_config_count << '\n';
    std::cout << "max_abs_arrival_time_diff=" << comparison.max_abs_arrival_time_diff << '\n';
    std::cout << "max_abs_alpha_diff=" << comparison.max_abs_alpha_diff << '\n';
    std::cout << "max_abs_theta_prime_diff=" << comparison.max_abs_theta_prime_diff << '\n';
    std::cout << "matches_baseline=" << (comparison.matches_baseline ? 1 : 0) << '\n';
}

const RunRecord* find_record(
    const std::vector<RunRecord>& records,
    const std::string& case_name,
    const std::string& config_name
) {
    for (const auto& record : records) {
        if (record.case_name == case_name && record.config_name == config_name) {
            return &record;
        }
    }
    return nullptr;
}

double record_ms(
    const std::vector<RunRecord>& records,
    const std::string& case_name,
    const std::string& config_name
) {
    const auto* record = find_record(records, case_name, config_name);
    return record == nullptr ? -1.0 : record->result.profile.total_ms;
}

int record_solution_count(
    const std::vector<RunRecord>& records,
    const std::string& case_name,
    const std::string& config_name
) {
    const auto* record = find_record(records, case_name, config_name);
    return record == nullptr ? -1 : static_cast<int>(record->result.solutions.size());
}

void print_pair_summary(const SelectedCase& selected, const std::vector<RunRecord>& records) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::cout << "Problem2EfficiencyAblationPairSummary\n";
    std::cout << "from_planet=" << planet_params::planet_name(selected.attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(selected.attempt.to_planet) << '\n';
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "baseline_total_ms=" << record_ms(records, selected.name, "baseline_full") << '\n';
    std::cout << "fast_no_adaptive_total_ms=" << record_ms(records, selected.name, "fast_no_adaptive") << '\n';
    std::cout << "depth2_total_ms=" << record_ms(records, selected.name, "adaptive_depth2") << '\n';
    std::cout << "depth4_total_ms=" << record_ms(records, selected.name, "adaptive_depth4") << '\n';
    std::cout << "depth6_total_ms=" << record_ms(records, selected.name, "adaptive_depth6") << '\n';
    std::cout << "sample32_total_ms=" << record_ms(records, selected.name, "sample32_full_adaptive") << '\n';
    std::cout << "cheap_scout_total_ms=" << record_ms(records, selected.name, "cheap_scout_16_no_adaptive") << '\n';
    std::cout << "baseline_solution_count=" << record_solution_count(records, selected.name, "baseline_full") << '\n';
    std::cout << "fast_no_adaptive_solution_count="
              << record_solution_count(records, selected.name, "fast_no_adaptive") << '\n';
    std::cout << "depth2_solution_count=" << record_solution_count(records, selected.name, "adaptive_depth2") << '\n';
    std::cout << "depth4_solution_count=" << record_solution_count(records, selected.name, "adaptive_depth4") << '\n';
    std::cout << "depth6_solution_count=" << record_solution_count(records, selected.name, "adaptive_depth6") << '\n';
    std::cout << "sample32_solution_count="
              << record_solution_count(records, selected.name, "sample32_full_adaptive") << '\n';
    std::cout << "cheap_scout_solution_count="
              << record_solution_count(records, selected.name, "cheap_scout_16_no_adaptive") << '\n';
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
        std::cout << "problem2_efficiency_ablation_skipped_missing_table\n";
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

    problem2::Problem2GravityAssistSolverOptions bfs_problem2_options{};
    bfs_problem2_options.theta_sample_count = 64;
    bfs_problem2_options.topology_adaptive_enabled = true;
    bfs_problem2_options.topology_max_depth = 10;
    bfs_problem2_options.max_transfer_revolution = 1;
    bfs_problem2_options.max_target_revolution = 1;

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

    const auto search = bfs::run_fixed_launch_theta_beam_search_with_table(
        loader,
        launch_planet,
        terminal_planet,
        launch_time,
        initial_theta,
        allowed_first_targets,
        allowed_transfer_planets,
        initial_options,
        bfs_problem2_options,
        search_options);
    if (!search.ok) {
        std::cout << "problem2_efficiency_ablation_failed_to_collect_bfs_attempts\n";
        std::cout << "error_message=" << search.error_message << '\n';
        return 1;
    }

    const auto& attempts = search.expansion_attempt_profiles;
    std::vector<SelectedCase> selected;
    std::set<std::string> selected_names;

    const auto fastest = std::min_element(attempts.begin(), attempts.end(), [](const auto& lhs, const auto& rhs) {
        const bool lhs_normal = lhs.problem2_solver_profile.total_ms >= 10.0 &&
            lhs.problem2_solver_profile.total_ms <= 30.0;
        const bool rhs_normal = rhs.problem2_solver_profile.total_ms >= 10.0 &&
            rhs.problem2_solver_profile.total_ms <= 30.0;
        if (lhs_normal != rhs_normal) {
            return lhs_normal;
        }
        return lhs.problem2_solver_profile.total_ms < rhs.problem2_solver_profile.total_ms;
    });
    if (fastest != attempts.end()) {
        add_case_if_present(&selected, &selected_names, "fast_normal", &*fastest);
    }

    const auto slowest_venus_mars = std::max_element(
        attempts.begin(), attempts.end(), [](const auto& lhs, const auto& rhs) {
            const bool lhs_vm = same_pair(lhs, spaceship_cpp::planet_params::PlanetId::Venus,
                                          spaceship_cpp::planet_params::PlanetId::Mars);
            const bool rhs_vm = same_pair(rhs, spaceship_cpp::planet_params::PlanetId::Venus,
                                          spaceship_cpp::planet_params::PlanetId::Mars);
            if (lhs_vm != rhs_vm) {
                return !lhs_vm;
            }
            return lhs.problem2_solver_profile.total_ms < rhs.problem2_solver_profile.total_ms;
        });
    if (slowest_venus_mars != attempts.end() &&
        same_pair(*slowest_venus_mars, planet_params::PlanetId::Venus, planet_params::PlanetId::Mars)) {
        add_case_if_present(&selected, &selected_names, "slow_venus_to_mars", &*slowest_venus_mars);
    }

    add_case_if_present(
        &selected,
        &selected_names,
        "zero_raw_mars_to_mercury",
        find_first(attempts, planet_params::PlanetId::Mars, planet_params::PlanetId::Mercury, true, false));
    add_case_if_present(
        &selected,
        &selected_names,
        "zero_raw_venus_to_mercury",
        find_first(attempts, planet_params::PlanetId::Venus, planet_params::PlanetId::Mercury, true, false));

    add_case_if_present(
        &selected,
        &selected_names,
        "useful_earth_to_mars",
        find_first(attempts, planet_params::PlanetId::Earth, planet_params::PlanetId::Mars, false, true));
    add_case_if_present(
        &selected,
        &selected_names,
        "useful_mars_to_earth",
        find_first(attempts, planet_params::PlanetId::Mars, planet_params::PlanetId::Earth, false, true));
    add_case_if_present(
        &selected,
        &selected_names,
        "useful_earth_to_mercury",
        find_first(attempts, planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury, false, true));

    std::cout << "Problem2EfficiencyAblationSelectionSummary\n";
    std::cout << "attempt_count=" << attempts.size() << '\n';
    std::cout << "selected_case_count=" << selected.size() << '\n';
    for (const auto& selected_case : selected) {
        std::cout << "selected_case=" << selected_case.name
                  << " pair=" << pair_key(selected_case.attempt.from_planet, selected_case.attempt.to_planet)
                  << " baseline_attempt_ms=" << selected_case.attempt.problem2_solver_profile.total_ms
                  << " raw_edge_count=" << selected_case.attempt.raw_edge_count
                  << " accepted_edge_count=" << selected_case.attempt.accepted_edge_count << '\n';
    }

    const auto configs = make_configs();
    std::vector<RunRecord> records;
    for (const auto& selected_case : selected) {
        for (const auto& config : configs) {
            const double local_theta = bfs::global_periapsis_angle_to_problem2_local(
                selected_case.attempt.from_planet,
                selected_case.attempt.parent_incoming_theta);
            auto result = problem2::solve_problem2_gravity_assist_with_table(
                loader,
                selected_case.attempt.from_planet,
                selected_case.attempt.to_planet,
                selected_case.attempt.parent_time,
                selected_case.attempt.parent_incoming_e,
                local_theta,
                config.options);
            print_case_result(selected_case, config, result);
            records.push_back(RunRecord{
                selected_case.name,
                config.name,
                selected_case.attempt.from_planet,
                selected_case.attempt.to_planet,
                std::move(result)});
        }
    }

    for (const auto& selected_case : selected) {
        const auto* baseline = find_record(records, selected_case.name, "baseline_full");
        if (baseline == nullptr) {
            continue;
        }
        for (const auto& config : configs) {
            const auto* candidate = find_record(records, selected_case.name, config.name);
            if (candidate == nullptr) {
                continue;
            }
            print_comparison(selected_case, config.name, baseline->result, candidate->result);
        }
        print_pair_summary(selected_case, records);
    }

    std::cout << "Problem2EfficiencyAblationSummary\n";
    std::cout << "selected_case_count=" << selected.size() << '\n';
    std::cout << "config_count=" << configs.size() << '\n';
    std::cout << "ablation_ok=1\n";
    return 0;
}
