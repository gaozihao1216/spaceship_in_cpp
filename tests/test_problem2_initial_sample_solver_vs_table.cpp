#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct SelectedCase {
    std::string name;
    spaceship_cpp::bfs::BfsExpansionAttemptProfile attempt;
};

struct SampleBranchKey {
    int sample_index = 0;
    int k = 0;
    int q = 0;
    double alpha = 0.0;
    double time = 0.0;
};

struct MethodResult {
    bool ok = true;
    int sample_count = 0;
    int valid_branch_count_total = 0;
    double total_ms = 0.0;
    int nearest_query_count = 0;
    double nearest_query_total_ms = 0.0;
    double table_load_ms = 0.0;
    double table_offset_build_ms = 0.0;
    double node_read_ms = 0.0;
    double seed_generation_ms = 0.0;
    double refine_ms = 0.0;
    double derivative_attach_ms = 0.0;
    std::vector<SampleBranchKey> branches;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') return raw;
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
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
        if (!same_pair(attempt, from, to)) continue;
        if (require_zero_raw && !attempt.zero_raw) continue;
        if (require_accepted && attempt.accepted_edge_count <= 0) continue;
        return &attempt;
    }
    return nullptr;
}

void add_case(
    std::vector<SelectedCase>* cases,
    const std::string& name,
    const spaceship_cpp::bfs::BfsExpansionAttemptProfile* attempt
) {
    if (attempt != nullptr) cases->push_back(SelectedCase{name, *attempt});
}

void add_branch_keys(
    int sample_index,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    std::vector<SampleBranchKey>* keys
) {
    for (const auto& branch : branches) {
        if (!branch.valid) continue;
        keys->push_back(SampleBranchKey{
            sample_index,
            branch.transfer_revolution,
            branch.target_revolution,
            branch.encounter_global_angle,
            branch.time_of_flight_seconds});
    }
}

MethodResult run_table_method(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double encounter_time,
    int theta_sample_count
) {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    MethodResult out{};
    const auto departure_state = planet_params::planet_state_at_time(from, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(to, encounter_time);
    problem1::Problem1NearestNodeQueryOptions options{};
    options.fallback_direct_solve = true;
    const auto start = Clock::now();
    for (int i = 0; i < theta_sample_count; ++i) {
        const double theta_prime = common::kTwoPi * static_cast<double>(i) / static_cast<double>(theta_sample_count);
        const double theta_A = common::normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
        const auto query = problem1::query_problem1_from_2deg_nearest_node(
            loader, from, to, departure_state.varphi, target_state.varphi, theta_A, 1, 1, options);
        out.sample_count += 1;
        out.nearest_query_count += 1;
        out.nearest_query_total_ms += query.profile.total_ms;
        out.table_load_ms += query.profile.table_load_ms;
        out.table_offset_build_ms += query.profile.table_offset_build_ms;
        out.node_read_ms += query.profile.table_node_read_ms;
        out.seed_generation_ms += query.profile.seed_generation_ms;
        out.refine_ms += query.profile.refine_ms;
        out.derivative_attach_ms += query.profile.derivative_attach_ms;
        out.valid_branch_count_total += static_cast<int>(query.branches.size());
        add_branch_keys(i, query.branches, &out.branches);
    }
    out.total_ms = elapsed_ms(start, Clock::now());
    return out;
}

MethodResult run_direct_method(
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double encounter_time,
    int theta_sample_count
) {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    MethodResult out{};
    const auto departure_state = planet_params::planet_state_at_time(from, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(to, encounter_time);
    const auto start = Clock::now();
    for (int i = 0; i < theta_sample_count; ++i) {
        const double theta_prime = common::kTwoPi * static_cast<double>(i) / static_cast<double>(theta_sample_count);
        const double theta_A = common::normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
        const auto branches = problem1::solve_problem1_from_departure_anomalies(
            from, to, departure_state.varphi, target_state.varphi, theta_A, 1, 1);
        out.sample_count += 1;
        out.valid_branch_count_total += static_cast<int>(branches.size());
        add_branch_keys(i, branches, &out.branches);
    }
    out.total_ms = elapsed_ms(start, Clock::now());
    return out;
}

int count_matches(const std::vector<SampleBranchKey>& table, const std::vector<SampleBranchKey>& direct) {
    std::vector<bool> used(direct.size(), false);
    int matches = 0;
    for (const auto& lhs : table) {
        int best = -1;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < direct.size(); ++i) {
            if (used[i]) continue;
            const auto& rhs = direct[i];
            if (lhs.sample_index != rhs.sample_index || lhs.k != rhs.k || lhs.q != rhs.q) continue;
            const double score = std::abs(lhs.time - rhs.time) + std::abs(lhs.alpha - rhs.alpha) * 1e6;
            if (score < best_score) {
                best_score = score;
                best = static_cast<int>(i);
            }
        }
        if (best >= 0 && best_score <= 1e-3 + 1e-8 * 1e6) {
            used[static_cast<std::size_t>(best)] = true;
            matches += 1;
        }
    }
    return matches;
}

void print_method(
    const SelectedCase& selected,
    const std::string& method,
    const MethodResult& result
) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::cout << "Problem2InitialSampleMethodCase\n";
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(selected.attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(selected.attempt.to_planet) << '\n';
    std::cout << "method=" << method << '\n';
    std::cout << "theta_sample_count=64\n";
    std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "sample_count=" << result.sample_count << '\n';
    std::cout << "valid_branch_count_total=" << result.valid_branch_count_total << '\n';
    std::cout << "initial_sampling_ms=" << result.total_ms << '\n';
    std::cout << "nearest_query_total_ms=" << result.nearest_query_total_ms << '\n';
    std::cout << "nearest_query_count=" << result.nearest_query_count << '\n';
    std::cout << "mean_sample_ms=" << result.total_ms / static_cast<double>(std::max(1, result.sample_count)) << '\n';
    std::cout << "table_load_ms_mean=" << result.table_load_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
    std::cout << "table_offset_build_ms_mean=" << result.table_offset_build_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
    std::cout << "node_read_ms_mean=" << result.node_read_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
    std::cout << "seed_generation_ms_mean=" << result.seed_generation_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
    std::cout << "refine_ms_mean=" << result.refine_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
    std::cout << "derivative_attach_ms_mean=" << result.derivative_attach_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
    std::cout << "total_query_ms_mean=" << result.nearest_query_total_ms / static_cast<double>(std::max(1, result.nearest_query_count)) << '\n';
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
        std::cout << "problem2_initial_sample_solver_vs_table_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);

    const double launch_time = 0.17 * planet_params::planet_orbital_period(planet_params::PlanetId::Earth) + 1262332.782636;
    const std::vector<planet_params::PlanetId> allowed_first_targets{
        planet_params::PlanetId::Mercury, planet_params::PlanetId::Venus, planet_params::PlanetId::Mars};
    const std::vector<planet_params::PlanetId> allowed_transfer_planets{
        planet_params::PlanetId::Mercury, planet_params::PlanetId::Venus,
        planet_params::PlanetId::Earth, planet_params::PlanetId::Mars};
    bfs::InitialLaunchExpansionOptions initial_options{};
    initial_options.max_transfer_revolution = 1;
    initial_options.max_target_revolution = 1;
    initial_options.max_launch_v_inf = 7000.0;
    problem2::Problem2GravityAssistSolverOptions p2{};
    p2.theta_sample_count = 64;
    p2.topology_adaptive_enabled = true;
    p2.topology_max_depth = 10;
    p2.max_transfer_revolution = 1;
    p2.max_target_revolution = 1;
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
        launch_time, 6.169738905800, allowed_first_targets, allowed_transfer_planets,
        initial_options, p2, search_options);
    if (!search.ok) return 1;
    const auto& attempts = search.expansion_attempt_profiles;
    std::vector<SelectedCase> selected;
    const auto slowest_vm = std::max_element(attempts.begin(), attempts.end(), [](const auto& a, const auto& b) {
        const bool av = same_pair(a, spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars);
        const bool bv = same_pair(b, spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars);
        if (av != bv) return !av;
        return a.problem2_solver_profile.total_ms < b.problem2_solver_profile.total_ms;
    });
    if (slowest_vm != attempts.end() && same_pair(*slowest_vm, planet_params::PlanetId::Venus, planet_params::PlanetId::Mars)) {
        add_case(&selected, "slow_venus_to_mars", &*slowest_vm);
    }
    add_case(&selected, "zero_raw_mars_to_mercury",
             find_first(attempts, planet_params::PlanetId::Mars, planet_params::PlanetId::Mercury, true, false));
    add_case(&selected, "zero_raw_venus_to_mercury",
             find_first(attempts, planet_params::PlanetId::Venus, planet_params::PlanetId::Mercury, true, false));
    add_case(&selected, "useful_earth_to_mars",
             find_first(attempts, planet_params::PlanetId::Earth, planet_params::PlanetId::Mars, false, true));
    add_case(&selected, "useful_mars_to_earth",
             find_first(attempts, planet_params::PlanetId::Mars, planet_params::PlanetId::Earth, false, true));
    add_case(&selected, "useful_earth_to_mercury",
             find_first(attempts, planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury, false, true));

    for (const auto& selected_case : selected) {
        const auto table = run_table_method(
            loader, selected_case.attempt.from_planet, selected_case.attempt.to_planet,
            selected_case.attempt.parent_time, 64);
        const auto direct = run_direct_method(
            selected_case.attempt.from_planet, selected_case.attempt.to_planet,
            selected_case.attempt.parent_time, 64);
        print_method(selected_case, "table_nearest", table);
        print_method(selected_case, "direct_problem1_solve", direct);
        const int matched = count_matches(table.branches, direct.branches);
        std::cout << "Problem2InitialSampleMethodComparison\n";
        std::cout << "case_name=" << selected_case.name << '\n';
        std::cout << "table_initial_sampling_ms=" << table.total_ms << '\n';
        std::cout << "direct_initial_sampling_ms=" << direct.total_ms << '\n';
        std::cout << "table_valid_branch_count=" << table.valid_branch_count_total << '\n';
        std::cout << "direct_valid_branch_count=" << direct.valid_branch_count_total << '\n';
        std::cout << "matched_branch_count=" << matched << '\n';
        std::cout << "missing_from_direct_count=" << std::max(0, table.valid_branch_count_total - matched) << '\n';
        std::cout << "new_in_direct_count=" << std::max(0, direct.valid_branch_count_total - matched) << '\n';
        std::cout << "direct_faster_than_table=" << (direct.total_ms < table.total_ms ? 1 : 0) << '\n';
        std::cout << "comparison_ok=1\n";
        if (direct.total_ms / 64.0 < table.nearest_query_total_ms / 64.0) {
            std::cout << "Problem2InitialSamplingTableNotWorthItWarning\n";
            std::cout << "case_name=" << selected_case.name << '\n';
            std::cout << "table_mean_query_ms=" << table.nearest_query_total_ms / 64.0 << '\n';
            std::cout << "direct_mean_solve_ms=" << direct.total_ms / 64.0 << '\n';
            std::cout << "recommendation=consider_direct_problem1_for_initial_sampling\n";
        }
    }
    std::cout << "initial_sample_solver_vs_table_ok=1\n";
    return 0;
}
