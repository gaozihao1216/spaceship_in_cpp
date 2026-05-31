#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"
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
#include <string>
#include <vector>

namespace {

struct SelectedCase {
    std::string name;
    spaceship_cpp::bfs::BfsExpansionAttemptProfile attempt;
};

struct MatchStats {
    int matched_count = 0;
    int missing_from_direct_count = 0;
    int new_in_direct_count = 0;
    double max_abs_arrival_time_diff = 0.0;
    double max_abs_alpha_diff = 0.0;
    double max_abs_theta_prime_diff = 0.0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') return raw;
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
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
    bool require_accepted
) {
    for (const auto& attempt : attempts) {
        if (!same_pair(attempt, from, to)) continue;
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
    if (attempt != nullptr) {
        cases->push_back(SelectedCase{name, *attempt});
    }
}

double angle_diff(double lhs, double rhs) {
    constexpr double two_pi = 6.283185307179586476925286766559;
    double diff = std::fmod(lhs - rhs + 3.1415926535897932384626433832795, two_pi);
    if (diff < 0.0) diff += two_pi;
    return std::abs(diff - 3.1415926535897932384626433832795);
}

MatchStats compare_solutions(
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& table_solutions,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& direct_solutions
) {
    MatchStats stats{};
    std::vector<bool> used(direct_solutions.size(), false);
    for (const auto& lhs : table_solutions) {
        int best = -1;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < direct_solutions.size(); ++i) {
            if (used[i]) continue;
            const auto& rhs = direct_solutions[i];
            if (lhs.transfer_revolution != rhs.transfer_revolution ||
                lhs.target_revolution != rhs.target_revolution) {
                continue;
            }
            const double theta_diff = angle_diff(lhs.theta_prime, rhs.theta_prime);
            const double alpha_diff = angle_diff(lhs.alpha, rhs.alpha);
            const double time_diff = std::abs(lhs.target_time_seconds - rhs.target_time_seconds);
            const double score = theta_diff * 1e9 + alpha_diff * 1e9 + time_diff;
            if (score < best_score) {
                best_score = score;
                best = static_cast<int>(i);
            }
        }
        if (best >= 0) {
            const auto& rhs = direct_solutions[static_cast<std::size_t>(best)];
            const double theta_diff = angle_diff(lhs.theta_prime, rhs.theta_prime);
            const double alpha_diff = angle_diff(lhs.alpha, rhs.alpha);
            const double time_diff = std::abs(lhs.target_time_seconds - rhs.target_time_seconds);
            if (theta_diff <= 1e-6 && alpha_diff <= 1e-6 && time_diff <= 1e3) {
                used[static_cast<std::size_t>(best)] = true;
                stats.matched_count += 1;
                stats.max_abs_theta_prime_diff = std::max(stats.max_abs_theta_prime_diff, theta_diff);
                stats.max_abs_alpha_diff = std::max(stats.max_abs_alpha_diff, alpha_diff);
                stats.max_abs_arrival_time_diff = std::max(stats.max_abs_arrival_time_diff, time_diff);
            }
        }
    }
    stats.missing_from_direct_count = static_cast<int>(table_solutions.size()) - stats.matched_count;
    stats.new_in_direct_count = static_cast<int>(direct_solutions.size()) - stats.matched_count;
    return stats;
}

void print_case(
    const SelectedCase& selected,
    const std::string& config_name,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& result
) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::cout << "Problem2Problem1DirectSeedModeCase\n";
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(selected.attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(selected.attempt.to_planet) << '\n';
    std::cout << "config_name=" << config_name << '\n';
    std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "solution_count=" << result.solutions.size() << '\n';
    std::cout << "total_ms=" << result.profile.total_ms << '\n';
    std::cout << "initial_sampling_ms=" << result.profile.initial_sampling_ms << '\n';
    std::cout << "topology_adaptive_sampling_ms=" << result.profile.topology_adaptive_sampling_ms << '\n';
    std::cout << "nearest_query_count=" << result.profile.problem1_nearest_node_query_count << '\n';
    std::cout << "nearest_query_total_ms=" << result.profile.nearest_query_total_ms << '\n';
    std::cout << "direct_problem1_seed_solve_count=" << result.profile.direct_problem1_seed_solve_count << '\n';
    std::cout << "direct_problem1_seed_solve_ms=" << result.profile.direct_problem1_seed_solve_ms << '\n';
    std::cout << "direct_problem1_seed_derivative_attach_ms=" << result.profile.direct_problem1_seed_derivative_attach_ms << '\n';
    std::cout << "direct_problem1_seed_branch_count=" << result.profile.direct_problem1_seed_branch_count << '\n';
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
        std::cout << "problem2_problem1_direct_seed_mode_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);

    const double launch_time =
        0.17 * planet_params::planet_orbital_period(planet_params::PlanetId::Earth) + 1262332.782636;
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
    p2.topology_max_depth = 6;
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
    if (!search.ok) {
        std::cout << "problem2_problem1_direct_seed_mode_search_failed\n";
        return 1;
    }

    const auto& attempts = search.expansion_attempt_profiles;
    std::vector<SelectedCase> selected;
    add_case(&selected, "useful_earth_to_mars",
             find_first(attempts, planet_params::PlanetId::Earth, planet_params::PlanetId::Mars, true));
    add_case(&selected, "useful_mars_to_earth",
             find_first(attempts, planet_params::PlanetId::Mars, planet_params::PlanetId::Earth, true));
    add_case(&selected, "useful_earth_to_mercury",
             find_first(attempts, planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury, true));
    add_case(&selected, "slow_venus_to_mars",
             find_first(attempts, planet_params::PlanetId::Venus, planet_params::PlanetId::Mars, false));
    add_case(&selected, "zero_mars_to_mercury",
             find_first(attempts, planet_params::PlanetId::Mars, planet_params::PlanetId::Mercury, false));
    add_case(&selected, "zero_venus_to_mercury",
             find_first(attempts, planet_params::PlanetId::Venus, planet_params::PlanetId::Mercury, false));

    double table_total_ms = 0.0;
    double direct_total_ms = 0.0;
    int table_solution_total = 0;
    int direct_solution_total = 0;
    int missing_from_direct_total = 0;
    int new_in_direct_total = 0;

    for (const auto& selected_case : selected) {
        problem2::Problem2GravityAssistSolverOptions table_options = p2;
        table_options.use_problem1_solve_for_problem2_seed = false;
        const auto table_result = problem2::solve_problem2_gravity_assist_with_table(
            loader,
            selected_case.attempt.from_planet,
            selected_case.attempt.to_planet,
            selected_case.attempt.parent_time,
            selected_case.attempt.parent_incoming_e,
            selected_case.attempt.parent_incoming_theta,
            table_options);

        problem2::Problem2GravityAssistSolverOptions direct_options = p2;
        direct_options.use_problem1_solve_for_problem2_seed = true;
        const auto direct_result = problem2::solve_problem2_gravity_assist_with_table(
            loader,
            selected_case.attempt.from_planet,
            selected_case.attempt.to_planet,
            selected_case.attempt.parent_time,
            selected_case.attempt.parent_incoming_e,
            selected_case.attempt.parent_incoming_theta,
            direct_options);

        print_case(selected_case, "table_seed", table_result);
        print_case(selected_case, "problem1_solve_seed", direct_result);

        const auto cmp = compare_solutions(table_result.solutions, direct_result.solutions);
        std::cout << "Problem2Problem1DirectSeedModeComparison\n";
        std::cout << "case_name=" << selected_case.name << '\n';
        std::cout << "table_solution_count=" << table_result.solutions.size() << '\n';
        std::cout << "direct_solution_count=" << direct_result.solutions.size() << '\n';
        std::cout << "matched_count=" << cmp.matched_count << '\n';
        std::cout << "missing_from_direct_count=" << cmp.missing_from_direct_count << '\n';
        std::cout << "new_in_direct_count=" << cmp.new_in_direct_count << '\n';
        std::cout << "max_abs_arrival_time_diff=" << cmp.max_abs_arrival_time_diff << '\n';
        std::cout << "max_abs_alpha_diff=" << cmp.max_abs_alpha_diff << '\n';
        std::cout << "max_abs_theta_prime_diff=" << cmp.max_abs_theta_prime_diff << '\n';
        std::cout << "matches_table_seed=" << (cmp.missing_from_direct_count == 0 && cmp.new_in_direct_count == 0 ? 1 : 0) << '\n';

        table_total_ms += table_result.profile.total_ms;
        direct_total_ms += direct_result.profile.total_ms;
        table_solution_total += static_cast<int>(table_result.solutions.size());
        direct_solution_total += static_cast<int>(direct_result.solutions.size());
        missing_from_direct_total += cmp.missing_from_direct_count;
        new_in_direct_total += cmp.new_in_direct_count;
    }

    const int case_count = static_cast<int>(selected.size());
    std::cout << "Problem2Problem1DirectSeedModeSummary\n";
    std::cout << "case_count=" << case_count << '\n';
    std::cout << "table_mean_ms=" << table_total_ms / static_cast<double>(std::max(1, case_count)) << '\n';
    std::cout << "direct_mean_ms=" << direct_total_ms / static_cast<double>(std::max(1, case_count)) << '\n';
    std::cout << "table_solution_total=" << table_solution_total << '\n';
    std::cout << "direct_solution_total=" << direct_solution_total << '\n';
    std::cout << "missing_from_direct_total=" << missing_from_direct_total << '\n';
    std::cout << "new_in_direct_total=" << new_in_direct_total << '\n';
    std::cout << "problem2_problem1_direct_seed_mode_ok=1\n";
    return 0;
}
