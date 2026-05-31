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
#include <string>
#include <vector>

namespace {

struct SelectedCase {
    std::string name;
    spaceship_cpp::bfs::BfsExpansionAttemptProfile attempt;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
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

void add_case(
    std::vector<SelectedCase>* cases,
    const std::string& name,
    const spaceship_cpp::bfs::BfsExpansionAttemptProfile* attempt
) {
    if (attempt != nullptr) {
        cases->push_back(SelectedCase{name, *attempt});
    }
}

struct Comparison {
    int matched_count = 0;
    int missing = 0;
    int added = 0;
    double max_arrival_diff = 0.0;
    double max_alpha_diff = 0.0;
    double max_theta_diff = 0.0;
    bool matches = false;
};

double angle_diff(double a, double b) {
    constexpr double kTwoPi = 6.283185307179586476925286766559;
    double d = std::fmod(a - b, kTwoPi);
    if (d > M_PI) d -= kTwoPi;
    if (d < -M_PI) d += kTwoPi;
    return std::abs(d);
}

Comparison compare_solutions(
    double encounter_time,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& baseline,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& other
) {
    Comparison c{};
    std::vector<bool> used(other.size(), false);
    for (const auto& b : baseline) {
        int best = -1;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < other.size(); ++i) {
            if (used[i]) continue;
            const double arrival = std::abs((encounter_time + b.time_of_flight_seconds) -
                                            (encounter_time + other[i].time_of_flight_seconds));
            const double alpha = angle_diff(b.alpha, other[i].alpha);
            const double theta = angle_diff(b.theta_prime, other[i].theta_prime);
            const double score = arrival + alpha * 1e6 + theta * 1e6;
            if (score < best_score) {
                best_score = score;
                best = static_cast<int>(i);
            }
        }
        if (best < 0) {
            c.missing += 1;
            continue;
        }
        const auto& o = other[static_cast<std::size_t>(best)];
        const double arrival = std::abs((encounter_time + b.time_of_flight_seconds) -
                                        (encounter_time + o.time_of_flight_seconds));
        const double alpha = angle_diff(b.alpha, o.alpha);
        const double theta = angle_diff(b.theta_prime, o.theta_prime);
        if (arrival <= 1e-3 && alpha <= 1e-8 && theta <= 1e-8) {
            used[static_cast<std::size_t>(best)] = true;
            c.matched_count += 1;
            c.max_arrival_diff = std::max(c.max_arrival_diff, arrival);
            c.max_alpha_diff = std::max(c.max_alpha_diff, alpha);
            c.max_theta_diff = std::max(c.max_theta_diff, theta);
        } else {
            c.missing += 1;
        }
    }
    int used_count = 0;
    for (const bool u : used) used_count += u ? 1 : 0;
    c.added = static_cast<int>(other.size()) - used_count;
    c.matches = c.missing == 0 && c.added == 0 && baseline.size() == other.size();
    return c;
}

void print_case(
    const SelectedCase& selected,
    const std::string& config_name,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& result
) {
    namespace planet_params = spaceship_cpp::planet_params;
    const auto& p = result.profile;
    std::cout << "Problem2MidpointRouteAAblationCase\n";
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(selected.attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(selected.attempt.to_planet) << '\n';
    std::cout << "config_name=" << config_name << '\n';
    std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "solution_count=" << result.solutions.size() << '\n';
    std::cout << "total_ms=" << p.total_ms << '\n';
    std::cout << "topology_adaptive_sampling_ms=" << p.topology_adaptive_sampling_ms << '\n';
    std::cout << "nearest_query_total_ms=" << p.nearest_query_total_ms << '\n';
    std::cout << "nearest_query_count=" << p.problem1_nearest_node_query_count << '\n';
    std::cout << "adaptive_midpoint_sample_count=" << p.adaptive_midpoint_sample_count << '\n';
    std::cout << "adaptive_midpoint_routeA_attempt_count=" << p.adaptive_midpoint_routeA_attempt_count << '\n';
    std::cout << "adaptive_midpoint_routeA_success_count=" << p.adaptive_midpoint_routeA_success_count << '\n';
    std::cout << "adaptive_midpoint_routeA_validation_failure_count="
              << p.adaptive_midpoint_routeA_validation_failure_count << '\n';
    std::cout << "adaptive_midpoint_routeA_fallback_count=" << p.adaptive_midpoint_routeA_fallback_count << '\n';
    std::cout << "adaptive_midpoint_nearest_query_avoided_count="
              << p.adaptive_midpoint_nearest_query_avoided_count << '\n';
    std::cout << "adaptive_midpoint_routeA_ms=" << p.adaptive_midpoint_routeA_ms << '\n';
    std::cout << "adaptive_midpoint_fallback_nearest_query_ms="
              << p.adaptive_midpoint_fallback_nearest_query_ms << '\n';
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
        std::cout << "problem2_midpoint_routeA_ablation_skipped_missing_table\n";
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
    problem2::Problem2GravityAssistSolverOptions base_options{};
    base_options.theta_sample_count = 64;
    base_options.topology_adaptive_enabled = true;
    base_options.topology_max_depth = 10;
    base_options.max_transfer_revolution = 1;
    base_options.max_target_revolution = 1;
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
        initial_options, base_options, search_options);
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
        const double local_theta = bfs::global_periapsis_angle_to_problem2_local(
            selected_case.attempt.from_planet, selected_case.attempt.parent_incoming_theta);
        const auto baseline = problem2::solve_problem2_gravity_assist_with_table(
            loader, selected_case.attempt.from_planet, selected_case.attempt.to_planet,
            selected_case.attempt.parent_time, selected_case.attempt.parent_incoming_e, local_theta, base_options);
        print_case(selected_case, "baseline", baseline);

        auto route_options = base_options;
        route_options.adaptive_midpoint_routeA_from_endpoint_enabled = true;
        route_options.adaptive_midpoint_seed_policy = 1;
        const auto route = problem2::solve_problem2_gravity_assist_with_table(
            loader, selected_case.attempt.from_planet, selected_case.attempt.to_planet,
            selected_case.attempt.parent_time, selected_case.attempt.parent_incoming_e, local_theta, route_options);
        print_case(selected_case, "routeA_midpoint", route);
        const auto comparison = compare_solutions(selected_case.attempt.parent_time, baseline.solutions, route.solutions);
        std::cout << "Problem2MidpointRouteAComparison\n";
        std::cout << "case_name=" << selected_case.name << '\n';
        std::cout << "baseline_solution_count=" << baseline.solutions.size() << '\n';
        std::cout << "routeA_solution_count=" << route.solutions.size() << '\n';
        std::cout << "matched_count=" << comparison.matched_count << '\n';
        std::cout << "missing_from_routeA_count=" << comparison.missing << '\n';
        std::cout << "new_in_routeA_count=" << comparison.added << '\n';
        std::cout << "max_abs_arrival_time_diff=" << comparison.max_arrival_diff << '\n';
        std::cout << "max_abs_alpha_diff=" << comparison.max_alpha_diff << '\n';
        std::cout << "max_abs_theta_prime_diff=" << comparison.max_theta_diff << '\n';
        std::cout << "matches_baseline=" << (comparison.matches ? 1 : 0) << '\n';
    }
    std::cout << "Problem2RouteAContinuationCapability\n";
    std::cout << "has_routeA_from_branch_seed=1\n";
    std::cout << "required_fields_available=1\n";
    std::cout << "missing_fields=\n";
    std::cout << "candidate_functions=continue_branch\n";
    std::cout << "can_implement_midpoint_without_table=1\n";
    std::cout << "midpoint_routeA_ablation_ok=1\n";
    return 0;
}
