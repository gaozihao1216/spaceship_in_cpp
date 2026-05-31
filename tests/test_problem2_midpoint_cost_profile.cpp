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
#include <string>
#include <vector>

namespace {

struct SelectedCase {
    std::string name;
    spaceship_cpp::bfs::BfsExpansionAttemptProfile attempt;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') return raw;
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

bool same_pair(
    const spaceship_cpp::bfs::BfsExpansionAttemptProfile& a,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to
) {
    return a.from_planet == from && a.to_planet == to;
}

const spaceship_cpp::bfs::BfsExpansionAttemptProfile* find_first(
    const std::vector<spaceship_cpp::bfs::BfsExpansionAttemptProfile>& attempts,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    bool zero_raw,
    bool accepted
) {
    for (const auto& a : attempts) {
        if (!same_pair(a, from, to)) continue;
        if (zero_raw && !a.zero_raw) continue;
        if (accepted && a.accepted_edge_count <= 0) continue;
        return &a;
    }
    return nullptr;
}

void add_case(std::vector<SelectedCase>* cases, const std::string& name,
              const spaceship_cpp::bfs::BfsExpansionAttemptProfile* attempt) {
    if (attempt != nullptr) cases->push_back({name, *attempt});
}

void print_query_details(
    const std::string& case_name,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& profile
) {
    for (const auto& q : profile.theta_sample_traces) {
        std::cout << "Problem1NearestNodeQueryDetailedProfile\n";
        std::cout << "case_name=" << case_name << '\n';
        std::cout << "phase=" << (q.adaptive ? "adaptive" : "initial") << '\n';
        std::cout << "theta_prime=" << q.theta_prime << '\n';
        std::cout << "total_ms=" << q.total_ms << '\n';
        std::cout << "table_offset_build_ms=" << q.table_offset_build_ms << '\n';
        std::cout << "table_node_read_ms=" << q.table_node_read_ms << '\n';
        std::cout << "seed_generation_ms=" << q.seed_generation_ms << '\n';
        std::cout << "refine_ms=" << q.refine_ms << '\n';
        std::cout << "derivative_attach_ms=" << q.derivative_attach_ms << '\n';
        std::cout << "dedup_ms=" << q.dedup_ms << '\n';
        std::cout << "fallback_direct_solve_ms=" << q.fallback_direct_solve_ms << '\n';
        std::cout << "loaded_branch_count=" << q.loaded_branch_count << '\n';
        std::cout << "seed_branch_count=" << q.seed_branch_count << '\n';
        std::cout << "refine_attempt_count=" << q.refine_attempt_count << '\n';
        std::cout << "refine_success_count=" << q.refine_success_count << '\n';
        std::cout << "derivative_attach_attempt_count=" << q.derivative_attach_attempt_count << '\n';
        std::cout << "fallback_direct_solve_count=" << q.fallback_direct_solve_count << '\n';
    }
}

void print_routeA_details(
    const std::string& case_name,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& profile
) {
    namespace pp = spaceship_cpp::planet_params;
    for (const auto& r : profile.routeA_midpoint_traces) {
        std::cout << "Problem2RouteAMidpointDetailedProfile\n";
        std::cout << "case_name=" << case_name << '\n';
        std::cout << "from_planet=" << pp::planet_name(from) << '\n';
        std::cout << "to_planet=" << pp::planet_name(to) << '\n';
        std::cout << "theta_mid=" << r.theta_mid << '\n';
        std::cout << "seed_policy=" << r.seed_policy << '\n';
        std::cout << "seed_endpoint=" << (r.seed_endpoint_left ? "left" : "right") << '\n';
        std::cout << "seed_valid_branch_count=" << r.seed_valid_branch_count << '\n';
        std::cout << "routeA_attempt_count_for_midpoint=" << r.routeA_attempt_count_for_midpoint << '\n';
        std::cout << "routeA_success_count_for_midpoint=" << r.routeA_success_count_for_midpoint << '\n';
        std::cout << "routeA_validation_failure_count_for_midpoint="
                  << r.routeA_validation_failure_count_for_midpoint << '\n';
        std::cout << "fallback_used=" << (r.fallback_used ? 1 : 0) << '\n';
        std::cout << "routeA_total_ms_for_midpoint=" << r.routeA_total_ms_for_midpoint << '\n';
        std::cout << "routeA_mean_ms_per_branch="
                  << (r.routeA_attempt_count_for_midpoint > 0
                          ? r.routeA_total_ms_for_midpoint / r.routeA_attempt_count_for_midpoint
                          : 0.0)
                  << '\n';
        std::cout << "validation_ms_for_midpoint=" << r.validation_ms_for_midpoint << '\n';
        std::cout << "fallback_nearest_query_ms_for_midpoint="
                  << r.fallback_nearest_query_ms_for_midpoint << '\n';
    }
}

double sum_adaptive_query_ms(const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& p) {
    double total = 0.0;
    for (const auto& q : p.theta_sample_traces) {
        if (q.adaptive) total += q.total_ms;
    }
    return total;
}

int adaptive_query_count(const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& p) {
    int total = 0;
    for (const auto& q : p.theta_sample_traces) {
        total += q.adaptive ? 1 : 0;
    }
    return total;
}

double sum_routeA_trace_ms(const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& p) {
    double total = 0.0;
    for (const auto& r : p.routeA_midpoint_traces) total += r.routeA_total_ms_for_midpoint;
    return total;
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace pp = spaceship_cpp::planet_params;
    namespace p1 = spaceship_cpp::problem1;
    namespace p2 = spaceship_cpp::problem2;
    namespace tr = spaceship_cpp::trajectory;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_midpoint_cost_profile_skipped_missing_table\n";
        return 0;
    }
    const auto loader = p1::Problem1RootTable2DegLoader::open(table_dir);

    const double launch_time = 0.17 * pp::planet_orbital_period(pp::PlanetId::Earth) + 1262332.782636;
    const std::vector<pp::PlanetId> first{pp::PlanetId::Mercury, pp::PlanetId::Venus, pp::PlanetId::Mars};
    const std::vector<pp::PlanetId> transfer{
        pp::PlanetId::Mercury, pp::PlanetId::Venus, pp::PlanetId::Earth, pp::PlanetId::Mars};
    bfs::InitialLaunchExpansionOptions initial{};
    initial.max_transfer_revolution = 1;
    initial.max_target_revolution = 1;
    initial.max_launch_v_inf = 7000.0;
    p2::Problem2GravityAssistSolverOptions base{};
    base.theta_sample_count = 64;
    base.topology_adaptive_enabled = true;
    base.topology_max_depth = 10;
    base.max_transfer_revolution = 1;
    base.max_target_revolution = 1;
    bfs::FixedLaunchThetaSearchOptions search_options{};
    search_options.max_depth = 5;
    search_options.enable_beam_pruning = false;
    search_options.beam_width = 3;
    search_options.flyby_physical_options.enabled = true;
    search_options.flyby_physical_options.mode = tr::FlybyPhysicalFilterMode::Enforce;
    search_options.flyby_physical_options.min_flyby_altitude_m = 300000.0;
    search_options.max_launch_v_inf = 7000.0;
    search_options.continue_after_reaching_terminal = true;
    const auto search = bfs::run_fixed_launch_theta_beam_search_with_table(
        loader, pp::PlanetId::Earth, pp::PlanetId::Mercury, launch_time, 6.169738905800,
        first, transfer, initial, base, search_options);
    if (!search.ok) return 1;

    const auto& attempts = search.expansion_attempt_profiles;
    std::vector<SelectedCase> selected;
    const auto slowest_vm = std::max_element(attempts.begin(), attempts.end(), [](const auto& a, const auto& b) {
        const bool av = same_pair(a, spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars);
        const bool bv = same_pair(b, spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars);
        if (av != bv) return !av;
        return a.problem2_solver_profile.total_ms < b.problem2_solver_profile.total_ms;
    });
    if (slowest_vm != attempts.end() && same_pair(*slowest_vm, pp::PlanetId::Venus, pp::PlanetId::Mars)) {
        add_case(&selected, "slow_venus_to_mars", &*slowest_vm);
    }
    add_case(&selected, "useful_mars_to_earth",
             find_first(attempts, pp::PlanetId::Mars, pp::PlanetId::Earth, false, true));
    add_case(&selected, "useful_earth_to_mars",
             find_first(attempts, pp::PlanetId::Earth, pp::PlanetId::Mars, false, true));
    add_case(&selected, "zero_raw_mars_to_mercury",
             find_first(attempts, pp::PlanetId::Mars, pp::PlanetId::Mercury, true, false));
    add_case(&selected, "zero_raw_venus_to_mercury",
             find_first(attempts, pp::PlanetId::Venus, pp::PlanetId::Mercury, true, false));
    add_case(&selected, "useful_earth_to_mercury",
             find_first(attempts, pp::PlanetId::Earth, pp::PlanetId::Mercury, false, true));

    for (const auto& c : selected) {
        const double local_theta = bfs::global_periapsis_angle_to_problem2_local(
            c.attempt.from_planet, c.attempt.parent_incoming_theta);
        auto baseline_options = base;
        baseline_options.adaptive_interval_trace_enabled = true;
        const auto baseline = p2::solve_problem2_gravity_assist_with_table(
            loader, c.attempt.from_planet, c.attempt.to_planet,
            c.attempt.parent_time, c.attempt.parent_incoming_e, local_theta, baseline_options);
        print_query_details(c.name, baseline.profile);

        auto route_options = baseline_options;
        route_options.adaptive_midpoint_routeA_from_endpoint_enabled = true;
        route_options.adaptive_midpoint_seed_policy = 1;
        const auto route = p2::solve_problem2_gravity_assist_with_table(
            loader, c.attempt.from_planet, c.attempt.to_planet,
            c.attempt.parent_time, c.attempt.parent_incoming_e, local_theta, route_options);
        print_routeA_details(c.name, c.attempt.from_planet, c.attempt.to_planet, route.profile);

        const int baseline_count = adaptive_query_count(baseline.profile);
        const double baseline_ms = sum_adaptive_query_ms(baseline.profile);
        const double routeA_ms = sum_routeA_trace_ms(route.profile);
        const double fallback_ms = route.profile.adaptive_midpoint_fallback_nearest_query_ms;
        std::cout << "Problem2MidpointCostComparison\n";
        std::cout << "case_name=" << c.name << '\n';
        std::cout << "baseline_midpoint_nearest_query_count=" << baseline_count << '\n';
        std::cout << "baseline_midpoint_nearest_query_total_ms=" << baseline_ms << '\n';
        std::cout << "baseline_midpoint_nearest_query_mean_ms="
                  << (baseline_count > 0 ? baseline_ms / baseline_count : 0.0) << '\n';
        std::cout << "routeA_midpoint_attempt_count=" << route.profile.adaptive_midpoint_routeA_attempt_count << '\n';
        std::cout << "routeA_midpoint_total_ms=" << routeA_ms << '\n';
        std::cout << "routeA_midpoint_mean_ms="
                  << (route.profile.adaptive_midpoint_routeA_attempt_count > 0
                          ? routeA_ms / route.profile.adaptive_midpoint_routeA_attempt_count
                          : 0.0)
                  << '\n';
        std::cout << "routeA_success_count=" << route.profile.adaptive_midpoint_routeA_success_count << '\n';
        std::cout << "routeA_fallback_count=" << route.profile.adaptive_midpoint_routeA_fallback_count << '\n';
        std::cout << "nearest_query_avoided_count="
                  << route.profile.adaptive_midpoint_nearest_query_avoided_count << '\n';
        std::cout << "net_time_saved_ms=" << baseline_ms - routeA_ms - fallback_ms << '\n';
        std::cout << "routeA_overhead_ms=" << routeA_ms << '\n';
    }
    std::cout << "midpoint_cost_profile_ok=1\n";
    return 0;
}
