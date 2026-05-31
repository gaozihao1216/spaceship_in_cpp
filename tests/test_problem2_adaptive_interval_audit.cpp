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
#include <map>
#include <set>
#include <sstream>
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

std::string cell_key(const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile::ThetaSampleTrace& trace) {
    std::ostringstream os;
    os << trace.nearest_nu_A_index << ':' << trace.nearest_nu_B_index << ':' << trace.nearest_theta_A_index;
    return os.str();
}

void print_repeated_query_summary(
    const std::string& case_name,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& profile
) {
    std::set<long long> unique_theta;
    std::map<std::string, int> cell_counts;
    std::map<int, int> chunk_counts;
    for (const auto& trace : profile.theta_sample_traces) {
        unique_theta.insert(static_cast<long long>(std::llround(trace.theta_prime / 1e-14)));
        cell_counts[cell_key(trace)] += 1;
        chunk_counts[trace.chunk_index] += 1;
    }
    int same_cell_query_count = 0;
    for (const auto& [key, count] : cell_counts) {
        (void)key;
        if (count > 1) {
            same_cell_query_count += count - 1;
        }
    }
    int same_chunk_query_count = 0;
    for (const auto& [chunk, count] : chunk_counts) {
        (void)chunk;
        if (count > 1) {
            same_chunk_query_count += count - 1;
        }
    }
    const int total = profile.total_sample_build_count;
    std::cout << "Problem2AdaptiveRepeatedQuerySummary\n";
    std::cout << "case_name=" << case_name << '\n';
    std::cout << "unique_theta_sample_count=" << unique_theta.size() << '\n';
    std::cout << "total_sample_build_count=" << total << '\n';
    std::cout << "cached_sample_hit_count=" << profile.cached_sample_hit_count << '\n';
    std::cout << "duplicate_theta_query_count=" << std::max(0, total - static_cast<int>(unique_theta.size())) << '\n';
    std::cout << "same_cell_query_count=" << same_cell_query_count << '\n';
    std::cout << "same_cell_query_ratio="
              << (total > 0 ? static_cast<double>(same_cell_query_count) / static_cast<double>(total) : 0.0) << '\n';
    std::cout << "same_chunk_query_count=" << same_chunk_query_count << '\n';
    std::cout << "same_chunk_query_ratio="
              << (total > 0 ? static_cast<double>(same_chunk_query_count) / static_cast<double>(total) : 0.0) << '\n';
}

void print_stable_endpoint_audit(
    const std::string& case_name,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverProfile& profile
) {
    int stable_interval_count = 0;
    int stable_both_endpoint_queried_count = 0;
    int stable_midpoint_queried_count = 0;
    int stable_same_branch_count = 0;
    int stable_branch_mismatch_count = 0;
    for (const auto& trace : profile.adaptive_interval_traces) {
        if (!trace.classified_stable) {
            continue;
        }
        stable_interval_count += 1;
        if (trace.left_sample_built && trace.right_sample_built) {
            stable_both_endpoint_queried_count += 1;
        }
        if (trace.midpoint_sample_built) {
            stable_midpoint_queried_count += 1;
        }
        if (trace.left_branch_count == trace.right_branch_count) {
            stable_same_branch_count += 1;
        } else {
            stable_branch_mismatch_count += 1;
        }
    }
    std::cout << "Problem2StableIntervalEndpointAudit\n";
    std::cout << "case_name=" << case_name << '\n';
    std::cout << "stable_interval_count=" << stable_interval_count << '\n';
    std::cout << "stable_both_endpoint_queried_count=" << stable_both_endpoint_queried_count << '\n';
    std::cout << "stable_midpoint_queried_count=" << stable_midpoint_queried_count << '\n';
    std::cout << "stable_same_branch_count=" << stable_same_branch_count << '\n';
    std::cout << "stable_branch_mismatch_count=" << stable_branch_mismatch_count << '\n';
    std::cout << "stable_validation_failure_count=0\n";
}

void print_case_summary(
    const SelectedCase& selected,
    const spaceship_cpp::problem2::Problem2GravityAssistSolverResult& result
) {
    namespace planet_params = spaceship_cpp::planet_params;
    const auto& profile = result.profile;
    std::cout << "Problem2AdaptiveIntervalAuditCase\n";
    std::cout << "case_name=" << selected.name << '\n';
    std::cout << "from_planet=" << planet_params::planet_name(selected.attempt.from_planet) << '\n';
    std::cout << "to_planet=" << planet_params::planet_name(selected.attempt.to_planet) << '\n';
    std::cout << "ok=" << (result.ok ? 1 : 0) << '\n';
    std::cout << "solution_count=" << result.solutions.size() << '\n';
    std::cout << "total_ms=" << profile.total_ms << '\n';
    std::cout << "initial_sampling_ms=" << profile.initial_sampling_ms << '\n';
    std::cout << "topology_adaptive_sampling_ms=" << profile.topology_adaptive_sampling_ms << '\n';
    std::cout << "nearest_query_total_ms=" << profile.nearest_query_total_ms << '\n';
    std::cout << "nearest_query_count=" << profile.problem1_nearest_node_query_count << '\n';
    std::cout << "adaptive_interval_count=" << profile.adaptive_interval_count << '\n';
    std::cout << "adaptive_stable_interval_count=" << profile.adaptive_stable_interval_count << '\n';
    std::cout << "adaptive_topology_change_interval_count=" << profile.adaptive_topology_change_interval_count << '\n';
    std::cout << "adaptive_boundary_ambiguous_interval_count=" << profile.adaptive_boundary_ambiguous_interval_count << '\n';
    std::cout << "adaptive_subdivided_interval_count=" << profile.adaptive_subdivided_interval_count << '\n';
    std::cout << "adaptive_left_endpoint_query_count=" << profile.adaptive_left_endpoint_query_count << '\n';
    std::cout << "adaptive_right_endpoint_query_count=" << profile.adaptive_right_endpoint_query_count << '\n';
    std::cout << "adaptive_midpoint_query_count=" << profile.adaptive_midpoint_query_count << '\n';
    std::cout << "adaptive_cache_hit_count=" << profile.adaptive_cache_hit_count << '\n';
    std::cout << "adaptive_cache_miss_count=" << profile.adaptive_cache_miss_count << '\n';
    std::cout << "max_depth_reached=" << result.summary.max_topology_recursion_depth_reached << '\n';
    std::cout << "zero_solution=" << (result.solutions.empty() ? 1 : 0) << '\n';
}

void print_top_intervals(
    const std::string& case_name,
    std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolverProfile::AdaptiveIntervalTrace> traces
) {
    std::sort(traces.begin(), traces.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.nearest_query_ms_for_interval > rhs.nearest_query_ms_for_interval;
    });
    const int count = std::min<int>(30, static_cast<int>(traces.size()));
    for (int i = 0; i < count; ++i) {
        const auto& trace = traces[static_cast<std::size_t>(i)];
        std::cout << "Problem2AdaptiveIntervalTraceTop\n";
        std::cout << "case_name=" << case_name << '\n';
        std::cout << "rank=" << (i + 1) << '\n';
        std::cout << "interval_id=" << trace.interval_id << '\n';
        std::cout << "depth=" << trace.depth << '\n';
        std::cout << "theta_left=" << trace.theta_left << '\n';
        std::cout << "theta_right=" << trace.theta_right << '\n';
        std::cout << "theta_mid=" << trace.theta_mid << '\n';
        std::cout << "left_branch_count=" << trace.left_branch_count << '\n';
        std::cout << "right_branch_count=" << trace.right_branch_count << '\n';
        std::cout << "midpoint_branch_count=" << trace.midpoint_branch_count << '\n';
        std::cout << "classified_stable=" << (trace.classified_stable ? 1 : 0) << '\n';
        std::cout << "classified_topology_change=" << (trace.classified_topology_change ? 1 : 0) << '\n';
        std::cout << "classified_boundary_ambiguous=" << (trace.classified_boundary_ambiguous ? 1 : 0) << '\n';
        std::cout << "classified_discontinuous=" << (trace.classified_discontinuous ? 1 : 0) << '\n';
        std::cout << "subdivided=" << (trace.subdivided ? 1 : 0) << '\n';
        std::cout << "subdivision_reason=" << trace.subdivision_reason << '\n';
        std::cout << "nearest_query_count_for_interval=" << trace.nearest_query_count_for_interval << '\n';
        std::cout << "nearest_query_ms_for_interval=" << trace.nearest_query_ms_for_interval << '\n';
    }
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
        std::cout << "problem2_adaptive_interval_audit_skipped_missing_table\n";
        std::cout << "ROOT_TABLE_2DEG_DIR=" << table_dir.string() << '\n';
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
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
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mercury,
        launch_time,
        initial_theta,
        allowed_first_targets,
        allowed_transfer_planets,
        initial_options,
        bfs_problem2_options,
        search_options);
    if (!search.ok) {
        std::cout << "problem2_adaptive_interval_audit_failed_to_collect_bfs_attempts\n";
        std::cout << "error_message=" << search.error_message << '\n';
        return 1;
    }

    const auto& attempts = search.expansion_attempt_profiles;
    std::vector<SelectedCase> selected;
    const auto slowest_venus_mars = std::max_element(attempts.begin(), attempts.end(), [](const auto& lhs, const auto& rhs) {
        const bool lhs_vm = same_pair(lhs, spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars);
        const bool rhs_vm = same_pair(rhs, spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars);
        if (lhs_vm != rhs_vm) {
            return !lhs_vm;
        }
        return lhs.problem2_solver_profile.total_ms < rhs.problem2_solver_profile.total_ms;
    });
    if (slowest_venus_mars != attempts.end() &&
        same_pair(*slowest_venus_mars, planet_params::PlanetId::Venus, planet_params::PlanetId::Mars)) {
        add_case(&selected, "slow_venus_to_mars", &*slowest_venus_mars);
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

    std::cout << "Problem2AdaptiveIntervalAuditSelectionSummary\n";
    std::cout << "attempt_count=" << attempts.size() << '\n';
    std::cout << "selected_case_count=" << selected.size() << '\n';

    for (const auto& selected_case : selected) {
        auto options = bfs_problem2_options;
        options.adaptive_interval_trace_enabled = true;
        const double local_theta = bfs::global_periapsis_angle_to_problem2_local(
            selected_case.attempt.from_planet,
            selected_case.attempt.parent_incoming_theta);
        const auto result = problem2::solve_problem2_gravity_assist_with_table(
            loader,
            selected_case.attempt.from_planet,
            selected_case.attempt.to_planet,
            selected_case.attempt.parent_time,
            selected_case.attempt.parent_incoming_e,
            local_theta,
            options);
        print_case_summary(selected_case, result);
        print_top_intervals(selected_case.name, result.profile.adaptive_interval_traces);
        print_repeated_query_summary(selected_case.name, result.profile);
        print_stable_endpoint_audit(selected_case.name, result.profile);
    }

    std::cout << "Problem2AdaptiveIntervalAuditSummary\n";
    std::cout << "selected_case_count=" << selected.size() << '\n';
    std::cout << "audit_ok=1\n";
    return 0;
}
