/*
 * 文件作用：实现从 leg0 种子出发的自由路径 BFS（Step 2）。
 */
#include "spaceship_cpp/bfs/free_path_bfs.hpp"

#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"
#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include <algorithm>
#include <chrono>
#include <deque>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::config::global_config;
using spaceship_cpp::config::make_problem2_flyby_solve_config;
using spaceship_cpp::planet_params::get_solar_system_physical_params;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem1::Problem1Candidate;
using spaceship_cpp::problem1::Problem1SolveInput;
using spaceship_cpp::problem2::Problem2FlybySolveConfig;
using spaceship_cpp::problem2::Problem2FlybySolveInput;
using spaceship_cpp::problem2::Problem2FlybySolution;
using spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult;
using spaceship_cpp::trajectory::FlybyPhysicalFilterMode;

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double problem1_residual_to_seconds(double residual_scale_free) {
    const double mu = get_solar_system_physical_params().GM_sun;
    return residual_scale_free / std::sqrt(mu);
}

struct SearchNode {
    TrajectorySearchState state{};
    int parent_index = -1;
    TrajectorySearchEdge edge{};
    int legs_completed = 0;
    std::vector<PlanetId> planet_sequence;
};

bool path_contains_planet(const std::vector<PlanetId>& path, PlanetId planet) {
    return std::find(path.begin(), path.end(), planet) != path.end();
}

bool validate_free_path_config(
    const TrajectorySearchGlobalConfig& config,
    const Leg0ThetaSeed& seed,
    std::string& error_message
) {
    if (!is_finite(config.mission.launch_time_seconds_since_j2000)) {
        error_message = "invalid_launch_time";
        return false;
    }
    if (!(config.max_search_legs >= 1)) {
        error_message = "max_search_legs_must_be_at_least_one";
        return false;
    }
    if (!is_finite(config.constraints.max_total_time_seconds) ||
        !(config.constraints.max_total_time_seconds > 0.0)) {
        error_message = "invalid_max_total_time";
        return false;
    }
    if (!path_contains_planet(config.mission.visit_planets, seed.first_leg_target_planet)) {
        error_message = "seed_first_leg_target_not_in_visit_planets";
        return false;
    }
    return true;
}

bool passes_flyby_turn_angle_pruning(
    PlanetId flyby_planet,
    double flyby_time_seconds,
    double incoming_e,
    double incoming_theta_global,
    const Problem2FlybySolution& solution,
    const trajectory::FlybyPhysicalFeasibilityOptions& options
) {
    if (!options.enabled || options.mode == FlybyPhysicalFilterMode::Disabled) {
        return true;
    }

    const auto feasibility = trajectory::evaluate_flyby_physical_feasibility(
        flyby_planet,
        flyby_time_seconds,
        incoming_e,
        incoming_theta_global,
        solution.outgoing_eccentricity,
        solution.outgoing_theta_prime_global,
        options);
    if (!feasibility.valid) {
        return options.mode != FlybyPhysicalFilterMode::Enforce;
    }
    if (options.mode == FlybyPhysicalFilterMode::ObserveOnly) {
        return true;
    }
    return !feasibility.rejected_by_turn_angle;
}

Problem2FlybySolveConfig make_flyby_solve_config(
    PlanetId flyby_planet,
    PlanetId target_planet,
    double flyby_time_seconds,
    const Leg0Problem1Options& problem1
) {
    const auto& defaults = global_config();
    return make_problem2_flyby_solve_config(
        flyby_planet,
        target_planet,
        flyby_time_seconds,
        defaults.problem2_theta_prime_scan,
        defaults.problem2_route_a_newton,
        make_problem1_solve_defaults(problem1));
}

TrajectorySearchEdge make_edge_from_problem1_candidate(
    PlanetId from_planet,
    PlanetId to_planet,
    const Problem1Candidate& candidate
) {
    const auto& residual = candidate.residual_result;
    TrajectorySearchEdge edge{};
    edge.valid = true;
    edge.from_planet = from_planet;
    edge.to_planet = to_planet;
    edge.departure_time = candidate.launch_time_seconds_since_j2000;
    edge.arrival_time = candidate.arrival_time_seconds_since_j2000;
    edge.transfer_time_seconds = candidate.time_of_flight_seconds;
    edge.outgoing_e = residual.transfer_e;
    edge.outgoing_p = residual.transfer_p;
    edge.outgoing_theta = residual.transfer_perihelion_angle_used;
    edge.theta_prime = global_periapsis_angle_to_problem2_local(
        from_planet,
        residual.transfer_perihelion_angle_used);
    edge.alpha = candidate.encounter_global_angle;
    edge.transfer_revolution = candidate.transfer_revolution;
    edge.target_revolution = candidate.target_revolution;
    edge.slingshot_residual = 0.0;
    edge.problem1_residual_seconds =
        problem1_residual_to_seconds(residual.residual);
    return edge;
}

TrajectorySearchEdge make_edge_from_flyby_solution(
    PlanetId from_planet,
    PlanetId to_planet,
    double departure_time,
    const Problem2FlybySolution& solution
) {
    TrajectorySearchEdge edge{};
    edge.valid = true;
    edge.from_planet = from_planet;
    edge.to_planet = to_planet;
    edge.departure_time = departure_time;
    edge.arrival_time = solution.arrival_time_seconds_since_j2000;
    edge.transfer_time_seconds = solution.time_of_flight_seconds;
    edge.outgoing_e = solution.outgoing_eccentricity;
    edge.outgoing_p = solution.outgoing_semi_latus_rectum;
    edge.outgoing_theta = solution.outgoing_theta_prime_global;
    edge.theta_prime = solution.outgoing_theta_prime_local;
    edge.alpha = solution.encounter_global_angle;
    edge.transfer_revolution = solution.transfer_revolution;
    edge.target_revolution = solution.target_revolution;
    edge.slingshot_residual = solution.flyby_constraint_G;
    edge.problem1_residual_seconds = 0.0;
    return edge;
}

SearchNode make_node_after_leg0(
    PlanetId target_planet,
    const Problem1Candidate& candidate,
    double launch_v_inf
) {
    const auto& residual = candidate.residual_result;
    SearchNode node{};
    node.legs_completed = 1;
    node.planet_sequence = {PlanetId::Earth, target_planet};
    node.state.valid = true;
    node.state.current_planet = target_planet;
    node.state.current_time = candidate.arrival_time_seconds_since_j2000;
    node.state.incoming_e = residual.transfer_e;
    node.state.incoming_theta = residual.transfer_perihelion_angle_used;
    node.state.depth = 1;
    node.state.accumulated_time_seconds = candidate.time_of_flight_seconds;
    node.state.launch_v_inf = launch_v_inf;
    node.state.accumulated_score = 0.0;
    node.edge = make_edge_from_problem1_candidate(PlanetId::Earth, target_planet, candidate);
    return node;
}

SearchNode make_node_after_flyby_leg(
    PlanetId target_planet,
    const SearchNode& parent,
    const Problem2FlybySolution& solution
) {
    SearchNode node{};
    node.legs_completed = parent.legs_completed + 1;
    node.planet_sequence = parent.planet_sequence;
    node.planet_sequence.push_back(target_planet);
    node.state.valid = true;
    node.state.current_planet = target_planet;
    node.state.current_time = solution.arrival_time_seconds_since_j2000;
    node.state.incoming_e = solution.outgoing_eccentricity;
    node.state.incoming_theta = solution.outgoing_theta_prime_global;
    node.state.depth = parent.state.depth + 1;
    node.state.accumulated_time_seconds =
        parent.state.accumulated_time_seconds + solution.time_of_flight_seconds;
    node.state.launch_v_inf = parent.state.launch_v_inf;
    node.state.accumulated_score = 0.0;
    node.edge = make_edge_from_flyby_solution(
        parent.state.current_planet,
        target_planet,
        parent.state.current_time,
        solution);
    return node;
}

std::vector<SearchNode> expand_leg0_from_seed(
    const TrajectorySearchGlobalConfig& config,
    const Leg0ThetaSeed& seed,
    FreePathBfsStats& stats
) {
    std::vector<SearchNode> children;
    const Leg0PruningLimits pruning = make_leg0_pruning_limits(config.constraints);
    const Problem1SolveInput solve_input = make_problem1_solve_input(
        PlanetId::Earth,
        seed.first_leg_target_planet,
        config.mission.launch_time_seconds_since_j2000,
        seed.transfer_theta_global,
        make_problem1_solve_defaults(config.problem1));
    const Clock::time_point leg0_start = Clock::now();
    const std::vector<Problem1Candidate> candidates = problem1::solve_problem1(solve_input);
    stats.leg0_solve_ms += elapsed_ms(leg0_start, Clock::now());
    stats.leg0_candidates = candidates.size();

    for (const Problem1Candidate& candidate : candidates) {
        Leg0CandidateFilterStats filter_stats{};
        if (!leg0_candidate_passes_pruning(
                candidate,
                config.mission.launch_time_seconds_since_j2000,
                pruning,
                &filter_stats)) {
            stats.dead_by_time_limit += filter_stats.dead_by_time_limit;
            stats.dead_by_invalid_v_inf += filter_stats.dead_by_invalid_v_inf;
            stats.dead_by_launch_v_inf_limit += filter_stats.dead_by_launch_v_inf_limit;
            continue;
        }

        const auto& residual = candidate.residual_result;
        const double launch_v_inf = trajectory::relative_speed_to_planet(
            PlanetId::Earth,
            config.mission.launch_time_seconds_since_j2000,
            residual.transfer_e,
            residual.transfer_perihelion_angle_used);
        children.push_back(make_node_after_leg0(
            seed.first_leg_target_planet,
            candidate,
            launch_v_inf));
    }

    if (children.empty()) {
        ++stats.dead_by_leg0_no_candidate;
    }
    return children;
}

std::vector<SearchNode> expand_flyby_to_target(
    const TrajectorySearchGlobalConfig& config,
    const SearchNode& parent,
    PlanetId target_planet,
    FreePathBfsStats& stats
) {
    std::vector<SearchNode> children;
    const PlanetId flyby_planet = parent.state.current_planet;
    const double flyby_time = parent.state.current_time;

    Problem2FlybySolveConfig solve_config =
        make_flyby_solve_config(flyby_planet, target_planet, flyby_time, config.problem1);
    solve_config.collect_g_search_profile = config.collect_g_search_profile;

    Problem2FlybySolveInput scan_input{};
    scan_input.flyby_planet = flyby_planet;
    scan_input.target_planet = target_planet;
    scan_input.flyby_time_seconds_since_j2000 = flyby_time;
    const Clock::time_point scan_start = Clock::now();
    const Problem2ThetaPrimeInitialScanResult scan =
        problem2::run_problem2_flyby_theta_prime_initial_scan(scan_input, solve_config);
    const double scan_ms = elapsed_ms(scan_start, Clock::now());
    stats.p2_scan_ms += scan_ms;
    if (!scan.ok) {
        ++stats.dead_by_no_p2_solution;
        record_p2_expansion_branching(stats, 0);
        if (config.collect_p2_expansion_timings) {
            stats.p2_expansion_timings.push_back(P2ExpansionTimingRecord{
                .flyby_planet = flyby_planet,
                .target_planet = target_planet,
                .flyby_time_seconds = flyby_time,
                .incoming_e = parent.state.incoming_e,
                .incoming_theta = parent.state.incoming_theta,
                .scan_ms = scan_ms,
                .solve_ms = 0.0,
                .scan_ok = false,
                .solve_ok = false,
                .solution_count = 0,
            });
        }
        return children;
    }

    Problem2FlybySolveInput solve_input{};
    solve_input.flyby_planet = flyby_planet;
    solve_input.target_planet = target_planet;
    solve_input.flyby_time_seconds_since_j2000 = flyby_time;
    solve_input.incoming_eccentricity = parent.state.incoming_e;
    solve_input.incoming_theta = parent.state.incoming_theta;
    solve_input.incoming_theta_is_local = false;

    ++stats.p2_solve_calls;
    const Clock::time_point p2_start = Clock::now();
    const auto solved = problem2::solve_problem2_flyby_with_scan(solve_input, solve_config, scan);
    const double solve_ms = elapsed_ms(p2_start, Clock::now());
    stats.p2_solve_ms += solve_ms;
    if (config.collect_g_search_profile && solved.has_g_search_profile) {
        problem2::merge_problem2_flyby_G_search_profile(
            stats.g_search_profile,
            solved.g_search_profile);
        ++stats.g_search_profile_samples;
    }
    if (!solved.ok || solved.solutions.empty()) {
        ++stats.dead_by_no_p2_solution;
        record_p2_expansion_branching(stats, 0);
        if (config.collect_p2_expansion_timings) {
            stats.p2_expansion_timings.push_back(P2ExpansionTimingRecord{
                .flyby_planet = flyby_planet,
                .target_planet = target_planet,
                .flyby_time_seconds = flyby_time,
                .incoming_e = parent.state.incoming_e,
                .incoming_theta = parent.state.incoming_theta,
                .scan_ms = scan_ms,
                .solve_ms = solve_ms,
                .scan_ok = true,
                .solve_ok = false,
                .solution_count = 0,
            });
        }
        return children;
    }

    for (const Problem2FlybySolution& solution : solved.solutions) {
        if (!is_finite(solution.time_of_flight_seconds) || !(solution.time_of_flight_seconds > 0.0)) {
            continue;
        }
        const double total_time =
            parent.state.accumulated_time_seconds + solution.time_of_flight_seconds;
        if (total_time > config.constraints.max_total_time_seconds) {
            ++stats.dead_by_time_limit;
            continue;
        }
        if (!passes_flyby_turn_angle_pruning(
                flyby_planet,
                flyby_time,
                parent.state.incoming_e,
                parent.state.incoming_theta,
                solution,
                config.constraints.flyby_physical_filter)) {
            ++stats.dead_by_turn_angle;
            continue;
        }

        children.push_back(make_node_after_flyby_leg(target_planet, parent, solution));
    }
    record_p2_expansion_branching(stats, children.size());
    if (config.collect_p2_expansion_timings) {
        stats.p2_expansion_timings.push_back(P2ExpansionTimingRecord{
            .flyby_planet = flyby_planet,
            .target_planet = target_planet,
            .flyby_time_seconds = flyby_time,
            .incoming_e = parent.state.incoming_e,
            .incoming_theta = parent.state.incoming_theta,
            .scan_ms = scan_ms,
            .solve_ms = solve_ms,
            .scan_ok = true,
            .solve_ok = true,
            .solution_count = children.size(),
        });
    }
    return children;
}

std::vector<SearchNode> expand_free_path_node(
    const TrajectorySearchGlobalConfig& config,
    const SearchNode& parent,
    FreePathBfsStats& stats
) {
    std::vector<SearchNode> children;
    for (const PlanetId target_planet : config.mission.visit_planets) {
        const std::vector<SearchNode> leg_children =
            expand_flyby_to_target(config, parent, target_planet, stats);
        for (SearchNode child : leg_children) {
            children.push_back(std::move(child));
        }
    }
    return children;
}

std::vector<TrajectorySearchEdge> reconstruct_path(
    const std::vector<SearchNode>& nodes,
    int leaf_index
) {
    std::vector<TrajectorySearchEdge> edges;
    int index = leaf_index;
    while (index >= 0) {
        const SearchNode& node = nodes[static_cast<std::size_t>(index)];
        if (node.edge.valid) {
            edges.push_back(node.edge);
        }
        index = node.parent_index;
    }
    std::reverse(edges.begin(), edges.end());
    return edges;
}

std::optional<FreePathBfsSolution> make_solution_from_leaf(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<SearchNode>& nodes,
    int leaf_index,
    FreePathBfsStats& stats
) {
    const SearchNode& node = nodes[static_cast<std::size_t>(leaf_index)];
    const double arrival_v_inf = trajectory::compute_arrival_v_inf_at_planet(
        node.state.current_planet,
        node.state.current_time,
        node.state.incoming_e,
        node.state.incoming_theta);
    if (!is_finite(arrival_v_inf) || !(arrival_v_inf > 0.0) ||
        !is_finite(node.state.launch_v_inf) || !(node.state.launch_v_inf > 0.0)) {
        ++stats.dead_by_invalid_v_inf;
        return std::nullopt;
    }
    if (node.state.accumulated_time_seconds > config.constraints.max_total_time_seconds) {
        ++stats.dead_by_time_limit;
        return std::nullopt;
    }

    FreePathBfsSolution solution{};
    solution.planet_sequence = node.planet_sequence;
    solution.launch_v_inf = node.state.launch_v_inf;
    solution.arrival_v_inf = arrival_v_inf;
    solution.score = node.state.launch_v_inf + arrival_v_inf;
    solution.total_time_seconds = node.state.accumulated_time_seconds;
    solution.edges = reconstruct_path(nodes, leaf_index);
    return solution;
}

}  // namespace

FreePathBfsResult search_free_path_from_leg0_seed(
    const TrajectorySearchGlobalConfig& config,
    const Leg0ThetaSeed& seed
) {
    FreePathBfsResult result{};
    result.seed = seed;

    std::string error_message;
    if (!validate_free_path_config(config, seed, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }

    FreePathBfsStats stats{};
    std::vector<SearchNode> nodes;
    std::deque<int> queue;

    const std::vector<SearchNode> leg0_children = expand_leg0_from_seed(config, seed, stats);
    for (SearchNode child : leg0_children) {
        child.parent_index = -1;
        const int index = static_cast<int>(nodes.size());
        nodes.push_back(std::move(child));
        queue.push_back(index);
        ++stats.enqueued_nodes;
    }

    const int max_search_legs = config.max_search_legs;

    while (!queue.empty()) {
        const int node_index = queue.front();
        queue.pop_front();
        ++stats.expanded_nodes;

        const SearchNode& node = nodes[static_cast<std::size_t>(node_index)];
        if (node.legs_completed >= max_search_legs) {
            ++stats.leaf_nodes;
            const auto solution = make_solution_from_leaf(config, nodes, node_index, stats);
            if (solution.has_value()) {
                result.solutions.push_back(*solution);
                ++stats.recorded_solutions;
            }
            continue;
        }

        const std::vector<SearchNode> children =
            expand_free_path_node(config, node, stats);
        if (children.empty() && node.legs_completed >= 1) {
            ++stats.leaf_nodes;
            const auto solution = make_solution_from_leaf(config, nodes, node_index, stats);
            if (solution.has_value()) {
                result.solutions.push_back(*solution);
                ++stats.recorded_solutions;
            }
            continue;
        }
        for (SearchNode child : children) {
            child.parent_index = node_index;
            const int child_index = static_cast<int>(nodes.size());
            nodes.push_back(std::move(child));
            queue.push_back(child_index);
            ++stats.enqueued_nodes;
        }
    }

    result.ok = true;
    result.stats = stats;
    if (result.solutions.empty()) {
        result.error_message = "no_feasible_free_path_from_seed";
    }
    return result;
}

Step2FreePathSearchResult run_step2_free_path_search(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<Leg0ThetaSeed>& seeds
) {
    Step2FreePathSearchResult result{};
    if (seeds.empty()) {
        result.error_message = "step2_requires_at_least_one_seed";
        return result;
    }

    result.by_seed.reserve(seeds.size());
    result.seeds_processed = seeds;
    for (const Leg0ThetaSeed& seed : seeds) {
        FreePathBfsResult seed_result = search_free_path_from_leg0_seed(config, seed);
        if (!seed_result.ok) {
            result.error_message = seed_result.error_message;
            return result;
        }
        result.by_seed.push_back(std::move(seed_result));
    }

    result.ok = true;
    return result;
}

Step2FreePathSearchResult run_step2_free_path_search(
    const TrajectorySearchGlobalConfig& config
) {
    const Leg0MultiTargetThetaResult leg0 = find_leg0_feasible_theta_for_first_leg_targets(config);
    Step2FreePathSearchResult result{};
    if (!leg0.ok) {
        result.error_message = leg0.error_message;
        return result;
    }

    const std::vector<Leg0ThetaSeed> seeds = discretize_leg0_theta_seeds(config, leg0);
    return run_step2_free_path_search(config, seeds);
}

void record_p2_expansion_branching(FreePathBfsStats& stats, std::size_t child_count) {
    ++stats.p2_expansion_attempts;
    ++stats.p2_expansion_child_count_hist[child_count];
}

void merge_free_path_bfs_stats(FreePathBfsStats& into, const FreePathBfsStats& from) {
    into.leg0_candidates += from.leg0_candidates;
    into.expanded_nodes += from.expanded_nodes;
    into.enqueued_nodes += from.enqueued_nodes;
    into.leaf_nodes += from.leaf_nodes;
    into.recorded_solutions += from.recorded_solutions;
    into.dead_by_time_limit += from.dead_by_time_limit;
    into.dead_by_no_p2_solution += from.dead_by_no_p2_solution;
    into.dead_by_turn_angle += from.dead_by_turn_angle;
    into.dead_by_invalid_v_inf += from.dead_by_invalid_v_inf;
    into.dead_by_launch_v_inf_limit += from.dead_by_launch_v_inf_limit;
    into.dead_by_leg0_no_candidate += from.dead_by_leg0_no_candidate;
    into.p2_solve_calls += from.p2_solve_calls;
    into.p2_expansion_attempts += from.p2_expansion_attempts;
    into.leg0_solve_ms += from.leg0_solve_ms;
    into.p2_scan_ms += from.p2_scan_ms;
    into.p2_solve_ms += from.p2_solve_ms;
    into.g_search_profile_samples += from.g_search_profile_samples;
    if (from.g_search_profile_samples > 0U) {
        problem2::merge_problem2_flyby_G_search_profile(
            into.g_search_profile,
            from.g_search_profile);
    }
    if (!from.p2_expansion_timings.empty()) {
        into.p2_expansion_timings.insert(
            into.p2_expansion_timings.end(),
            from.p2_expansion_timings.begin(),
            from.p2_expansion_timings.end());
    }
    for (const auto& [child_count, frequency] : from.p2_expansion_child_count_hist) {
        into.p2_expansion_child_count_hist[child_count] += frequency;
    }
}

FreePathBfsStats aggregate_free_path_bfs_stats(
    const std::vector<FreePathBfsResult>& by_seed
) {
    FreePathBfsStats aggregate{};
    for (const FreePathBfsResult& seed_result : by_seed) {
        merge_free_path_bfs_stats(aggregate, seed_result.stats);
    }
    return aggregate;
}

std::vector<P2ExpansionTimingRecord> aggregate_p2_expansion_timings(
    const std::vector<FreePathBfsResult>& by_seed
) {
    std::vector<P2ExpansionTimingRecord> timings;
    for (const FreePathBfsResult& seed_result : by_seed) {
        timings.insert(
            timings.end(),
            seed_result.stats.p2_expansion_timings.begin(),
            seed_result.stats.p2_expansion_timings.end());
    }
    return timings;
}

void print_p2_expansion_branching_histogram(
    std::ostream& output,
    const FreePathBfsStats& stats
) {
    output << "P2 expansion attempts=" << stats.p2_expansion_attempts
           << " solve_calls=" << stats.p2_solve_calls << '\n';
    if (stats.p2_expansion_attempts == 0U) {
        output << "  (no P2 expansions)\n";
        return;
    }

    std::size_t max_child_count = 0;
    for (const auto& [child_count, _] : stats.p2_expansion_child_count_hist) {
        max_child_count = std::max(max_child_count, child_count);
    }

    for (std::size_t child_count = 0; child_count <= max_child_count; ++child_count) {
        const auto found = stats.p2_expansion_child_count_hist.find(child_count);
        const std::size_t frequency =
            found == stats.p2_expansion_child_count_hist.end() ? 0U : found->second;
        if (frequency == 0U) {
            continue;
        }
        const double probability =
            static_cast<double>(frequency) /
            static_cast<double>(stats.p2_expansion_attempts);
        output << "  child_count=" << child_count
               << " frequency=" << frequency
               << " probability=" << probability << '\n';
    }
}

void print_step2_timing_breakdown(
    std::ostream& output,
    const FreePathBfsStats& stats,
    double total_wall_ms
) {
    const double oracle_ms = stats.leg0_solve_ms + stats.p2_scan_ms + stats.p2_solve_ms;
    const double bfs_overhead_ms = std::max(0.0, total_wall_ms - oracle_ms);

    output << std::fixed << std::setprecision(1);
    output << "Step 2 timing breakdown (ms):\n";
    output << "  total_wall=" << total_wall_ms << '\n';
    output << "  leg0_solve (P1 per seed)=" << stats.leg0_solve_ms;
    if (total_wall_ms > 0.0) {
        output << " (" << (100.0 * stats.leg0_solve_ms / total_wall_ms) << "%)";
    }
    output << '\n';
    output << "  p2_scan (theta' initial scan)=" << stats.p2_scan_ms;
    if (total_wall_ms > 0.0) {
        output << " (" << (100.0 * stats.p2_scan_ms / total_wall_ms) << "%)";
    }
    output << ", p2_expansion_attempts=" << stats.p2_expansion_attempts << '\n';
    output << "  p2_solve (flyby with scan)=" << stats.p2_solve_ms;
    if (total_wall_ms > 0.0) {
        output << " (" << (100.0 * stats.p2_solve_ms / total_wall_ms) << "%)";
    }
    output << ", calls=" << stats.p2_solve_calls << '\n';
    if (stats.p2_solve_calls > 0U) {
        output << std::setprecision(2);
        output << "  mean_ms_per_p2_solve="
               << (stats.p2_solve_ms / static_cast<double>(stats.p2_solve_calls)) << '\n';
    }
    output << std::setprecision(1);
    output << "  bfs_overhead (rest)=" << bfs_overhead_ms;
    if (total_wall_ms > 0.0) {
        output << " (" << (100.0 * bfs_overhead_ms / total_wall_ms) << "%)";
    }
    output << '\n';
    output << "Step 2 BFS scale:\n";
    output << "  expanded_nodes=" << stats.expanded_nodes
           << ", enqueued_nodes=" << stats.enqueued_nodes
           << ", p2_expansion_attempts=" << stats.p2_expansion_attempts << '\n';
    if (stats.p2_solve_calls > 0U && stats.expanded_nodes > 0U) {
        output << std::setprecision(2);
        output << "  p2_solve_calls_per_expanded_node="
               << (static_cast<double>(stats.p2_solve_calls) /
                   static_cast<double>(stats.expanded_nodes))
               << '\n';
        output << "  p2_expansion_attempts_per_expanded_node="
               << (static_cast<double>(stats.p2_expansion_attempts) /
                   static_cast<double>(stats.expanded_nodes))
               << " (visit_planets per node)\n";
    }
}

void print_step2_g_search_profile_breakdown(
    std::ostream& output,
    const FreePathBfsStats& stats
) {
    if (stats.g_search_profile_samples == 0U) {
        output << "G-search profile: (not collected; set collect_g_search_profile=true)\n";
        return;
    }

    const auto& p = stats.g_search_profile;
    const double n = static_cast<double>(stats.g_search_profile_samples);
    const double mean_solve_ms =
        stats.p2_solve_ms / static_cast<double>(stats.p2_solve_calls);
    const double mean_route_a = p.route_a_ms / n;
    const double mean_case_c = p.case_c_middle_ms / n;
    const double mean_incoming = p.incoming_cache_ms / n;
    const double mean_enrich = p.enrich_ms / n;
    const double mean_framework =
        mean_solve_ms - mean_route_a - mean_case_c - mean_incoming - mean_enrich;

    output << std::fixed << std::setprecision(2);
    output << "G-search profile (aggregated over " << stats.g_search_profile_samples
           << " P2 solves, mean_ms_per_p2_solve=" << mean_solve_ms << "):\n";
    output << "  time_ms/solve: route_a=" << mean_route_a
           << " (" << (100.0 * mean_route_a / std::max(mean_solve_ms, 1e-9)) << "%)"
           << ", case_c_middle=" << mean_case_c
           << " (" << (100.0 * mean_case_c / std::max(mean_solve_ms, 1e-9)) << "%)"
           << ", framework=" << mean_framework
           << " (" << (100.0 * mean_framework / std::max(mean_solve_ms, 1e-9)) << "%)"
           << ", incoming_F=" << mean_incoming
           << ", enrich=" << mean_enrich << '\n';
    output << std::setprecision(1);
    output << "  counts/solve: interval_visits=" << (p.interval_visits / n)
           << " equal=" << (p.interval_equal / n)
           << " case_b=" << (p.interval_case_b / n)
           << " case_c=" << (p.interval_case_c / n)
           << " discarded=" << (p.interval_discarded / n) << '\n';
    output << "  calls/solve: route_a=" << (p.route_a_calls / n)
           << " (iters=" << (p.route_a_iterations / n) << ")"
           << ", case_c_middle=" << (p.case_c_middle_calls / n)
           << ", case_b_probe=" << (p.case_b_probe_calls / n)
           << ", g_newton=" << (p.g_newton_calls / n)
           << " (iters=" << (p.g_newton_iterations / n) << ")\n";
}

}  // namespace spaceship_cpp::bfs
