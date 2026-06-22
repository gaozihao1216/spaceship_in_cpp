/*
 * 文件作用：实现固定行星序列上的 BFS 轨迹搜索。
 */
#include "spaceship_cpp/bfs/fixed_sequence_bfs.hpp"

#include "spaceship_cpp/bfs/trajectory_search_config.hpp"
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
#include <cmath>
#include <deque>
#include <limits>
#include <optional>
#include <utility>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::config::global_config;
using spaceship_cpp::config::make_problem1_solve_input;
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

double problem1_residual_to_seconds(double residual_scale_free) {
    const double mu = get_solar_system_physical_params().GM_sun;
    return residual_scale_free / std::sqrt(mu);
}

struct SearchNode {
    TrajectorySearchState state{};
    int parent_index = -1;
    TrajectorySearchEdge edge{};
    int legs_completed = 0;
};

bool validate_config(const FixedSequenceBfsConfig& config, std::string& error_message) {
    if (!is_finite(config.launch_time_seconds_since_j2000) ||
        !is_finite(config.launch_transfer_theta_global) ||
        !is_finite(config.max_total_time_seconds) ||
        !(config.max_total_time_seconds > 0.0)) {
        error_message = "invalid_launch_time_theta_or_max_total_time";
        return false;
    }
    if (is_finite(config.max_launch_v_inf_mps) && !(config.max_launch_v_inf_mps > 0.0)) {
        error_message = "invalid_max_launch_v_inf";
        return false;
    }
    if (config.planet_sequence.size() < 2U) {
        error_message = "planet_sequence_must_have_at_least_two_planets";
        return false;
    }
    if (config.planet_sequence.front() != PlanetId::Earth) {
        error_message = "planet_sequence_must_start_at_earth";
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

Problem1SolveInput make_leg0_problem1_input(
    const FixedSequenceBfsConfig& config,
    PlanetId target_planet
) {
    return make_problem1_solve_input(
        PlanetId::Earth,
        target_planet,
        config.launch_time_seconds_since_j2000,
        config.launch_transfer_theta_global,
        make_problem1_solve_defaults(config.problem1));
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

std::vector<SearchNode> expand_leg0(
    const FixedSequenceBfsConfig& config,
    FixedSequenceBfsStats& stats
) {
    std::vector<SearchNode> children;
    const PlanetId target_planet = config.planet_sequence[1];
    const Problem1SolveInput solve_input = make_leg0_problem1_input(config, target_planet);
    const std::vector<Problem1Candidate> candidates = problem1::solve_problem1(solve_input);
    stats.leg0_candidates = candidates.size();

    for (const Problem1Candidate& candidate : candidates) {
        Leg0CandidateFilterStats filter_stats{};
        if (!leg0_candidate_passes_pruning(
                candidate,
                config.launch_time_seconds_since_j2000,
                Leg0PruningLimits{
                    .max_total_time_seconds = config.max_total_time_seconds,
                    .max_launch_v_inf_mps = config.max_launch_v_inf_mps,
                },
                &filter_stats)) {
            stats.dead_by_time_limit += filter_stats.dead_by_time_limit;
            stats.dead_by_invalid_v_inf += filter_stats.dead_by_invalid_v_inf;
            stats.dead_by_launch_v_inf_limit += filter_stats.dead_by_launch_v_inf_limit;
            continue;
        }

        const auto& residual = candidate.residual_result;
        const double launch_v_inf = trajectory::relative_speed_to_planet(
            PlanetId::Earth,
            config.launch_time_seconds_since_j2000,
            residual.transfer_e,
            residual.transfer_perihelion_angle_used);

        children.push_back(make_node_after_leg0(target_planet, candidate, launch_v_inf));
    }
    return children;
}

std::vector<SearchNode> expand_flyby_leg(
    const FixedSequenceBfsConfig& config,
    const SearchNode& parent,
    PlanetId target_planet,
    FixedSequenceBfsStats& stats
) {
    std::vector<SearchNode> children;
    const PlanetId flyby_planet = parent.state.current_planet;
    const double flyby_time = parent.state.current_time;

    const Problem2FlybySolveConfig solve_config =
        make_flyby_solve_config(flyby_planet, target_planet, flyby_time, config.problem1);

    Problem2FlybySolveInput scan_input{};
    scan_input.flyby_planet = flyby_planet;
    scan_input.target_planet = target_planet;
    scan_input.flyby_time_seconds_since_j2000 = flyby_time;
    const Problem2ThetaPrimeInitialScanResult scan =
        problem2::run_problem2_flyby_theta_prime_initial_scan(scan_input, solve_config);
    if (!scan.ok) {
        ++stats.dead_by_no_p2_solution;
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
    const auto solved = problem2::solve_problem2_flyby_with_scan(solve_input, solve_config, scan);
    if (!solved.ok || solved.solutions.empty()) {
        ++stats.dead_by_no_p2_solution;
        return children;
    }

    for (const Problem2FlybySolution& solution : solved.solutions) {
        if (!is_finite(solution.time_of_flight_seconds) || !(solution.time_of_flight_seconds > 0.0)) {
            continue;
        }
        const double total_time =
            parent.state.accumulated_time_seconds + solution.time_of_flight_seconds;
        if (total_time > config.max_total_time_seconds) {
            ++stats.dead_by_time_limit;
            continue;
        }
        if (!passes_flyby_turn_angle_pruning(
                flyby_planet,
                flyby_time,
                parent.state.incoming_e,
                parent.state.incoming_theta,
                solution,
                config.flyby_physical_filter)) {
            ++stats.dead_by_turn_angle;
            continue;
        }

        children.push_back(make_node_after_flyby_leg(target_planet, parent, solution));
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

}  // namespace

FixedSequenceBfsResult search_fixed_sequence_bfs(const FixedSequenceBfsConfig& config) {
    FixedSequenceBfsResult result{};
    std::string error_message;
    if (!validate_config(config, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }

    FixedSequenceBfsStats stats{};
    std::vector<SearchNode> nodes;
    std::deque<int> queue;

    const std::vector<SearchNode> leg0_children = expand_leg0(config, stats);
    for (SearchNode child : leg0_children) {
        child.parent_index = -1;
        const int index = static_cast<int>(nodes.size());
        nodes.push_back(std::move(child));
        queue.push_back(index);
        ++stats.enqueued_nodes;
    }

    const int final_leg_index = static_cast<int>(config.planet_sequence.size()) - 1;

    while (!queue.empty()) {
        const int node_index = queue.front();
        queue.pop_front();
        ++stats.expanded_nodes;

        const SearchNode& node = nodes[static_cast<std::size_t>(node_index)];
        if (node.legs_completed >= final_leg_index) {
            ++stats.leaf_nodes;
            const double arrival_v_inf = trajectory::compute_arrival_v_inf_at_planet(
                node.state.current_planet,
                node.state.current_time,
                node.state.incoming_e,
                node.state.incoming_theta);
            if (!is_finite(arrival_v_inf) || !(arrival_v_inf > 0.0) ||
                !is_finite(node.state.launch_v_inf) || !(node.state.launch_v_inf > 0.0)) {
                ++stats.dead_by_invalid_v_inf;
                continue;
            }
            if (node.state.accumulated_time_seconds > config.max_total_time_seconds) {
                ++stats.dead_by_time_limit;
                continue;
            }

            const double score = node.state.launch_v_inf + arrival_v_inf;
            if (score < result.best_score) {
                result.found_solution = true;
                result.best_score = score;
                result.best_launch_v_inf = node.state.launch_v_inf;
                result.best_arrival_v_inf = arrival_v_inf;
                result.best_total_time_seconds = node.state.accumulated_time_seconds;
                result.best_edges = reconstruct_path(nodes, node_index);
            }
            continue;
        }

        const PlanetId target_planet =
            config.planet_sequence[static_cast<std::size_t>(node.legs_completed + 1)];
        const std::vector<SearchNode> children =
            expand_flyby_leg(config, node, target_planet, stats);
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
    if (!result.found_solution) {
        result.error_message = "no_feasible_fixed_sequence_path";
    }
    return result;
}

}  // namespace spaceship_cpp::bfs
