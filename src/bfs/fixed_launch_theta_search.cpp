#include "spaceship_cpp/bfs/fixed_launch_theta_search.hpp"

#include "spaceship_cpp/bfs/gravity_assist_step_expansion.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>
#include <utility>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double time_score(double time_weight_m_per_s_per_day, double seconds) {
    return time_weight_m_per_s_per_day * (seconds / 86400.0);
}

bool strict_residual_source(const TrajectorySearchEdge& edge) {
    return edge.residual_source == problem2::Problem2ResidualSource::Strict;
}

bool node_less_for_beam(const FixedLaunchThetaSearchNode& lhs, const FixedLaunchThetaSearchNode& rhs) {
    if (lhs.partial_score != rhs.partial_score) {
        return lhs.partial_score < rhs.partial_score;
    }
    if (strict_residual_source(lhs.incoming_edge) != strict_residual_source(rhs.incoming_edge)) {
        return strict_residual_source(lhs.incoming_edge);
    }
    const double lhs_problem1 = std::abs(lhs.incoming_edge.problem1_residual_seconds);
    const double rhs_problem1 = std::abs(rhs.incoming_edge.problem1_residual_seconds);
    if (lhs_problem1 != rhs_problem1) {
        return lhs_problem1 < rhs_problem1;
    }
    return lhs.incoming_edge.arrival_time < rhs.incoming_edge.arrival_time;
}

bool terminal_less(
    const FixedLaunchThetaTerminalSolution& lhs,
    const FixedLaunchThetaTerminalSolution& rhs
) {
    if (lhs.score != rhs.score) {
        return lhs.score < rhs.score;
    }
    if (lhs.total_delta_v != rhs.total_delta_v) {
        return lhs.total_delta_v < rhs.total_delta_v;
    }
    return lhs.total_flight_time_seconds < rhs.total_flight_time_seconds;
}

bool edge_passes_filters(
    const TrajectorySearchEdge& edge,
    const TrajectorySearchState& next_state,
    const FixedLaunchThetaSearchOptions& options
) {
    return edge.valid &&
        next_state.valid &&
        edge.arrival_time > edge.departure_time &&
        edge.transfer_time_seconds > 0.0 &&
        edge.outgoing_p > 0.0 &&
        std::abs(edge.slingshot_residual) <= options.max_abs_slingshot_residual &&
        std::abs(edge.problem1_residual_seconds) <= options.max_abs_problem1_residual_seconds &&
        edge.arrival_time <= options.max_arrival_time &&
        is_finite(edge.arrival_time) &&
        is_finite(next_state.incoming_e) &&
        is_finite(next_state.incoming_theta);
}

std::string invalid_planet_state_message(
    const FixedLaunchThetaSearchNode& node,
    int node_index,
    const char* reason
) {
    std::ostringstream os;
    os << "BfsInvalidPlanetStateDiagnostic"
       << " depth=" << node.state.depth
       << " node_index=" << node_index
       << " current_planet_raw=" << planet_params::planet_id_raw_value(node.state.current_planet)
       << " parent_index=" << node.parent_index
       << " incoming_edge_from_raw=" << planet_params::planet_id_raw_value(node.incoming_edge.from_planet)
       << " incoming_edge_to_raw=" << planet_params::planet_id_raw_value(node.incoming_edge.to_planet)
       << " reason=" << reason;
    return os.str();
}

bool node_has_valid_planets(
    const FixedLaunchThetaSearchNode& node,
    int node_index,
    std::string* error_message
) {
    if (!planet_params::is_valid_planet_id(node.state.current_planet)) {
        if (error_message) {
            *error_message = invalid_planet_state_message(node, node_index, "invalid_current_planet");
        }
        return false;
    }
    if (node.incoming_edge.valid &&
        (!planet_params::is_valid_planet_id(node.incoming_edge.from_planet) ||
         !planet_params::is_valid_planet_id(node.incoming_edge.to_planet))) {
        if (error_message) {
            *error_message = invalid_planet_state_message(node, node_index, "invalid_incoming_edge_planet");
        }
        return false;
    }
    return true;
}

bool edge_passes_flyby_physical_filter(
    const TrajectorySearchState& parent_state,
    const TrajectorySearchEdge& edge,
    const FixedLaunchThetaSearchOptions& options
) {
    if (!options.flyby_physical_options.enabled ||
        options.flyby_physical_options.mode == trajectory::FlybyPhysicalFilterMode::Disabled ||
        options.flyby_physical_options.mode == trajectory::FlybyPhysicalFilterMode::ObserveOnly) {
        return true;
    }
    const auto flyby = trajectory::evaluate_flyby_physical_feasibility(
        parent_state.current_planet,
        parent_state.current_time,
        parent_state.incoming_e,
        parent_state.incoming_theta,
        edge.outgoing_e,
        edge.outgoing_theta,
        options.flyby_physical_options);
    return flyby.valid && flyby.feasible;
}

std::string expansion_exception_message(
    int depth,
    int parent_index,
    const FixedLaunchThetaSearchNode& parent,
    planet_params::PlanetId target_planet,
    const char* phase,
    const std::exception& ex
) {
    std::ostringstream os;
    os << "BfsOpenSearchExceptionDiagnostic"
       << " phase=" << phase
       << " depth=" << depth
       << " parent_index=" << parent_index
       << " parent_current_planet_raw=" << planet_params::planet_id_raw_value(parent.state.current_planet)
       << " target_planet_raw=" << planet_params::planet_id_raw_value(target_planet)
       << " parent_incoming_edge_from_raw=" << planet_params::planet_id_raw_value(parent.incoming_edge.from_planet)
       << " parent_incoming_edge_to_raw=" << planet_params::planet_id_raw_value(parent.incoming_edge.to_planet)
       << " message=" << ex.what();
    return os.str();
}

std::string terminal_exception_message(
    int node_index,
    const FixedLaunchThetaSearchNode& node,
    const char* phase,
    const std::exception& ex
) {
    std::ostringstream os;
    os << "BfsOpenSearchExceptionDiagnostic"
       << " phase=" << phase
       << " depth=" << node.state.depth
       << " node_index=" << node_index
       << " current_planet_raw=" << planet_params::planet_id_raw_value(node.state.current_planet)
       << " incoming_edge_from_raw=" << planet_params::planet_id_raw_value(node.incoming_edge.from_planet)
       << " incoming_edge_to_raw=" << planet_params::planet_id_raw_value(node.incoming_edge.to_planet)
       << " message=" << ex.what();
    return os.str();
}

FixedLaunchThetaTerminalSolution make_terminal_solution(
    const FixedLaunchThetaSearchNode& node,
    int node_index,
    const FixedLaunchThetaSearchOptions& options
) {
    FixedLaunchThetaTerminalSolution terminal{};
    terminal.node_index = node_index;
    terminal.launch_v_inf = node.state.launch_v_inf;
    terminal.arrival_v_inf = compute_arrival_v_inf_from_state(node.state);
    terminal.total_delta_v = terminal.launch_v_inf + terminal.arrival_v_inf;
    terminal.total_flight_time_seconds = node.state.accumulated_time_seconds;
    terminal.score = terminal.total_delta_v +
        time_score(options.time_weight_m_per_s_per_day, terminal.total_flight_time_seconds);
    terminal.depth = node.state.depth;
    terminal.valid = is_finite(terminal.launch_v_inf) &&
        is_finite(terminal.arrival_v_inf) &&
        is_finite(terminal.total_delta_v) &&
        is_finite(terminal.score);
    return terminal;
}

std::vector<int> select_frontier(
    const std::vector<FixedLaunchThetaSearchNode>& nodes,
    std::vector<int> candidate_indices,
    int beam_width
) {
    std::sort(candidate_indices.begin(), candidate_indices.end(), [&](int lhs, int rhs) {
        return node_less_for_beam(
            nodes[static_cast<std::size_t>(lhs)],
            nodes[static_cast<std::size_t>(rhs)]);
    });
    if (static_cast<int>(candidate_indices.size()) > beam_width) {
        candidate_indices.resize(static_cast<std::size_t>(beam_width));
    }
    return candidate_indices;
}

std::vector<int> select_frontier_optional_beam(
    const std::vector<FixedLaunchThetaSearchNode>& nodes,
    std::vector<int> candidate_indices,
    const FixedLaunchThetaSearchOptions& options
) {
    if (!options.enable_beam_pruning) {
        return candidate_indices;
    }
    return select_frontier(nodes, std::move(candidate_indices), options.beam_width);
}

}  // namespace

FixedLaunchThetaSearchResult run_fixed_launch_theta_beam_search_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    planet_params::PlanetId terminal_planet,
    double launch_time,
    double initial_theta,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const std::vector<planet_params::PlanetId>& allowed_transfer_planets,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedLaunchThetaSearchOptions& search_options
) {
    FixedLaunchThetaSearchResult result{};
    result.launch_time = launch_time;
    result.initial_theta = initial_theta;
    const auto search_start = Clock::now();
    auto finalize_success = [&]() {
        if (!result.terminal_solutions.empty()) {
            result.best_terminal_solution = *std::min_element(
                result.terminal_solutions.begin(), result.terminal_solutions.end(), terminal_less);
        }
        result.ok = true;
        return result;
    };
    auto profile_limits_exceeded = [&]() {
        if (search_options.profile_max_node_count > 0 &&
            static_cast<int>(result.nodes.size()) >= search_options.profile_max_node_count) {
            result.stopped_by_node_limit = true;
        }
        if (search_options.profile_max_wall_time_ms > 0.0 &&
            elapsed_ms(search_start, Clock::now()) >= search_options.profile_max_wall_time_ms) {
            result.stopped_by_wall_time_limit = true;
        }
        return result.stopped_by_node_limit || result.stopped_by_wall_time_limit;
    };
    if (!planet_params::is_valid_planet_id(launch_planet) ||
        !planet_params::is_valid_planet_id(terminal_planet)) {
        result.error_message = "invalid_launch_or_terminal_planet";
        return result;
    }
    for (const auto planet : allowed_first_targets) {
        if (!planet_params::is_valid_planet_id(planet)) {
            result.error_message = "invalid_allowed_first_target";
            return result;
        }
    }
    for (const auto planet : allowed_transfer_planets) {
        if (!planet_params::is_valid_planet_id(planet)) {
            result.error_message = "invalid_allowed_transfer_planet";
            return result;
        }
    }
    if (!is_finite(launch_time) || !is_finite(initial_theta)) {
        result.error_message = "invalid_launch_time_or_initial_theta";
        return result;
    }
    if (search_options.max_depth < 1 || search_options.beam_width < 1) {
        result.error_message = "invalid_search_depth_or_beam_width";
        return result;
    }

    FixedLaunchThetaSearchNode root{};
    root.valid = true;
    root.state.valid = true;
    root.state.current_planet = launch_planet;
    root.state.current_time = launch_time;
    root.state.depth = 0;
    root.parent_index = -1;
    result.nodes.push_back(root);

    auto launch_options = initial_launch_options;
    launch_options.max_launch_v_inf = std::min(launch_options.max_launch_v_inf, search_options.max_launch_v_inf);
    launch_options.time_weight_m_per_s_per_day = search_options.time_weight_m_per_s_per_day;

    FixedLaunchThetaSearchLayerSummary initial_summary{};
    initial_summary.depth = 1;
    initial_summary.input_state_count = 1;
    initial_summary.attempted_expansion_count = static_cast<int>(allowed_first_targets.size());
    const auto initial_start = Clock::now();
    InitialLaunchExpansionResult initial{};
    try {
        initial = expand_initial_launch_with_problem1_table(
            loader, launch_planet, launch_time, initial_theta, allowed_first_targets, launch_options);
    } catch (const std::exception& ex) {
        std::ostringstream os;
        os << "BfsOpenSearchExceptionDiagnostic"
           << " phase=initial_launch_expansion"
           << " depth=1"
           << " launch_planet_raw=" << planet_params::planet_id_raw_value(launch_planet)
           << " terminal_planet_raw=" << planet_params::planet_id_raw_value(terminal_planet)
           << " message=" << ex.what();
        result.error_message = os.str();
        return result;
    }
    if (!initial.ok) {
        result.error_message = initial.error_message.empty() ? "initial_launch_expansion_failed" : initial.error_message;
        return result;
    }
    result.initial_candidate_count = initial.accepted_candidate_count;
    result.launch_v_inf_pruned_count = initial.launch_v_inf_pruned_count;
    result.initial_target_error_messages = initial.target_error_messages;
    initial_summary.generated_edge_count = initial.raw_branch_count;

    int initial_terminal_solution_count_at_depth = 0;
    std::vector<int> initial_frontier_candidates;
    for (const auto& candidate : initial.candidates) {
        if (!candidate.valid || !candidate.next_state.valid) {
            continue;
        }
        FixedLaunchThetaSearchNode node{};
        node.valid = true;
        node.state = candidate.next_state;
        node.state.accumulated_score = candidate.next_state.launch_v_inf +
            time_score(search_options.time_weight_m_per_s_per_day, candidate.next_state.accumulated_time_seconds);
        node.parent_index = 0;
        node.incoming_edge = candidate.edge;
        node.partial_score = node.state.accumulated_score;
        const int node_index = static_cast<int>(result.nodes.size());
        std::string invalid_node_message;
        if (!node_has_valid_planets(node, node_index, &invalid_node_message)) {
            result.error_message = invalid_node_message;
            return result;
        }
        result.nodes.push_back(node);
        initial_summary.accepted_edge_count += 1;

        const bool is_terminal = node.state.current_planet == terminal_planet;
        if (is_terminal) {
            FixedLaunchThetaTerminalSolution terminal{};
            try {
                terminal = make_terminal_solution(node, node_index, search_options);
            } catch (const std::exception& ex) {
                result.error_message = terminal_exception_message(
                    node_index, node, "initial_terminal_solution", ex);
                return result;
            }
            if (terminal.valid) {
                result.terminal_solutions.push_back(terminal);
                initial_terminal_solution_count_at_depth += 1;
            }
        }
        if (!is_terminal || search_options.continue_after_reaching_terminal) {
            initial_frontier_candidates.push_back(node_index);
        }
    }
    const int initial_output_before_beam = static_cast<int>(initial_frontier_candidates.size());
    result.frontier_node_indices = select_frontier_optional_beam(
        result.nodes, std::move(initial_frontier_candidates), search_options);
    initial_summary.output_state_count = static_cast<int>(result.frontier_node_indices.size());
    initial_summary.terminal_solution_count_after_layer = static_cast<int>(result.terminal_solutions.size());
    initial_summary.layer_wall_time_ms = elapsed_ms(initial_start, Clock::now());
    result.layer_summaries.push_back(initial_summary);
    OpenSearchDepthWidthStats initial_width{};
    initial_width.depth = 1;
    initial_width.input_state_count = 1;
    initial_width.attempted_expansion_count = static_cast<int>(allowed_first_targets.size());
    initial_width.raw_edge_count_total = initial.raw_branch_count;
    initial_width.accepted_edge_count_total = initial_summary.accepted_edge_count;
    initial_width.terminal_solution_count_at_depth = initial_terminal_solution_count_at_depth;
    initial_width.output_state_count_before_beam = initial_output_before_beam;
    initial_width.output_state_count_after_beam = initial_summary.output_state_count;
    initial_width.beam_pruned_count =
        std::max(0, initial_width.output_state_count_before_beam - initial_width.output_state_count_after_beam);
    initial_width.layer_wall_time_ms = initial_summary.layer_wall_time_ms;
    result.depth_width_stats.push_back(initial_width);
    if (profile_limits_exceeded()) {
        return finalize_success();
    }

    for (int depth = 2; depth <= search_options.max_depth; ++depth) {
        if (result.frontier_node_indices.empty()) {
            break;
        }
        FixedLaunchThetaSearchLayerSummary summary{};
        summary.depth = depth;
        summary.input_state_count = static_cast<int>(result.frontier_node_indices.size());
        const auto layer_start = Clock::now();
        OpenSearchDepthWidthStats width{};
        width.depth = depth;
        width.input_state_count = summary.input_state_count;

        std::vector<int> next_frontier_candidates;
        for (const int parent_index : result.frontier_node_indices) {
            const auto parent = result.nodes[static_cast<std::size_t>(parent_index)];
            std::string invalid_parent_message;
            if (!node_has_valid_planets(parent, parent_index, &invalid_parent_message)) {
                result.error_message = invalid_parent_message;
                return result;
            }
            for (const auto target_planet : allowed_transfer_planets) {
                summary.attempted_expansion_count += 1;
                width.attempted_expansion_count += 1;
                TrajectorySearchExpansionResult expanded{};
                BfsExpansionAttemptProfile attempt_profile{};
                attempt_profile.depth = depth;
                attempt_profile.parent_index = parent_index;
                attempt_profile.from_planet = parent.state.current_planet;
                attempt_profile.to_planet = target_planet;
                attempt_profile.parent_time = parent.state.current_time;
                attempt_profile.parent_incoming_e = parent.state.incoming_e;
                attempt_profile.parent_incoming_theta = parent.state.incoming_theta;
                const auto expansion_start = Clock::now();
                try {
                    expanded = expand_trajectory_state_by_gravity_assist(
                        loader, parent.state, target_planet, problem2_options);
                } catch (const std::exception& ex) {
                    result.error_message = expansion_exception_message(
                        depth, parent_index, parent, target_planet, "problem2_expansion", ex);
                    return result;
                }
                attempt_profile.outer_expansion_wall_ms = elapsed_ms(expansion_start, Clock::now());
                attempt_profile.problem2_solver_profile = expanded.problem2_solver_profile;
                attempt_profile.raw_edge_count = static_cast<int>(expanded.edges.size());
                if (!expanded.ok) {
                    attempt_profile.zero_raw = true;
                    attempt_profile.zero_accepted = true;
                    result.expansion_attempt_profiles.push_back(attempt_profile);
                    continue;
                }
                if (expanded.edges.empty()) {
                    width.zero_raw_edge_solve_count += 1;
                }
                width.raw_edge_count_total += static_cast<int>(expanded.edges.size());
                const std::size_t pair_count = std::min(expanded.edges.size(), expanded.next_states.size());
                summary.generated_edge_count += static_cast<int>(pair_count);
                const int accepted_edge_count_before_solve = width.accepted_edge_count_total;
                const int flyby_reject_count_before_solve = width.reject_flyby_physical_infeasible_count;
                for (std::size_t i = 0; i < pair_count; ++i) {
                    const auto& edge = expanded.edges[i];
                    auto next_state = expanded.next_states[i];
                    if (!edge_passes_filters(edge, next_state, search_options)) {
                        continue;
                    }
                    bool flyby_passes = false;
                    try {
                        flyby_passes = edge_passes_flyby_physical_filter(parent.state, edge, search_options);
                    } catch (const std::exception& ex) {
                        result.error_message = expansion_exception_message(
                            depth, parent_index, parent, target_planet, "flyby_physical_filter", ex);
                        return result;
                    }
                    if (!flyby_passes) {
                        width.reject_flyby_physical_infeasible_count += 1;
                        continue;
                    }
                    next_state.depth = parent.state.depth + 1;
                    next_state.launch_v_inf = parent.state.launch_v_inf;
                    next_state.accumulated_time_seconds =
                        parent.state.accumulated_time_seconds + edge.transfer_time_seconds;
                    next_state.accumulated_score = parent.state.launch_v_inf +
                        time_score(search_options.time_weight_m_per_s_per_day, next_state.accumulated_time_seconds);

                    FixedLaunchThetaSearchNode node{};
                    node.valid = true;
                    node.state = next_state;
                    node.parent_index = parent_index;
                    node.incoming_edge = edge;
                    node.partial_score = next_state.accumulated_score;
                    const int node_index = static_cast<int>(result.nodes.size());
                    std::string invalid_node_message;
                    if (!node_has_valid_planets(node, node_index, &invalid_node_message)) {
                        result.error_message = invalid_node_message;
                        return result;
                    }
                    result.nodes.push_back(node);
                    summary.accepted_edge_count += 1;
                    width.accepted_edge_count_total += 1;

                    const bool is_terminal = next_state.current_planet == terminal_planet;
                    if (is_terminal) {
                        FixedLaunchThetaTerminalSolution terminal{};
                        try {
                            terminal = make_terminal_solution(node, node_index, search_options);
                        } catch (const std::exception& ex) {
                            result.error_message = terminal_exception_message(
                                node_index, node, "problem2_terminal_solution", ex);
                            return result;
                        }
                        if (terminal.valid) {
                            result.terminal_solutions.push_back(terminal);
                            width.terminal_solution_count_at_depth += 1;
                        }
                    }
                    if (!is_terminal || search_options.continue_after_reaching_terminal) {
                        next_frontier_candidates.push_back(node_index);
                    }
                }
                if (width.accepted_edge_count_total == accepted_edge_count_before_solve) {
                    width.zero_accepted_edge_solve_count += 1;
                }
                attempt_profile.accepted_edge_count =
                    width.accepted_edge_count_total - accepted_edge_count_before_solve;
                attempt_profile.flyby_reject_count =
                    width.reject_flyby_physical_infeasible_count - flyby_reject_count_before_solve;
                attempt_profile.zero_raw = attempt_profile.raw_edge_count == 0;
                attempt_profile.zero_accepted = attempt_profile.accepted_edge_count == 0;
                result.expansion_attempt_profiles.push_back(attempt_profile);
            }
        }

        width.output_state_count_before_beam = static_cast<int>(next_frontier_candidates.size());
        result.frontier_node_indices = select_frontier_optional_beam(
            result.nodes, std::move(next_frontier_candidates), search_options);
        summary.output_state_count = static_cast<int>(result.frontier_node_indices.size());
        summary.terminal_solution_count_after_layer = static_cast<int>(result.terminal_solutions.size());
        summary.layer_wall_time_ms = elapsed_ms(layer_start, Clock::now());
        result.layer_summaries.push_back(summary);
        width.output_state_count_after_beam = summary.output_state_count;
        width.beam_pruned_count =
            std::max(0, width.output_state_count_before_beam - width.output_state_count_after_beam);
        width.layer_wall_time_ms = summary.layer_wall_time_ms;
        result.depth_width_stats.push_back(width);
        if (profile_limits_exceeded()) {
            return finalize_success();
        }
    }

    return finalize_success();
}

std::vector<TrajectorySearchEdge> reconstruct_fixed_launch_theta_path(
    const FixedLaunchThetaSearchResult& result,
    int node_index
) {
    std::vector<TrajectorySearchEdge> reversed;
    while (node_index > 0 && node_index < static_cast<int>(result.nodes.size())) {
        const auto& node = result.nodes[static_cast<std::size_t>(node_index)];
        if (!node.valid || !node.incoming_edge.valid) {
            break;
        }
        reversed.push_back(node.incoming_edge);
        node_index = node.parent_index;
    }
    std::reverse(reversed.begin(), reversed.end());
    return reversed;
}

}  // namespace spaceship_cpp::bfs
