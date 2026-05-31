#include "spaceship_cpp/bfs/fixed_sequence_theta_search.hpp"

#include "spaceship_cpp/bfs/gravity_assist_step_expansion.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/trajectory/flyby_physics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>

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

FixedLaunchThetaSearchOptions as_fixed_launch_options(
    const FixedSequenceThetaSearchOptions& options,
    int max_depth
) {
    FixedLaunchThetaSearchOptions converted{};
    converted.max_depth = max_depth;
    converted.beam_width = options.beam_width;
    converted.max_launch_v_inf = options.max_launch_v_inf;
    converted.time_weight_m_per_s_per_day = options.time_weight_m_per_s_per_day;
    converted.max_abs_slingshot_residual = options.max_abs_slingshot_residual;
    converted.max_abs_problem1_residual_seconds = options.max_abs_problem1_residual_seconds;
    converted.max_arrival_time = options.max_arrival_time;
    converted.flyby_physical_options = options.flyby_physical_options;
    converted.continue_after_reaching_terminal = false;
    return converted;
}

Problem2BranchRejectReason classify_problem2_edge_filter_failure(
    const TrajectorySearchEdge& edge,
    const TrajectorySearchState& next_state,
    const FixedLaunchThetaSearchOptions& options
) {
    if (!edge.valid) {
        return Problem2BranchRejectReason::InvalidEdge;
    }
    if (!next_state.valid) {
        return Problem2BranchRejectReason::InvalidNextState;
    }
    if (!(edge.arrival_time > edge.departure_time)) {
        return Problem2BranchRejectReason::NonIncreasingTime;
    }
    if (!(edge.transfer_time_seconds > 0.0)) {
        return Problem2BranchRejectReason::NonPositiveTransferTime;
    }
    if (!(edge.outgoing_p > 0.0)) {
        return Problem2BranchRejectReason::NonPositiveOutgoingP;
    }
    if (!(std::abs(edge.slingshot_residual) <= options.max_abs_slingshot_residual)) {
        return Problem2BranchRejectReason::SlingshotResidualTooLarge;
    }
    if (!(std::abs(edge.problem1_residual_seconds) <= options.max_abs_problem1_residual_seconds)) {
        return Problem2BranchRejectReason::Problem1ResidualTooLarge;
    }
    if (!(edge.arrival_time <= options.max_arrival_time)) {
        return Problem2BranchRejectReason::ArrivalTimeTooLarge;
    }
    if (!is_finite(edge.arrival_time)) {
        return Problem2BranchRejectReason::NonFiniteArrivalTime;
    }
    if (!is_finite(next_state.incoming_e) || !is_finite(next_state.incoming_theta)) {
        return Problem2BranchRejectReason::NonFiniteIncomingOrbit;
    }
    return Problem2BranchRejectReason::None;
}

bool edge_passes_filters(
    const TrajectorySearchEdge& edge,
    const TrajectorySearchState& next_state,
    const FixedLaunchThetaSearchOptions& options
) {
    return classify_problem2_edge_filter_failure(edge, next_state, options) == Problem2BranchRejectReason::None;
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

void accumulate_reject_reason(
    Problem2SolveBranchStats* stats,
    Problem2BranchRejectReason reason
) {
    switch (reason) {
    case Problem2BranchRejectReason::None:
        break;
    case Problem2BranchRejectReason::InvalidEdge:
        stats->reject_invalid_edge_count += 1;
        break;
    case Problem2BranchRejectReason::InvalidNextState:
        stats->reject_invalid_next_state_count += 1;
        break;
    case Problem2BranchRejectReason::NonFiniteArrivalTime:
        stats->reject_nonfinite_arrival_time_count += 1;
        break;
    case Problem2BranchRejectReason::NonFiniteIncomingOrbit:
        stats->reject_nonfinite_incoming_orbit_count += 1;
        break;
    case Problem2BranchRejectReason::NonIncreasingTime:
        stats->reject_non_increasing_time_count += 1;
        break;
    case Problem2BranchRejectReason::NonPositiveTransferTime:
        stats->reject_non_positive_transfer_time_count += 1;
        break;
    case Problem2BranchRejectReason::NonPositiveOutgoingP:
        stats->reject_non_positive_outgoing_p_count += 1;
        break;
    case Problem2BranchRejectReason::SlingshotResidualTooLarge:
        stats->reject_slingshot_residual_count += 1;
        break;
    case Problem2BranchRejectReason::Problem1ResidualTooLarge:
        stats->reject_problem1_residual_count += 1;
        break;
    case Problem2BranchRejectReason::ArrivalTimeTooLarge:
        stats->reject_arrival_time_count += 1;
        break;
    case Problem2BranchRejectReason::FlybyEvaluatorInvalid:
        stats->reject_flyby_evaluator_invalid_count += 1;
        break;
    case Problem2BranchRejectReason::FlybyPhysicalInfeasible:
        stats->reject_flyby_physical_infeasible_count += 1;
        break;
    }
}

void update_solve_ranges(
    Problem2SolveBranchStats* stats,
    const TrajectorySearchEdge& edge
) {
    if (is_finite(edge.arrival_time)) {
        stats->min_arrival_time = std::min(stats->min_arrival_time, edge.arrival_time);
        stats->max_arrival_time = std::max(stats->max_arrival_time, edge.arrival_time);
    }
    if (is_finite(edge.transfer_time_seconds)) {
        stats->min_transfer_time_seconds =
            std::min(stats->min_transfer_time_seconds, edge.transfer_time_seconds);
        stats->max_transfer_time_seconds =
            std::max(stats->max_transfer_time_seconds, edge.transfer_time_seconds);
    }
    const double abs_slingshot = std::abs(edge.slingshot_residual);
    if (is_finite(abs_slingshot)) {
        stats->min_abs_slingshot_residual =
            std::min(stats->min_abs_slingshot_residual, abs_slingshot);
        stats->max_abs_slingshot_residual =
            std::max(stats->max_abs_slingshot_residual, abs_slingshot);
    }
    const double abs_problem1 = std::abs(edge.problem1_residual_seconds);
    if (is_finite(abs_problem1)) {
        stats->min_abs_problem1_residual_seconds =
            std::min(stats->min_abs_problem1_residual_seconds, abs_problem1);
        stats->max_abs_problem1_residual_seconds =
            std::max(stats->max_abs_problem1_residual_seconds, abs_problem1);
    }
}

void update_flyby_ranges(
    Problem2SolveBranchStats* stats,
    const trajectory::FlybyPhysicalFeasibilityResult& flyby,
    const trajectory::FlybyPhysicalFeasibilityOptions& options
) {
    stats->flyby_evaluated_count += 1;
    if (!flyby.valid) {
        stats->flyby_invalid_count += 1;
        stats->would_reject_flyby_evaluator_invalid_count += 1;
        return;
    }
    stats->flyby_valid_count += 1;
    if (flyby.feasible) {
        stats->flyby_feasible_count += 1;
    } else {
        stats->flyby_infeasible_count += 1;
        stats->would_reject_flyby_physical_infeasible_count += 1;
    }
    if (flyby.rejected_by_vinf_mismatch) {
        stats->rejected_by_vinf_mismatch_count += 1;
    }
    if (flyby.rejected_by_turn_angle) {
        stats->rejected_by_turn_angle_count += 1;
    }
    if (flyby.rejected_by_periapsis_radius) {
        stats->rejected_by_periapsis_radius_count += 1;
    }
    const double v_inf = 0.5 * (flyby.v_infinity_in + flyby.v_infinity_out);
    const double mismatch_tolerance = options.relative_speed_tolerance * std::max(1.0, v_inf);
    if (is_finite(flyby.v_infinity_in)) {
        stats->min_v_inf_in = std::min(stats->min_v_inf_in, flyby.v_infinity_in);
        stats->max_v_inf_in = std::max(stats->max_v_inf_in, flyby.v_infinity_in);
    }
    if (is_finite(flyby.v_infinity_out)) {
        stats->min_v_inf_out = std::min(stats->min_v_inf_out, flyby.v_infinity_out);
        stats->max_v_inf_out = std::max(stats->max_v_inf_out, flyby.v_infinity_out);
    }
    if (is_finite(flyby.v_infinity_mismatch)) {
        stats->min_v_inf_mismatch = std::min(stats->min_v_inf_mismatch, flyby.v_infinity_mismatch);
        stats->max_v_inf_mismatch = std::max(stats->max_v_inf_mismatch, flyby.v_infinity_mismatch);
    }
    if (is_finite(mismatch_tolerance)) {
        stats->min_v_inf_mismatch_tolerance =
            std::min(stats->min_v_inf_mismatch_tolerance, mismatch_tolerance);
        stats->max_v_inf_mismatch_tolerance =
            std::max(stats->max_v_inf_mismatch_tolerance, mismatch_tolerance);
    }
    if (is_finite(flyby.required_periapsis_radius_m)) {
        stats->min_flyby_required_periapsis_radius_m =
            std::min(stats->min_flyby_required_periapsis_radius_m, flyby.required_periapsis_radius_m);
        stats->max_flyby_required_periapsis_radius_m =
            std::max(stats->max_flyby_required_periapsis_radius_m, flyby.required_periapsis_radius_m);
    }
    if (is_finite(flyby.turn_angle_rad)) {
        stats->min_flyby_turn_angle_rad =
            std::min(stats->min_flyby_turn_angle_rad, flyby.turn_angle_rad);
        stats->max_flyby_turn_angle_rad =
            std::max(stats->max_flyby_turn_angle_rad, flyby.turn_angle_rad);
    }
    if (is_finite(flyby.max_turn_angle_rad)) {
        stats->min_flyby_max_turn_angle_rad =
            std::min(stats->min_flyby_max_turn_angle_rad, flyby.max_turn_angle_rad);
        stats->max_flyby_max_turn_angle_rad =
            std::max(stats->max_flyby_max_turn_angle_rad, flyby.max_turn_angle_rad);
    }
    if (is_finite(flyby.min_allowed_periapsis_radius_m)) {
        stats->min_flyby_allowed_periapsis_radius_m =
            std::min(stats->min_flyby_allowed_periapsis_radius_m, flyby.min_allowed_periapsis_radius_m);
        stats->max_flyby_allowed_periapsis_radius_m =
            std::max(stats->max_flyby_allowed_periapsis_radius_m, flyby.min_allowed_periapsis_radius_m);
    }
}

Problem2LayerBranchStats summarize_layer_branch_stats(
    const std::vector<Problem2SolveBranchStats>& solve_stats,
    std::size_t begin_index,
    std::size_t end_index,
    int depth,
    planet_params::PlanetId from_planet,
    planet_params::PlanetId to_planet,
    int input_state_count,
    int output_state_count_after_beam,
    double layer_wall_time_ms
) {
    Problem2LayerBranchStats layer{};
    layer.depth = depth;
    layer.from_planet = from_planet;
    layer.to_planet = to_planet;
    layer.input_state_count = input_state_count;
    layer.problem2_solve_count = static_cast<int>(end_index - begin_index);
    layer.output_state_count_after_beam = output_state_count_after_beam;
    layer.layer_wall_time_ms = layer_wall_time_ms;
    if (begin_index == end_index) {
        return layer;
    }

    layer.raw_edge_count_min = std::numeric_limits<int>::max();
    layer.accepted_edge_count_min = std::numeric_limits<int>::max();
    for (std::size_t i = begin_index; i < end_index; ++i) {
        const auto& stats = solve_stats[i];
        layer.raw_edge_count_total += stats.raw_edge_count;
        layer.accepted_edge_count_total += stats.accepted_edge_count;
        layer.raw_edge_count_min = std::min(layer.raw_edge_count_min, stats.raw_edge_count);
        layer.raw_edge_count_max = std::max(layer.raw_edge_count_max, stats.raw_edge_count);
        layer.accepted_edge_count_min = std::min(layer.accepted_edge_count_min, stats.accepted_edge_count);
        layer.accepted_edge_count_max = std::max(layer.accepted_edge_count_max, stats.accepted_edge_count);
    }
    layer.raw_edge_count_mean =
        static_cast<double>(layer.raw_edge_count_total) / static_cast<double>(layer.problem2_solve_count);
    layer.accepted_edge_count_mean =
        static_cast<double>(layer.accepted_edge_count_total) / static_cast<double>(layer.problem2_solve_count);
    layer.beam_pruned_count =
        std::max(0, layer.accepted_edge_count_total - layer.output_state_count_after_beam);
    return layer;
}

}  // namespace

FixedLaunchThetaSearchResult run_fixed_sequence_theta_search_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    double initial_theta,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceThetaSearchOptions& sequence_options
) {
    FixedLaunchThetaSearchResult result{};
    result.launch_time = launch_time;
    result.initial_theta = initial_theta;
    if (!is_finite(launch_time) || !is_finite(initial_theta)) {
        result.error_message = "invalid_launch_time_or_initial_theta";
        return result;
    }
    if (sequence.size() < 2) {
        result.error_message = "fixed_sequence_too_short";
        return result;
    }
    if (sequence_options.beam_width < 1) {
        result.error_message = "invalid_fixed_sequence_beam_width";
        return result;
    }

    const auto search_options = as_fixed_launch_options(
        sequence_options,
        static_cast<int>(sequence.size()) - 1);
    const auto launch_planet = sequence.front();
    const auto terminal_planet = sequence.back();

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
    initial_summary.attempted_expansion_count = 1;
    const auto initial_start = Clock::now();
    const std::vector<planet_params::PlanetId> allowed_first_targets{sequence[1]};
    const auto initial = expand_initial_launch_with_problem1_table(
        loader, launch_planet, launch_time, initial_theta, allowed_first_targets, launch_options);
    if (!initial.ok) {
        result.error_message = initial.error_message.empty() ? "initial_launch_expansion_failed" : initial.error_message;
        return result;
    }
    result.initial_candidate_count = initial.accepted_candidate_count;
    result.launch_v_inf_pruned_count = initial.launch_v_inf_pruned_count;
    initial_summary.generated_edge_count = initial.raw_branch_count;

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
        result.nodes.push_back(node);
        initial_summary.accepted_edge_count += 1;

        const bool is_final_depth = sequence.size() == 2;
        const bool is_terminal = is_final_depth && node.state.current_planet == terminal_planet;
        if (is_terminal) {
            const auto terminal = make_terminal_solution(node, node_index, search_options);
            if (terminal.valid) {
                result.terminal_solutions.push_back(terminal);
            }
        } else if (!is_final_depth) {
            initial_frontier_candidates.push_back(node_index);
        }
    }
    result.frontier_node_indices = select_frontier(result.nodes, initial_frontier_candidates, search_options.beam_width);
    initial_summary.output_state_count = static_cast<int>(result.frontier_node_indices.size());
    initial_summary.terminal_solution_count_after_layer = static_cast<int>(result.terminal_solutions.size());
    initial_summary.layer_wall_time_ms = elapsed_ms(initial_start, Clock::now());
    result.layer_summaries.push_back(initial_summary);

    for (std::size_t depth = 2; depth < sequence.size(); ++depth) {
        if (result.frontier_node_indices.empty()) {
            break;
        }
        FixedLaunchThetaSearchLayerSummary summary{};
        summary.depth = static_cast<int>(depth);
        summary.input_state_count = static_cast<int>(result.frontier_node_indices.size());
        const auto layer_start = Clock::now();
        const std::size_t layer_stats_begin = result.problem2_solve_branch_stats.size();

        std::vector<int> next_frontier_candidates;
        const auto target_planet = sequence[depth];
        const bool is_final_depth = depth + 1 == sequence.size();
        planet_params::PlanetId layer_from_planet = result.nodes[
            static_cast<std::size_t>(result.frontier_node_indices.front())].state.current_planet;
        for (const int parent_index : result.frontier_node_indices) {
            const auto parent = result.nodes[static_cast<std::size_t>(parent_index)];
            summary.attempted_expansion_count += 1;
            Problem2SolveBranchStats solve_stats{};
            solve_stats.valid = true;
            solve_stats.depth = static_cast<int>(depth);
            solve_stats.parent_node_index = parent_index;
            solve_stats.from_planet = parent.state.current_planet;
            solve_stats.to_planet = target_planet;
            solve_stats.parent_time = parent.state.current_time;
            solve_stats.parent_incoming_e = parent.state.incoming_e;
            solve_stats.parent_incoming_theta = parent.state.incoming_theta;

            const auto solve_start = Clock::now();
            const auto expanded = expand_trajectory_state_by_gravity_assist(
                loader, parent.state, target_planet, problem2_options);
            solve_stats.solve_wall_time_ms = elapsed_ms(solve_start, Clock::now());
            if (!expanded.ok) {
                result.problem2_solve_branch_stats.push_back(solve_stats);
                continue;
            }
            solve_stats.raw_edge_count = static_cast<int>(expanded.edges.size());
            const std::size_t pair_count = std::min(expanded.edges.size(), expanded.next_states.size());
            if (expanded.edges.size() > expanded.next_states.size()) {
                solve_stats.reject_invalid_next_state_count +=
                    static_cast<int>(expanded.edges.size() - expanded.next_states.size());
            }
            summary.generated_edge_count += static_cast<int>(pair_count);
            for (std::size_t i = 0; i < pair_count; ++i) {
                const auto& edge = expanded.edges[i];
                auto next_state = expanded.next_states[i];
                update_solve_ranges(&solve_stats, edge);
                const auto reject_reason =
                    classify_problem2_edge_filter_failure(edge, next_state, search_options);
                if (reject_reason != Problem2BranchRejectReason::None) {
                    accumulate_reject_reason(&solve_stats, reject_reason);
                    continue;
                }
                const bool evaluate_flyby = search_options.flyby_physical_options.enabled &&
                    search_options.flyby_physical_options.mode !=
                        trajectory::FlybyPhysicalFilterMode::Disabled;
                if (evaluate_flyby) {
                    const auto flyby = trajectory::evaluate_flyby_physical_feasibility(
                        parent.state.current_planet,
                        parent.state.current_time,
                        parent.state.incoming_e,
                        parent.state.incoming_theta,
                        edge.outgoing_e,
                        edge.outgoing_theta,
                        search_options.flyby_physical_options);
                    update_flyby_ranges(&solve_stats, flyby, search_options.flyby_physical_options);
                    if (search_options.flyby_physical_options.mode ==
                        trajectory::FlybyPhysicalFilterMode::Enforce) {
                        if (!flyby.valid) {
                            accumulate_reject_reason(
                                &solve_stats,
                                Problem2BranchRejectReason::FlybyEvaluatorInvalid);
                            continue;
                        }
                        if (!flyby.feasible) {
                            accumulate_reject_reason(
                                &solve_stats,
                                Problem2BranchRejectReason::FlybyPhysicalInfeasible);
                            continue;
                        }
                    }
                }
                solve_stats.accepted_edge_count += 1;
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
                result.nodes.push_back(node);
                summary.accepted_edge_count += 1;

                const bool is_terminal = is_final_depth && next_state.current_planet == terminal_planet;
                if (is_terminal) {
                    const auto terminal = make_terminal_solution(node, node_index, search_options);
                    if (terminal.valid) {
                        result.terminal_solutions.push_back(terminal);
                    }
                } else if (!is_final_depth) {
                    next_frontier_candidates.push_back(node_index);
                }
            }
            result.problem2_solve_branch_stats.push_back(solve_stats);
        }

        result.frontier_node_indices = select_frontier(result.nodes, next_frontier_candidates, search_options.beam_width);
        summary.output_state_count = static_cast<int>(result.frontier_node_indices.size());
        summary.terminal_solution_count_after_layer = static_cast<int>(result.terminal_solutions.size());
        summary.layer_wall_time_ms = elapsed_ms(layer_start, Clock::now());
        result.layer_summaries.push_back(summary);
        result.problem2_layer_branch_stats.push_back(summarize_layer_branch_stats(
            result.problem2_solve_branch_stats,
            layer_stats_begin,
            result.problem2_solve_branch_stats.size(),
            static_cast<int>(depth),
            layer_from_planet,
            target_planet,
            summary.input_state_count,
            summary.output_state_count,
            summary.layer_wall_time_ms));
    }

    if (!result.terminal_solutions.empty()) {
        result.best_terminal_solution = *std::min_element(
            result.terminal_solutions.begin(), result.terminal_solutions.end(), terminal_less);
    }
    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
