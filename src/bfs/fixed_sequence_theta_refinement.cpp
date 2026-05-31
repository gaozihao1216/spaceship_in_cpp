#include "spaceship_cpp/bfs/fixed_sequence_theta_refinement.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::normalize_angle_0_2pi;

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

FixedSequenceThetaSearchOptions make_sequence_options(
    const FixedSequenceThetaRefinementOptions& options
) {
    FixedSequenceThetaSearchOptions sequence_options{};
    sequence_options.beam_width = options.beam_width;
    sequence_options.max_launch_v_inf = options.max_launch_v_inf;
    sequence_options.time_weight_m_per_s_per_day = options.time_weight_m_per_s_per_day;
    sequence_options.max_abs_slingshot_residual = options.max_abs_slingshot_residual;
    sequence_options.max_abs_problem1_residual_seconds = options.max_abs_problem1_residual_seconds;
    sequence_options.max_arrival_time = options.max_arrival_time;
    return sequence_options;
}

std::vector<double> sample_theta_interval(
    double theta_center,
    double half_width,
    int samples_per_round,
    bool include_center
) {
    std::vector<double> samples;
    if (samples_per_round <= 0) {
        return samples;
    }
    samples.reserve(static_cast<std::size_t>(samples_per_round));
    if (samples_per_round == 1) {
        samples.push_back(normalize_angle_0_2pi(theta_center));
        return samples;
    }

    const double left = theta_center - half_width;
    const double step = (2.0 * half_width) / static_cast<double>(samples_per_round - 1);
    for (int i = 0; i < samples_per_round; ++i) {
        samples.push_back(normalize_angle_0_2pi(left + step * static_cast<double>(i)));
    }

    if (include_center) {
        const double normalized_center = normalize_angle_0_2pi(theta_center);
        auto center_it = std::min_element(samples.begin(), samples.end(), [&](double lhs, double rhs) {
            return std::abs(normalize_angle_0_2pi(lhs - normalized_center)) <
                std::abs(normalize_angle_0_2pi(rhs - normalized_center));
        });
        if (center_it != samples.end()) {
            *center_it = normalized_center;
        }
    }
    return samples;
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

void copy_terminal_to_round_summary(
    FixedSequenceThetaRefinementRoundSummary* summary,
    double theta,
    const FixedLaunchThetaTerminalSolution& terminal
) {
    summary->best_valid = true;
    summary->best_theta = theta;
    summary->best_score = terminal.score;
    summary->best_total_delta_v = terminal.total_delta_v;
    summary->best_launch_v_inf = terminal.launch_v_inf;
    summary->best_arrival_v_inf = terminal.arrival_v_inf;
    summary->best_total_flight_time_seconds = terminal.total_flight_time_seconds;
}

void copy_solution_to_result(
    FixedSequenceThetaRefinementResult* result,
    FixedLaunchTimeThetaScanSolution solution
) {
    result->best_valid = true;
    result->best_theta = solution.theta;
    result->best_score = solution.terminal_solution.score;
    result->best_total_delta_v = solution.terminal_solution.total_delta_v;
    result->best_launch_v_inf = solution.terminal_solution.launch_v_inf;
    result->best_arrival_v_inf = solution.terminal_solution.arrival_v_inf;
    result->best_total_flight_time_seconds = solution.terminal_solution.total_flight_time_seconds;
    result->best_solution = std::move(solution);
}

}  // namespace

FixedSequenceThetaRefinementResult refine_fixed_sequence_theta_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceThetaRefinementOptions& refinement_options
) {
    FixedSequenceThetaRefinementResult result{};
    result.sequence = sequence;
    result.launch_time = launch_time;

    if (!is_finite(launch_time) || sequence.size() < 2 ||
        !is_finite(refinement_options.theta_center) ||
        !is_finite(refinement_options.initial_half_width) ||
        refinement_options.initial_half_width <= 0.0 ||
        refinement_options.samples_per_round <= 0 ||
        refinement_options.refinement_rounds <= 0 ||
        !is_finite(refinement_options.shrink_factor) ||
        refinement_options.shrink_factor <= 0.0 ||
        refinement_options.beam_width < 1) {
        result.error_message = "invalid_fixed_sequence_theta_refinement_options";
        return result;
    }

    const auto sequence_options = make_sequence_options(refinement_options);
    double theta_center = normalize_angle_0_2pi(refinement_options.theta_center);
    double half_width = refinement_options.initial_half_width;

    for (int round_index = 0; round_index < refinement_options.refinement_rounds; ++round_index) {
        FixedSequenceThetaRefinementRoundSummary summary{};
        summary.round_index = round_index;
        summary.theta_center = theta_center;
        summary.half_width = half_width;
        summary.theta_left = theta_center - half_width;
        summary.theta_right = theta_center + half_width;

        const auto theta_samples = sample_theta_interval(
            theta_center,
            half_width,
            refinement_options.samples_per_round,
            refinement_options.include_center);
        summary.theta_sample_count = static_cast<int>(theta_samples.size());

        const auto round_start = Clock::now();
        FixedLaunchTimeThetaScanSolution round_best_solution{};
        for (const double theta : theta_samples) {
            const auto search = run_fixed_sequence_theta_search_with_table(
                loader,
                launch_time,
                theta,
                sequence,
                initial_launch_options,
                problem2_options,
                sequence_options);
            if (!search.ok || search.terminal_solutions.empty() || !search.best_terminal_solution.valid) {
                continue;
            }
            summary.terminal_theta_count += 1;
            const auto& terminal = search.best_terminal_solution;
            if (!summary.best_valid || terminal_less(terminal, round_best_solution.terminal_solution)) {
                copy_terminal_to_round_summary(&summary, theta, terminal);

                round_best_solution.valid = true;
                round_best_solution.theta = theta;
                round_best_solution.estimated_min_launch_v_inf = std::numeric_limits<double>::infinity();
                round_best_solution.terminal_solution = terminal;
                round_best_solution.path = reconstruct_fixed_launch_theta_path(search, terminal.node_index);
                round_best_solution.sequence = sequence;
            }
        }
        summary.wall_time_ms = elapsed_ms(round_start, Clock::now());
        result.round_summaries.push_back(summary);

        if (!summary.best_valid) {
            result.ok = true;
            return result;
        }

        if (!result.best_valid || round_best_solution.terminal_solution.score < result.best_score) {
            copy_solution_to_result(&result, std::move(round_best_solution));
        }

        theta_center = normalize_angle_0_2pi(summary.best_theta);
        half_width *= refinement_options.shrink_factor;
    }

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
