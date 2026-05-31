#include "spaceship_cpp/bfs/fixed_sequence_launch_time_scan.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <chrono>
#include <utility>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

FixedSequenceThetaRefinementOptions make_refinement_options(
    const FixedSequenceLaunchTimeScanOptions& options
) {
    FixedSequenceThetaRefinementOptions refinement_options{};
    refinement_options.theta_center = options.theta_center;
    refinement_options.initial_half_width = options.theta_half_width;
    refinement_options.samples_per_round = options.theta_samples_per_round;
    refinement_options.refinement_rounds = options.theta_refinement_rounds;
    refinement_options.shrink_factor = options.theta_shrink_factor;
    refinement_options.include_center = true;
    refinement_options.beam_width = options.beam_width;
    refinement_options.max_launch_v_inf = options.max_launch_v_inf;
    refinement_options.time_weight_m_per_s_per_day = options.time_weight_m_per_s_per_day;
    return refinement_options;
}

FixedSequenceThetaRefinementResult run_single_theta_as_refinement_result(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceLaunchTimeScanOptions& options
) {
    FixedSequenceThetaRefinementResult result{};
    result.sequence = sequence;
    result.launch_time = launch_time;

    FixedSequenceThetaSearchOptions sequence_options{};
    sequence_options.beam_width = options.beam_width;
    sequence_options.max_launch_v_inf = options.max_launch_v_inf;
    sequence_options.time_weight_m_per_s_per_day = options.time_weight_m_per_s_per_day;

    const auto search = run_fixed_sequence_theta_search_with_table(
        loader,
        launch_time,
        options.theta_center,
        sequence,
        initial_launch_options,
        problem2_options,
        sequence_options);
    if (!search.ok) {
        result.error_message = search.error_message;
        return result;
    }

    FixedSequenceThetaRefinementRoundSummary summary{};
    summary.round_index = 0;
    summary.theta_left = options.theta_center;
    summary.theta_right = options.theta_center;
    summary.theta_center = options.theta_center;
    summary.half_width = 0.0;
    summary.theta_sample_count = 1;
    summary.terminal_theta_count = search.best_terminal_solution.valid ? 1 : 0;
    if (search.best_terminal_solution.valid) {
        const auto& terminal = search.best_terminal_solution;
        summary.best_valid = true;
        summary.best_theta = options.theta_center;
        summary.best_score = terminal.score;
        summary.best_total_delta_v = terminal.total_delta_v;
        summary.best_launch_v_inf = terminal.launch_v_inf;
        summary.best_arrival_v_inf = terminal.arrival_v_inf;
        summary.best_total_flight_time_seconds = terminal.total_flight_time_seconds;

        result.best_valid = true;
        result.best_theta = options.theta_center;
        result.best_score = terminal.score;
        result.best_total_delta_v = terminal.total_delta_v;
        result.best_launch_v_inf = terminal.launch_v_inf;
        result.best_arrival_v_inf = terminal.arrival_v_inf;
        result.best_total_flight_time_seconds = terminal.total_flight_time_seconds;
    }
    result.round_summaries.push_back(summary);
    result.ok = true;
    return result;
}

void copy_best_to_sample(
    FixedSequenceLaunchTimeScanSampleSummary* sample,
    const FixedSequenceThetaRefinementResult& refinement
) {
    sample->best_valid = refinement.best_valid;
    sample->best_theta = refinement.best_theta;
    sample->best_total_delta_v = refinement.best_total_delta_v;
    sample->best_launch_v_inf = refinement.best_launch_v_inf;
    sample->best_arrival_v_inf = refinement.best_arrival_v_inf;
    sample->best_total_flight_time_seconds = refinement.best_total_flight_time_seconds;
}

void copy_best_to_result(
    FixedSequenceLaunchTimeScanResult* result,
    const FixedSequenceLaunchTimeScanSampleSummary& sample,
    FixedSequenceThetaRefinementResult refinement
) {
    result->best_valid = true;
    result->best_launch_time = sample.launch_time;
    result->best_theta = sample.best_theta;
    result->best_total_delta_v = sample.best_total_delta_v;
    result->best_launch_v_inf = sample.best_launch_v_inf;
    result->best_arrival_v_inf = sample.best_arrival_v_inf;
    result->best_total_flight_time_seconds = sample.best_total_flight_time_seconds;
    result->best_refinement_result = std::move(refinement);
}

}  // namespace

FixedSequenceLaunchTimeScanResult scan_fixed_sequence_launch_time_locally_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceLaunchTimeScanOptions& options
) {
    FixedSequenceLaunchTimeScanResult result{};
    result.sequence = sequence;

    if (sequence.size() < 2 ||
        !is_finite(options.launch_time_center) ||
        !is_finite(options.initial_half_width_seconds) ||
        options.initial_half_width_seconds < 0.0 ||
        options.launch_time_sample_count <= 0 ||
        !is_finite(options.theta_center) ||
        !is_finite(options.theta_half_width) ||
        options.theta_half_width < 0.0 ||
        options.theta_samples_per_round <= 0 ||
        options.theta_refinement_rounds <= 0 ||
        !is_finite(options.theta_shrink_factor) ||
        options.theta_shrink_factor <= 0.0 ||
        options.beam_width < 1) {
        result.error_message = "invalid_fixed_sequence_launch_time_scan_options";
        return result;
    }

    const double left = options.launch_time_center - options.initial_half_width_seconds;
    const double step = options.launch_time_sample_count == 1
        ? 0.0
        : (2.0 * options.initial_half_width_seconds) /
            static_cast<double>(options.launch_time_sample_count - 1);

    for (int i = 0; i < options.launch_time_sample_count; ++i) {
        FixedSequenceLaunchTimeScanSampleSummary sample{};
        sample.valid = true;
        sample.launch_time = options.launch_time_sample_count == 1
            ? options.launch_time_center
            : left + step * static_cast<double>(i);
        sample.launch_time_offset_seconds = sample.launch_time - options.launch_time_center;

        const auto sample_start = Clock::now();
        auto refinement = options.refine_theta_at_each_time
            ? refine_fixed_sequence_theta_with_table(
                  loader,
                  sample.launch_time,
                  sequence,
                  initial_launch_options,
                  problem2_options,
                  make_refinement_options(options))
            : run_single_theta_as_refinement_result(
                  loader,
                  sample.launch_time,
                  sequence,
                  initial_launch_options,
                  problem2_options,
                  options);
        sample.wall_time_ms = elapsed_ms(sample_start, Clock::now());

        if (refinement.ok && refinement.best_valid) {
            copy_best_to_sample(&sample, refinement);
            if (!result.best_valid || sample.best_total_delta_v < result.best_total_delta_v) {
                copy_best_to_result(&result, sample, std::move(refinement));
            }
        }
        result.samples.push_back(sample);
    }

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
