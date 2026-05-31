#include "spaceship_cpp/bfs/fixed_sequence_theta_scan.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <map>
#include <sstream>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;

bool scan_solution_less(
    const FixedLaunchTimeThetaScanSolution& lhs,
    const FixedLaunchTimeThetaScanSolution& rhs
) {
    if (lhs.terminal_solution.score != rhs.terminal_solution.score) {
        return lhs.terminal_solution.score < rhs.terminal_solution.score;
    }
    return lhs.theta < rhs.theta;
}

std::string sequence_key(const std::vector<planet_params::PlanetId>& sequence) {
    std::ostringstream os;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        if (i > 0) {
            os << '>';
        }
        os << static_cast<int>(sequence[i]);
    }
    return os.str();
}

}  // namespace

FixedLaunchTimeThetaScanResult scan_fixed_sequence_over_theta_candidates_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    double launch_time,
    const std::vector<ThetaCandidate>& theta_candidates,
    const std::vector<planet_params::PlanetId>& sequence,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedSequenceThetaSearchOptions& sequence_options,
    const FixedLaunchTimeThetaScanOptions& scan_options
) {
    FixedLaunchTimeThetaScanResult result{};
    result.launch_time = launch_time;
    result.theta_candidate_count = static_cast<int>(theta_candidates.size());
    if (!is_finite(launch_time) || sequence.size() < 2 || scan_options.max_theta_to_scan <= 0 ||
        scan_options.top_solution_count < 0 || scan_options.top_sequence_count < 0) {
        result.error_message = "invalid_fixed_sequence_theta_scan_options";
        return result;
    }

    const auto selected = select_fixed_launch_time_theta_scan_candidates(theta_candidates, scan_options);
    result.theta_scanned_count = static_cast<int>(selected.size());

    std::vector<FixedLaunchTimeThetaScanSolution> all_solutions;
    std::map<std::string, FixedLaunchTimeThetaScanSequenceSummary> sequence_summaries;

    for (const auto& candidate : selected) {
        const auto search = run_fixed_sequence_theta_search_with_table(
            loader,
            launch_time,
            candidate.theta,
            sequence,
            initial_launch_options,
            problem2_options,
            sequence_options);

        FixedLaunchTimeThetaScanThetaSummary summary{};
        summary.valid = search.ok;
        summary.theta = candidate.theta;
        summary.estimated_min_launch_v_inf = candidate.estimated_min_launch_v_inf;
        summary.initial_candidate_count = search.initial_candidate_count;
        summary.terminal_solution_count = static_cast<int>(search.terminal_solutions.size());
        summary.node_count = static_cast<int>(search.nodes.size());
        result.total_terminal_solution_count += summary.terminal_solution_count;

        if (search.best_terminal_solution.valid) {
            summary.has_terminal = true;
            summary.best_score = search.best_terminal_solution.score;
            summary.best_total_delta_v = search.best_terminal_solution.total_delta_v;
            summary.best_launch_v_inf = search.best_terminal_solution.launch_v_inf;
            summary.best_arrival_v_inf = search.best_terminal_solution.arrival_v_inf;
            summary.best_total_flight_time_seconds = search.best_terminal_solution.total_flight_time_seconds;
            summary.best_node_index = search.best_terminal_solution.node_index;
            summary.best_sequence = sequence;
            result.terminal_theta_count += 1;

            FixedLaunchTimeThetaScanSolution solution{};
            solution.valid = true;
            solution.theta = candidate.theta;
            solution.estimated_min_launch_v_inf = candidate.estimated_min_launch_v_inf;
            solution.terminal_solution = search.best_terminal_solution;
            solution.path = reconstruct_fixed_launch_theta_path(search, search.best_terminal_solution.node_index);
            solution.sequence = sequence;
            all_solutions.push_back(solution);

            auto& sequence_summary = sequence_summaries[sequence_key(solution.sequence)];
            if (!sequence_summary.valid) {
                sequence_summary.valid = true;
                sequence_summary.sequence = solution.sequence;
            }
            sequence_summary.occurrence_count += 1;
            if (solution.terminal_solution.score < sequence_summary.best_score) {
                sequence_summary.best_theta = solution.theta;
                sequence_summary.best_score = solution.terminal_solution.score;
                sequence_summary.best_total_delta_v = solution.terminal_solution.total_delta_v;
                sequence_summary.best_launch_v_inf = solution.terminal_solution.launch_v_inf;
                sequence_summary.best_arrival_v_inf = solution.terminal_solution.arrival_v_inf;
                sequence_summary.best_total_flight_time_seconds =
                    solution.terminal_solution.total_flight_time_seconds;
            }
        }
        result.theta_summaries.push_back(summary);
    }

    std::sort(all_solutions.begin(), all_solutions.end(), scan_solution_less);
    if (scan_options.top_solution_count > 0 &&
        static_cast<int>(all_solutions.size()) > scan_options.top_solution_count) {
        all_solutions.resize(static_cast<std::size_t>(scan_options.top_solution_count));
    }
    result.top_solutions = all_solutions;
    if (!result.top_solutions.empty()) {
        result.best_solution = result.top_solutions.front();
    }

    for (const auto& [key, summary] : sequence_summaries) {
        (void)key;
        result.top_sequences.push_back(summary);
    }
    if (scan_options.top_sequence_count > 0 &&
        static_cast<int>(result.top_sequences.size()) > scan_options.top_sequence_count) {
        result.top_sequences.resize(static_cast<std::size_t>(scan_options.top_sequence_count));
    }

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
