#include "spaceship_cpp/bfs/fixed_launch_time_theta_scan.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

double wrapped_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

bool theta_candidate_less(const ThetaCandidate& lhs, const ThetaCandidate& rhs) {
    if (lhs.estimated_min_launch_v_inf != rhs.estimated_min_launch_v_inf) {
        return lhs.estimated_min_launch_v_inf < rhs.estimated_min_launch_v_inf;
    }
    return lhs.theta < rhs.theta;
}

bool theta_less(const ThetaCandidate& lhs, const ThetaCandidate& rhs) {
    return lhs.theta < rhs.theta;
}

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

bool contains_near_theta(
    const std::vector<ThetaCandidate>& candidates,
    double theta,
    const FixedLaunchTimeThetaScanOptions& options
) {
    if (!options.deduplicate_nearby_theta) {
        return false;
    }
    for (const auto& existing : candidates) {
        if (wrapped_distance(theta, existing.theta) <= options.duplicate_theta_tolerance) {
            return true;
        }
    }
    return false;
}

bool add_unique_candidate(
    std::vector<ThetaCandidate>* selected,
    const ThetaCandidate& candidate,
    const FixedLaunchTimeThetaScanOptions& options
) {
    if (!candidate.valid || !is_finite(candidate.theta)) {
        return false;
    }
    if (contains_near_theta(*selected, candidate.theta, options)) {
        return false;
    }
    selected->push_back(candidate);
    return true;
}

std::vector<ThetaCandidate> valid_unique_candidates_by_theta(
    const std::vector<ThetaCandidate>& candidates,
    const FixedLaunchTimeThetaScanOptions& options
) {
    auto sorted = candidates;
    std::sort(sorted.begin(), sorted.end(), theta_less);
    std::vector<ThetaCandidate> unique;
    for (const auto& candidate : sorted) {
        add_unique_candidate(&unique, candidate, options);
    }
    return unique;
}

void append_lowest_vinf_candidates(
    std::vector<ThetaCandidate>* selected,
    std::vector<ThetaCandidate> candidates,
    int count,
    const FixedLaunchTimeThetaScanOptions& options
) {
    if (count <= 0) {
        return;
    }
    std::sort(candidates.begin(), candidates.end(), theta_candidate_less);
    for (const auto& candidate : candidates) {
        if (static_cast<int>(selected->size()) >= options.max_theta_to_scan) {
            break;
        }
        add_unique_candidate(selected, candidate, options);
        if (static_cast<int>(selected->size()) >= count) {
            break;
        }
    }
}

void append_stratified_theta_candidates(
    std::vector<ThetaCandidate>* selected,
    const std::vector<ThetaCandidate>& candidates_by_theta,
    int count,
    const FixedLaunchTimeThetaScanOptions& options
) {
    if (count <= 0 || candidates_by_theta.empty()) {
        return;
    }
    const int initial_selected_count = static_cast<int>(selected->size());
    const int remaining_capacity = options.max_theta_to_scan - initial_selected_count;
    const int target_count = std::min(count, remaining_capacity);
    if (target_count <= 0) {
        return;
    }

    std::vector<ThetaCandidate> remaining;
    for (const auto& candidate : candidates_by_theta) {
        if (!contains_near_theta(*selected, candidate.theta, options)) {
            remaining.push_back(candidate);
        }
    }
    if (remaining.empty()) {
        return;
    }
    if (target_count >= static_cast<int>(remaining.size())) {
        for (const auto& candidate : remaining) {
            if (static_cast<int>(selected->size()) >= options.max_theta_to_scan) {
                break;
            }
            add_unique_candidate(selected, candidate, options);
        }
        return;
    }

    for (int i = 0; i < target_count; ++i) {
        const double position = target_count == 1
            ? 0.5 * static_cast<double>(remaining.size() - 1)
            : static_cast<double>(i) * static_cast<double>(remaining.size() - 1) /
                static_cast<double>(target_count - 1);
        const std::size_t index = static_cast<std::size_t>(std::llround(position));
        add_unique_candidate(selected, remaining[index], options);
    }

    for (const auto& candidate : remaining) {
        if (static_cast<int>(selected->size()) - initial_selected_count >= target_count ||
            static_cast<int>(selected->size()) >= options.max_theta_to_scan) {
            break;
        }
        add_unique_candidate(selected, candidate, options);
    }
}

}  // namespace

const char* theta_scan_selection_mode_name(ThetaScanSelectionMode mode) {
    switch (mode) {
    case ThetaScanSelectionMode::LowestEstimatedVInf:
        return "LowestEstimatedVInf";
    case ThetaScanSelectionMode::StratifiedByTheta:
        return "StratifiedByTheta";
    case ThetaScanSelectionMode::Hybrid:
        return "Hybrid";
    }
    return "Unknown";
}

std::vector<ThetaCandidate> select_fixed_launch_time_theta_scan_candidates(
    const std::vector<ThetaCandidate>& theta_candidates,
    const FixedLaunchTimeThetaScanOptions& options
) {
    const auto candidates_by_theta = valid_unique_candidates_by_theta(theta_candidates, options);
    std::vector<ThetaCandidate> selected;
    if (options.max_theta_to_scan <= 0) {
        return selected;
    }

    switch (options.theta_selection_mode) {
    case ThetaScanSelectionMode::LowestEstimatedVInf:
        append_lowest_vinf_candidates(&selected, candidates_by_theta, options.max_theta_to_scan, options);
        break;
    case ThetaScanSelectionMode::StratifiedByTheta:
        append_stratified_theta_candidates(&selected, candidates_by_theta, options.max_theta_to_scan, options);
        break;
    case ThetaScanSelectionMode::Hybrid: {
        const int lowest_count = std::min(options.hybrid_lowest_vinf_count, options.max_theta_to_scan);
        append_lowest_vinf_candidates(&selected, candidates_by_theta, lowest_count, options);
        const int stratified_count = std::min(
            options.hybrid_stratified_count,
            options.max_theta_to_scan - static_cast<int>(selected.size()));
        append_stratified_theta_candidates(&selected, candidates_by_theta, stratified_count, options);
        break;
    }
    }

    if (static_cast<int>(selected.size()) > options.max_theta_to_scan) {
        selected.resize(static_cast<std::size_t>(options.max_theta_to_scan));
    }
    return selected;
}

std::vector<planet_params::PlanetId> sequence_from_path(
    planet_params::PlanetId launch_planet,
    const std::vector<TrajectorySearchEdge>& path
) {
    std::vector<planet_params::PlanetId> sequence;
    sequence.reserve(path.size() + 1);
    sequence.push_back(launch_planet);
    for (const auto& edge : path) {
        sequence.push_back(edge.to_planet);
    }
    return sequence;
}

FixedLaunchTimeThetaScanResult scan_fixed_launch_time_theta_candidates_with_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    planet_params::PlanetId terminal_planet,
    double launch_time,
    const std::vector<ThetaCandidate>& theta_candidates,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const std::vector<planet_params::PlanetId>& allowed_transfer_planets,
    const InitialLaunchExpansionOptions& initial_launch_options,
    const problem2::Problem2GravityAssistSolverOptions& problem2_options,
    const FixedLaunchThetaSearchOptions& search_options,
    const FixedLaunchTimeThetaScanOptions& scan_options
) {
    FixedLaunchTimeThetaScanResult result{};
    result.launch_time = launch_time;
    result.theta_candidate_count = static_cast<int>(theta_candidates.size());
    if (!is_finite(launch_time) || scan_options.max_theta_to_scan <= 0 ||
        scan_options.top_solution_count < 0 || scan_options.top_sequence_count < 0) {
        result.error_message = "invalid_fixed_launch_time_theta_scan_options";
        return result;
    }

    const auto selected = select_fixed_launch_time_theta_scan_candidates(theta_candidates, scan_options);
    result.theta_scanned_count = static_cast<int>(selected.size());

    std::vector<FixedLaunchTimeThetaScanSolution> all_solutions;
    std::map<std::string, FixedLaunchTimeThetaScanSequenceSummary> sequence_summaries;

    for (const auto& candidate : selected) {
        const auto search = run_fixed_launch_theta_beam_search_with_table(
            loader,
            launch_planet,
            terminal_planet,
            launch_time,
            candidate.theta,
            allowed_first_targets,
            allowed_transfer_planets,
            initial_launch_options,
            problem2_options,
            search_options);

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
            result.terminal_theta_count += 1;

            FixedLaunchTimeThetaScanSolution solution{};
            solution.valid = true;
            solution.theta = candidate.theta;
            solution.estimated_min_launch_v_inf = candidate.estimated_min_launch_v_inf;
            solution.terminal_solution = search.best_terminal_solution;
            solution.path = reconstruct_fixed_launch_theta_path(search, search.best_terminal_solution.node_index);
            solution.sequence = sequence_from_path(launch_planet, solution.path);
            summary.best_sequence = solution.sequence;
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
    std::sort(result.top_sequences.begin(), result.top_sequences.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.best_score != rhs.best_score) {
            return lhs.best_score < rhs.best_score;
        }
        return lhs.occurrence_count > rhs.occurrence_count;
    });
    if (scan_options.top_sequence_count > 0 &&
        static_cast<int>(result.top_sequences.size()) > scan_options.top_sequence_count) {
        result.top_sequences.resize(static_cast<std::size_t>(scan_options.top_sequence_count));
    }

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
