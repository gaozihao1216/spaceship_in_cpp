#include "spaceship_cpp/bfs/theta_candidate_selector.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

double wrapped_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

bool theta_inside_interval(double theta, const ThetaLaunchFeasibilityInterval& interval, double tolerance) {
    double unwrapped = theta;
    if (unwrapped + tolerance < interval.theta_left) {
        unwrapped += kTwoPi;
    }
    return unwrapped + tolerance >= interval.theta_left && unwrapped <= interval.theta_right + tolerance;
}

bool overlaps_valid_intervals(
    const ThetaLaunchFeasibilityInterval& interval,
    const std::vector<ThetaLaunchFeasibilityInterval>& valid_intervals
) {
    for (const auto& valid : valid_intervals) {
        const double left = std::max(interval.theta_left, valid.theta_left);
        const double right = std::min(interval.theta_right, valid.theta_right);
        if (right > left) {
            return true;
        }
    }
    return false;
}

double estimate_v_inf_from_nearest_sample(
    const ThetaLaunchFeasibilityScoutResult& scout,
    double theta
) {
    double best_distance = std::numeric_limits<double>::infinity();
    double best_v_inf = std::numeric_limits<double>::infinity();
    for (const auto& sample : scout.samples) {
        if (!sample.valid || !is_finite(sample.min_launch_v_inf)) {
            continue;
        }
        const double distance = wrapped_distance(theta, sample.theta);
        if (distance < best_distance) {
            best_distance = distance;
            best_v_inf = sample.min_launch_v_inf;
        }
    }
    return best_v_inf;
}

void add_or_merge_candidate(
    std::vector<ThetaCandidate>* candidates,
    ThetaCandidate candidate,
    double duplicate_theta_tolerance
) {
    if (!is_finite(candidate.theta)) {
        return;
    }
    candidate.valid = true;
    candidate.theta = normalize_angle_0_2pi(candidate.theta);
    for (auto& existing : *candidates) {
        if (wrapped_distance(existing.theta, candidate.theta) > duplicate_theta_tolerance) {
            continue;
        }
        if (candidate.estimated_min_launch_v_inf < existing.estimated_min_launch_v_inf) {
            const bool was_boundary = existing.is_interval_boundary;
            const bool was_midpoint = existing.is_interval_midpoint;
            const bool was_minimum = existing.is_interval_minimum;
            const bool was_valid = existing.from_valid_interval;
            const bool was_near = existing.from_near_valid_interval;
            existing = candidate;
            existing.is_interval_boundary = existing.is_interval_boundary || was_boundary;
            existing.is_interval_midpoint = existing.is_interval_midpoint || was_midpoint;
            existing.is_interval_minimum = existing.is_interval_minimum || was_minimum;
            existing.from_valid_interval = existing.from_valid_interval || was_valid;
            existing.from_near_valid_interval = existing.from_near_valid_interval || was_near;
        } else {
            existing.is_interval_boundary = existing.is_interval_boundary || candidate.is_interval_boundary;
            existing.is_interval_midpoint = existing.is_interval_midpoint || candidate.is_interval_midpoint;
            existing.is_interval_minimum = existing.is_interval_minimum || candidate.is_interval_minimum;
            existing.from_valid_interval = existing.from_valid_interval || candidate.from_valid_interval;
            existing.from_near_valid_interval = existing.from_near_valid_interval || candidate.from_near_valid_interval;
        }
        return;
    }
    candidates->push_back(candidate);
}

ThetaCandidate make_candidate(
    const ThetaLaunchFeasibilityScoutResult& scout,
    double theta,
    int interval_index,
    bool from_valid,
    bool from_near,
    bool boundary,
    bool midpoint,
    bool minimum,
    double explicit_estimate
) {
    ThetaCandidate candidate{};
    candidate.theta = theta;
    candidate.estimated_min_launch_v_inf =
        is_finite(explicit_estimate) ? explicit_estimate : estimate_v_inf_from_nearest_sample(scout, theta);
    candidate.from_valid_interval = from_valid;
    candidate.from_near_valid_interval = from_near;
    candidate.is_interval_boundary = boundary;
    candidate.is_interval_midpoint = midpoint;
    candidate.is_interval_minimum = minimum;
    candidate.interval_index = interval_index;
    return candidate;
}

int target_count_for_interval(
    const std::vector<ThetaLaunchFeasibilityInterval>& intervals,
    std::size_t index,
    const ThetaCandidateSelectorOptions& options
) {
    if (intervals.size() == 1) {
        return std::min(options.max_theta_candidates, options.max_samples_per_valid_interval);
    }
    double total_width = 0.0;
    for (const auto& interval : intervals) {
        total_width += std::max(0.0, interval.theta_right - interval.theta_left);
    }
    if (!(total_width > 0.0)) {
        return options.min_samples_per_valid_interval;
    }
    const double width = std::max(0.0, intervals[index].theta_right - intervals[index].theta_left);
    const int proportional = static_cast<int>(std::llround(
        static_cast<double>(options.max_theta_candidates) * width / total_width));
    return std::clamp(
        proportional,
        options.min_samples_per_valid_interval,
        options.max_samples_per_valid_interval);
}

bool is_priority_candidate(const ThetaCandidate& candidate) {
    return candidate.is_interval_boundary || candidate.is_interval_midpoint || candidate.is_interval_minimum;
}

void sort_candidates(std::vector<ThetaCandidate>* candidates) {
    std::sort(candidates->begin(), candidates->end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.estimated_min_launch_v_inf != rhs.estimated_min_launch_v_inf) {
            return lhs.estimated_min_launch_v_inf < rhs.estimated_min_launch_v_inf;
        }
        return lhs.theta < rhs.theta;
    });
}

void truncate_candidates(std::vector<ThetaCandidate>* candidates, int max_theta_candidates) {
    if (static_cast<int>(candidates->size()) <= max_theta_candidates) {
        return;
    }
    std::vector<ThetaCandidate> priority;
    std::vector<ThetaCandidate> ordinary;
    for (const auto& candidate : *candidates) {
        if (is_priority_candidate(candidate)) {
            priority.push_back(candidate);
        } else {
            ordinary.push_back(candidate);
        }
    }
    sort_candidates(&priority);
    sort_candidates(&ordinary);

    std::vector<ThetaCandidate> kept;
    kept.reserve(static_cast<std::size_t>(max_theta_candidates));
    for (const auto& candidate : priority) {
        if (static_cast<int>(kept.size()) >= max_theta_candidates) {
            break;
        }
        kept.push_back(candidate);
    }
    for (const auto& candidate : ordinary) {
        if (static_cast<int>(kept.size()) >= max_theta_candidates) {
            break;
        }
        kept.push_back(candidate);
    }
    sort_candidates(&kept);
    *candidates = std::move(kept);
}

}  // namespace

ThetaCandidateSelectionResult select_theta_candidates_from_launch_scout(
    const ThetaLaunchFeasibilityScoutResult& scout,
    const ThetaCandidateSelectorOptions& options
) {
    ThetaCandidateSelectionResult result{};
    result.valid_interval_count = static_cast<int>(scout.valid_intervals.size());
    result.near_valid_interval_count = static_cast<int>(scout.near_valid_intervals.size());
    if (!scout.ok) {
        result.error_message = "invalid_launch_scout_result";
        return result;
    }
    if (options.max_theta_candidates <= 0 ||
        options.min_samples_per_valid_interval <= 0 ||
        options.max_samples_per_valid_interval < options.min_samples_per_valid_interval ||
        options.duplicate_theta_tolerance < 0.0) {
        result.error_message = "invalid_theta_candidate_selector_options";
        return result;
    }

    std::vector<ThetaCandidate> candidates;
    for (std::size_t i = 0; i < scout.valid_intervals.size(); ++i) {
        const auto& interval = scout.valid_intervals[i];
        if (!interval.valid || !(interval.theta_right > interval.theta_left)) {
            continue;
        }
        const int interval_index = static_cast<int>(i);
        if (options.include_interval_boundaries) {
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, interval.theta_left, interval_index, true, false, true, false, false,
                               std::numeric_limits<double>::infinity()),
                options.duplicate_theta_tolerance);
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, interval.theta_right, interval_index, true, false, true, false, false,
                               std::numeric_limits<double>::infinity()),
                options.duplicate_theta_tolerance);
        }
        if (options.include_interval_midpoint) {
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, 0.5 * (interval.theta_left + interval.theta_right), interval_index,
                               true, false, false, true, false, std::numeric_limits<double>::infinity()),
                options.duplicate_theta_tolerance);
        }
        if (options.include_interval_minimum) {
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, interval.theta_at_min, interval_index, true, false, false, false, true,
                               interval.min_v_inf_inside),
                options.duplicate_theta_tolerance);
        }

        const int target_count = target_count_for_interval(scout.valid_intervals, i, options);
        const int interior_count = std::max(0, target_count - 2);
        for (int j = 1; j <= interior_count; ++j) {
            const double fraction = static_cast<double>(j) / static_cast<double>(interior_count + 1);
            const double theta = interval.theta_left + fraction * (interval.theta_right - interval.theta_left);
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, theta, interval_index, true, false, false, false, false,
                               std::numeric_limits<double>::infinity()),
                options.duplicate_theta_tolerance);
        }
    }

    if (options.include_near_valid_intervals && options.max_near_valid_extra_candidates > 0) {
        int added_near = 0;
        for (std::size_t i = 0; i < scout.near_valid_intervals.size(); ++i) {
            if (added_near >= options.max_near_valid_extra_candidates) {
                break;
            }
            const auto& interval = scout.near_valid_intervals[i];
            if (!interval.valid || overlaps_valid_intervals(interval, scout.valid_intervals)) {
                continue;
            }
            const int interval_index = static_cast<int>(i);
            const double midpoint = 0.5 * (interval.theta_left + interval.theta_right);
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, interval.theta_at_min, interval_index, false, true, false, false, true,
                               interval.min_v_inf_inside),
                options.duplicate_theta_tolerance);
            added_near += 1;
            if (added_near >= options.max_near_valid_extra_candidates) {
                break;
            }
            add_or_merge_candidate(
                &candidates,
                make_candidate(scout, midpoint, interval_index, false, true, false, true, false,
                               std::numeric_limits<double>::infinity()),
                options.duplicate_theta_tolerance);
            added_near += 1;
        }
    }

    sort_candidates(&candidates);
    truncate_candidates(&candidates, options.max_theta_candidates);
    result.candidates = std::move(candidates);
    result.candidate_count = static_cast<int>(result.candidates.size());
    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
