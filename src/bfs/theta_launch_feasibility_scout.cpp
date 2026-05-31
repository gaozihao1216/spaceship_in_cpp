#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"

#include "spaceship_cpp/bfs/initial_launch_expansion.hpp"
#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

struct CyclicRun {
    int start_index = 0;
    int count = 0;
};

ThetaLaunchFeasibilitySample evaluate_theta(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    double theta,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const ThetaLaunchFeasibilityScoutOptions& options
) {
    InitialLaunchExpansionOptions launch_options{};
    launch_options.max_transfer_revolution = options.max_transfer_revolution;
    launch_options.max_target_revolution = options.max_target_revolution;
    launch_options.problem1_residual_tolerance_seconds = options.problem1_residual_tolerance_seconds;
    launch_options.problem1_max_newton_iterations = options.problem1_max_newton_iterations;
    launch_options.max_launch_v_inf = std::numeric_limits<double>::infinity();
    launch_options.time_weight_m_per_s_per_day = 0.0;

    const auto launch = expand_initial_launch_with_problem1_table(
        loader, launch_planet, launch_time, normalize_angle_0_2pi(theta), allowed_first_targets, launch_options);

    ThetaLaunchFeasibilitySample sample{};
    sample.theta = normalize_angle_0_2pi(theta);
    if (!launch.ok) {
        return sample;
    }
    sample.raw_branch_count = launch.raw_branch_count;
    sample.accepted_candidate_count = launch.accepted_candidate_count;
    if (launch.candidates.empty()) {
        return sample;
    }

    for (const auto& candidate : launch.candidates) {
        if (!candidate.valid || !is_finite(candidate.launch_v_inf) || candidate.launch_v_inf < 0.0) {
            continue;
        }
        if (candidate.launch_v_inf < sample.min_launch_v_inf) {
            sample.second_min_launch_v_inf = sample.min_launch_v_inf;
            sample.min_launch_v_inf = candidate.launch_v_inf;
            sample.best_target_planet = candidate.target_planet;
            sample.best_transfer_revolution = candidate.transfer_revolution;
            sample.best_target_revolution = candidate.target_revolution;
            sample.best_arrival_time = candidate.arrival_time;
            sample.best_transfer_time_seconds = candidate.transfer_time_seconds;
        } else if (candidate.launch_v_inf < sample.second_min_launch_v_inf) {
            sample.second_min_launch_v_inf = candidate.launch_v_inf;
        }
    }
    sample.valid = is_finite(sample.min_launch_v_inf);
    return sample;
}

std::vector<CyclicRun> collect_cyclic_runs(
    const std::vector<ThetaLaunchFeasibilitySample>& samples,
    const std::function<bool(const ThetaLaunchFeasibilitySample&)>& predicate
) {
    const int n = static_cast<int>(samples.size());
    std::vector<CyclicRun> runs;
    if (n == 0) {
        return runs;
    }

    int false_index = -1;
    int true_count = 0;
    for (int i = 0; i < n; ++i) {
        if (predicate(samples[static_cast<std::size_t>(i)])) {
            true_count += 1;
        } else if (false_index < 0) {
            false_index = i;
        }
    }
    if (true_count == 0) {
        return runs;
    }
    if (true_count == n) {
        runs.push_back(CyclicRun{0, n});
        return runs;
    }

    int offset = 1;
    while (offset <= n) {
        const int index = (false_index + offset) % n;
        if (!predicate(samples[static_cast<std::size_t>(index)])) {
            offset += 1;
            continue;
        }
        CyclicRun run{};
        run.start_index = index;
        while (offset <= n) {
            const int current = (false_index + offset) % n;
            if (!predicate(samples[static_cast<std::size_t>(current)])) {
                break;
            }
            run.count += 1;
            offset += 1;
        }
        runs.push_back(run);
    }
    return runs;
}

ThetaLaunchFeasibilityInterval make_interval_from_run(
    const std::vector<ThetaLaunchFeasibilitySample>& samples,
    const CyclicRun& run,
    double theta_step
) {
    ThetaLaunchFeasibilityInterval interval{};
    interval.valid = run.count > 0;
    interval.sample_count_inside = run.count;
    const int n = static_cast<int>(samples.size());
    interval.theta_left = static_cast<double>(run.start_index) * theta_step;
    interval.theta_right = interval.theta_left + static_cast<double>(run.count) * theta_step;
    for (int k = 0; k < run.count; ++k) {
        const int index = (run.start_index + k) % n;
        const auto& sample = samples[static_cast<std::size_t>(index)];
        if (sample.min_launch_v_inf < interval.min_v_inf_inside) {
            interval.min_v_inf_inside = sample.min_launch_v_inf;
            interval.theta_at_min = sample.theta;
            if (index < run.start_index) {
                interval.theta_at_min += kTwoPi;
            }
        }
    }
    return interval;
}

std::vector<ThetaLaunchFeasibilityInterval> collect_cyclic_intervals(
    const std::vector<ThetaLaunchFeasibilitySample>& samples,
    double theta_step,
    const std::function<bool(const ThetaLaunchFeasibilitySample&)>& predicate
) {
    std::vector<ThetaLaunchFeasibilityInterval> intervals;
    const auto runs = collect_cyclic_runs(samples, predicate);
    intervals.reserve(runs.size());
    for (const auto& run : runs) {
        intervals.push_back(make_interval_from_run(samples, run, theta_step));
    }
    return intervals;
}

double g_at_theta(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    double theta,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const ThetaLaunchFeasibilityScoutOptions& options
) {
    const auto sample = evaluate_theta(loader, launch_planet, launch_time, theta, allowed_first_targets, options);
    if (!sample.valid) {
        return std::numeric_limits<double>::infinity();
    }
    return sample.min_launch_v_inf - options.max_launch_v_inf;
}

bool refine_boundary(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const ThetaLaunchFeasibilityScoutOptions& options,
    double invalid_theta,
    double valid_theta,
    double* refined_theta
) {
    double invalid_g = g_at_theta(loader, launch_planet, launch_time, invalid_theta, allowed_first_targets, options);
    double valid_g = g_at_theta(loader, launch_planet, launch_time, valid_theta, allowed_first_targets, options);
    if (!(invalid_g > 0.0 && valid_g <= 0.0)) {
        return false;
    }
    double left = invalid_theta;
    double right = valid_theta;
    for (int i = 0; i < options.max_boundary_refine_iterations; ++i) {
        if (std::abs(right - left) <= options.theta_boundary_tolerance) {
            break;
        }
        const double mid = 0.5 * (left + right);
        const double mid_g = g_at_theta(loader, launch_planet, launch_time, mid, allowed_first_targets, options);
        if (mid_g <= 0.0) {
            right = mid;
            valid_g = mid_g;
        } else {
            left = mid;
            invalid_g = mid_g;
        }
        (void)valid_g;
        (void)invalid_g;
    }
    *refined_theta = 0.5 * (left + right);
    return true;
}

void refine_valid_intervals(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const ThetaLaunchFeasibilityScoutOptions& options,
    const std::vector<ThetaLaunchFeasibilitySample>& samples,
    double theta_step,
    std::vector<ThetaLaunchFeasibilityInterval>* intervals
) {
    if (!options.refine_interval_boundaries || intervals->empty()) {
        return;
    }
    const auto runs = collect_cyclic_runs(samples, [&](const auto& sample) {
        return sample.min_launch_v_inf <= options.max_launch_v_inf;
    });
    if (runs.size() != intervals->size()) {
        return;
    }
    for (std::size_t i = 0; i < runs.size(); ++i) {
        const auto& run = runs[i];
        if (run.count == static_cast<int>(samples.size())) {
            continue;
        }
        const double left_valid = static_cast<double>(run.start_index) * theta_step;
        const double left_invalid = left_valid - theta_step;
        double refined_left = 0.0;
        if (refine_boundary(
                loader, launch_planet, launch_time, allowed_first_targets, options,
                left_invalid, left_valid, &refined_left)) {
            (*intervals)[i].theta_left = refined_left;
        }

        const double right_valid = left_valid + static_cast<double>(run.count - 1) * theta_step;
        const double right_invalid = left_valid + static_cast<double>(run.count) * theta_step;
        double refined_right = 0.0;
        if (refine_boundary(
                loader, launch_planet, launch_time, allowed_first_targets, options,
                right_invalid, right_valid, &refined_right)) {
            (*intervals)[i].theta_right = refined_right;
        }
    }
}

}  // namespace

ThetaLaunchFeasibilityScoutResult scout_theta_launch_feasibility_with_problem1_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const ThetaLaunchFeasibilityScoutOptions& options
) {
    ThetaLaunchFeasibilityScoutResult result{};
    result.theta_scout_count = options.theta_scout_count;
    if (options.theta_scout_count < 16) {
        result.error_message = "theta_scout_count_must_be_at_least_16";
        return result;
    }
    if (!is_finite(launch_time) || !is_finite(options.max_launch_v_inf) ||
        !is_finite(options.near_v_inf_buffer) || options.near_v_inf_buffer < 0.0) {
        result.error_message = "invalid_scout_options";
        return result;
    }

    const double theta_step = kTwoPi / static_cast<double>(options.theta_scout_count);
    result.samples.reserve(static_cast<std::size_t>(options.theta_scout_count));
    for (int i = 0; i < options.theta_scout_count; ++i) {
        const double theta = theta_step * static_cast<double>(i);
        auto sample = evaluate_theta(loader, launch_planet, launch_time, theta, allowed_first_targets, options);
        if (!sample.valid) {
            result.no_candidate_sample_count += 1;
        }
        if (sample.min_launch_v_inf <= options.max_launch_v_inf) {
            result.valid_sample_count += 1;
        }
        if (sample.min_launch_v_inf <= options.max_launch_v_inf + options.near_v_inf_buffer) {
            result.near_valid_sample_count += 1;
        }
        if (sample.min_launch_v_inf < result.global_min_launch_v_inf) {
            result.global_min_launch_v_inf = sample.min_launch_v_inf;
            result.theta_at_global_min = sample.theta;
            result.best_sample = sample;
        }
        result.samples.push_back(sample);
    }

    result.valid_intervals = collect_cyclic_intervals(result.samples, theta_step, [&](const auto& sample) {
        return sample.min_launch_v_inf <= options.max_launch_v_inf;
    });
    result.near_valid_intervals = collect_cyclic_intervals(result.samples, theta_step, [&](const auto& sample) {
        return sample.min_launch_v_inf <= options.max_launch_v_inf + options.near_v_inf_buffer;
    });
    refine_valid_intervals(
        loader, launch_planet, launch_time, allowed_first_targets, options,
        result.samples, theta_step, &result.valid_intervals);

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
