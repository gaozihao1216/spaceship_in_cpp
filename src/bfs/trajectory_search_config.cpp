/*
 * 文件作用：实现轨迹搜索管线的全局参量与子配置工厂。
 */
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <vector>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::planet_params::PlanetId;

bool theta_in_half_open_interval(double theta, double start_rad, double end_rad) {
    return theta >= start_rad && theta < end_rad;
}

std::vector<double> discretize_feasible_theta_interval_impl(
    const FeasibleThetaInterval& interval,
    double spacing_d_rad
) {
    std::vector<double> thetas;
    const double start_rad = interval.start_rad;
    const double end_rad = interval.end_rad;
    const double width = end_rad - start_rad;
    if (!(width > 0.0) || !(spacing_d_rad > 0.0)) {
        return thetas;
    }

    const double mid = 0.5 * (start_rad + end_rad);

    for (int layer = 1; ; ++layer) {
        const double offset = static_cast<double>(2 * layer - 1) * spacing_d_rad * 0.5;
        const double left = mid - offset;
        const double right = mid + offset;
        bool added = false;
        if (theta_in_half_open_interval(left, start_rad, end_rad)) {
            thetas.push_back(left);
            added = true;
        }
        if (theta_in_half_open_interval(right, start_rad, end_rad)) {
            thetas.push_back(right);
            added = true;
        }
        if (!added) {
            break;
        }
    }

    if (thetas.empty() && theta_in_half_open_interval(mid, start_rad, end_rad)) {
        thetas.push_back(mid);
    }
    return thetas;
}

struct IntervalCandidateQueue {
    PlanetId target_planet{};
    std::vector<double> candidates{};
    std::size_t cursor = 0;
};

double sum_interval_widths_for_target(const Leg0TargetThetaFeasibility& entry) {
    double total_width = 0.0;
    for (const FeasibleThetaInterval& interval : entry.result.feasible_intervals) {
        const double width = interval.end_rad - interval.start_rad;
        if (width > 0.0) {
            total_width += width;
        }
    }
    return total_width;
}

std::vector<IntervalCandidateQueue> build_interval_candidate_queues(
    const Leg0TargetThetaFeasibility& entry,
    double spacing_d_rad
) {
    std::vector<IntervalCandidateQueue> queues;
    if (!(spacing_d_rad > 0.0)) {
        return queues;
    }

    for (const FeasibleThetaInterval& interval : entry.result.feasible_intervals) {
        std::vector<double> candidates =
            discretize_feasible_theta_interval_impl(interval, spacing_d_rad);
        if (candidates.empty()) {
            continue;
        }
        queues.push_back(IntervalCandidateQueue{
            .target_planet = entry.target_planet,
            .candidates = std::move(candidates),
            .cursor = 0,
        });
    }
    return queues;
}

std::vector<Leg0ThetaSeed> discretize_leg0_theta_seeds_for_target(
    const Leg0TargetThetaFeasibility& entry,
    int seed_budget_per_target
) {
    std::vector<Leg0ThetaSeed> seeds;
    if (!entry.result.ok || seed_budget_per_target <= 0) {
        return seeds;
    }

    const double total_width = sum_interval_widths_for_target(entry);
    if (!(total_width > 0.0)) {
        return seeds;
    }

    const double spacing_d_rad =
        total_width / static_cast<double>(seed_budget_per_target);

    std::vector<IntervalCandidateQueue> queues =
        build_interval_candidate_queues(entry, spacing_d_rad);
    if (queues.empty()) {
        return seeds;
    }

    seeds.reserve(static_cast<std::size_t>(seed_budget_per_target));
    while (static_cast<int>(seeds.size()) < seed_budget_per_target) {
        bool progressed = false;
        for (IntervalCandidateQueue& queue : queues) {
            if (static_cast<int>(seeds.size()) >= seed_budget_per_target) {
                break;
            }
            if (queue.cursor >= queue.candidates.size()) {
                continue;
            }
            seeds.push_back(Leg0ThetaSeed{
                .first_leg_target_planet = queue.target_planet,
                .transfer_theta_global = queue.candidates[queue.cursor],
            });
            ++queue.cursor;
            progressed = true;
        }
        if (!progressed) {
            break;
        }
    }
    return seeds;
}

}  // namespace

std::vector<double> discretize_leg0_theta_samples_for_target(
    const Leg0TargetThetaFeasibility& entry,
    int sample_budget
) {
    const std::vector<Leg0ThetaSeed> seeds =
        discretize_leg0_theta_seeds_for_target(entry, sample_budget);
    std::vector<double> samples;
    samples.reserve(seeds.size());
    for (const Leg0ThetaSeed& seed : seeds) {
        samples.push_back(seed.transfer_theta_global);
    }
    return samples;
}

std::vector<double> discretize_feasible_theta_interval(
    const FeasibleThetaInterval& interval,
    double spacing_d_rad
) {
    return discretize_feasible_theta_interval_impl(interval, spacing_d_rad);
}

std::vector<Leg0ThetaSeed> discretize_leg0_theta_seeds(
    const Leg0MultiTargetThetaResult& leg0_result,
    const TrajectorySearchDiscretization& discretization
) {
    std::vector<Leg0ThetaSeed> seeds;
    if (!leg0_result.ok || discretization.leg0_theta_seeds_per_first_leg_target <= 0) {
        return seeds;
    }

    const int seed_budget_per_target = discretization.leg0_theta_seeds_per_first_leg_target;
    for (const Leg0TargetThetaFeasibility& entry : leg0_result.by_target) {
        const std::vector<Leg0ThetaSeed> target_seeds =
            discretize_leg0_theta_seeds_for_target(entry, seed_budget_per_target);
        seeds.insert(seeds.end(), target_seeds.begin(), target_seeds.end());
    }
    return seeds;
}

std::vector<Leg0ThetaSeed> discretize_leg0_theta_seeds(
    const TrajectorySearchGlobalConfig& config,
    const Leg0MultiTargetThetaResult& leg0_result
) {
    return discretize_leg0_theta_seeds(leg0_result, config.discretization);
}

std::vector<double> discretize_leg0_theta_samples_for_target(
    const Leg0MultiTargetThetaResult& leg0_result,
    PlanetId first_leg_target_planet,
    int sample_budget
) {
    if (!leg0_result.ok || sample_budget <= 0) {
        return {};
    }
    for (const Leg0TargetThetaFeasibility& entry : leg0_result.by_target) {
        if (entry.target_planet == first_leg_target_planet) {
            return discretize_leg0_theta_samples_for_target(entry, sample_budget);
        }
    }
    return {};
}

TrajectorySearchDiscretization make_trajectory_search_discretization(
    const config::TrajectorySearchDefaults& defaults
) {
    return TrajectorySearchDiscretization{
        .leg0_theta_coarse_scan_count = defaults.leg0_theta_coarse_scan_count,
        .leg0_theta_fine_scan_count = defaults.leg0_theta_fine_scan_count,
        .leg0_theta_seeds_per_first_leg_target = defaults.leg0_theta_seeds_per_first_leg_target,
        .top_k_sequences = defaults.top_k_sequences,
    };
}

TrajectorySearchGlobalConfig default_trajectory_search_global_config() {
    TrajectorySearchGlobalConfig config{};
    const auto& search_defaults = config::global_config().trajectory_search;
    config.discretization = make_trajectory_search_discretization(search_defaults);
    config.max_search_legs = search_defaults.max_search_legs;
    config.problem1.max_transfer_revolution = search_defaults.max_transfer_revolution;
    config.problem1.max_target_revolution = search_defaults.max_target_revolution;
    config.constraints.flyby_physical_filter.enabled = true;
    config.constraints.flyby_physical_filter.mode =
        trajectory::FlybyPhysicalFilterMode::Enforce;
    return config;
}

Leg0PruningLimits make_leg0_pruning_limits(const TrajectorySearchConstraints& constraints) {
    return Leg0PruningLimits{
        .max_total_time_seconds = constraints.max_total_time_seconds,
        .max_launch_v_inf_mps = constraints.max_launch_v_inf_mps,
    };
}

config::Problem1SolveDefaults make_problem1_solve_defaults(
    const Leg0Problem1Options& problem1
) {
    config::Problem1SolveDefaults defaults = config::global_config().problem1_solve;
    defaults.max_transfer_revolution = problem1.max_transfer_revolution;
    defaults.max_target_revolution = problem1.max_target_revolution;
    if (problem1.phi_scan_count_override >= 0) {
        defaults.phi_scan_count = problem1.phi_scan_count_override;
    }
    return defaults;
}

config::Problem1SolveDefaults make_problem1_solve_defaults(
    const TrajectorySearchGlobalConfig& config
) {
    return make_problem1_solve_defaults(config.problem1);
}

Leg0ThetaScanConfig make_leg0_theta_scan_config(
    const TrajectorySearchGlobalConfig& config,
    planet_params::PlanetId first_leg_target_planet
) {
    Leg0ThetaScanConfig scan_config{};
    scan_config.launch_time_seconds_since_j2000 = config.mission.launch_time_seconds_since_j2000;
    scan_config.target_planet = first_leg_target_planet;
    scan_config.pruning = make_leg0_pruning_limits(config.constraints);
    scan_config.problem1 = config.problem1;
    scan_config.theta_start_rad = 0.0;
    scan_config.theta_end_rad = kTwoPi;
    scan_config.theta_scan_count = config.discretization.leg0_theta_coarse_scan_count;
    return scan_config;
}

FixedSequenceBfsConfig make_fixed_sequence_bfs_config(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<planet_params::PlanetId>& planet_sequence,
    double launch_transfer_theta_global
) {
    FixedSequenceBfsConfig bfs_config{};
    bfs_config.launch_time_seconds_since_j2000 = config.mission.launch_time_seconds_since_j2000;
    bfs_config.launch_transfer_theta_global = launch_transfer_theta_global;
    bfs_config.planet_sequence = planet_sequence;
    bfs_config.max_total_time_seconds = config.constraints.max_total_time_seconds;
    bfs_config.max_launch_v_inf_mps = config.constraints.max_launch_v_inf_mps;
    bfs_config.problem1 = config.problem1;
    bfs_config.flyby_physical_filter = config.constraints.flyby_physical_filter;
    return bfs_config;
}

Leg0MultiTargetThetaResult find_leg0_feasible_theta_for_first_leg_targets(
    const TrajectorySearchGlobalConfig& config
) {
    Leg0MultiTargetThetaResult result{};
    if (config.mission.first_leg_target_planets.empty()) {
        result.error_message = "first_leg_target_planets_must_not_be_empty";
        return result;
    }

    result.by_target.reserve(config.mission.first_leg_target_planets.size());
    for (const planet_params::PlanetId target_planet : config.mission.first_leg_target_planets) {
        Leg0TargetThetaFeasibility entry{};
        entry.target_planet = target_planet;
        entry.result = find_leg0_feasible_theta_intervals(
            make_leg0_theta_scan_config(config, target_planet));
        if (!entry.result.ok) {
            result.error_message = entry.result.error_message;
            return result;
        }
        result.by_target.push_back(std::move(entry));
    }

    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
