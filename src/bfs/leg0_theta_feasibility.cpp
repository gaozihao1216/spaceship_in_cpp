/*
 * 文件作用：实现固定发射时刻下 leg0 可行 θ 区间扫描。
 */
#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/trajectory/orbit_velocity.hpp"

#include <cmath>

namespace spaceship_cpp::bfs {
namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::config::global_config;
using spaceship_cpp::config::make_problem1_solve_input;
using spaceship_cpp::config::Problem1SolveDefaults;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem1::Problem1Candidate;
using spaceship_cpp::problem1::Problem1SolveInput;

bool validate_theta_scan_config(const Leg0ThetaScanConfig& config, std::string& error_message) {
    if (!is_finite(config.launch_time_seconds_since_j2000)) {
        error_message = "invalid_launch_time";
        return false;
    }
    if (!is_finite(config.theta_start_rad) || !is_finite(config.theta_end_rad)) {
        error_message = "invalid_theta_range";
        return false;
    }
    if (!(config.theta_end_rad > config.theta_start_rad)) {
        error_message = "theta_end_must_be_greater_than_theta_start";
        return false;
    }
    if (config.theta_scan_count < 2) {
        error_message = "theta_scan_count_must_be_at_least_two";
        return false;
    }
    if (!is_finite(config.pruning.max_total_time_seconds) ||
        !(config.pruning.max_total_time_seconds > 0.0)) {
        error_message = "invalid_max_total_time";
        return false;
    }
    if (is_finite(config.pruning.max_launch_v_inf_mps) &&
        !(config.pruning.max_launch_v_inf_mps > 0.0)) {
        error_message = "invalid_max_launch_v_inf";
        return false;
    }
    return true;
}

Problem1SolveInput make_leg0_solve_input(
    double launch_time_seconds_since_j2000,
    PlanetId target_planet,
    double transfer_theta_global,
    const Leg0Problem1Options& problem1
) {
    config::Problem1SolveDefaults defaults = config::global_config().problem1_solve;
    defaults.max_transfer_revolution = problem1.max_transfer_revolution;
    defaults.max_target_revolution = problem1.max_target_revolution;
    if (problem1.phi_scan_count_override >= 0) {
        defaults.phi_scan_count = problem1.phi_scan_count_override;
    }
    return make_problem1_solve_input(
        PlanetId::Earth,
        target_planet,
        launch_time_seconds_since_j2000,
        transfer_theta_global,
        defaults);
}

std::vector<FeasibleThetaInterval> merge_feasible_theta_samples(
    const std::vector<double>& sample_thetas,
    const std::vector<bool>& sample_feasible,
    double sample_step_rad
) {
    std::vector<FeasibleThetaInterval> intervals;
    if (sample_thetas.empty() || sample_thetas.size() != sample_feasible.size()) {
        return intervals;
    }

    std::size_t index = 0;
    while (index < sample_feasible.size()) {
        if (!sample_feasible[index]) {
            ++index;
            continue;
        }

        const double interval_start = sample_thetas[index];
        std::size_t last_feasible_index = index;
        ++index;
        while (index < sample_feasible.size() && sample_feasible[index]) {
            last_feasible_index = index;
            ++index;
        }

        intervals.push_back(FeasibleThetaInterval{
            .start_rad = interval_start,
            .end_rad = sample_thetas[last_feasible_index] + sample_step_rad,
        });
    }
    return intervals;
}

}  // namespace

bool leg0_candidate_passes_pruning(
    const Problem1Candidate& candidate,
    double launch_time_seconds_since_j2000,
    const Leg0PruningLimits& pruning,
    Leg0CandidateFilterStats* stats
) {
    if (stats != nullptr) {
        ++stats->raw_candidates;
    }

    const auto& residual = candidate.residual_result;
    // (-e, θ+π) 与 (e, θ) 物理等价；leg0 只保留 transfer_e_raw > 0 的正分支表示。
    if (!is_finite(residual.transfer_e_raw) || !(residual.transfer_e_raw > 0.0)) {
        if (stats != nullptr) {
            ++stats->dead_by_non_positive_transfer_e_raw;
        }
        return false;
    }

    if (!is_finite(candidate.time_of_flight_seconds) || !(candidate.time_of_flight_seconds > 0.0) ||
        !is_finite(candidate.arrival_time_seconds_since_j2000)) {
        return false;
    }
    if (candidate.time_of_flight_seconds > pruning.max_total_time_seconds) {
        if (stats != nullptr) {
            ++stats->dead_by_time_limit;
        }
        return false;
    }

    const double launch_v_inf = trajectory::relative_speed_to_planet(
        PlanetId::Earth,
        launch_time_seconds_since_j2000,
        residual.transfer_e,
        residual.transfer_perihelion_angle_used);
    if (!is_finite(launch_v_inf) || !(launch_v_inf > 0.0)) {
        if (stats != nullptr) {
            ++stats->dead_by_invalid_v_inf;
        }
        return false;
    }
    if (is_finite(pruning.max_launch_v_inf_mps) &&
        launch_v_inf > pruning.max_launch_v_inf_mps) {
        if (stats != nullptr) {
            ++stats->dead_by_launch_v_inf_limit;
        }
        return false;
    }

    if (stats != nullptr) {
        ++stats->passed;
    }
    return true;
}

bool leg0_has_feasible_candidate_at_theta(
    double launch_time_seconds_since_j2000,
    PlanetId target_planet,
    double transfer_theta_global,
    const Leg0PruningLimits& pruning,
    const Leg0Problem1Options& problem1,
    Leg0CandidateFilterStats* stats
) {
    const Problem1SolveInput solve_input = make_leg0_solve_input(
        launch_time_seconds_since_j2000,
        target_planet,
        transfer_theta_global,
        problem1);
    const std::vector<Problem1Candidate> candidates = problem1::solve_problem1(solve_input);

    for (const Problem1Candidate& candidate : candidates) {
        if (leg0_candidate_passes_pruning(
                candidate,
                launch_time_seconds_since_j2000,
                pruning,
                stats)) {
            return true;
        }
    }
    return false;
}

Leg0ThetaFeasibilityResult find_leg0_feasible_theta_intervals(const Leg0ThetaScanConfig& config) {
    Leg0ThetaFeasibilityResult result{};
    std::string error_message;
    if (!validate_theta_scan_config(config, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }

    const double span = config.theta_end_rad - config.theta_start_rad;
    const double sample_step = span / static_cast<double>(config.theta_scan_count);

    std::vector<double> sample_thetas;
    std::vector<bool> sample_feasible;
    sample_thetas.reserve(static_cast<std::size_t>(config.theta_scan_count));
    sample_feasible.reserve(static_cast<std::size_t>(config.theta_scan_count));

    for (int sample_index = 0; sample_index < config.theta_scan_count; ++sample_index) {
        const double theta =
            config.theta_start_rad + static_cast<double>(sample_index) * sample_step;
        sample_thetas.push_back(theta);

        const bool feasible = leg0_has_feasible_candidate_at_theta(
            config.launch_time_seconds_since_j2000,
            config.target_planet,
            theta,
            config.pruning,
            config.problem1,
            nullptr);
        sample_feasible.push_back(feasible);
        ++result.theta_samples_tested;
        if (feasible) {
            ++result.theta_samples_feasible;
        }
    }

    result.feasible_intervals =
        merge_feasible_theta_samples(sample_thetas, sample_feasible, sample_step);
    result.ok = true;
    return result;
}

}  // namespace spaceship_cpp::bfs
