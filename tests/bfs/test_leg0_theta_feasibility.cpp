/*
 * 文件作用：验证固定发射时刻下 leg0 可行 θ 区间扫描。
 */
#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

using spaceship_cpp::bfs::FeasibleThetaInterval;
using spaceship_cpp::bfs::Leg0CandidateFilterStats;
using spaceship_cpp::bfs::Leg0PruningLimits;
using spaceship_cpp::bfs::Leg0ThetaScanConfig;
using spaceship_cpp::bfs::find_leg0_feasible_theta_intervals;
using spaceship_cpp::bfs::leg0_candidate_passes_pruning;
using spaceship_cpp::bfs::leg0_has_feasible_candidate_at_theta;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::config::make_problem1_solve_input;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem1::Problem1Candidate;

constexpr double kDayInSeconds = 86400.0;
constexpr double kFortyYearsInSeconds = 40.0 * 365.25 * kDayInSeconds;

bool theta_in_interval(double theta, const FeasibleThetaInterval& interval) {
    return theta >= interval.start_rad && theta < interval.end_rad;
}

bool theta_in_any_interval(
    double theta,
    const std::vector<FeasibleThetaInterval>& intervals
) {
    for (const FeasibleThetaInterval& interval : intervals) {
        if (theta_in_interval(theta, interval)) {
            return true;
        }
    }
    return false;
}

Leg0PruningLimits make_default_leg0_pruning() {
    Leg0PruningLimits pruning{};
    pruning.max_total_time_seconds = kFortyYearsInSeconds;
    pruning.max_launch_v_inf_mps = 7500.0;
    return pruning;
}

}  // namespace

int main() {
    {
        Leg0ThetaScanConfig config{};
        config.launch_time_seconds_since_j2000 = 0.0;
        config.target_planet = PlanetId::Mars;
        config.pruning = make_default_leg0_pruning();
        config.theta_start_rad = 0.0;
        config.theta_end_rad = kTwoPi;
        config.theta_scan_count = 180;

        const auto result = find_leg0_feasible_theta_intervals(config);
        assert(result.ok);
        assert(result.theta_samples_tested == 180U);

        const bool theta_half_feasible = leg0_has_feasible_candidate_at_theta(
            config.launch_time_seconds_since_j2000,
            config.target_planet,
            0.5,
            config.pruning,
            config.problem1,
            nullptr);
        assert(theta_half_feasible == theta_in_any_interval(0.5, result.feasible_intervals));

        std::cout << "Earth->Mars leg0 feasible theta: samples_feasible="
                  << result.theta_samples_feasible << '/' << result.theta_samples_tested
                  << " intervals=" << result.feasible_intervals.size() << '\n';
        for (const FeasibleThetaInterval& interval : result.feasible_intervals) {
            std::cout << "  [" << interval.start_rad << ", " << interval.end_rad << ") rad\n";
        }
    }

    {
        Leg0ThetaScanConfig config{};
        config.launch_time_seconds_since_j2000 = 0.0;
        config.target_planet = PlanetId::Mars;
        config.pruning = make_default_leg0_pruning();
        config.pruning.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
        config.theta_start_rad = 0.0;
        config.theta_end_rad = kTwoPi;
        config.theta_scan_count = 180;

        const auto capped = find_leg0_feasible_theta_intervals(
            Leg0ThetaScanConfig{
                .launch_time_seconds_since_j2000 = 0.0,
                .target_planet = PlanetId::Mars,
                .pruning = make_default_leg0_pruning(),
                .theta_start_rad = 0.0,
                .theta_end_rad = kTwoPi,
                .theta_scan_count = 180,
            });
        const auto uncapped = find_leg0_feasible_theta_intervals(config);
        assert(capped.ok && uncapped.ok);
        assert(uncapped.theta_samples_feasible >= capped.theta_samples_feasible);

        const bool theta_half_uncapped = leg0_has_feasible_candidate_at_theta(
            0.0,
            PlanetId::Mars,
            0.5,
            config.pruning,
            config.problem1,
            nullptr);
        assert(theta_half_uncapped);
        std::cout << "Uncapped launch v_inf: theta=0.5 feasible, samples_feasible="
                  << uncapped.theta_samples_feasible << '\n';
    }

    {
        Leg0PruningLimits pruning = make_default_leg0_pruning();
        pruning.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

        const auto solve_input = make_problem1_solve_input(
            PlanetId::Earth,
            PlanetId::Mars,
            0.0,
            0.5,
            spaceship_cpp::config::global_config().problem1_solve);
        const std::vector<Problem1Candidate> candidates =
            spaceship_cpp::problem1::solve_problem1(solve_input);

        bool saw_negative_e_raw = false;
        for (const Problem1Candidate& candidate : candidates) {
            if (candidate.residual_result.transfer_e_raw < 0.0) {
                saw_negative_e_raw = true;
                Leg0CandidateFilterStats stats{};
                assert(!leg0_candidate_passes_pruning(candidate, 0.0, pruning, &stats));
                assert(stats.dead_by_non_positive_transfer_e_raw >= 1U);
            }
            if (leg0_candidate_passes_pruning(candidate, 0.0, pruning, nullptr)) {
                assert(candidate.residual_result.transfer_e_raw > 0.0);
            }
        }
        std::cout << "transfer_e_raw>0 filter: candidates=" << candidates.size()
                  << " saw_negative_e_raw=" << saw_negative_e_raw << '\n';
    }

    {
        Leg0ThetaScanConfig config{};
        config.theta_scan_count = 1;
        const auto result = find_leg0_feasible_theta_intervals(config);
        assert(!result.ok);
        assert(result.error_message == "theta_scan_count_must_be_at_least_two");
    }

    return 0;
}
