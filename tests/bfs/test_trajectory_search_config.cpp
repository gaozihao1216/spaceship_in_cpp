/*
 * 文件作用：验证轨迹搜索全局参量与子配置工厂。
 */
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"
#include "spaceship_cpp/config/global_config.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::discretize_feasible_theta_interval;
using spaceship_cpp::bfs::discretize_leg0_theta_seeds;
using spaceship_cpp::bfs::find_leg0_feasible_theta_for_first_leg_targets;
using spaceship_cpp::bfs::FeasibleThetaInterval;
using spaceship_cpp::bfs::kDefaultMaxTotalTimeSeconds;
using spaceship_cpp::bfs::Leg0MultiTargetThetaResult;
using spaceship_cpp::bfs::Leg0TargetThetaFeasibility;
using spaceship_cpp::bfs::make_fixed_sequence_bfs_config;
using spaceship_cpp::bfs::make_leg0_theta_scan_config;
using spaceship_cpp::bfs::make_problem1_solve_defaults;
using spaceship_cpp::bfs::make_trajectory_search_discretization;
using spaceship_cpp::config::global_config;
using spaceship_cpp::planet_params::PlanetId;

bool near_equal(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-12;
}

}  // namespace

int main() {
    const auto& search_defaults = global_config().trajectory_search;

    {
        const auto defaults = default_trajectory_search_global_config();
        assert(defaults.discretization.leg0_theta_coarse_scan_count ==
            search_defaults.leg0_theta_coarse_scan_count);
        assert(defaults.discretization.leg0_theta_seeds_per_first_leg_target ==
            search_defaults.leg0_theta_seeds_per_first_leg_target);
        assert(defaults.max_search_legs == search_defaults.max_search_legs);
        assert(defaults.problem1.max_transfer_revolution ==
            search_defaults.max_transfer_revolution);
        assert(defaults.problem1.max_target_revolution ==
            search_defaults.max_target_revolution);
        assert(defaults.problem1.max_transfer_revolution == 1);
        assert(defaults.problem1.max_target_revolution == 1);
    }

    {
        constexpr double d = 0.2;
        const FeasibleThetaInterval interval{.start_rad = 1.0, .end_rad = 3.0};
        const auto samples = discretize_feasible_theta_interval(interval, d);
        assert(samples.size() >= 4U);
        assert(near_equal(samples[0], 2.0 - 0.5 * d));
        assert(near_equal(samples[1], 2.0 + 0.5 * d));
        assert(near_equal(samples[2], 2.0 - 1.5 * d));
        assert(near_equal(samples[3], 2.0 + 1.5 * d));
    }

    {
        Leg0MultiTargetThetaResult leg0{};
        leg0.ok = true;
        leg0.by_target = {
            Leg0TargetThetaFeasibility{
                .target_planet = PlanetId::Mars,
                .result{
                    .ok = true,
                    .feasible_intervals = {
                        FeasibleThetaInterval{.start_rad = 0.0, .end_rad = 1.0},
                        FeasibleThetaInterval{.start_rad = 2.0, .end_rad = 3.0},
                    },
                },
            },
            Leg0TargetThetaFeasibility{
                .target_planet = PlanetId::Venus,
                .result{
                    .ok = true,
                    .feasible_intervals = {
                        FeasibleThetaInterval{.start_rad = 4.0, .end_rad = 6.0},
                    },
                },
            },
        };

        const auto discretization = make_trajectory_search_discretization(search_defaults);
        const int seed_budget = discretization.leg0_theta_seeds_per_first_leg_target;
        const auto seeds = discretize_leg0_theta_seeds(leg0, discretization);
        assert(seeds.size() == static_cast<std::size_t>(seed_budget * 2));
        std::size_t mars_count = 0;
        std::size_t venus_count = 0;
        for (const auto& seed : seeds) {
            if (seed.first_leg_target_planet == PlanetId::Mars) {
                ++mars_count;
            } else if (seed.first_leg_target_planet == PlanetId::Venus) {
                ++venus_count;
            }
        }
        assert(mars_count == static_cast<std::size_t>(seed_budget));
        assert(venus_count == static_cast<std::size_t>(seed_budget));
    }

    auto global = default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.discretization.leg0_theta_coarse_scan_count = 90;

    const auto p1_defaults = make_problem1_solve_defaults(global);
    assert(p1_defaults.max_transfer_revolution == 1);
    assert(p1_defaults.max_target_revolution == 1);

    const auto mars_scan = make_leg0_theta_scan_config(global, PlanetId::Mars);
    assert(mars_scan.target_planet == PlanetId::Mars);
    assert(mars_scan.pruning.max_total_time_seconds == kDefaultMaxTotalTimeSeconds);
    assert(mars_scan.problem1.max_transfer_revolution == 1);

    const auto multi = find_leg0_feasible_theta_for_first_leg_targets(global);
    assert(multi.ok);
    assert(multi.by_target.size() == 3U);

    std::cout << "Global leg0 theta scan (default k=q=1, 90 samples/target):\n";
    for (const auto& entry : multi.by_target) {
        std::cout << "  target=" << static_cast<int>(entry.target_planet)
                  << " feasible_samples=" << entry.result.theta_samples_feasible
                  << " intervals=" << entry.result.feasible_intervals.size() << '\n';
    }

    const auto seeds = discretize_leg0_theta_seeds(global, multi);
    const int seed_budget = global.discretization.leg0_theta_seeds_per_first_leg_target;
    assert(seeds.size() <= static_cast<std::size_t>(seed_budget) * multi.by_target.size());
    std::cout << "Step2 theta seeds: count=" << seeds.size() << " (per-target budget="
              << seed_budget << ")\n";
    for (const auto& seed : seeds) {
        std::cout << "  target=" << static_cast<int>(seed.first_leg_target_planet)
                  << " theta=" << seed.transfer_theta_global << '\n';
    }

    const auto bfs_config = make_fixed_sequence_bfs_config(
        global,
        {PlanetId::Earth, PlanetId::Mars},
        0.5);
    assert(bfs_config.problem1.max_transfer_revolution == 1);
    assert(bfs_config.max_total_time_seconds == kDefaultMaxTotalTimeSeconds);

    return 0;
}
