/*
 * 文件作用：验证 Step 2 自由路径 BFS。
 */
#include "spaceship_cpp/bfs/free_path_bfs.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

constexpr double kDayInSeconds = 86400.0;

using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::Leg0ThetaSeed;
using spaceship_cpp::bfs::run_step2_free_path_search;
using spaceship_cpp::bfs::search_free_path_from_leg0_seed;
using spaceship_cpp::planet_params::PlanetId;

}  // namespace

int main() {
    auto global = default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 700.0 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    {
        const Leg0ThetaSeed seed{
            .first_leg_target_planet = PlanetId::Mars,
            .transfer_theta_global = 0.5,
        };
        // 单跳 Earth→Mars 只需 1 段转移（含 leg0）。
        global.max_search_legs = 1;

        const auto result = search_free_path_from_leg0_seed(global, seed);
        assert(result.ok);
        assert(!result.solutions.empty());
        assert(result.solutions.front().planet_sequence.size() == 2U);
        assert(result.solutions.front().planet_sequence.front() == PlanetId::Earth);
        assert(result.solutions.front().planet_sequence.back() == PlanetId::Mars);
        assert(std::isfinite(result.solutions.front().score));
        std::cout << "Single seed Earth->Mars: solutions=" << result.solutions.size()
                  << " score=" << result.solutions.front().score
                  << " expanded=" << result.stats.expanded_nodes << '\n';
    }

    {
        global.max_search_legs = 1;
        const Leg0ThetaSeed seed{
            .first_leg_target_planet = PlanetId::Mars,
            .transfer_theta_global = 0.5,
        };
        const auto step2 = run_step2_free_path_search(global, {seed});
        assert(step2.ok);
        assert(step2.by_seed.size() == 1U);
        assert(!step2.by_seed.front().solutions.empty());
        std::cout << "Step2 single seed: solutions="
                  << step2.by_seed.front().solutions.size() << '\n';
    }

    return 0;
}
