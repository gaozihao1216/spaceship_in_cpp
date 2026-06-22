/*
 * 文件作用：验证固定行星序列 BFS 能产出可行解并满足 TOF 与评分语义。
 */
#include "spaceship_cpp/bfs/fixed_sequence_bfs.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

using spaceship_cpp::bfs::FixedSequenceBfsConfig;
using spaceship_cpp::bfs::search_fixed_sequence_bfs;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::trajectory::FlybyPhysicalFeasibilityOptions;
using spaceship_cpp::trajectory::FlybyPhysicalFilterMode;

constexpr double kDayInSeconds = 86400.0;

FlybyPhysicalFeasibilityOptions make_default_bfs_filter() {
    FlybyPhysicalFeasibilityOptions options{};
    options.enabled = true;
    options.mode = FlybyPhysicalFilterMode::Enforce;
    return options;
}

}  // namespace

int main() {
    {
        FixedSequenceBfsConfig config{};
        config.launch_time_seconds_since_j2000 = 0.0;
        config.launch_transfer_theta_global = 0.5;
        config.planet_sequence = {PlanetId::Earth, PlanetId::Mars};
        config.max_total_time_seconds = 700.0 * kDayInSeconds;
        config.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
        config.flyby_physical_filter = make_default_bfs_filter();

        const auto result = search_fixed_sequence_bfs(config);
        if (!result.ok || !result.found_solution) {
            std::cerr << "Earth->Mars BFS failed: ok=" << result.ok
                      << " err=" << result.error_message
                      << " leg0=" << result.stats.leg0_candidates
                      << " dead_time=" << result.stats.dead_by_time_limit << '\n';
            return 1;
        }
        assert(std::isfinite(result.best_score));
        assert(result.best_score > 0.0);
        assert(std::isfinite(result.best_launch_v_inf));
        assert(std::isfinite(result.best_arrival_v_inf));
        assert(result.best_score == result.best_launch_v_inf + result.best_arrival_v_inf);
        assert(result.best_total_time_seconds <= config.max_total_time_seconds);
        assert(result.best_edges.size() == 1U);
        assert(result.best_edges.front().from_planet == PlanetId::Earth);
        assert(result.best_edges.front().to_planet == PlanetId::Mars);
        assert(result.stats.leg0_candidates > 0U);
        assert(result.stats.leaf_nodes > 0U);

        std::cout << "Earth->Mars BFS: score=" << result.best_score
                  << " tof_days=" << (result.best_total_time_seconds / kDayInSeconds)
                  << " leg0_candidates=" << result.stats.leg0_candidates
                  << " expanded=" << result.stats.expanded_nodes << '\n';
    }

    {
        FixedSequenceBfsConfig config{};
        config.launch_time_seconds_since_j2000 = 0.0;
        config.launch_transfer_theta_global = 0.5;
        config.planet_sequence = {PlanetId::Earth, PlanetId::Mars};
        config.max_total_time_seconds = 700.0 * kDayInSeconds;
        // 默认 max_launch_v_inf_mps = 7500；theta=0.5 的 Earth→Mars 发射 v_∞ 高于此上限。
        const auto result = search_fixed_sequence_bfs(config);
        assert(result.ok);
        assert(!result.found_solution);
        assert(result.stats.leg0_candidates >= 1U);
        assert(result.stats.dead_by_launch_v_inf_limit >= 1U);
        assert(result.stats.enqueued_nodes == 0U);
        std::cout << "Earth->Mars launch_v_inf cap: dead_launch_cap="
                  << result.stats.dead_by_launch_v_inf_limit << '\n';
    }

    {
        FixedSequenceBfsConfig config{};
        config.launch_time_seconds_since_j2000 = 0.0;
        config.launch_transfer_theta_global = 0.5;
        config.planet_sequence = {PlanetId::Earth, PlanetId::Mars, PlanetId::Jupiter};
        config.max_total_time_seconds = 3500.0 * kDayInSeconds;
        config.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
        config.flyby_physical_filter = make_default_bfs_filter();

        const auto result = search_fixed_sequence_bfs(config);
        assert(result.ok);
        if (result.found_solution) {
            assert(result.best_edges.size() == 2U);
            assert(result.best_edges.front().from_planet == PlanetId::Earth);
            assert(result.best_edges.back().to_planet == PlanetId::Jupiter);
            std::cout << "Earth->Mars->Jupiter BFS: score=" << result.best_score
                      << " tof_days=" << (result.best_total_time_seconds / kDayInSeconds)
                      << " expanded=" << result.stats.expanded_nodes
                      << " p2_calls=" << result.stats.p2_solve_calls << '\n';
        } else {
            std::cout << "Earth->Mars->Jupiter BFS: no solution within T_max"
                      << " leg0=" << result.stats.leg0_candidates
                      << " dead_time=" << result.stats.dead_by_time_limit
                      << " dead_launch_cap=" << result.stats.dead_by_launch_v_inf_limit
                      << " dead_p2=" << result.stats.dead_by_no_p2_solution
                      << " dead_turn=" << result.stats.dead_by_turn_angle << '\n';
        }
    }

    {
        FixedSequenceBfsConfig config{};
        config.planet_sequence = {PlanetId::Mars, PlanetId::Jupiter};
        const auto result = search_fixed_sequence_bfs(config);
        assert(!result.ok);
        assert(result.error_message == "planet_sequence_must_start_at_earth");
    }

    return 0;
}
