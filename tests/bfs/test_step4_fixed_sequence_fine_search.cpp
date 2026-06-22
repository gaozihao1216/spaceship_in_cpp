/*
 * 文件作用：验证 Step 4 固定序列 θ 细扫 + BFS 精搜。
 */
#include "spaceship_cpp/bfs/step4_fixed_sequence_fine_search.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

constexpr double kDayInSeconds = 86400.0;

using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::find_leg0_feasible_theta_for_first_leg_targets;
using spaceship_cpp::bfs::Leg0ThetaSeed;
using spaceship_cpp::bfs::RankedPlanetSequence;
using spaceship_cpp::bfs::run_step2_free_path_search;
using spaceship_cpp::bfs::run_step3_select_top_k_sequences;
using spaceship_cpp::bfs::run_step4_fixed_sequence_fine_search;
using spaceship_cpp::bfs::search_fixed_sequence_over_theta_samples;
using spaceship_cpp::bfs::Step3TopKSequencesResult;
using spaceship_cpp::planet_params::PlanetId;

}  // namespace

int main() {
    auto global = default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 700.0 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
    global.discretization.leg0_theta_fine_scan_count = 5;

    {
        const std::vector<double> theta_samples{0.3, 0.5, 0.7};
        const auto result = search_fixed_sequence_over_theta_samples(
            global,
            {PlanetId::Earth, PlanetId::Mars},
            theta_samples);
        assert(result.theta_samples_tried == 3U);
        if (result.found_solution) {
            assert(std::isfinite(result.best_score));
            assert(result.best_score == result.best_bfs.best_score);
            assert(result.best_bfs.best_score ==
                result.best_bfs.best_launch_v_inf + result.best_bfs.best_arrival_v_inf);
            std::cout << "Direct theta sweep Earth->Mars: score=" << result.best_score
                      << " theta=" << result.best_theta_global
                      << " hits=" << result.theta_samples_with_solution << '\n';
        } else {
            std::cout << "Direct theta sweep Earth->Mars: no solution\n";
        }
    }

    {
        const auto leg0 = find_leg0_feasible_theta_for_first_leg_targets(global);
        assert(leg0.ok);

        Step3TopKSequencesResult step3{};
        step3.ok = true;
        step3.sequences.push_back(RankedPlanetSequence{
            .planet_sequence = {PlanetId::Earth, PlanetId::Mars},
            .best_seed = Leg0ThetaSeed{
                .first_leg_target_planet = PlanetId::Mars,
                .transfer_theta_global = 0.5,
            },
        });

        const auto step4 = run_step4_fixed_sequence_fine_search(global, step3, leg0);
        assert(step4.ok);
        assert(step4.by_sequence.size() == 1U);
        assert(step4.stats.total_bfs_calls == step4.by_sequence.front().theta_samples_tried);
        if (step4.by_sequence.front().found_solution) {
            assert(step4.best_sequence_index == 0);
            assert(step4.global_best_score == step4.by_sequence.front().best_score);
            std::cout << "Step4 with leg0 intervals: score=" << step4.global_best_score
                      << " theta=" << step4.by_sequence.front().best_theta_global
                      << " samples=" << step4.by_sequence.front().theta_samples_tried << '\n';
        }
    }

    {
        global.max_search_legs = 1;
        global.discretization.leg0_theta_coarse_scan_count = 30;
        global.discretization.leg0_theta_seeds_per_first_leg_target = 2;
        global.discretization.top_k_sequences = 1;
        global.discretization.leg0_theta_fine_scan_count = 4;

        const Leg0ThetaSeed seed{
            .first_leg_target_planet = PlanetId::Mars,
            .transfer_theta_global = 0.5,
        };
        const auto step2 = run_step2_free_path_search(global, {seed});
        assert(step2.ok);
        const auto step3 = run_step3_select_top_k_sequences(global, step2);
        assert(step3.ok);
        const auto step4 = run_step4_fixed_sequence_fine_search(global, step3);
        assert(step4.ok);
        if (step4.stats.sequences_with_solution > 0U) {
            assert(step4.best_sequence_index >= 0);
            assert(std::isfinite(step4.global_best_score));
            std::cout << "Pipeline Step2->4: global_score=" << step4.global_best_score
                      << " bfs_calls=" << step4.stats.total_bfs_calls << '\n';
        } else {
            std::cout << "Pipeline Step2->4: no fine-search solution\n";
        }
    }

    return 0;
}
