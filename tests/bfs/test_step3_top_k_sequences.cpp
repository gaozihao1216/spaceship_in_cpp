/*
 * 文件作用：验证 Step 3 行星序列 Top-K 筛选。
 */
#include "spaceship_cpp/bfs/step3_top_k_sequences.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

constexpr double kDayInSeconds = 86400.0;

using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::FreePathBfsResult;
using spaceship_cpp::bfs::FreePathBfsSolution;
using spaceship_cpp::bfs::Leg0ThetaSeed;
using spaceship_cpp::bfs::run_step2_free_path_search;
using spaceship_cpp::bfs::run_step3_select_top_k_sequences;
using spaceship_cpp::bfs::select_top_k_planet_sequences;
using spaceship_cpp::bfs::Step2FreePathSearchResult;
using spaceship_cpp::planet_params::PlanetId;

FreePathBfsSolution make_solution(
    std::vector<PlanetId> sequence,
    double score,
    double launch_v_inf,
    double arrival_v_inf
) {
    FreePathBfsSolution solution{};
    solution.planet_sequence = std::move(sequence);
    solution.score = score;
    solution.launch_v_inf = launch_v_inf;
    solution.arrival_v_inf = arrival_v_inf;
    solution.total_time_seconds = 100.0 * kDayInSeconds;
    return solution;
}

Step2FreePathSearchResult make_synthetic_step2() {
    Step2FreePathSearchResult step2{};
    step2.ok = true;

    FreePathBfsResult seed_a{};
    seed_a.ok = true;
    seed_a.seed = Leg0ThetaSeed{
        .first_leg_target_planet = PlanetId::Mars,
        .transfer_theta_global = 1.0,
    };
    seed_a.solutions = {
        make_solution({PlanetId::Earth, PlanetId::Mars}, 300.0, 100.0, 200.0),
        make_solution({PlanetId::Earth, PlanetId::Venus, PlanetId::Mars}, 250.0, 120.0, 130.0),
        make_solution({PlanetId::Earth, PlanetId::Venus, PlanetId::Mars}, 240.0, 110.0, 130.0),
    };

    FreePathBfsResult seed_b{};
    seed_b.ok = true;
    seed_b.seed = Leg0ThetaSeed{
        .first_leg_target_planet = PlanetId::Venus,
        .transfer_theta_global = 2.0,
    };
    seed_b.solutions = {
        make_solution({PlanetId::Earth, PlanetId::Venus}, 200.0, 80.0, 120.0),
        make_solution({PlanetId::Earth, PlanetId::Mars, PlanetId::Mercury}, 260.0, 140.0, 120.0),
    };

    step2.by_seed = {seed_a, seed_b};
    return step2;
}

}  // namespace

int main() {
    {
        const auto step2 = make_synthetic_step2();
        const auto step3 = select_top_k_planet_sequences(step2, 3);
        assert(step3.ok);
        assert(step3.stats.solutions_pooled == 5U);
        assert(step3.stats.unique_sequences == 4U);
        assert(step3.sequences.size() == 3U);

        assert(step3.sequences[0].planet_sequence ==
            (std::vector<PlanetId>{PlanetId::Earth, PlanetId::Venus}));
        assert(step3.sequences[0].best_score == 200.0);
        assert(step3.sequences[0].best_seed.first_leg_target_planet == PlanetId::Venus);

        assert(step3.sequences[1].planet_sequence ==
            (std::vector<PlanetId>{PlanetId::Earth, PlanetId::Venus, PlanetId::Mars}));
        assert(step3.sequences[1].best_score == 240.0);
        assert(step3.sequences[1].best_seed.first_leg_target_planet == PlanetId::Mars);

        assert(step3.sequences[2].planet_sequence ==
            (std::vector<PlanetId>{PlanetId::Earth, PlanetId::Mars, PlanetId::Mercury}));
        assert(step3.sequences[2].best_score == 260.0);
        assert(step3.sequences[2].best_seed.first_leg_target_planet == PlanetId::Venus);

        std::cout << "Synthetic Top-3:\n";
        for (const auto& entry : step3.sequences) {
            std::cout << "  score=" << entry.best_score
                      << " len=" << entry.planet_sequence.size() << '\n';
        }
    }

    {
        const auto step2 = make_synthetic_step2();
        const auto step3 = select_top_k_planet_sequences(step2, 10);
        assert(step3.ok);
        assert(step3.sequences.size() == 4U);
        assert(step3.stats.sequences_returned == 4U);
    }

    {
        Step2FreePathSearchResult step2{};
        step2.ok = false;
        step2.error_message = "failed";
        const auto step3 = select_top_k_planet_sequences(step2, 3);
        assert(!step3.ok);
        assert(step3.error_message == "failed");
    }

    {
        auto global = default_trajectory_search_global_config();
        global.mission.launch_time_seconds_since_j2000 = 0.0;
        global.constraints.max_total_time_seconds = 700.0 * kDayInSeconds;
        global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
        global.max_search_legs = 1;
        global.discretization.top_k_sequences = 2;

        const Leg0ThetaSeed seed{
            .first_leg_target_planet = PlanetId::Mars,
            .transfer_theta_global = 0.5,
        };
        const auto step2 = run_step2_free_path_search(global, {seed});
        assert(step2.ok);

        const auto step3 = run_step3_select_top_k_sequences(global, step2);
        assert(step3.ok);
        assert(step3.sequences.size() == 1U);
        assert(step3.sequences.front().planet_sequence.size() == 2U);
        assert(step3.sequences.front().planet_sequence.front() == PlanetId::Earth);
        assert(step3.sequences.front().planet_sequence.back() == PlanetId::Mars);
        assert(std::isfinite(step3.sequences.front().best_score));
        std::cout << "Integrated Step3: count=" << step3.sequences.size()
                  << " score=" << step3.sequences.front().best_score << '\n';
    }

    return 0;
}
