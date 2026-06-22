/*
 * 文件作用：诊断 Earth 出发能否到达 Mercury，以及 max_search_legs 与死胡同叶节点记解的影响。
 */
#include "spaceship_cpp/bfs/fixed_sequence_bfs.hpp"
#include "spaceship_cpp/bfs/free_path_bfs.hpp"
#include "spaceship_cpp/bfs/leg0_theta_feasibility.hpp"
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"

#include <iostream>
#include <limits>

namespace {

using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::FixedSequenceBfsConfig;
using spaceship_cpp::bfs::Leg0ThetaSeed;
using spaceship_cpp::bfs::search_fixed_sequence_bfs;
using spaceship_cpp::bfs::search_free_path_from_leg0_seed;
using spaceship_cpp::bfs::leg0_has_feasible_candidate_at_theta;
using spaceship_cpp::bfs::make_leg0_pruning_limits;
using spaceship_cpp::planet_params::PlanetId;

void print_sequence(const std::vector<PlanetId>& seq) {
    for (PlanetId p : seq) {
        std::cout << static_cast<int>(p) << ' ';
    }
}

}  // namespace

int main() {
    auto global = default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
    const auto pruning = make_leg0_pruning_limits(global.constraints);

    std::cout << "=== leg0 Earth->Mercury feasibility ===\n";
    for (double theta : {4.1, 4.8, 1.4}) {
        const bool ok = leg0_has_feasible_candidate_at_theta(
            0.0, PlanetId::Mercury, theta, pruning, global.problem1, nullptr);
        std::cout << "  Mercury theta=" << theta << " feasible=" << ok << '\n';
    }

    std::cout << "\n=== fixed sequence Earth->Mercury (1 leg) ===\n";
    {
        FixedSequenceBfsConfig cfg{};
        cfg.launch_time_seconds_since_j2000 = 0.0;
        cfg.launch_transfer_theta_global = 4.1;
        cfg.planet_sequence = {PlanetId::Earth, PlanetId::Mercury};
        cfg.max_total_time_seconds = global.constraints.max_total_time_seconds;
        cfg.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();
        cfg.problem1 = global.problem1;
        const auto r = search_fixed_sequence_bfs(cfg);
        std::cout << "  theta=4.1 found=" << r.found_solution
                  << " score=" << r.best_score << '\n';
    }

    std::cout << "\n=== free path Mercury seed, varying max_search_legs ===\n";
    const Leg0ThetaSeed mercury_seed{
        .first_leg_target_planet = PlanetId::Mercury,
        .transfer_theta_global = 4.1,
    };
    for (int legs : {1, 2, 3, 6}) {
        auto g = global;
        g.max_search_legs = legs;
        const auto r = search_free_path_from_leg0_seed(g, mercury_seed);
        std::cout << "  max_search_legs=" << legs
                  << " solutions=" << r.solutions.size()
                  << " expanded=" << r.stats.expanded_nodes
                  << " dead_p2=" << r.stats.dead_by_no_p2_solution
                  << " dead_turn=" << r.stats.dead_by_turn_angle << '\n';
        if (!r.solutions.empty()) {
            std::cout << "    seq: ";
            print_sequence(r.solutions.front().planet_sequence);
            std::cout << " score=" << r.solutions.front().score << '\n';
        }
    }

    std::cout << "\n=== free path Mercury seed legs=3, flyby filter Disabled ===\n";
    {
        auto g = global;
        g.max_search_legs = 3;
        g.constraints.flyby_physical_filter.mode =
            spaceship_cpp::trajectory::FlybyPhysicalFilterMode::Disabled;
        const auto r = search_free_path_from_leg0_seed(g, mercury_seed);
        std::cout << "  solutions=" << r.solutions.size()
                  << " expanded=" << r.stats.expanded_nodes
                  << " p2_solve=" << r.stats.p2_solve_calls << '\n';
    }

    return 0;
}
