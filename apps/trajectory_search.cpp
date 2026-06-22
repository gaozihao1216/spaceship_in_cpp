/*
 * 文件作用：Earth → Mercury 轨迹搜索命令行入口。
 * 主要工作：调用 search_best_trajectory，并列出 Step 2 中所有到达 Mercury 的候选序列。
 */
#include "spaceship_cpp/bfs/bfs.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

constexpr double kDayInSeconds = 86400.0;

using spaceship_cpp::bfs::default_trajectory_search_global_config;
using spaceship_cpp::bfs::discretize_leg0_theta_seeds;
using spaceship_cpp::bfs::find_leg0_feasible_theta_for_first_leg_targets;
using spaceship_cpp::bfs::run_step2_free_path_search;
using spaceship_cpp::bfs::search_best_trajectory;
using spaceship_cpp::bfs::select_top_k_planet_sequences_reaching;
using spaceship_cpp::bfs::RankedPlanetSequence;
using spaceship_cpp::bfs::TrajectorySearchInput;
using spaceship_cpp::bfs::TransferLegDescriptor;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::planet_params::planet_name;

using Clock = std::chrono::steady_clock;

void print_leg(const TransferLegDescriptor& leg, int index) {
    std::cout << std::fixed;
    std::cout << "  leg " << index << ": "
              << planet_name(leg.from_planet) << " -> " << planet_name(leg.to_planet) << '\n';
    std::cout << std::setprecision(1)
              << "    t_dep=" << leg.departure_time_seconds_since_j2000 / kDayInSeconds << " d"
              << ", t_arr=" << leg.arrival_time_seconds_since_j2000 / kDayInSeconds << " d"
              << ", tof=" << leg.time_of_flight_seconds / kDayInSeconds << " d\n";
    std::cout << std::setprecision(6)
              << "    e=" << leg.eccentricity
              << ", p=" << leg.semi_latus_rectum_au << " AU"
              << ", theta=" << leg.perihelion_angle_global_rad << " rad"
              << ", k=" << leg.transfer_revolution
              << ", q=" << leg.target_revolution << '\n';
}

void print_sequence(const std::vector<PlanetId>& sequence) {
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        if (i > 0) {
            std::cout << " -> ";
        }
        std::cout << planet_name(sequence[i]);
    }
}

}  // namespace

int main() {
    TrajectorySearchInput input{};
    input.launch_time_seconds_since_j2000 = 0.0;
    input.departure_planet = PlanetId::Earth;
    input.destination_planet = PlanetId::Mercury;

    auto config = default_trajectory_search_global_config();
    config.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    std::cout << "=== Earth -> Mercury trajectory search ===\n";
    std::cout << "launch_time=" << input.launch_time_seconds_since_j2000 << " s (J2000)"
              << ", max_search_legs=" << config.max_search_legs
              << ", top_k=" << config.discretization.top_k_sequences << '\n';

    const Clock::time_point t0 = Clock::now();
    const auto leg0 = find_leg0_feasible_theta_for_first_leg_targets(config);
    const auto seeds = discretize_leg0_theta_seeds(config, leg0);
    const auto step2 = run_step2_free_path_search(config, seeds);
    const Clock::time_point t1 = Clock::now();

    if (!step2.ok) {
        std::cerr << "Step 2 failed: " << step2.error_message << '\n';
        return 1;
    }

    constexpr int kListTop = 10;
    const auto mercury_candidates = select_top_k_planet_sequences_reaching(
        step2, kListTop, PlanetId::Mercury);

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "\n--- Step 2 candidates reaching Mercury (coarse, top "
              << kListTop << ", " << std::chrono::duration<double>(t1 - t0).count() << " s) ---\n";
    std::cout << "  unique_sequences=" << mercury_candidates.stats.unique_sequences << '\n';
    for (std::size_t i = 0; i < mercury_candidates.sequences.size(); ++i) {
        const auto& seq = mercury_candidates.sequences[i];
        std::cout << "  [" << i << "] score=" << seq.best_score
                  << ", legs=" << (seq.planet_sequence.size() - 1) << ", ";
        print_sequence(seq.planet_sequence);
        std::cout << '\n';
    }

    const Clock::time_point t2 = Clock::now();
    const auto result = search_best_trajectory(input, config);
    const Clock::time_point t3 = Clock::now();

    std::cout << "\n--- Best after Step 4 fine search ("
              << std::chrono::duration<double>(t3 - t2).count() << " s) ---\n";

    if (!result.ok) {
        std::cerr << "search failed: " << result.error_message << '\n';
        return 1;
    }
    if (!result.found_solution) {
        std::cerr << "no solution: " << result.error_message << '\n';
        return 1;
    }

    std::cout << "visit_sequence: ";
    print_sequence(result.visit_sequence);
    std::cout << "\nscore=" << result.score << " m/s"
              << ", legs=" << result.legs.size()
              << ", v_inf_launch=" << result.launch_v_inf_mps
              << ", v_inf_arrival=" << result.arrival_v_inf_mps
              << ", total_tof=" << result.total_time_seconds / kDayInSeconds << " d\n";
    std::cout << std::setprecision(4)
              << "leg0_theta=" << result.leg0_theta_global_rad << " rad\n";
    if (!result.error_message.empty()) {
        std::cout << "note: " << result.error_message << '\n';
    }

    std::cout << "\ntransfer legs:\n";
    for (std::size_t i = 0; i < result.legs.size(); ++i) {
        print_leg(result.legs[i], static_cast<int>(i) + 1);
    }

    if (!mercury_candidates.sequences.empty()) {
        const auto& best_coarse = mercury_candidates.sequences.front();
        const RankedPlanetSequence* direct_coarse = nullptr;
        for (const auto& seq : mercury_candidates.sequences) {
            if (seq.planet_sequence.size() == 2U) {
                direct_coarse = &seq;
                break;
            }
        }

        std::cout << std::setprecision(1)
                  << "\n--- vs direct Earth -> Mercury (coarse Step 2) ---\n";
        if (direct_coarse != nullptr) {
            std::cout << "  direct score=" << direct_coarse->best_score << '\n';
        } else {
            std::cout << "  direct path not in top coarse list\n";
        }
        std::cout << "  best coarse reaching Mercury score=" << best_coarse.best_score
                  << ", legs=" << (best_coarse.planet_sequence.size() - 1) << ", ";
        print_sequence(best_coarse.planet_sequence);
        std::cout << "\n  refined best score=" << result.score << '\n';
        if (direct_coarse != nullptr && best_coarse.best_score < direct_coarse->best_score) {
            std::cout << "  multi-leg beats direct by "
                      << (direct_coarse->best_score - best_coarse.best_score) << " m/s (coarse)\n";
        } else if (direct_coarse != nullptr) {
            std::cout << "  direct still best or equal at coarse stage\n";
        }
    }

    return 0;
}
