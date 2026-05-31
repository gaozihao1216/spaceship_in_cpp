#include "spaceship_cpp/bfs/fixed_sequence_theta_refinement.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace {

constexpr double kBaselineTotalDeltaV = 15183.26387674;

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

std::string sequence_string(const std::vector<spaceship_cpp::planet_params::PlanetId>& sequence) {
    namespace planet_params = spaceship_cpp::planet_params;
    std::ostringstream os;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        if (i > 0) {
            os << "->";
        }
        os << planet_params::planet_name(sequence[i]);
    }
    return os.str();
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "bfs_fixed_sequence_theta_refinement_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet);
    const std::vector<planet_params::PlanetId> sequence{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };

    bfs::InitialLaunchExpansionOptions initial_options{};
    initial_options.max_transfer_revolution = 1;
    initial_options.max_target_revolution = 1;
    initial_options.max_launch_v_inf = 7000.0;
    initial_options.time_weight_m_per_s_per_day = 0.0;

    problem2::Problem2GravityAssistSolverOptions problem2_options{};
    problem2_options.theta_sample_count = 64;
    problem2_options.topology_adaptive_enabled = true;

    bfs::FixedSequenceThetaRefinementOptions refinement_options{};
    refinement_options.theta_center = 2.552022761644;
    refinement_options.initial_half_width = 0.05;
    refinement_options.samples_per_round = 31;
    refinement_options.refinement_rounds = 3;
    refinement_options.shrink_factor = 0.35;
    refinement_options.include_center = true;
    refinement_options.beam_width = 10;
    refinement_options.max_launch_v_inf = 7000.0;
    refinement_options.time_weight_m_per_s_per_day = 0.0;

    const auto result = bfs::refine_fixed_sequence_theta_with_table(
        loader,
        launch_time,
        sequence,
        initial_options,
        problem2_options,
        refinement_options);
    assert(result.ok);
    assert(!result.round_summaries.empty());
    if (result.best_valid) {
        assert(std::isfinite(result.best_total_delta_v));
        assert(result.best_launch_v_inf <= 7000.0);
        assert(std::isfinite(result.best_arrival_v_inf));
        assert(result.best_solution.sequence == sequence);
    }

    std::cout << "BfsFixedSequenceThetaRefinementSummary\n";
    std::cout << "sequence=" << sequence_string(sequence) << '\n';
    std::cout << "round_count=" << result.round_summaries.size() << '\n';
    std::cout << "best_valid=" << (result.best_valid ? 1 : 0) << '\n';
    std::cout << "best_theta=" << result.best_theta << '\n';
    std::cout << "best_score=" << result.best_score << '\n';
    std::cout << "best_total_delta_v=" << result.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << result.best_total_flight_time_seconds << '\n';
    std::cout << "baseline_total_delta_v=" << kBaselineTotalDeltaV << '\n';
    std::cout << "improvement_vs_baseline="
              << (result.best_valid ? kBaselineTotalDeltaV - result.best_total_delta_v
                                    : -std::numeric_limits<double>::infinity())
              << '\n';
    std::cout << "refinement_ok=1\n";

    for (const auto& summary : result.round_summaries) {
        std::cout << "BfsFixedSequenceThetaRefinementRoundSummary\n";
        std::cout << "round_index=" << summary.round_index << '\n';
        std::cout << "theta_left=" << summary.theta_left << '\n';
        std::cout << "theta_right=" << summary.theta_right << '\n';
        std::cout << "theta_center=" << summary.theta_center << '\n';
        std::cout << "half_width=" << summary.half_width << '\n';
        std::cout << "theta_sample_count=" << summary.theta_sample_count << '\n';
        std::cout << "terminal_theta_count=" << summary.terminal_theta_count << '\n';
        std::cout << "best_valid=" << (summary.best_valid ? 1 : 0) << '\n';
        std::cout << "best_theta=" << summary.best_theta << '\n';
        std::cout << "best_score=" << summary.best_score << '\n';
        std::cout << "best_total_delta_v=" << summary.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << summary.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << summary.best_arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds=" << summary.best_total_flight_time_seconds << '\n';
        std::cout << "wall_time_ms=" << summary.wall_time_ms << '\n';
    }

    return 0;
}
