#include "spaceship_cpp/bfs/fixed_sequence_launch_time_scan.hpp"
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

constexpr double kBaselineTotalDeltaV = 14984.43332290;

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
        std::cout << "bfs_fixed_sequence_launch_time_scan_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time_center = 0.17 * planet_params::planet_orbital_period(launch_planet);
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

    bfs::FixedSequenceLaunchTimeScanOptions scan_options{};
    scan_options.launch_time_center = launch_time_center;
    scan_options.initial_half_width_seconds = 0.01 * planet_params::planet_orbital_period(launch_planet);
    scan_options.launch_time_sample_count = 5;
    scan_options.refine_theta_at_each_time = true;
    scan_options.theta_center = 2.579506094977;
    scan_options.theta_half_width = 0.02;
    scan_options.theta_samples_per_round = 15;
    scan_options.theta_refinement_rounds = 2;
    scan_options.theta_shrink_factor = 0.4;
    scan_options.beam_width = 10;
    scan_options.max_launch_v_inf = 7000.0;
    scan_options.time_weight_m_per_s_per_day = 0.0;

    const auto result = bfs::scan_fixed_sequence_launch_time_locally_with_table(
        loader,
        sequence,
        initial_options,
        problem2_options,
        scan_options);
    assert(result.ok);
    assert(static_cast<int>(result.samples.size()) == scan_options.launch_time_sample_count);
    if (result.best_valid) {
        assert(std::isfinite(result.best_total_delta_v));
        assert(result.best_launch_v_inf <= 7000.0);
        assert(std::isfinite(result.best_arrival_v_inf));
    }

    std::cout << "BfsFixedSequenceLaunchTimeScanSummary\n";
    std::cout << "sequence=" << sequence_string(sequence) << '\n';
    std::cout << "sample_count=" << result.samples.size() << '\n';
    std::cout << "best_valid=" << (result.best_valid ? 1 : 0) << '\n';
    std::cout << "best_launch_time=" << result.best_launch_time << '\n';
    std::cout << "best_launch_time_offset_seconds="
              << (result.best_valid ? result.best_launch_time - scan_options.launch_time_center
                                    : std::numeric_limits<double>::infinity())
              << '\n';
    std::cout << "best_theta=" << result.best_theta << '\n';
    std::cout << "best_total_delta_v=" << result.best_total_delta_v << '\n';
    std::cout << "best_launch_v_inf=" << result.best_launch_v_inf << '\n';
    std::cout << "best_arrival_v_inf=" << result.best_arrival_v_inf << '\n';
    std::cout << "best_total_flight_time_seconds=" << result.best_total_flight_time_seconds << '\n';
    std::cout << "baseline_total_delta_v=" << kBaselineTotalDeltaV << '\n';
    std::cout << "improvement_vs_baseline="
              << (result.best_valid ? kBaselineTotalDeltaV - result.best_total_delta_v
                                    : -std::numeric_limits<double>::infinity())
              << '\n';
    std::cout << "launch_time_scan_ok=1\n";

    for (std::size_t i = 0; i < result.samples.size(); ++i) {
        const auto& sample = result.samples[i];
        std::cout << "BfsFixedSequenceLaunchTimeScanSample\n";
        std::cout << "index=" << i << '\n';
        std::cout << "launch_time=" << sample.launch_time << '\n';
        std::cout << "launch_time_offset_seconds=" << sample.launch_time_offset_seconds << '\n';
        std::cout << "best_valid=" << (sample.best_valid ? 1 : 0) << '\n';
        std::cout << "best_theta=" << sample.best_theta << '\n';
        std::cout << "best_total_delta_v=" << sample.best_total_delta_v << '\n';
        std::cout << "best_launch_v_inf=" << sample.best_launch_v_inf << '\n';
        std::cout << "best_arrival_v_inf=" << sample.best_arrival_v_inf << '\n';
        std::cout << "best_total_flight_time_seconds=" << sample.best_total_flight_time_seconds << '\n';
        std::cout << "wall_time_ms=" << sample.wall_time_ms << '\n';
    }

    return 0;
}
