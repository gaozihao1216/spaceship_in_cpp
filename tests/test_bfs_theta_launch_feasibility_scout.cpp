#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

void print_intervals(
    const char* label,
    const std::vector<spaceship_cpp::bfs::ThetaLaunchFeasibilityInterval>& intervals
) {
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        const auto& interval = intervals[i];
        std::cout << label << '\n';
        std::cout << "index=" << i << '\n';
        std::cout << "theta_left=" << interval.theta_left << '\n';
        std::cout << "theta_right=" << interval.theta_right << '\n';
        std::cout << "width=" << interval.theta_right - interval.theta_left << '\n';
        std::cout << "sample_count_inside=" << interval.sample_count_inside << '\n';
        std::cout << "min_v_inf_inside=" << interval.min_v_inf_inside << '\n';
        std::cout << "theta_at_min=" << interval.theta_at_min << '\n';
    }
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "bfs_theta_launch_feasibility_scout_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet);
    const std::vector<planet_params::PlanetId> allowed_first_targets{
        planet_params::PlanetId::Mercury,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mars,
    };

    bfs::ThetaLaunchFeasibilityScoutOptions dense_options{};
    dense_options.theta_scout_count = 720;
    dense_options.max_launch_v_inf = 7000.0;
    dense_options.near_v_inf_buffer = 1000.0;
    dense_options.max_transfer_revolution = 1;
    dense_options.max_target_revolution = 1;
    dense_options.refine_interval_boundaries = true;

    const auto dense = bfs::scout_theta_launch_feasibility_with_problem1_table(
        loader, launch_planet, launch_time, allowed_first_targets, dense_options);
    assert(dense.ok);
    assert(dense.samples.size() == 720);
    assert(std::isfinite(dense.global_min_launch_v_inf));
    assert(std::isfinite(dense.theta_at_global_min));
    for (const auto& interval : dense.valid_intervals) {
        assert(interval.valid);
        assert(interval.theta_right > interval.theta_left);
        assert(interval.min_v_inf_inside <= 7000.0);
    }
    for (const auto& interval : dense.near_valid_intervals) {
        assert(interval.valid);
        assert(interval.theta_right > interval.theta_left);
        assert(interval.min_v_inf_inside <= 8000.0);
    }

    std::cout << "BfsThetaLaunchFeasibilityScoutSummary\n";
    std::cout << "theta_scout_count=" << dense.theta_scout_count << '\n';
    std::cout << "valid_sample_count=" << dense.valid_sample_count << '\n';
    std::cout << "near_valid_sample_count=" << dense.near_valid_sample_count << '\n';
    std::cout << "no_candidate_sample_count=" << dense.no_candidate_sample_count << '\n';
    std::cout << "valid_interval_count=" << dense.valid_intervals.size() << '\n';
    std::cout << "near_valid_interval_count=" << dense.near_valid_intervals.size() << '\n';
    std::cout << "global_min_launch_v_inf=" << dense.global_min_launch_v_inf << '\n';
    std::cout << "theta_at_global_min=" << dense.theta_at_global_min << '\n';
    std::cout << "best_target_planet="
              << planet_params::planet_name(dense.best_sample.best_target_planet) << '\n';
    std::cout << "scout_ok=1\n";

    print_intervals("BfsThetaLaunchFeasibilityValidInterval", dense.valid_intervals);
    print_intervals("BfsThetaLaunchFeasibilityNearValidInterval", dense.near_valid_intervals);

    auto coarse_options = dense_options;
    coarse_options.theta_scout_count = 180;
    const auto coarse = bfs::scout_theta_launch_feasibility_with_problem1_table(
        loader, launch_planet, launch_time, allowed_first_targets, coarse_options);
    assert(coarse.ok);
    assert(coarse.samples.size() == 180);

    std::cout << "BfsThetaLaunchFeasibilityScoutDensityComparison\n";
    std::cout << "coarse_count=180\n";
    std::cout << "dense_count=720\n";
    std::cout << "coarse_global_min_launch_v_inf=" << coarse.global_min_launch_v_inf << '\n';
    std::cout << "dense_global_min_launch_v_inf=" << dense.global_min_launch_v_inf << '\n';
    std::cout << "coarse_theta_at_global_min=" << coarse.theta_at_global_min << '\n';
    std::cout << "dense_theta_at_global_min=" << dense.theta_at_global_min << '\n';
    std::cout << "coarse_valid_sample_count=" << coarse.valid_sample_count << '\n';
    std::cout << "dense_valid_sample_count=" << dense.valid_sample_count << '\n';
    std::cout << "density_comparison_ok=1\n";
    return 0;
}
