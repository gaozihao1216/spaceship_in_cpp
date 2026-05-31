#include "spaceship_cpp/bfs/theta_candidate_selector.hpp"
#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
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

double wrapped_distance(double lhs, double rhs) {
    return std::abs(spaceship_cpp::common::normalize_angle_minus_pi_pi(lhs - rhs));
}

bool theta_inside_interval(
    double theta,
    const spaceship_cpp::bfs::ThetaLaunchFeasibilityInterval& interval,
    double tolerance
) {
    double unwrapped = theta;
    if (unwrapped + tolerance < interval.theta_left) {
        unwrapped += spaceship_cpp::common::kTwoPi;
    }
    return unwrapped + tolerance >= interval.theta_left && unwrapped <= interval.theta_right + tolerance;
}

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "bfs_theta_candidate_selector_skipped_missing_table\n";
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

    bfs::ThetaLaunchFeasibilityScoutOptions scout_options{};
    scout_options.theta_scout_count = 720;
    scout_options.max_launch_v_inf = 7000.0;
    scout_options.near_v_inf_buffer = 1000.0;
    scout_options.max_transfer_revolution = 1;
    scout_options.max_target_revolution = 1;
    scout_options.refine_interval_boundaries = true;
    const auto scout = bfs::scout_theta_launch_feasibility_with_problem1_table(
        loader, launch_planet, launch_time, allowed_first_targets, scout_options);
    assert(scout.ok);
    assert(!scout.valid_intervals.empty());

    bfs::ThetaCandidateSelectorOptions selector_options{};
    selector_options.max_theta_candidates = 50;
    selector_options.min_samples_per_valid_interval = 5;
    selector_options.max_samples_per_valid_interval = 30;
    selector_options.include_near_valid_intervals = false;
    const auto selection = bfs::select_theta_candidates_from_launch_scout(scout, selector_options);
    assert(selection.ok);
    assert(selection.candidate_count > 0);
    assert(selection.candidate_count <= 50);

    double min_candidate_theta = std::numeric_limits<double>::infinity();
    double max_candidate_theta = 0.0;
    bool contains_global_min = false;
    const double scout_theta_step = spaceship_cpp::common::kTwoPi / static_cast<double>(scout.theta_scout_count);
    for (const auto& candidate : selection.candidates) {
        assert(candidate.valid);
        assert(std::isfinite(candidate.theta));
        assert(candidate.theta >= 0.0);
        assert(candidate.theta < spaceship_cpp::common::kTwoPi);
        min_candidate_theta = std::min(min_candidate_theta, candidate.theta);
        max_candidate_theta = std::max(max_candidate_theta, candidate.theta);
        if (wrapped_distance(candidate.theta, scout.theta_at_global_min) <= std::max(scout_theta_step, 1e-3)) {
            contains_global_min = true;
        }
        if (scout.valid_intervals.size() == 1) {
            assert(theta_inside_interval(candidate.theta, scout.valid_intervals.front(), 1e-4));
        }
    }
    assert(contains_global_min);

    const auto& best = selection.candidates.front();
    std::cout << "BfsThetaCandidateSelectorSummary\n";
    std::cout << "valid_interval_count=" << selection.valid_interval_count << '\n';
    std::cout << "near_valid_interval_count=" << selection.near_valid_interval_count << '\n';
    std::cout << "candidate_count=" << selection.candidate_count << '\n';
    std::cout << "min_candidate_theta=" << min_candidate_theta << '\n';
    std::cout << "max_candidate_theta=" << max_candidate_theta << '\n';
    std::cout << "best_candidate_theta=" << best.theta << '\n';
    std::cout << "best_candidate_estimated_v_inf=" << best.estimated_min_launch_v_inf << '\n';
    std::cout << "selector_ok=1\n";

    for (std::size_t i = 0; i < selection.candidates.size(); ++i) {
        const auto& candidate = selection.candidates[i];
        std::cout << "BfsThetaCandidate\n";
        std::cout << "index=" << i << '\n';
        std::cout << "theta=" << candidate.theta << '\n';
        std::cout << "estimated_min_launch_v_inf=" << candidate.estimated_min_launch_v_inf << '\n';
        std::cout << "from_valid_interval=" << (candidate.from_valid_interval ? 1 : 0) << '\n';
        std::cout << "from_near_valid_interval=" << (candidate.from_near_valid_interval ? 1 : 0) << '\n';
        std::cout << "is_interval_boundary=" << (candidate.is_interval_boundary ? 1 : 0) << '\n';
        std::cout << "is_interval_midpoint=" << (candidate.is_interval_midpoint ? 1 : 0) << '\n';
        std::cout << "is_interval_minimum=" << (candidate.is_interval_minimum ? 1 : 0) << '\n';
        std::cout << "interval_index=" << candidate.interval_index << '\n';
    }

    auto near_options = selector_options;
    near_options.include_near_valid_intervals = true;
    near_options.max_near_valid_extra_candidates = 10;
    near_options.max_theta_candidates = 60;
    const auto near_selection = bfs::select_theta_candidates_from_launch_scout(scout, near_options);
    assert(near_selection.ok);
    assert(near_selection.candidate_count <= 60);
    const int near_valid_candidate_count = static_cast<int>(std::count_if(
        near_selection.candidates.begin(), near_selection.candidates.end(), [](const auto& candidate) {
            return candidate.from_near_valid_interval && !candidate.from_valid_interval;
        }));

    std::cout << "BfsThetaCandidateSelectorNearValidSummary\n";
    std::cout << "candidate_count=" << near_selection.candidate_count << '\n';
    std::cout << "near_valid_candidate_count=" << near_valid_candidate_count << '\n';
    std::cout << "near_valid_selector_ok=1\n";
    return 0;
}
