#include "spaceship_cpp/bfs/initial_launch_expansion.hpp"
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

}  // namespace

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "bfs_initial_launch_expansion_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto launch_planet = planet_params::PlanetId::Earth;
    const double launch_time = 0.17 * planet_params::planet_orbital_period(launch_planet);
    const double initial_theta = 0.4;
    const std::vector<planet_params::PlanetId> allowed_first_targets{
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };

    bfs::InitialLaunchExpansionOptions options{};
    options.max_transfer_revolution = 1;
    options.max_target_revolution = 1;
    options.max_launch_v_inf = std::numeric_limits<double>::infinity();
    options.time_weight_m_per_s_per_day = 0.0;

    const auto result = bfs::expand_initial_launch_with_problem1_table(
        loader, launch_planet, launch_time, initial_theta, allowed_first_targets, options);
    assert(result.ok);
    assert(result.problem1_query_count == static_cast<int>(allowed_first_targets.size()));
    assert(result.raw_branch_count >= result.accepted_candidate_count);
    assert(result.periapsis_raw_selected_count + result.periapsis_flipped_selected_count ==
           result.accepted_candidate_count);
    if (result.accepted_candidate_count > 0) {
        for (const auto& candidate : result.candidates) {
            assert(candidate.valid);
            assert(std::isfinite(candidate.launch_v_inf));
            assert(candidate.launch_v_inf >= 0.0);
            assert(candidate.next_state.depth == 1);
            assert(candidate.next_state.launch_v_inf == candidate.launch_v_inf);
            assert(candidate.next_state.current_planet == candidate.target_planet);
            assert(candidate.edge.arrival_time > candidate.edge.departure_time);
            assert(candidate.edge.transfer_time_seconds > 0.0);
            const double arrival_v_inf = bfs::compute_arrival_v_inf_from_state(candidate.next_state);
            assert(std::isfinite(arrival_v_inf));
            assert(arrival_v_inf >= 0.0);
        }
    }

    std::cout << "BfsInitialLaunchExpansionSummary\n";
    std::cout << "target_planet_count=" << result.target_planet_count << '\n';
    std::cout << "problem1_query_count=" << result.problem1_query_count << '\n';
    std::cout << "raw_branch_count=" << result.raw_branch_count << '\n';
    std::cout << "accepted_candidate_count=" << result.accepted_candidate_count << '\n';
    std::cout << "launch_v_inf_pruned_count=" << result.launch_v_inf_pruned_count << '\n';
    std::cout << "min_launch_v_inf=" << result.min_launch_v_inf << '\n';
    std::cout << "max_launch_v_inf=" << result.max_launch_v_inf << '\n';
    std::cout << "initial_launch_expansion_ok=1\n";
    std::cout << "BfsInitialLaunchPeriapsisResolutionSummary\n";
    std::cout << "raw_selected_count=" << result.periapsis_raw_selected_count << '\n';
    std::cout << "flipped_selected_count=" << result.periapsis_flipped_selected_count << '\n';
    std::cout << "max_transfer_p_resolution_error=" << result.max_transfer_p_resolution_error << '\n';
    std::cout << "resolution_ok=1\n";

    bool pruning_ok = true;
    double threshold = 0.0;
    bfs::InitialLaunchExpansionResult filtered{};
    if (!result.candidates.empty()) {
        std::vector<double> speeds;
        speeds.reserve(result.candidates.size());
        for (const auto& candidate : result.candidates) {
            speeds.push_back(candidate.launch_v_inf);
        }
        std::sort(speeds.begin(), speeds.end());
        threshold = speeds[speeds.size() / 2];
        auto prune_options = options;
        prune_options.max_launch_v_inf = threshold;
        filtered = bfs::expand_initial_launch_with_problem1_table(
            loader, launch_planet, launch_time, initial_theta, allowed_first_targets, prune_options);
        pruning_ok = filtered.ok &&
            filtered.accepted_candidate_count <= result.accepted_candidate_count &&
            filtered.launch_v_inf_pruned_count > 0;
        assert(pruning_ok);
    }

    std::cout << "BfsInitialLaunchVinfPruningSummary\n";
    std::cout << "unfiltered_count=" << result.accepted_candidate_count << '\n';
    std::cout << "threshold=" << threshold << '\n';
    std::cout << "filtered_count=" << filtered.accepted_candidate_count << '\n';
    std::cout << "pruned_count=" << filtered.launch_v_inf_pruned_count << '\n';
    std::cout << "vinf_pruning_ok=" << (pruning_ok ? 1 : 0) << '\n';
    return 0;
}
