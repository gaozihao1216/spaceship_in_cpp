#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_transition.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

int strict_count(const spaceship_cpp::problem2::Problem2GravityAssistTransitionResult& result) {
    return static_cast<int>(std::count_if(result.candidates.begin(), result.candidates.end(), [](const auto& candidate) {
        return candidate.residual_source == spaceship_cpp::problem2::Problem2ResidualSource::Strict;
    }));
}

int relaxed_count(const spaceship_cpp::problem2::Problem2GravityAssistTransitionResult& result) {
    return static_cast<int>(std::count_if(result.candidates.begin(), result.candidates.end(), [](const auto& candidate) {
        return candidate.residual_source ==
               spaceship_cpp::problem2::Problem2ResidualSource::BoundaryAmbiguousRoundoff;
    }));
}

int topology_change_origin_count(const spaceship_cpp::problem2::Problem2GravityAssistTransitionResult& result) {
    return static_cast<int>(std::count_if(result.candidates.begin(), result.candidates.end(), [](const auto& candidate) {
        return candidate.origin_was_topology_change;
    }));
}

double max_abs_slingshot_residual(const spaceship_cpp::problem2::Problem2GravityAssistTransitionResult& result) {
    double max_abs = 0.0;
    for (const auto& candidate : result.candidates) {
        max_abs = std::max(max_abs, std::abs(candidate.slingshot_residual));
    }
    return max_abs;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_gravity_assist_transition_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto encounter_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(encounter_planet);

    problem2::Problem2GravityAssistSolverOptions options{};
    options.theta_sample_count = 64;
    assert(options.topology_adaptive_enabled);

    const double scenario0_encounter_time = 0.17 * earth_period;
    const auto scenario0 = problem2::generate_problem2_gravity_assist_transitions(
        loader, encounter_planet, target_planet, scenario0_encounter_time, 0.3, 0.4, options);
    assert(scenario0.ok);
    assert(scenario0.candidates.size() == 16);
    assert(scenario0.solver_summary.dedup_success_count == 16);
    assert(scenario0.solver_summary.table_fallback_count == 0);
    assert(max_abs_slingshot_residual(scenario0) <= 1e-8);

    double min_arrival_time = std::numeric_limits<double>::infinity();
    double max_arrival_time = 0.0;
    for (const auto& candidate : scenario0.candidates) {
        assert(candidate.valid);
        assert(candidate.arrival_time > candidate.encounter_time);
        assert(candidate.outgoing_semi_latus_rectum > 0.0);
        min_arrival_time = std::min(min_arrival_time, candidate.arrival_time);
        max_arrival_time = std::max(max_arrival_time, candidate.arrival_time);
    }

    std::cout << "Problem2GravityAssistTransitionTestSummary\n";
    std::cout << "scenario_id=0\n";
    std::cout << "candidate_count=" << scenario0.candidates.size() << '\n';
    std::cout << "strict_count=" << strict_count(scenario0) << '\n';
    std::cout << "relaxed_boundary_count=" << relaxed_count(scenario0) << '\n';
    std::cout << "topology_change_origin_count=" << topology_change_origin_count(scenario0) << '\n';
    std::cout << "min_arrival_time=" << min_arrival_time << '\n';
    std::cout << "max_arrival_time=" << max_arrival_time << '\n';
    std::cout << "max_abs_slingshot_residual=" << max_abs_slingshot_residual(scenario0) << '\n';
    std::cout << "max_abs_problem1_residual_seconds="
              << scenario0.solver_summary.max_abs_problem1_residual_seconds_at_root << '\n';
    std::cout << "transition_test_ok=1\n";

    const double scenario2_encounter_time = 0.25 * earth_period;
    const auto scenario2 = problem2::generate_problem2_gravity_assist_transitions(
        loader, encounter_planet, target_planet, scenario2_encounter_time, 0.4, 1.0, options);
    const int scenario2_topology_count = topology_change_origin_count(scenario2);
    assert(scenario2.ok);
    assert(scenario2.candidates.size() == 16);
    assert(scenario2_topology_count >= 4);
    assert(max_abs_slingshot_residual(scenario2) <= 1e-8);

    std::cout << "Problem2GravityAssistTransitionScenario2Summary\n";
    std::cout << "candidate_count=" << scenario2.candidates.size() << '\n';
    std::cout << "topology_change_origin_count=" << scenario2_topology_count << '\n';
    std::cout << "strict_count=" << strict_count(scenario2) << '\n';
    std::cout << "relaxed_boundary_count=" << relaxed_count(scenario2) << '\n';
    std::cout << "max_abs_slingshot_residual=" << max_abs_slingshot_residual(scenario2) << '\n';
    std::cout << "scenario2_transition_ok=1\n";

    return 0;
}
