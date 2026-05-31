#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_expansion.hpp"

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

int strict_count(const std::vector<spaceship_cpp::problem2::Problem2GravityAssistExpandedState>& states) {
    return static_cast<int>(std::count_if(states.begin(), states.end(), [](const auto& state) {
        return state.residual_source == spaceship_cpp::problem2::Problem2ResidualSource::Strict;
    }));
}

int relaxed_count(const std::vector<spaceship_cpp::problem2::Problem2GravityAssistExpandedState>& states) {
    return static_cast<int>(std::count_if(states.begin(), states.end(), [](const auto& state) {
        return state.residual_source == spaceship_cpp::problem2::Problem2ResidualSource::BoundaryAmbiguousRoundoff;
    }));
}

int topology_change_origin_count(
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistExpandedState>& states
) {
    return static_cast<int>(std::count_if(states.begin(), states.end(), [](const auto& state) {
        return state.origin_was_topology_change;
    }));
}

double max_abs_slingshot_residual(
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistExpandedState>& states
) {
    double max_abs = 0.0;
    for (const auto& state : states) {
        max_abs = std::max(max_abs, std::abs(state.slingshot_residual));
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
        std::cout << "problem2_gravity_assist_expansion_nearest_node_table_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto earth = planet_params::PlanetId::Earth;
    const auto mars = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(earth);

    problem2::Problem2GravityAssistSolverOptions options{};
    options.theta_sample_count = 64;
    assert(options.topology_adaptive_enabled);

    problem2::Problem2GravityAssistExpansionInput scenario0{};
    scenario0.encounter_planet = earth;
    scenario0.target_planet = mars;
    scenario0.encounter_time = 0.17 * earth_period;
    scenario0.incoming_e = 0.3;
    scenario0.incoming_theta = 0.4;

    const auto states0 = problem2::expand_one_gravity_assist_step_with_table(loader, scenario0, options);
    assert(states0.size() == 16);

    double min_arrival_time = std::numeric_limits<double>::infinity();
    double max_arrival_time = 0.0;
    for (const auto& state : states0) {
        assert(state.valid);
        assert(state.next_planet == mars);
        assert(state.arrival_time > state.departure_time);
        assert(state.transfer_time_seconds > 0.0);
        assert(state.outgoing_p > 0.0);
        min_arrival_time = std::min(min_arrival_time, state.arrival_time);
        max_arrival_time = std::max(max_arrival_time, state.arrival_time);
    }
    assert(max_abs_slingshot_residual(states0) <= 1e-8);

    std::cout << "Problem2GravityAssistExpansionSmokeSummary\n";
    std::cout << "scenario_id=0\n";
    std::cout << "expanded_state_count=" << states0.size() << '\n';
    std::cout << "strict_count=" << strict_count(states0) << '\n';
    std::cout << "relaxed_boundary_count=" << relaxed_count(states0) << '\n';
    std::cout << "topology_change_origin_count=" << topology_change_origin_count(states0) << '\n';
    std::cout << "min_arrival_time=" << min_arrival_time << '\n';
    std::cout << "max_arrival_time=" << max_arrival_time << '\n';
    std::cout << "max_abs_slingshot_residual=" << max_abs_slingshot_residual(states0) << '\n';
    std::cout << "expansion_smoke_ok=1\n";

    problem2::Problem2GravityAssistExpansionInput scenario2{};
    scenario2.encounter_planet = earth;
    scenario2.target_planet = mars;
    scenario2.encounter_time = 0.25 * earth_period;
    scenario2.incoming_e = 0.4;
    scenario2.incoming_theta = 1.0;

    const auto states2 = problem2::expand_one_gravity_assist_step_with_table(loader, scenario2, options);
    const int scenario2_topology_count = topology_change_origin_count(states2);
    assert(states2.size() == 16);
    assert(scenario2_topology_count >= 4);
    assert(max_abs_slingshot_residual(states2) <= 1e-8);

    std::cout << "Problem2GravityAssistExpansionScenario2Summary\n";
    std::cout << "expanded_state_count=" << states2.size() << '\n';
    std::cout << "topology_change_origin_count=" << scenario2_topology_count << '\n';
    std::cout << "strict_count=" << strict_count(states2) << '\n';
    std::cout << "relaxed_boundary_count=" << relaxed_count(states2) << '\n';
    std::cout << "max_abs_slingshot_residual=" << max_abs_slingshot_residual(states2) << '\n';
    std::cout << "scenario2_expansion_ok=1\n";

    return 0;
}
