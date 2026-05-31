#include "spaceship_cpp/bfs/gravity_assist_step_expansion.hpp"
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

int strict_count(const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& edges) {
    return static_cast<int>(std::count_if(edges.begin(), edges.end(), [](const auto& edge) {
        return edge.residual_source == spaceship_cpp::problem2::Problem2ResidualSource::Strict;
    }));
}

int relaxed_count(const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& edges) {
    return static_cast<int>(std::count_if(edges.begin(), edges.end(), [](const auto& edge) {
        return edge.residual_source == spaceship_cpp::problem2::Problem2ResidualSource::BoundaryAmbiguousRoundoff;
    }));
}

int topology_change_origin_count(const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& edges) {
    return static_cast<int>(std::count_if(edges.begin(), edges.end(), [](const auto& edge) {
        return edge.origin_was_topology_change;
    }));
}

double max_abs_slingshot_residual(const std::vector<spaceship_cpp::bfs::TrajectorySearchEdge>& edges) {
    double max_abs = 0.0;
    for (const auto& edge : edges) {
        max_abs = std::max(max_abs, std::abs(edge.slingshot_residual));
    }
    return max_abs;
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
        std::cout << "bfs_gravity_assist_step_expansion_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto earth = planet_params::PlanetId::Earth;
    const auto mars = planet_params::PlanetId::Mars;
    const double earth_period = planet_params::planet_orbital_period(earth);

    problem2::Problem2GravityAssistSolverOptions options{};
    options.theta_sample_count = 64;
    assert(options.topology_adaptive_enabled);

    bfs::TrajectorySearchState state0{};
    state0.valid = true;
    state0.current_planet = earth;
    state0.current_time = 0.17 * earth_period;
    state0.incoming_e = 0.3;
    state0.incoming_theta = 0.4;
    state0.depth = 0;
    state0.accumulated_time_seconds = 0.0;

    const auto expansion0 = bfs::expand_trajectory_state_by_gravity_assist(loader, state0, mars, options);
    assert(expansion0.ok);
    assert(expansion0.edges.size() == 16);
    assert(expansion0.next_states.size() == 16);

    double min_arrival_time = std::numeric_limits<double>::infinity();
    double max_arrival_time = 0.0;
    for (std::size_t i = 0; i < expansion0.edges.size(); ++i) {
        const auto& edge = expansion0.edges[i];
        const auto& next_state = expansion0.next_states[i];
        assert(edge.valid);
        assert(next_state.valid);
        assert(next_state.current_planet == mars);
        assert(next_state.current_time > state0.current_time);
        assert(next_state.depth == 1);
        min_arrival_time = std::min(min_arrival_time, edge.arrival_time);
        max_arrival_time = std::max(max_arrival_time, edge.arrival_time);
    }
    assert(max_abs_slingshot_residual(expansion0.edges) <= 1e-8);

    std::cout << "BfsGravityAssistStepExpansionSmokeSummary\n";
    std::cout << "scenario_id=0\n";
    std::cout << "edge_count=" << expansion0.edges.size() << '\n';
    std::cout << "next_state_count=" << expansion0.next_states.size() << '\n';
    std::cout << "strict_count=" << strict_count(expansion0.edges) << '\n';
    std::cout << "relaxed_boundary_count=" << relaxed_count(expansion0.edges) << '\n';
    std::cout << "topology_change_origin_count=" << topology_change_origin_count(expansion0.edges) << '\n';
    std::cout << "min_arrival_time=" << min_arrival_time << '\n';
    std::cout << "max_arrival_time=" << max_arrival_time << '\n';
    std::cout << "max_abs_slingshot_residual=" << max_abs_slingshot_residual(expansion0.edges) << '\n';
    std::cout << "smoke_ok=1\n";

    bfs::TrajectorySearchState state2{};
    state2.valid = true;
    state2.current_planet = earth;
    state2.current_time = 0.25 * earth_period;
    state2.incoming_e = 0.4;
    state2.incoming_theta = 1.0;
    state2.depth = 0;
    state2.accumulated_time_seconds = 0.0;

    const auto expansion2 = bfs::expand_trajectory_state_by_gravity_assist(loader, state2, mars, options);
    const int scenario2_topology_count = topology_change_origin_count(expansion2.edges);
    assert(expansion2.ok);
    assert(expansion2.edges.size() == 16);
    assert(expansion2.next_states.size() == 16);
    assert(scenario2_topology_count >= 4);
    assert(max_abs_slingshot_residual(expansion2.edges) <= 1e-8);

    std::cout << "BfsGravityAssistStepExpansionScenario2Summary\n";
    std::cout << "edge_count=" << expansion2.edges.size() << '\n';
    std::cout << "next_state_count=" << expansion2.next_states.size() << '\n';
    std::cout << "topology_change_origin_count=" << scenario2_topology_count << '\n';
    std::cout << "strict_count=" << strict_count(expansion2.edges) << '\n';
    std::cout << "relaxed_boundary_count=" << relaxed_count(expansion2.edges) << '\n';
    std::cout << "max_abs_slingshot_residual=" << max_abs_slingshot_residual(expansion2.edges) << '\n';
    std::cout << "scenario2_ok=1\n";

    return 0;
}
