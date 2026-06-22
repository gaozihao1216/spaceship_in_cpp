/*
 * 文件作用：测试 search_best_trajectory 公共 API（Earth -> Mercury）。
 */
#include "spaceship_cpp/bfs/bfs.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <cassert>
#include <limits>

int main() {
    namespace bfs = spaceship_cpp::bfs;
    namespace planet_params = spaceship_cpp::planet_params;

    bfs::TrajectorySearchInput input{};
    input.launch_time_seconds_since_j2000 = 0.0;
    input.departure_planet = planet_params::PlanetId::Earth;
    input.destination_planet = planet_params::PlanetId::Mercury;

    auto config = bfs::default_trajectory_search_global_config();
    config.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    const bfs::TrajectorySearchOutput result = bfs::search_best_trajectory(input, config);

    assert(result.ok);
    assert(result.found_solution);
    assert(result.departure_planet == planet_params::PlanetId::Earth);
    assert(result.destination_planet == planet_params::PlanetId::Mercury);
    assert(!result.visit_sequence.empty());
    assert(result.visit_sequence.front() == planet_params::PlanetId::Earth);
    assert(result.visit_sequence.back() == planet_params::PlanetId::Mercury);
    assert(!result.legs.empty());
    assert(result.legs.size() + 1 == result.visit_sequence.size());
    assert(result.legs.front().from_planet == planet_params::PlanetId::Earth);
    assert(result.legs.back().to_planet == planet_params::PlanetId::Mercury);
    assert(result.legs.front().eccentricity > 0.0);
    assert(result.legs.front().semi_latus_rectum_au > 0.0);
    assert(result.score > 0.0);

    return 0;
}
