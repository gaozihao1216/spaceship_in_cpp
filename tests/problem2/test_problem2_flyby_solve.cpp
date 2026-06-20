/*
 * 文件作用：测试 Problem 2 飞掠顶层求解封装。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

namespace {

using spaceship_cpp::planet_params::PlanetId;

std::optional<spaceship_cpp::problem2::Problem2OutgoingBranchSolution>
find_solution_with_derivatives(
    const spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult& scan,
    std::size_t& out_node_index
) {
    for (std::size_t node_index = 1; node_index + 1 < scan.nodes.size(); ++node_index) {
        for (const auto& solution : scan.nodes[node_index].solutions) {
            if (solution.has_dphi_dtheta_prime && solution.has_de_dtheta_prime) {
                out_node_index = node_index;
                return solution;
            }
        }
    }
    return std::nullopt;
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;
    namespace problem2 = spaceship_cpp::problem2;

    const auto& defaults = config::global_config();
    const auto solve_config = config::make_problem2_flyby_solve_config(
        PlanetId::Earth,
        PlanetId::Mars,
        0.0,
        defaults.problem2_theta_prime_scan,
        defaults.problem2_route_a_newton,
        defaults.problem1_solve);

    problem2::Problem2FlybySolveInput input{};
    input.flyby_planet = PlanetId::Earth;
    input.target_planet = PlanetId::Mars;
    input.flyby_time_seconds_since_j2000 = 0.0;

    const auto scan = problem2::run_problem2_flyby_theta_prime_initial_scan(input, solve_config);
    assert(scan.ok);

    std::size_t node_index = 0;
    const auto branch = find_solution_with_derivatives(scan, node_index);
    assert(branch.has_value());

    input.incoming_eccentricity = branch->outgoing_eccentricity;
    input.incoming_theta = scan.nodes[node_index].theta_prime_local;
    input.incoming_theta_is_local = true;

    const auto solved = problem2::solve_problem2_flyby_with_scan(input, solve_config, scan);
    assert(solved.ok);
    assert(!solved.solutions.empty());

    bool found_near_incoming_theta = false;
    for (const auto& solution : solved.solutions) {
        assert(std::isfinite(solution.outgoing_eccentricity));
        assert(std::isfinite(solution.outgoing_theta_prime_local));
        assert(std::isfinite(solution.outgoing_theta_prime_global));
        assert(std::isfinite(solution.encounter_global_angle));
        assert(std::isfinite(solution.time_of_flight_seconds));
        assert(solution.time_of_flight_seconds > 0.0);
        assert(std::isfinite(solution.arrival_time_seconds_since_j2000));
        assert(solution.arrival_time_seconds_since_j2000 ==
            solution.flyby_time_seconds_since_j2000 + solution.time_of_flight_seconds);
        assert(std::abs(solution.flyby_constraint_G) <=
            solve_config.g_search_config.near_zero_G_threshold + 1e-6);

        if (std::abs(solution.outgoing_theta_prime_local - input.incoming_theta) <= 1e-2) {
            found_near_incoming_theta = true;
        }
    }
    assert(found_near_incoming_theta);

    std::cout << "test_problem2_flyby_solve PASSED\n";
    return 0;
}
