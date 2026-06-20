/*
 * 文件作用：验证统一 G 搜索管线（端点 near-zero + 二次求根 + G-Newton）。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

namespace {

using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2FlybyGSearchConfig;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;

std::optional<Problem2OutgoingBranchSolution> find_solution_with_derivatives(
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

Problem2FlybyGSearchConfig make_test_G_search_config(
    const spaceship_cpp::problem2::Problem2ThetaPrimeScanConfig& scan_config
) {
    namespace config = spaceship_cpp::config;
    const auto& defaults = config::global_config();
    Problem2FlybyGSearchConfig search_config{};
    search_config.scan_config = scan_config;
    search_config.route_a_newton_options = config::make_problem2_route_a_newton_options(
        defaults.problem2_route_a_newton,
        defaults.problem1_solve);
    search_config.g_newton_options.max_iterations = 20;
    search_config.g_newton_options.G_tolerance = 1e-6;
    search_config.g_newton_options.theta_prime_tolerance = 1e-10;
    search_config.near_zero_G_threshold = 1e-3;
    search_config.solution_theta_prime_tolerance = 1e-7;
    return search_config;
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;
    namespace problem2 = spaceship_cpp::problem2;

    const auto& defaults = config::global_config();
    const auto scan_config = config::make_problem2_theta_prime_scan_config(
        PlanetId::Earth,
        PlanetId::Mars,
        0.0,
        defaults.problem2_theta_prime_scan,
        defaults.problem1_solve);
    const auto scan = problem2::run_problem2_theta_prime_initial_scan(scan_config);
    assert(scan.ok);

    std::size_t node_index = 0;
    const auto branch = find_solution_with_derivatives(scan, node_index);
    assert(branch.has_value());

    const double theta_prime_local = scan.nodes[node_index].theta_prime_local;
    const auto context = problem2::build_problem2_flyby_G_context(
        scan_config.flyby_planet,
        scan_config.flyby_time_seconds_since_j2000,
        branch->outgoing_eccentricity,
        theta_prime_local,
        true);
    assert(context.incoming_cache.valid);

    const auto G_at_branch = problem2::evaluate_flyby_constraint_G_at_branch(
        context.incoming_cache,
        *branch,
        theta_prime_local);
    assert(G_at_branch.has_value());
    assert(std::abs(*G_at_branch) < 1e-6);

    const auto search_config = make_test_G_search_config(scan_config);
    const auto search = problem2::search_flyby_constraint_G_zeros_from_initial_scan(
        context,
        search_config,
        scan);
    assert(search.ok);
    assert(!search.solutions.empty());

    bool found_near_zero = false;
    for (const auto& solution : search.solutions) {
        assert(std::abs(solution.G) <= search_config.near_zero_G_threshold ||
            std::abs(solution.G) <= search_config.g_newton_options.G_tolerance);
        if (std::abs(solution.theta_prime_local - theta_prime_local) <= 1e-2) {
            found_near_zero = true;
        }
    }
    assert(found_near_zero);

    if (node_index > 0U) {
        const auto& left_node = scan.nodes[node_index - 1];
        const auto& right_node = scan.nodes[node_index];
        const auto left_match = problem2::find_best_matching_outgoing_branch_index(
            *branch,
            left_node.solutions,
            scan_config.branch_phi_pairing_max_gap);
        const auto right_match = problem2::find_best_matching_outgoing_branch_index(
            *branch,
            right_node.solutions,
            scan_config.branch_phi_pairing_max_gap);
        if (left_match.has_value() && right_match.has_value()) {
            const problem2::Problem2FlybyGBranchIntervalInput interval{
                .theta_prime_left = left_node.theta_prime_local,
                .theta_prime_right = right_node.theta_prime_local,
                .left_branch = left_node.solutions[*left_match],
                .right_branch = right_node.solutions[*right_match],
            };
            const auto interval_solutions = problem2::process_flyby_constraint_G_branch_interval(
                context.incoming_cache,
                search_config,
                interval);
            assert(!interval_solutions.empty());
        }
    }

    {
        const auto newton = problem2::refine_flyby_constraint_G_zero_by_newton(
            context.incoming_cache,
            scan_config,
            search_config.route_a_newton_options,
            search_config.g_newton_options,
            theta_prime_local,
            *branch,
            *branch,
            theta_prime_local);
        assert(newton.ok);
        assert(std::abs(newton.final_G) <= search_config.g_newton_options.G_tolerance);
    }

    std::cout << "test_problem2_flyby_G_search PASSED\n";
    return 0;
}
