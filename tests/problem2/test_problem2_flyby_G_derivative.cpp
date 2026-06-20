/*
 * 文件作用：验证 dG/dθ' 解析链式求导与数值差分一致。
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

double central_G_derivative_along_branch(
    const spaceship_cpp::problem2::Problem2ThetaPrimeScanConfig& scan_config,
    const spaceship_cpp::problem2::FlybyConstraintIncomingCache& incoming_cache,
    const Problem2OutgoingBranchSolution& branch,
    double theta_prime_local,
    double step
) {
    namespace problem2 = spaceship_cpp::problem2;
    const auto plus_node = problem2::evaluate_problem2_theta_prime_node(scan_config, theta_prime_local + step);
    const auto minus_node = problem2::evaluate_problem2_theta_prime_node(scan_config, theta_prime_local - step);
    const auto plus_match = problem2::find_best_matching_outgoing_branch_index(
        branch,
        plus_node.solutions,
        scan_config.branch_phi_pairing_max_gap);
    const auto minus_match = problem2::find_best_matching_outgoing_branch_index(
        branch,
        minus_node.solutions,
        scan_config.branch_phi_pairing_max_gap);
    assert(plus_match.has_value() && minus_match.has_value());

    const auto G_plus = problem2::evaluate_flyby_constraint_G_at_branch(
        incoming_cache,
        plus_node.solutions[*plus_match],
        theta_prime_local + step);
    const auto G_minus = problem2::evaluate_flyby_constraint_G_at_branch(
        incoming_cache,
        minus_node.solutions[*minus_match],
        theta_prime_local - step);
    assert(G_plus.has_value() && G_minus.has_value());
    return (*G_plus - *G_minus) / (2.0 * step);
}

bool approx_equal(double a, double b, double rtol, double atol) {
    return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
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

    const auto analytical = problem2::evaluate_flyby_constraint_G_and_derivative(
        context.incoming_cache,
        branch->outgoing_eccentricity,
        theta_prime_local,
        branch->dphi_dtheta_prime,
        branch->de_dtheta_prime);
    assert(analytical.valid);

    const double coarse_step =
        scan.nodes[node_index + 1].theta_prime_local - scan.nodes[node_index].theta_prime_local;
    const double step = coarse_step / 80.0;
    const double numeric = central_G_derivative_along_branch(
        scan_config,
        context.incoming_cache,
        *branch,
        theta_prime_local,
        step);

    assert(approx_equal(analytical.dG_dtheta_prime, numeric, 5e-2, 1e-3));
    assert(std::abs(analytical.G) < 1e-6);

    std::cout << "test_problem2_flyby_G_derivative PASSED\n";
    return 0;
}
