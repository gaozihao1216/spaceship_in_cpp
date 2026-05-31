#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

struct Problem2SlingshotResidualBranchTestOnly {
    bool valid = false;
    std::string invalid_reason;
    double theta_prime = 0.0;
    double alpha = 0.0;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double problem1_residual_seconds = 0.0;
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double incoming_invariant = 0.0;
    double outgoing_invariant = 0.0;
    double slingshot_residual = 0.0;
};

struct TableBridgeResult {
    std::vector<Problem2SlingshotResidualBranchTestOnly> branches;
    bool fallback_used = false;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

bool residual_branch_less(
    const Problem2SlingshotResidualBranchTestOnly& lhs,
    const Problem2SlingshotResidualBranchTestOnly& rhs
) {
    return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.time_of_flight_seconds) <
           std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.time_of_flight_seconds);
}

std::vector<Problem2SlingshotResidualBranchTestOnly> convert_problem1_to_problem2_residuals(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_anomaly_phi,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem2 = spaceship_cpp::problem2;
    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);

    std::vector<Problem2SlingshotResidualBranchTestOnly> results;
    results.reserve(branches.size());
    for (const auto& branch : branches) {
        Problem2SlingshotResidualBranchTestOnly residual_branch{};
        residual_branch.theta_prime = theta_prime;
        residual_branch.alpha = branch.target_arrival_true_anomaly;
        residual_branch.transfer_revolution = branch.transfer_revolution;
        residual_branch.target_revolution = branch.target_revolution;
        residual_branch.time_of_flight_seconds = branch.time_of_flight_seconds;
        residual_branch.target_time_seconds = branch.target_time_seconds;
        residual_branch.problem1_residual_seconds = branch.residual_seconds;
        if (!branch.valid) {
            residual_branch.invalid_reason =
                branch.invalid_reason.empty() ? "problem1_branch_invalid" : branch.invalid_reason;
            results.push_back(residual_branch);
            continue;
        }
        const auto theta_alpha_residual = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            encounter_anomaly_phi,
            branch.target_arrival_true_anomaly,
            incoming_e,
            incoming_theta,
            theta_prime);
        if (!theta_alpha_residual.valid) {
            residual_branch.invalid_reason = theta_alpha_residual.invalid_reason;
            results.push_back(residual_branch);
            continue;
        }
        residual_branch.valid = true;
        residual_branch.outgoing_eccentricity = theta_alpha_residual.outgoing_eccentricity;
        residual_branch.outgoing_semi_latus_rectum = theta_alpha_residual.outgoing_semi_latus_rectum;
        residual_branch.incoming_invariant = theta_alpha_residual.incoming_invariant;
        residual_branch.outgoing_invariant = theta_alpha_residual.outgoing_invariant;
        residual_branch.slingshot_residual = theta_alpha_residual.slingshot_residual;
        results.push_back(residual_branch);
    }
    std::sort(results.begin(), results.end(), residual_branch_less);
    return results;
}

std::vector<Problem2SlingshotResidualBranchTestOnly>
evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    int max_transfer_revolution,
    int max_target_revolution
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        theta_A,
        max_transfer_revolution,
        max_target_revolution);
    return convert_problem1_to_problem2_residuals(
        departure_planet, target_planet, encounter_anomaly_phi, incoming_e, incoming_theta, theta_prime, branches);
}

TableBridgeResult evaluate_problem2_slingshot_residual_branches_for_theta_prime_nearest_node_table_test_only(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double encounter_anomaly_phi,
    double target_current_anomaly_beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    int max_transfer_revolution,
    int max_target_revolution
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto query = problem1::query_problem1_from_2deg_nearest_node(
        loader,
        departure_planet,
        target_planet,
        encounter_anomaly_phi,
        target_current_anomaly_beta,
        theta_A,
        max_transfer_revolution,
        max_target_revolution,
        options);
    TableBridgeResult result{};
    result.fallback_used = query.used_direct_solve_fallback;
    result.branches = convert_problem1_to_problem2_residuals(
        departure_planet, target_planet, encounter_anomaly_phi, incoming_e, incoming_theta, theta_prime, query.branches);
    return result;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_nearest_node_table_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    problem1::Problem1NearestNodeQueryOptions options{};
    options.residual_tolerance_seconds = 1e-2;
    options.max_newton_iterations = 80;
    options.fallback_direct_solve = true;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double incoming_e = 0.3;
    const double incoming_theta = 0.4;

    int total_direct = 0;
    int total_table = 0;
    int fallback_count = 0;
    double max_alpha_diff = 0.0;
    double max_time_diff = 0.0;
    double max_slingshot_diff = 0.0;
    bool ok = true;

    std::cout << std::setprecision(12) << std::scientific;
    for (const double theta_prime : std::array<double, 4>{0.1, 0.5, 1.0, 2.0}) {
        auto direct = evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
            departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
            theta_prime, 1, 1);
        auto table = evaluate_problem2_slingshot_residual_branches_for_theta_prime_nearest_node_table_test_only(
            loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
            theta_prime, 1, 1);
        total_direct += static_cast<int>(direct.size());
        total_table += static_cast<int>(table.branches.size());
        fallback_count += table.fallback_used ? 1 : 0;
        int valid_table = 0;
        const int pair_count = std::min<int>(direct.size(), table.branches.size());
        for (int i = 0; i < pair_count; ++i) {
            const auto& d = direct[static_cast<std::size_t>(i)];
            const auto& t = table.branches[static_cast<std::size_t>(i)];
            max_alpha_diff = std::max(max_alpha_diff, std::abs(normalize_angle_minus_pi_pi(d.alpha - t.alpha)));
            max_time_diff = std::max(max_time_diff, std::abs(d.time_of_flight_seconds - t.time_of_flight_seconds));
            if (d.valid && t.valid) {
                valid_table += 1;
                max_slingshot_diff =
                    std::max(max_slingshot_diff, std::abs(d.slingshot_residual - t.slingshot_residual));
            }
        }
        for (const auto& branch : table.branches) {
            if (branch.valid && std::isfinite(branch.slingshot_residual)) {
                valid_table += 0;
            } else if (branch.valid) {
                ok = false;
            }
        }
        if (table.branches.empty()) {
            ok = false;
        }
        std::cout << "Problem2ThetaPrimeNearestNodeResidualSummary\n";
        std::cout << "theta_prime=" << theta_prime << '\n';
        std::cout << "direct_branch_count=" << direct.size() << '\n';
        std::cout << "table_branch_count=" << table.branches.size() << '\n';
        std::cout << "table_fallback_used=" << (table.fallback_used ? 1 : 0) << '\n';
        std::cout << "valid_slingshot_branch_count=" << valid_table << '\n';
        std::cout << "max_alpha_diff=" << max_alpha_diff << '\n';
        std::cout << "max_time_diff_seconds=" << max_time_diff << '\n';
        std::cout << "max_slingshot_residual_diff=" << max_slingshot_diff << '\n';
    }

    std::cout << "Problem2ThetaPrimeNearestNodeTableSummary\n";
    std::cout << "theta_prime_sample_count=4\n";
    std::cout << "total_direct_branch_count=" << total_direct << '\n';
    std::cout << "total_table_branch_count=" << total_table << '\n';
    std::cout << "table_fallback_count=" << fallback_count << '\n';
    std::cout << "max_alpha_diff=" << max_alpha_diff << '\n';
    std::cout << "max_time_diff_seconds=" << max_time_diff << '\n';
    std::cout << "max_slingshot_residual_diff=" << max_slingshot_diff << '\n';
    std::cout << "problem2_theta_prime_table_bridge_ok=" << (ok ? 1 : 0) << '\n';
    return ok ? 0 : 1;
}
