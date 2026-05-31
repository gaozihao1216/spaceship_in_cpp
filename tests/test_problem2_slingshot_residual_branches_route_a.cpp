#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

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
    namespace problem2 = spaceship_cpp::problem2;

    (void)encounter_time;

    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);

    const std::vector<problem1::Problem1SolutionBranch> branches =
        problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            encounter_anomaly_phi,
            target_current_anomaly_beta,
            theta_A,
            max_transfer_revolution,
            max_target_revolution);

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

        const auto theta_alpha_residual =
            problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
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
    return results;
}

void print_theta_summary(
    double theta_prime,
    const std::vector<Problem2SlingshotResidualBranchTestOnly>& branches
) {
    int valid_count = 0;
    int invalid_geometry_count = 0;
    double min_residual = std::numeric_limits<double>::infinity();
    double max_residual = -std::numeric_limits<double>::infinity();

    for (const auto& branch : branches) {
        if (branch.valid) {
            valid_count += 1;
            min_residual = std::min(min_residual, branch.slingshot_residual);
            max_residual = std::max(max_residual, branch.slingshot_residual);
        } else {
            invalid_geometry_count += 1;
        }
    }

    std::cout << "Problem2ThetaPrimeResidualSummary\n";
    std::cout << "theta_prime=" << theta_prime << '\n';
    std::cout << "problem1_branch_count=" << branches.size() << '\n';
    std::cout << "valid_slingshot_branch_count=" << valid_count << '\n';
    std::cout << "invalid_geometry_count=" << invalid_geometry_count << '\n';
    std::cout << "min_slingshot_residual="
              << (valid_count > 0 ? min_residual : std::numeric_limits<double>::quiet_NaN()) << '\n';
    std::cout << "max_slingshot_residual="
              << (valid_count > 0 ? max_residual : std::numeric_limits<double>::quiet_NaN()) << '\n';

    const int to_print = std::min<int>(5, branches.size());
    for (int i = 0; i < to_print; ++i) {
        const auto& branch = branches[static_cast<std::size_t>(i)];
        std::cout << "Problem2SlingshotResidualBranch\n";
        std::cout << "theta_prime=" << branch.theta_prime << '\n';
        std::cout << "branch_index=" << i << '\n';
        std::cout << "k=" << branch.transfer_revolution << '\n';
        std::cout << "q=" << branch.target_revolution << '\n';
        std::cout << "alpha=" << branch.alpha << '\n';
        std::cout << "time_of_flight_seconds=" << branch.time_of_flight_seconds << '\n';
        std::cout << "problem1_residual_seconds=" << branch.problem1_residual_seconds << '\n';
        std::cout << "outgoing_eccentricity=" << branch.outgoing_eccentricity << '\n';
        std::cout << "outgoing_semi_latus_rectum=" << branch.outgoing_semi_latus_rectum << '\n';
        std::cout << "slingshot_residual=" << branch.slingshot_residual << '\n';
        std::cout << "valid=" << (branch.valid ? 1 : 0) << '\n';
        std::cout << "invalid_reason=" << branch.invalid_reason << '\n';
    }
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;

    std::cout << std::setprecision(6) << std::scientific;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time =
        0.17 * planet_params::planet_orbital_period(departure_planet);
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);

    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double incoming_e = 0.3;
    const double incoming_theta = 0.4;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;

    {
        const double theta_prime = 0.2;
        const auto branches =
            evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                theta_prime,
                max_transfer_revolution,
                max_target_revolution);
        assert(!branches.empty());
        int valid_problem1_branch_count = 0;
        int valid_slingshot_branch_count = 0;
        for (const auto& branch : branches) {
            if (std::isfinite(branch.time_of_flight_seconds)) {
                valid_problem1_branch_count += 1;
            }
            if (branch.valid) {
                valid_slingshot_branch_count += 1;
                assert(std::isfinite(branch.slingshot_residual));
            }
        }
        assert(valid_problem1_branch_count > 0);
        assert(valid_slingshot_branch_count > 0);
    }

    for (const double theta_prime : std::array<double, 4>{0.1, 0.5, 1.0, 2.0}) {
        const auto branches =
            evaluate_problem2_slingshot_residual_branches_for_theta_prime_route_a_test_only(
                departure_planet,
                target_planet,
                encounter_time,
                phi,
                beta,
                incoming_e,
                incoming_theta,
                theta_prime,
                max_transfer_revolution,
                max_target_revolution);
        print_theta_summary(theta_prime, branches);
    }

    std::cout << "problem2_slingshot_residual_branches_route_a_ok\n";
    return 0;
}
