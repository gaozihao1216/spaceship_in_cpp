/*
 * 文件作用：验证 Route A 线性外推穿越飞掠行星全局角 ω_J 时的 q 修正。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_route_a.hpp"

#include <cassert>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;

Problem2OutgoingBranchSolution make_test_branch(
    int target_revolution,
    double encounter_phi,
    double dphi_dtheta_prime
) {
    Problem2OutgoingBranchSolution branch{};
    branch.transfer_revolution = 0;
    branch.target_revolution = target_revolution;
    branch.encounter_global_angle = encounter_phi;
    branch.outgoing_eccentricity = 0.2;
    branch.outgoing_semi_latus_rectum = 1.0;
    branch.dphi_dtheta_prime = dphi_dtheta_prime;
    branch.de_dtheta_prime = 0.0;
    branch.has_dphi_dtheta_prime = true;
    branch.has_de_dtheta_prime = true;
    return branch;
}

}  // namespace

int main() {
    namespace problem2 = spaceship_cpp::problem2;

    {
        const auto branch = make_test_branch(1, -0.5, spaceship_cpp::common::kTwoPi * 20.0);
        const double delta_theta_prime = 0.25;
        const auto adjustment = problem2::adjust_target_revolution_for_route_a_linear_prediction(
            PlanetId::Earth,
            0.0,
            branch,
            delta_theta_prime,
            5);
        assert(adjustment.valid);

        const auto flyby_state =
            spaceship_cpp::planet_params::planet_state_at_time(PlanetId::Earth, 0.0);
        const double omega_J = flyby_state.theta_global;
        const double phi_pred = branch.encounter_global_angle +
            branch.dphi_dtheta_prime * delta_theta_prime;
        const int expected_delta =
            static_cast<int>(std::floor((phi_pred - omega_J) / spaceship_cpp::common::kTwoPi)) -
            static_cast<int>(std::floor((branch.encounter_global_angle - omega_J) /
                                        spaceship_cpp::common::kTwoPi));

        assert(adjustment.target_revolution_delta == expected_delta);
        assert(adjustment.adjusted_target_revolution ==
            std::clamp(1 + expected_delta, 0, 5));
        assert(adjustment.omega_J == omega_J);
    }

    {
        const auto branch = make_test_branch(2, 1.0, 0.1);
        const auto adjustment = problem2::adjust_target_revolution_for_route_a_linear_prediction(
            PlanetId::Earth,
            0.0,
            branch,
            0.5,
            2);
        assert(adjustment.valid);
        assert(adjustment.target_revolution_delta == 0);
        assert(adjustment.adjusted_target_revolution == 2);
    }

    {
        const auto branch = make_test_branch(0, 0.0, -spaceship_cpp::common::kTwoPi * 8.0);
        const auto adjustment = problem2::adjust_target_revolution_for_route_a_linear_prediction(
            PlanetId::Mars,
            0.0,
            branch,
            0.5,
            1);
        assert(adjustment.valid);
        assert(adjustment.target_revolution_delta <= 0);
        assert(adjustment.adjusted_target_revolution >= 0);
    }

    {
        Problem2OutgoingBranchSolution branch = make_test_branch(1, 0.2, 1.0);
        branch.has_dphi_dtheta_prime = false;
        const auto adjustment = problem2::adjust_target_revolution_for_route_a_linear_prediction(
            PlanetId::Earth,
            0.0,
            branch,
            0.3,
            2);
        assert(!adjustment.valid);
    }

    {
        const auto branch = make_test_branch(1, 0.3, spaceship_cpp::common::kTwoPi * 3.0);
        auto linear = problem2::estimate_problem2_route_a_solution_linear(0.0, branch, 0.4);
        assert(linear.valid);
        const auto adjustment = problem2::adjust_target_revolution_for_route_a_linear_prediction(
            PlanetId::Earth,
            0.0,
            branch,
            linear.delta_theta_prime,
            3);
        assert(adjustment.valid);
        problem2::apply_target_revolution_adjustment_to_route_a_linear_estimate(adjustment, linear);
        assert(linear.target_revolution_adjusted == (adjustment.target_revolution_delta != 0));
        assert(linear.adjusted_target_revolution == adjustment.adjusted_target_revolution);
    }

    std::cout << "test_problem2_route_a_q_adjustment PASSED\n";
    return 0;
}
