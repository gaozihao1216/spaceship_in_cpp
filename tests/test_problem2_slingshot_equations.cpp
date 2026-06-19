/*
 * 文件作用：测试 Problem 2 弹弓残差方程。
 * 主要工作：验证轨道半径、外向轨道几何和残差函数的基本有效性。
 */
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool approx_equal(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

}  // namespace

int main() {
    namespace problem2 = spaceship_cpp::problem2;

    {
        const auto invariant = problem2::evaluate_problem2_slingshot_invariant(
            1.0,
            0.1,
            0.3,
            0.4);
        assert(invariant.valid);
        assert(std::isfinite(invariant.invariant));
    }

    {
        const auto residual = problem2::evaluate_problem2_slingshot_residual(
            1.0,
            0.1,
            0.3,
            0.4,
            0.3,
            0.4);
        assert(residual.valid);
        assert(std::abs(residual.residual) < 1e-12);
    }

    {
        const double R_J = 1.0;
        const double e_J = 0.05;
        const double R_K = 1.5;
        const double e_K = 0.1;
        const double phi = 0.7;
        const double alpha = 1.8;
        const double theta_prime = 0.2;
        const auto geometry = problem2::solve_problem2_outgoing_orbit_from_two_points(
            R_J, e_J, R_K, e_K, phi, alpha, theta_prime);
        assert(geometry.valid);
        const double encounter_radius = problem2::problem2_orbit_radius(R_J, e_J, phi);
        const double target_radius = problem2::problem2_orbit_radius(R_K, e_K, alpha);
        const double outgoing_at_encounter =
            problem2::problem2_orbit_radius(geometry.semi_latus_rectum, geometry.eccentricity, phi - theta_prime);
        const double outgoing_at_target =
            problem2::problem2_orbit_radius(geometry.semi_latus_rectum, geometry.eccentricity, alpha - theta_prime);
        assert(approx_equal(outgoing_at_encounter, encounter_radius, 1e-10));
        assert(approx_equal(outgoing_at_target, target_radius, 1e-10));
    }

    {
        const auto geometry = problem2::solve_problem2_outgoing_orbit_from_two_points(
            1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        assert(!geometry.valid);
        assert(geometry.invalid_reason == "geometry_denominator_too_small");
    }

    {
        // 中文注释：旧简化公式只在特殊情形下与当前 invariant 等价，这里只保留为 sanity check。
        const auto residual = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
            1.0,
            0.05,
            1.5,
            0.1,
            0.7,
            1.8,
            0.3,
            0.4,
            0.2);
        assert(!residual.valid || std::isfinite(residual.slingshot_residual));
    }

    std::cout << "problem2_slingshot_equations_ok\n";
    return 0;
}
