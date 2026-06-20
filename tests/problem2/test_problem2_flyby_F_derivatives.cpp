/*
 * 文件作用：验证 F(e, theta_global) 解析偏导数与中心差分数值导数一致。
 */
#include "spaceship_cpp/problem2/problem2_flyby_constraint.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool approx_equal(double a, double b, double rtol, double atol) {
    return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

double central_derivative_e(
    double e,
    double theta_global,
    double phi,
    double e_J,
    double step
) {
    namespace problem2 = spaceship_cpp::problem2;
    const auto F_plus = problem2::evaluate_flyby_constraint_F(e + step, theta_global, phi, e_J);
    const auto F_minus = problem2::evaluate_flyby_constraint_F(e - step, theta_global, phi, e_J);
    assert(F_plus.valid && F_minus.valid);
    return (F_plus.value - F_minus.value) / (2.0 * step);
}

double central_derivative_theta(
    double e,
    double theta_global,
    double phi,
    double e_J,
    double step
) {
    namespace problem2 = spaceship_cpp::problem2;
    const auto F_plus = problem2::evaluate_flyby_constraint_F(e, theta_global + step, phi, e_J);
    const auto F_minus = problem2::evaluate_flyby_constraint_F(e, theta_global - step, phi, e_J);
    assert(F_plus.valid && F_minus.valid);
    return (F_plus.value - F_minus.value) / (2.0 * step);
}

}  // namespace

int main() {
    namespace problem2 = spaceship_cpp::problem2;

    const double phi = 1.17;
    const double e_J = 0.093;
    const double step = 1e-7;
    const double rtol = 1e-6;
    const double atol = 1e-8;

    const double test_cases[][2] = {
        {0.12, -0.85},
        {0.45, 0.33},
        {0.78, 1.95},
        {0.05, -2.4},
    };

    for (const auto& test_case : test_cases) {
        const double e = test_case[0];
        const double theta_global = test_case[1];

        const auto analytical = problem2::evaluate_flyby_constraint_F_partial_derivatives(
            e,
            theta_global,
            phi,
            e_J);
        assert(analytical.valid);

        const double dF_de_numeric = central_derivative_e(e, theta_global, phi, e_J, step);
        const double dF_dtheta_numeric = central_derivative_theta(e, theta_global, phi, e_J, step);

        assert(approx_equal(analytical.dF_de, dF_de_numeric, rtol, atol));
        assert(approx_equal(analytical.dF_dtheta_global, dF_dtheta_numeric, rtol, atol));
    }

    std::cout << "test_problem2_flyby_F_derivatives PASSED\n";
    return 0;
}
