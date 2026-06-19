/*
 * 文件作用：测试 Problem 1 轨道时间积分工具。
 * 主要工作：覆盖椭圆轨道时间积分的方向、周期和多圈分支行为。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/orbit_time_integral.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {

bool approx_equal(double a, double b, double eps = 1e-10) {
    return std::abs(a - b) <= eps;
}

double numerical_xi_derivative(double e, double xi, double h = 1e-6) {
    return (spaceship_cpp::problem1::orbit_F(e, xi + h) -
            spaceship_cpp::problem1::orbit_F(e, xi - h)) /
        (2.0 * h);
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace problem1 = spaceship_cpp::problem1;

    assert(approx_equal(problem1::orbit_F(0.0, 0.0), 0.0));
    assert(approx_equal(problem1::orbit_F(0.0, 1.2), 1.2));
    assert(approx_equal(problem1::orbit_F(0.0, -1.2), -1.2));
    assert(approx_equal(problem1::orbit_F(0.0, common::kTwoPi), common::kTwoPi));

    {
        const double e = 0.3;
        const double xi = 0.7;
        const double expected = 1.0 / std::pow(1.0 + e * std::cos(xi), 2.0);
        assert(approx_equal(problem1::orbit_F_xi_derivative(e, xi), expected, 1e-12));
    }

    {
        const double e = 1.0;
        const double xi = 0.5;
        const double expected = 1.0 / std::pow(1.0 + e * std::cos(xi), 2.0);
        assert(approx_equal(problem1::orbit_F_xi_derivative(e, xi), expected, 1e-12));
    }

    {
        const double e = 1.5;
        const double xi = 0.4;
        const double expected = 1.0 / std::pow(1.0 + e * std::cos(xi), 2.0);
        assert(approx_equal(problem1::orbit_F_xi_derivative(e, xi), expected, 1e-12));
    }

    {
        const double e = 0.3;
        const double xi = 0.7;
        const double derivative = numerical_xi_derivative(e, xi);
        const double integrand = problem1::orbit_F_integrand(e, xi);
        assert(approx_equal(derivative, integrand, 1e-7));
    }

    {
        const double e = 0.3;
        const double xi = 0.7;
        const double numerical = numerical_xi_derivative(e, xi);
        const double analytic = problem1::orbit_F_xi_derivative(e, xi);
        assert(approx_equal(numerical, analytic, 1e-7));
    }

    {
        const double e = 0.3;
        const double xi = 0.7;
        const double delta = 1.0 - e * e;
        const double expected = common::kTwoPi / (delta * std::sqrt(delta));
        const double actual = problem1::orbit_F(e, xi + common::kTwoPi) - problem1::orbit_F(e, xi);
        assert(approx_equal(actual, expected, 1e-10));
    }

    {
        const double xi = 0.8;
        const double tangent = std::tan(0.5 * xi);
        const double expected = 0.5 * tangent + tangent * tangent * tangent / 6.0;
        assert(approx_equal(problem1::orbit_F(1.0, xi), expected, 1e-12));
    }

    {
        const problem1::OrbitFResult result = problem1::try_orbit_F(1.0, common::kPi);
        assert(!result.ok);
        assert(std::isnan(result.value));
    }

    {
        const double e = 1.0;
        const double xi = 0.5;
        const double numerical = numerical_xi_derivative(e, xi);
        const double analytic = problem1::orbit_F_xi_derivative(e, xi);
        assert(approx_equal(numerical, analytic, 1e-7));
    }

    {
        const problem1::OrbitFResult result = problem1::try_orbit_F_xi_derivative(1.0, common::kPi);
        assert(!result.ok);
    }

    {
        const double e = 1.5;
        const double xi = 0.4;
        const double derivative = numerical_xi_derivative(e, xi);
        const double integrand = problem1::orbit_F_integrand(e, xi);
        assert(approx_equal(derivative, integrand, 1e-7));
    }

    {
        const double e = 1.5;
        const double xi = 0.4;
        const double numerical = numerical_xi_derivative(e, xi);
        const double analytic = problem1::orbit_F_xi_derivative(e, xi);
        assert(approx_equal(numerical, analytic, 1e-7));
    }

    {
        const double e = 1.5;
        const double xi_limit = std::acos(-1.0 / e);
        const problem1::OrbitFResult result = problem1::try_orbit_F(e, xi_limit);
        assert(!result.ok);
        assert(std::isnan(result.value));
    }

    {
        const double e = 1.5;
        const double xi_limit = std::acos(-1.0 / e);
        const problem1::OrbitFResult result = problem1::try_orbit_F_xi_derivative(e, xi_limit);
        assert(!result.ok);
    }

    {
        bool threw = false;
        try {
            (void)problem1::orbit_F(-0.1, 0.2);
        } catch (const std::domain_error&) {
            threw = true;
        }
        assert(threw);
    }

    return 0;
}
