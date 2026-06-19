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

// 中文说明：提供固定绝对容差的浮点比较，用于 orbit_F 及其导数数值断言。
bool approx_equal(double a, double b, double eps = 1e-10) {
    return std::abs(a - b) <= eps;
}

// 中文说明：用中心差分近似 orbit_F 对 xi 的导数，作为解析导数/被积函数的数值对照基准。
double numerical_xi_derivative(double e, double xi, double h = 1e-6) {
    return (spaceship_cpp::problem1::orbit_F(e, xi + h) -
            spaceship_cpp::problem1::orbit_F(e, xi - h)) /
        (2.0 * h);
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace problem1 = spaceship_cpp::problem1;

    // 中文说明：验证 e=0（圆轨道）时 orbit_F 退化为恒等映射 xi → xi。
    assert(approx_equal(problem1::orbit_F(0.0, 0.0), 0.0));
    assert(approx_equal(problem1::orbit_F(0.0, 1.2), 1.2));
    assert(approx_equal(problem1::orbit_F(0.0, -1.2), -1.2));
    assert(approx_equal(problem1::orbit_F(0.0, common::kTwoPi), common::kTwoPi));

    {
        // 中文说明：验证椭圆轨道 (e<1) 下 orbit_F_xi_derivative 与 1/(1+e cos xi)^2 解析式一致。
        const double e = 0.3;
        const double xi = 0.7;
        const double expected = 1.0 / std::pow(1.0 + e * std::cos(xi), 2.0);
        assert(approx_equal(problem1::orbit_F_xi_derivative(e, xi), expected, 1e-12));
    }

    {
        // 中文说明：验证抛物线轨道 (e=1) 下导数公式仍与标准被积函数形式一致。
        const double e = 1.0;
        const double xi = 0.5;
        const double expected = 1.0 / std::pow(1.0 + e * std::cos(xi), 2.0);
        assert(approx_equal(problem1::orbit_F_xi_derivative(e, xi), expected, 1e-12));
    }

    {
        // 中文说明：验证双曲线轨道 (e>1) 下导数公式在有效 xi 区间内仍正确。
        const double e = 1.5;
        const double xi = 0.4;
        const double expected = 1.0 / std::pow(1.0 + e * std::cos(xi), 2.0);
        assert(approx_equal(problem1::orbit_F_xi_derivative(e, xi), expected, 1e-12));
    }

    {
        // 中文说明：验证 orbit_F_integrand 与 orbit_F 对 xi 的数值导数一致（导数即被积函数）。
        const double e = 0.3;
        const double xi = 0.7;
        const double derivative = numerical_xi_derivative(e, xi);
        const double integrand = problem1::orbit_F_integrand(e, xi);
        assert(approx_equal(derivative, integrand, 1e-7));
    }

    {
        // 中文说明：验证椭圆轨道解析导数与中心差分数值导数在 e=0.3 处吻合。
        const double e = 0.3;
        const double xi = 0.7;
        const double numerical = numerical_xi_derivative(e, xi);
        const double analytic = problem1::orbit_F_xi_derivative(e, xi);
        assert(approx_equal(numerical, analytic, 1e-7));
    }

    {
        // 中文说明：验证椭圆轨道 orbit_F 沿 xi 增加 2pi 的增量等于一个轨道周期的标准公式。
        const double e = 0.3;
        const double xi = 0.7;
        const double delta = 1.0 - e * e;
        const double expected = common::kTwoPi / (delta * std::sqrt(delta));
        const double actual = problem1::orbit_F(e, xi + common::kTwoPi) - problem1::orbit_F(e, xi);
        assert(approx_equal(actual, expected, 1e-10));
    }

    {
        // 中文说明：验证抛物线轨道 e=1 时 orbit_F 使用 tan(xi/2) 级数展开的特殊形式。
        const double xi = 0.8;
        const double tangent = std::tan(0.5 * xi);
        const double expected = 0.5 * tangent + tangent * tangent * tangent / 6.0;
        assert(approx_equal(problem1::orbit_F(1.0, xi), expected, 1e-12));
    }

    {
        // 中文说明：验证抛物线在 xi=pi（奇点）处 try_orbit_F 返回失败而非 NaN 静默传播。
        const problem1::OrbitFResult result = problem1::try_orbit_F(1.0, common::kPi);
        assert(!result.ok);
        assert(std::isnan(result.value));
    }

    {
        // 中文说明：验证抛物线 e=1 时解析导数与数值导数在非奇点处一致。
        const double e = 1.0;
        const double xi = 0.5;
        const double numerical = numerical_xi_derivative(e, xi);
        const double analytic = problem1::orbit_F_xi_derivative(e, xi);
        assert(approx_equal(numerical, analytic, 1e-7));
    }

    {
        // 中文说明：验证抛物线在 xi=pi 处 try_orbit_F_xi_derivative 同样标记为无效。
        const problem1::OrbitFResult result = problem1::try_orbit_F_xi_derivative(1.0, common::kPi);
        assert(!result.ok);
    }

    {
        // 中文说明：验证双曲线 e=1.5 时被积函数与数值导数一致。
        const double e = 1.5;
        const double xi = 0.4;
        const double derivative = numerical_xi_derivative(e, xi);
        const double integrand = problem1::orbit_F_integrand(e, xi);
        assert(approx_equal(derivative, integrand, 1e-7));
    }

    {
        // 中文说明：验证双曲线解析导数与数值导数在有效 xi 处吻合。
        const double e = 1.5;
        const double xi = 0.4;
        const double numerical = numerical_xi_derivative(e, xi);
        const double analytic = problem1::orbit_F_xi_derivative(e, xi);
        assert(approx_equal(numerical, analytic, 1e-7));
    }

    {
        // 中文说明：验证双曲线在可达上界 xi=acos(-1/e) 处 try_orbit_F 因奇点返回失败。
        const double e = 1.5;
        const double xi_limit = std::acos(-1.0 / e);
        const problem1::OrbitFResult result = problem1::try_orbit_F(e, xi_limit);
        assert(!result.ok);
        assert(std::isnan(result.value));
    }

    {
        // 中文说明：验证双曲线在 xi 上界处 try_orbit_F_xi_derivative 同样标记无效。
        const double e = 1.5;
        const double xi_limit = std::acos(-1.0 / e);
        const problem1::OrbitFResult result = problem1::try_orbit_F_xi_derivative(e, xi_limit);
        assert(!result.ok);
    }

    {
        // 中文说明：验证负偏心率 e<0 时 orbit_F 抛出 domain_error 而非未定义行为。
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
