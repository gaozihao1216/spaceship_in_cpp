/*
 * 文件作用：验证 G 的 Hermite 二次零点预测使用求根公式（非极值公式）。
 */
#include "spaceship_cpp/problem2/problem2_flyby_constraint.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

bool approx_equal(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

double extremum_normalized_t(double G_left_derivative, double G_right_derivative) {
    const double denominator = G_right_derivative - G_left_derivative;
    if (std::abs(denominator) <= 1e-15) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return -G_left_derivative / denominator;
}

}  // namespace

int main() {
    namespace problem2 = spaceship_cpp::problem2;

    constexpr double kTol = 1e-10;

    // 线性模型：G(t) = 0.5 - t，根在 t = 0.5；极值公式在此退化。
    {
        const auto estimate = problem2::estimate_flyby_constraint_G_quadratic_root_on_theta_prime_interval(
            0.0,
            0.5,
            -1.0,
            1.0,
            -0.5,
            -1.0);
        assert(estimate.has_value());
        assert(estimate->valid);
        assert(estimate->has_root_in_interval);
        assert(approx_equal(estimate->selected_root.normalized_t, 0.5, kTol));
        assert(approx_equal(estimate->selected_root.theta_prime, 0.5, kTol));
        assert(approx_equal(estimate->selected_root.predicted_G, 0.0, kTol));

        const double t_extremum = extremum_normalized_t(-1.0, -1.0);
        assert(!std::isfinite(t_extremum) || !approx_equal(t_extremum, 0.5, 1e-3));
    }

    // 二次模型：G(t) = t^2 - 0.24，根在 t = sqrt(0.24)。
    {
        const double expected_root = std::sqrt(0.24);
        const auto estimate = problem2::estimate_flyby_constraint_G_quadratic_root_on_theta_prime_interval(
            1.0,
            -0.24,
            0.0,
            2.0,
            0.76,
            2.0);
        assert(estimate.has_value());
        assert(estimate->valid);
        assert(estimate->has_root_in_interval);
        assert(approx_equal(estimate->selected_root.normalized_t, expected_root, kTol));
        assert(approx_equal(estimate->selected_root.predicted_G, 0.0, 1e-9));

        const double t_extremum = extremum_normalized_t(0.0, 2.0);
        assert(approx_equal(t_extremum, 0.0, kTol));
        assert(!approx_equal(t_extremum, expected_root, 1e-2));
    }

    // 区间内无根：G(t) = t^2 + 1。
    {
        const auto estimate = problem2::estimate_flyby_constraint_G_quadratic_root_on_theta_prime_interval(
            0.0,
            1.0,
            0.0,
            1.0,
            2.0,
            2.0);
        assert(estimate.has_value());
        assert(estimate->valid);
        assert(!estimate->has_root_in_interval);
    }

    // 变号区间：二次模型在 (0,1) 内恰有一个根。
    {
        const auto estimate = problem2::estimate_flyby_constraint_G_quadratic_root_on_theta_prime_interval(
            0.0,
            0.2,
            -1.2,
            1.0,
            -0.1,
            0.6);
        assert(estimate.has_value());
        assert(estimate->valid);
        assert(estimate->has_root_in_interval);
        assert(estimate->selected_root.normalized_t > 0.0);
        assert(estimate->selected_root.normalized_t < 1.0);
        assert(approx_equal(estimate->selected_root.predicted_G, 0.0, 1e-9));
    }

    // 同一区间内两个根都在 (0,1)：应选取更接近线性根 t_lin 的那一个。
    {
        const double G_left = 0.24;
        const double G_right = -0.1;
        const double G_left_derivative = -1.0;
        const double G_right_derivative = 1.0;
        const double t_linear = -G_left / (G_right - G_left);

        const double interval_width = 1.0;
        const double a = 0.5 * interval_width * (G_right_derivative - G_left_derivative);
        const double b = interval_width * G_left_derivative;
        const double c = G_left;
        const double discriminant = b * b - 4.0 * a * c;
        assert(discriminant > 0.0);
        const double sqrt_discriminant = std::sqrt(discriminant);
        const double root_minus = (-b - sqrt_discriminant) / (2.0 * a);
        const double root_plus = (-b + sqrt_discriminant) / (2.0 * a);
        assert(root_minus > 0.0 && root_minus < 1.0);
        assert(root_plus > 0.0 && root_plus < 1.0);

        const double expected_root =
            std::abs(root_minus - t_linear) <= std::abs(root_plus - t_linear) ? root_minus : root_plus;

        const auto estimate = problem2::estimate_flyby_constraint_G_quadratic_root_on_theta_prime_interval(
            0.0,
            G_left,
            G_left_derivative,
            1.0,
            G_right,
            G_right_derivative);
        assert(estimate.has_value());
        assert(estimate->valid);
        assert(estimate->has_root_in_interval);
        assert(approx_equal(estimate->selected_root.normalized_t, expected_root, kTol));
        assert(approx_equal(estimate->selected_root.predicted_G, 0.0, 1e-9));
        assert(estimate->alternate_root.has_value());
    }

    std::cout << "test_problem2_G_quadratic_roots PASSED\n";
    return 0;
}
