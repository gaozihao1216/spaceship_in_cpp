/*
 * 文件作用：测试 Problem 1 二次 Hermite 极值估计。
 * 主要工作：验证区间内/外极值点判定与解析预测值。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cassert>
#include <cmath>

namespace {

bool approx_equal(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) <= eps;
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;

    auto make_success_residual_result = [](double residual) {
        return problem1::Problem1ResidualResult{
            problem1::Problem1ResidualStatus::Success,
            residual,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
        };
    };

    {
        // 中文说明：f(phi)=(phi-0.5)^2 在 [0,1] 上极小点在 0.5，预测残差为 0。
        const auto estimate = problem1::estimate_problem1_residual_quadratic_extremum_on_interval(
            0.0,
            0.25,
            -1.0,
            1.0,
            0.25,
            1.0);
        assert(estimate.has_value());
        assert(approx_equal(estimate->phi, 0.5));
        assert(approx_equal(estimate->predicted_residual, 0.0));
        assert(estimate->is_minimum);
    }

    {
        // 中文说明：f(phi)=-(phi-0.5)^2+0.25 在 [0,1] 上极大点在 0.5。
        const auto estimate = problem1::estimate_problem1_residual_quadratic_extremum_on_interval(
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            -1.0);
        assert(estimate.has_value());
        assert(approx_equal(estimate->phi, 0.5));
        assert(approx_equal(estimate->predicted_residual, 0.25));
        assert(!estimate->is_minimum);
    }

    {
        // 中文说明：线性函数无内点极值，应返回空。
        const auto estimate = problem1::estimate_problem1_residual_quadratic_extremum_on_interval(
            0.0,
            0.0,
            2.0,
            1.0,
            2.0,
            2.0);
        assert(!estimate.has_value());
    }

    {
        // 中文说明：极值落在右端点外（t*>1）时返回空；取 f'_L=-2, f'_R=-1 => t*=2。
        const auto estimate = problem1::estimate_problem1_residual_quadratic_extremum_on_interval(
            0.0,
            1.0,
            -2.0,
            1.0,
            0.0,
            -1.0);
        assert(!estimate.has_value());
    }

    {
        // 中文说明：极值落在左端点外（t*<0）时返回空；取 f'_L=1, f'_R=2 => t*=-1。
        const auto estimate = problem1::estimate_problem1_residual_quadratic_extremum_on_interval(
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            2.0);
        assert(!estimate.has_value());
    }

    {
        // 中文说明：真实 Earth→Mars 残差粗扫区间上，二次极值点应落在开区间内或返回空。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            0,
            0,
            240,
            1e-10,
            0.0,
            80,
        };
        const double coarse_step = spaceship_cpp::common::kTwoPi / 240.0;
        int interior_extrema_count = 0;

        for (int index = 1; index < input.phi_scan_count - 1; ++index) {
            const double left_phi = static_cast<double>(index - 1) * coarse_step;
            const double center_phi = static_cast<double>(index) * coarse_step;
            const double right_phi = static_cast<double>(index + 1) * coarse_step;

            const problem1::Problem1ResidualResult left_result = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    left_phi,
                    0,
                    0,
                });
            const problem1::Problem1ResidualResult center_result = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    center_phi,
                    0,
                    0,
                });
            const problem1::Problem1ResidualResult right_result = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    right_phi,
                    0,
                    0,
                });
            if (left_result.status != problem1::Problem1ResidualStatus::Success ||
                center_result.status != problem1::Problem1ResidualStatus::Success ||
                right_result.status != problem1::Problem1ResidualStatus::Success) {
                continue;
            }

            const double left_derivative =
                problem1::estimate_problem1_residual_derivative_at_right_endpoint(
                    left_phi,
                    left_result.residual,
                    center_phi,
                    center_result.residual);
            const double right_derivative =
                problem1::estimate_problem1_residual_derivative_at_left_endpoint(
                    center_phi,
                    center_result.residual,
                    right_phi,
                    right_result.residual);

            const auto estimate =
                problem1::estimate_problem1_residual_quadratic_extremum_on_interval(
                    left_phi,
                    left_result.residual,
                    left_derivative,
                    right_phi,
                    right_result.residual,
                    right_derivative);
            if (!estimate.has_value()) {
                continue;
            }

            assert(estimate->phi > left_phi);
            assert(estimate->phi < right_phi);
            assert(std::isfinite(estimate->predicted_residual));
            ++interior_extrema_count;
        }

        assert(interior_extrema_count > 0);
    }

    {
        // 中文说明：机制门控——抛物线切触时 ρ=0，应通过 detect；平坦同号区间 ρ≈1 应拒绝。
        const auto tangent_detect = problem1::detect_problem1_fold_interval_by_quadratic_extremum(
            0.0,
            make_success_residual_result(0.25),
            -1.0,
            1.0,
            make_success_residual_result(0.25),
            1.0);
        assert(tangent_detect.has_value());
        assert(approx_equal(tangent_detect->predicted_residual, 0.0));

        const auto flat_detect = problem1::detect_problem1_fold_interval_by_quadratic_extremum(
            0.0,
            make_success_residual_result(1.0),
            -1.0,
            1.0,
            make_success_residual_result(1.0),
            1.0);
        assert(!flat_detect.has_value());
    }

    return 0;
}
