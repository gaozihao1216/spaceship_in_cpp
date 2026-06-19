/*
 * 文件作用：测试 Problem 1 遇合角搜索基础函数。
 * 主要工作：验证导数差分、区间二分与三分法的可调用性与基本行为。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cassert>
#include <cmath>

int main() {
    namespace common = spaceship_cpp::common;
    namespace problem1 = spaceship_cpp::problem1;

    {
        // 中文说明：f(phi)=phi^2 在 phi=1 处导数应为 2；中心差分应接近精确值。
        const double center_phi = 1.0;
        const double previous_phi = 0.9;
        const double next_phi = 1.1;
        const double previous_residual = previous_phi * previous_phi;
        const double next_residual = next_phi * next_phi;
        const double derivative = problem1::estimate_problem1_residual_derivative_central(
            previous_phi,
            previous_residual,
            next_phi,
            next_residual);
        assert(std::isfinite(derivative));
        assert(std::abs(derivative - 2.0) <= 1e-12);

        const double left_derivative = problem1::estimate_problem1_residual_derivative_at_left_endpoint(
            center_phi,
            center_phi * center_phi,
            next_phi,
            next_residual);
        const double right_derivative = problem1::estimate_problem1_residual_derivative_at_right_endpoint(
            previous_phi,
            previous_residual,
            center_phi,
            center_phi * center_phi);
        assert(std::isfinite(left_derivative));
        assert(std::isfinite(right_derivative));
        assert(left_derivative > 2.0);
        assert(right_derivative < 2.0);
    }

    {
        // 中文说明：非法输入（phi 相同）时导数函数返回 NaN。
        const double nan = problem1::estimate_problem1_residual_derivative_at_left_endpoint(
            1.0,
            0.0,
            1.0,
            1.0);
        assert(std::isnan(nan));
    }

    {
        // 中文说明：Earth→Mars 粗扫后，对首个变号区间调用二分应得到有限残差根。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            0,
            0,
            120,
            1e-10,
            0.0,
            80,
        };

        double left_phi = 0.0;
        problem1::Problem1ResidualResult left_result;
        double right_phi = 0.0;
        problem1::Problem1ResidualResult right_result;
        bool found_bracket = false;

        for (int index = 1; index < input.phi_scan_count; ++index) {
            const double phi_left = static_cast<double>(index - 1) * common::kTwoPi /
                static_cast<double>(input.phi_scan_count);
            const double phi_right = static_cast<double>(index) * common::kTwoPi /
                static_cast<double>(input.phi_scan_count);
            const problem1::Problem1ResidualResult result_left = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    phi_left,
                    0,
                    0,
                });
            const problem1::Problem1ResidualResult result_right = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    phi_right,
                    0,
                    0,
                });
            if (result_left.status != problem1::Problem1ResidualStatus::Success ||
                result_right.status != problem1::Problem1ResidualStatus::Success ||
                !std::isfinite(result_left.residual) ||
                !std::isfinite(result_right.residual)) {
                continue;
            }
            if (result_left.residual * result_right.residual < 0.0) {
                left_phi = phi_left;
                left_result = result_left;
                right_phi = phi_right;
                right_result = result_right;
                found_bracket = true;
                break;
            }
        }

        assert(found_bracket);
        const auto bisection = problem1::bisect_problem1_residual_on_interval(
            input,
            left_phi,
            left_result,
            right_phi,
            right_result,
            0,
            0);
        assert(bisection.has_value());
        assert(bisection->success);
        assert(std::isfinite(bisection->phi));
        assert(bisection->phi > left_phi);
        assert(bisection->phi < right_phi);
        assert(bisection->residual_result.status == problem1::Problem1ResidualStatus::Success);
        assert(std::isfinite(bisection->residual_result.residual));
        assert(bisection->iterations_used > 0);
        assert(bisection->bracket_width <= right_phi - left_phi);
    }

    {
        // 中文说明：同号区间调用二分应返回空；三分法在指定极值类型下应返回区间内极值点。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            0,
            0,
            120,
            1e-10,
            0.0,
            80,
        };

        problem1::Problem1ResidualResult same_sign_left;
        problem1::Problem1ResidualResult same_sign_right;
        double same_sign_left_phi = 0.0;
        double same_sign_right_phi = 0.0;
        bool found_same_sign = false;

        for (int index = 1; index < input.phi_scan_count; ++index) {
            const double phi_left = static_cast<double>(index - 1) * common::kTwoPi /
                static_cast<double>(input.phi_scan_count);
            const double phi_right = static_cast<double>(index) * common::kTwoPi /
                static_cast<double>(input.phi_scan_count);
            const problem1::Problem1ResidualResult result_left = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    phi_left,
                    0,
                    0,
                });
            const problem1::Problem1ResidualResult result_right = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    input.departure_planet,
                    input.target_planet,
                    input.launch_time_seconds_since_j2000,
                    input.transfer_perihelion_angle,
                    phi_right,
                    0,
                    0,
                });
            if (result_left.status != problem1::Problem1ResidualStatus::Success ||
                result_right.status != problem1::Problem1ResidualStatus::Success ||
                !std::isfinite(result_left.residual) ||
                !std::isfinite(result_right.residual)) {
                continue;
            }
            if (result_left.residual * result_right.residual > 0.0) {
                same_sign_left_phi = phi_left;
                same_sign_left = result_left;
                same_sign_right_phi = phi_right;
                same_sign_right = result_right;
                found_same_sign = true;
                break;
            }
        }

        assert(found_same_sign);
        const auto rejected_bisection = problem1::bisect_problem1_residual_on_interval(
            input,
            same_sign_left_phi,
            same_sign_left,
            same_sign_right_phi,
            same_sign_right,
            0,
            0);
        assert(!rejected_bisection.has_value());

        const double left_derivative = problem1::estimate_problem1_residual_derivative_at_left_endpoint(
            same_sign_left_phi,
            same_sign_left.residual,
            same_sign_right_phi,
            same_sign_right.residual);
        const double right_derivative = problem1::estimate_problem1_residual_derivative_at_right_endpoint(
            same_sign_left_phi,
            same_sign_left.residual,
            same_sign_right_phi,
            same_sign_right.residual);
        const auto extremum_kind =
            left_derivative < 0.0 && right_derivative > 0.0
                ? problem1::Problem1ResidualExtremumKind::Minimum
                : problem1::Problem1ResidualExtremumKind::Maximum;
        const auto ternary = problem1::ternary_search_problem1_residual_extremum_on_interval(
            input,
            same_sign_left_phi,
            same_sign_left,
            same_sign_right_phi,
            same_sign_right,
            extremum_kind,
            0,
            0,
            24);
        assert(ternary.has_value());
        assert(ternary->success);
        assert(ternary->phi >= same_sign_left_phi);
        assert(ternary->phi <= same_sign_right_phi);
        assert(ternary->residual_result.status == problem1::Problem1ResidualStatus::Success);
        assert(std::isfinite(ternary->residual_result.residual));
        assert(ternary->iterations_used > 0);
    }

    return 0;
}
