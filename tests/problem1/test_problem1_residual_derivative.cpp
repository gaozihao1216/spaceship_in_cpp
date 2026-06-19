/*
 * 文件作用：测试 Problem 1 残差对遇合角 phi 的数值导数估计。
 * 主要工作：用细步长中心差分作参考，验证粗扫差分在含 e(phi) 耦合时的误差界。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cassert>
#include <cmath>
#include <optional>

namespace {

using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem1::Problem1ResidualInput;
using spaceship_cpp::problem1::Problem1ResidualResult;
using spaceship_cpp::problem1::Problem1ResidualStatus;

struct ResidualEvaluationContext {
    PlanetId departure_planet = PlanetId::Earth;
    PlanetId target_planet = PlanetId::Mars;
    double launch_time_seconds_since_j2000 = 0.0;
    double transfer_perihelion_angle = 0.5;
    int transfer_revolution = 0;
    int target_revolution = 0;
};

Problem1ResidualInput make_residual_input(const ResidualEvaluationContext& context, double encounter_global_angle) {
    return Problem1ResidualInput{
        context.departure_planet,
        context.target_planet,
        context.launch_time_seconds_since_j2000,
        context.transfer_perihelion_angle,
        encounter_global_angle,
        context.transfer_revolution,
        context.target_revolution,
    };
}

std::optional<Problem1ResidualResult> evaluate_residual_if_success(
    const ResidualEvaluationContext& context,
    double encounter_global_angle
) {
    const Problem1ResidualResult result =
        spaceship_cpp::problem1::evaluate_problem1_residual(make_residual_input(context, encounter_global_angle));
    if (result.status != Problem1ResidualStatus::Success || !std::isfinite(result.residual)) {
        return std::nullopt;
    }
    return result;
}

// 细步长中心差分：作为 residual(phi) 的数值导数参考，包含 e(phi)、p(phi)、deltaF(phi) 等全部耦合。
double numerical_residual_derivative_central(
    const ResidualEvaluationContext& context,
    double phi,
    double step
) {
    const auto minus_result = evaluate_residual_if_success(context, phi - step);
    const auto plus_result = evaluate_residual_if_success(context, phi + step);
    assert(minus_result.has_value());
    assert(plus_result.has_value());
    return (plus_result->residual - minus_result->residual) / (2.0 * step);
}

// 用三点中心差分估计 residual 对 phi 的二阶导，用于粗网格截断误差上界。
double numerical_residual_second_derivative_central(
    const ResidualEvaluationContext& context,
    double phi,
    double step
) {
    const auto minus_result = evaluate_residual_if_success(context, phi - step);
    const auto center_result = evaluate_residual_if_success(context, phi);
    const auto plus_result = evaluate_residual_if_success(context, phi + step);
    assert(minus_result.has_value());
    assert(center_result.has_value());
    assert(plus_result.has_value());
    return (plus_result->residual - 2.0 * center_result->residual + minus_result->residual) / (step * step);
}

double relative_error(double estimated, double reference) {
    return std::abs(estimated - reference) / std::max(std::abs(reference), 1.0);
}

// 中文说明：在单个 phi 点上，对比细步长参考导数与三种差分估计，并确认转移离心率随 phi 变化。
void verify_derivative_estimators_at_phi(
    const ResidualEvaluationContext& context,
    double phi,
    double coarse_step,
    double reference_step
) {
    const auto center_result = evaluate_residual_if_success(context, phi);
    const auto minus_coarse = evaluate_residual_if_success(context, phi - coarse_step);
    const auto plus_coarse = evaluate_residual_if_success(context, phi + coarse_step);
    const auto minus_fine = evaluate_residual_if_success(context, phi - reference_step);
    const auto plus_fine = evaluate_residual_if_success(context, phi + reference_step);
    assert(center_result.has_value());
    assert(minus_coarse.has_value());
    assert(plus_coarse.has_value());
    assert(minus_fine.has_value());
    assert(plus_fine.has_value());

    // 转移轨道离心率 e(phi) 会随遇合角变化；导数估计必须基于真实 residual，而非固定 e 的简化模型。
    const double eccentricity_span = std::max(
        std::abs(plus_coarse->transfer_e_raw - minus_coarse->transfer_e_raw),
        std::abs(plus_coarse->transfer_e - minus_coarse->transfer_e));
    assert(eccentricity_span > 1e-14);

    const double reference_derivative =
        numerical_residual_derivative_central(context, phi, reference_step);
    assert(std::isfinite(reference_derivative));

    const double fine_central_estimate =
        spaceship_cpp::problem1::estimate_problem1_residual_derivative_central(
            phi - reference_step,
            minus_fine->residual,
            phi + reference_step,
            plus_fine->residual);
    assert(std::isfinite(fine_central_estimate));
    assert(relative_error(fine_central_estimate, reference_derivative) <= 1e-9);

    const double coarse_central_estimate =
        spaceship_cpp::problem1::estimate_problem1_residual_derivative_central(
            phi - coarse_step,
            minus_coarse->residual,
            phi + coarse_step,
            plus_coarse->residual);
    assert(std::isfinite(coarse_central_estimate));

    const double half_step = 0.5 * coarse_step;
    const auto minus_half = evaluate_residual_if_success(context, phi - half_step);
    const auto plus_half = evaluate_residual_if_success(context, phi + half_step);
    assert(minus_half.has_value());
    assert(plus_half.has_value());
    const double half_step_derivative =
        spaceship_cpp::problem1::estimate_problem1_residual_derivative_central(
            phi - half_step,
            minus_half->residual,
            phi + half_step,
            plus_half->residual);

    const double abs_coarse_error = std::abs(coarse_central_estimate - reference_derivative);
    const double richardson_error_bound =
        std::abs(coarse_central_estimate - half_step_derivative) / 3.0 + 1e-9;
    double allowed_coarse_error = std::max(
        richardson_error_bound * 4.0,
        0.01 * std::max(std::abs(reference_derivative), 1.0));

    const double second_derivative =
        numerical_residual_second_derivative_central(context, phi, reference_step);
    if (std::isfinite(second_derivative)) {
        allowed_coarse_error = std::max(
            allowed_coarse_error,
            0.5 * std::abs(second_derivative) * coarse_step * coarse_step + 1e-9);
    }
    assert(abs_coarse_error <= allowed_coarse_error);

    const double left_endpoint_estimate =
        spaceship_cpp::problem1::estimate_problem1_residual_derivative_at_left_endpoint(
            phi - coarse_step,
            minus_coarse->residual,
            phi,
            center_result->residual);
    const double right_endpoint_estimate =
        spaceship_cpp::problem1::estimate_problem1_residual_derivative_at_right_endpoint(
            phi,
            center_result->residual,
            phi + coarse_step,
            plus_coarse->residual);
    assert(std::isfinite(left_endpoint_estimate));
    assert(std::isfinite(right_endpoint_estimate));

    const double left_reference =
        (center_result->residual - minus_coarse->residual) / coarse_step;
    const double right_reference =
        (plus_coarse->residual - center_result->residual) / coarse_step;
    assert(relative_error(left_endpoint_estimate, left_reference) <= 1e-12);
    assert(relative_error(right_endpoint_estimate, right_reference) <= 1e-12);
}

void verify_derivative_estimators_over_scan_grid(
    const ResidualEvaluationContext& context,
    int phi_scan_count,
    int minimum_valid_samples
) {
    const double coarse_step =
        spaceship_cpp::common::kTwoPi / static_cast<double>(phi_scan_count);
    const double reference_step = 1e-6;
    int valid_samples = 0;

    for (int index = 1; index < phi_scan_count - 1; ++index) {
        const double phi = static_cast<double>(index) * coarse_step;
        const auto center_result = evaluate_residual_if_success(context, phi);
        if (!center_result.has_value()) {
            continue;
        }
        verify_derivative_estimators_at_phi(context, phi, coarse_step, reference_step);
        ++valid_samples;
    }

    assert(valid_samples >= minimum_valid_samples);
}

}  // namespace

int main() {
    {
        // 中文说明：Earth→Mars、k=q=0 时，在 240 粗扫内点上导数估计应与细步长参考一致，且 e(phi) 非平凡变化。
        verify_derivative_estimators_over_scan_grid(
            ResidualEvaluationContext{},
            240,
            200);
    }

    {
        // 中文说明：k=q=1 多圈分支下，导数通道仍包含变分支后的 e(phi) 与 deltaF(phi) 耦合。
        verify_derivative_estimators_over_scan_grid(
            ResidualEvaluationContext{
                PlanetId::Earth,
                PlanetId::Mars,
                0.0,
                0.5,
                1,
                1,
            },
            180,
            120);
    }

    {
        // 中文说明：另一组发射/轨道方向参数，避免测试只覆盖单一几何点。
        verify_derivative_estimators_over_scan_grid(
            ResidualEvaluationContext{
                PlanetId::Earth,
                PlanetId::Mars,
                86400.0 * 120.0,
                1.1,
                0,
                0,
            },
            240,
            200);
    }

    return 0;
}
