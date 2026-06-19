/*
 * 文件作用：实现 Problem 1 的直接残差评估和求解器。
 * 主要工作：构造转移轨道残差，扫描遇合角根区间，并用二分法细化候选解。
 */
#include "spaceship_cpp/problem1/problem1.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace spaceship_cpp::problem1 {

// 由两端点半径和相对近日点角解转移轨道偏心率；
// 核心公式 e = (r2-r1)/(r1·cos(ξ1)-r2·cos(ξ2))。
double compute_transfer_e_from_two_points(double r1, double xi1, double r2, double xi2) {
    if (!spaceship_cpp::common::is_finite(r1) || !spaceship_cpp::common::is_finite(r2) ||
        !spaceship_cpp::common::is_finite(xi1) || !spaceship_cpp::common::is_finite(xi2)) {
        throw std::domain_error("compute_transfer_e_from_two_points requires finite inputs");
    }
    if (!(r1 > 0.0) || !(r2 > 0.0)) {
        throw std::domain_error("compute_transfer_e_from_two_points requires positive radii");
    }

    const double denominator = r1 * std::cos(xi1) - r2 * std::cos(xi2);
    if (spaceship_cpp::common::near_zero(denominator, spaceship_cpp::common::kDefaultEpsilon)) {
        throw std::domain_error("compute_transfer_e_from_two_points has singular denominator");
    }

    return (r2 - r1) / denominator;
}

// 由出发点半径、偏心率和相对角反推半通径 p = r1·(1+e·cos(ξ1))。
double compute_transfer_p_from_departure(double r1, double e_transfer, double xi1) {
    if (!spaceship_cpp::common::is_finite(r1) || !spaceship_cpp::common::is_finite(e_transfer) ||
        !spaceship_cpp::common::is_finite(xi1)) {
        throw std::domain_error("compute_transfer_p_from_departure requires finite inputs");
    }
    if (!(r1 > 0.0)) {
        throw std::domain_error("compute_transfer_p_from_departure requires positive radius");
    }

    const double factor = 1.0 + e_transfer * std::cos(xi1);
    if (!(spaceship_cpp::common::is_finite(factor) &&
          factor > spaceship_cpp::common::kDefaultEpsilon)) {
        throw std::domain_error("compute_transfer_p_from_departure requires positive orbit factor");
    }
    return r1 * factor;
}

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::near_zero;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::common::orbit_F;

// 构造一个仅含 status 字段、其余为 NaN 的残差结果；用于快速返回失败状态。
Problem1ResidualResult make_result(Problem1ResidualStatus status) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return Problem1ResidualResult{
        status,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
    };
}

// 检查 1+e·cos(angle) > 0；确保轨道半径公式分母为正。
bool positive_orbit_denominator(double e, double angle) {
    const double denominator = 1.0 + e * std::cos(angle);
    return is_finite(denominator) && denominator > kDefaultEpsilon;
}

// 校验 Problem1SolveInput 各字段的合法性；非法输入直接抛 invalid_argument。
void validate_problem1_solve_input(const Problem1SolveInput& input) {
    if (!is_finite(input.launch_time_seconds_since_j2000) ||
        !is_finite(input.transfer_perihelion_angle) ||
        input.max_transfer_revolution < 0 ||
        input.max_target_revolution < 0 ||
        input.phi_scan_count < 3 ||
        !is_finite(input.phi_tolerance) ||
        !(input.phi_tolerance > 0.0) ||
        !is_finite(input.residual_tolerance) ||
        input.residual_tolerance < 0.0 ||
        input.max_bisection_iterations <= 0 ||
        !is_finite(input.max_candidate_relative_residual) ||
        !(input.max_candidate_relative_residual > 0.0)) {
        throw std::invalid_argument("invalid Problem1SolveInput");
    }
}

// 发射时刻行星状态；与遇合角 phi 无关，可在 phi 扫描前一次性计算并复用。
struct Problem1LaunchState {
    double r1 = 0.0;
    double lambda1 = 0.0;
    double theta2_start = 0.0;
    double target_e = 0.0;
    double target_p = 0.0;
    double target_theta_0 = 0.0;
};

Problem1LaunchState compute_problem1_launch_state(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double launch_time_seconds_since_j2000
) {
    const planet_params::PlanetParams& target =
        planet_params::get_planet_params(target_planet);
    const planet_params::PlanetState state1 =
        planet_params::planet_state_at_time(departure_planet, launch_time_seconds_since_j2000);
    return Problem1LaunchState{
        state1.radius,
        state1.theta_global,
        planet_params::planet_true_anomaly_at_time(
            target_planet, launch_time_seconds_since_j2000),
        target.orbit.e,
        target.orbit.p,
        target.orbit.theta_0,
    };
}

Problem1ResidualResult evaluate_problem1_residual_with_launch_state(
    const Problem1ResidualInput& input,
    const Problem1LaunchState& launch
) {
    if (!is_finite(input.launch_time_seconds_since_j2000) ||
        !is_finite(input.transfer_perihelion_angle) ||
        !is_finite(input.encounter_global_angle) ||
        input.transfer_revolution < 0 ||
        input.target_revolution < 0) {
        return make_result(Problem1ResidualStatus::InvalidInput);
    }

    const double r1 = launch.r1;
    const double lambda1 = launch.lambda1;
    const double theta2_start = launch.theta2_start;

    const double u2_base = normalize_angle_0_2pi(
        input.encounter_global_angle - launch.target_theta_0);
    const double r2 = launch.target_p / (1.0 + launch.target_e * std::cos(u2_base));

    const double raw_phi0 = input.transfer_perihelion_angle;
    const double xi1_raw_base = normalize_angle_0_2pi(lambda1 - raw_phi0);
    const double xi2_raw_base = normalize_angle_0_2pi(input.encounter_global_angle - raw_phi0);

    double e_raw = std::numeric_limits<double>::quiet_NaN();
    double p_transfer = std::numeric_limits<double>::quiet_NaN();
    try {
        e_raw = compute_transfer_e_from_two_points(r1, xi1_raw_base, r2, xi2_raw_base);
        p_transfer = compute_transfer_p_from_departure(r1, e_raw, xi1_raw_base);
    } catch (const std::domain_error&) {
        return make_result(Problem1ResidualStatus::SingularGeometry);
    }

    double e_transfer = e_raw;
    double phi0_used = raw_phi0;
    if (e_transfer < 0.0) {
        e_transfer = -e_transfer;
        phi0_used = raw_phi0 + spaceship_cpp::common::kPi;
    }
    phi0_used = normalize_angle_0_2pi(phi0_used);

    double xi1_geometry = xi1_raw_base;
    double xi2_geometry = xi2_raw_base;
    if (e_raw < 0.0) {
        xi1_geometry = normalize_angle_0_2pi(lambda1 - phi0_used);
        xi2_geometry = normalize_angle_0_2pi(input.encounter_global_angle - phi0_used);
    }

    if (!is_finite(e_transfer) || !(e_transfer >= 0.0) ||
        !is_finite(p_transfer) || !(p_transfer > 0.0) ||
        !positive_orbit_denominator(e_transfer, xi1_geometry) ||
        !positive_orbit_denominator(e_transfer, xi2_geometry)) {
        return make_result(Problem1ResidualStatus::InvalidTransferOrbit);
    }

    double theta2_end = u2_base;
    while (theta2_end <= theta2_start) {
        theta2_end += kTwoPi;
    }
    theta2_end += static_cast<double>(input.target_revolution) * kTwoPi;

    double xi1 = std::numeric_limits<double>::quiet_NaN();
    double xi2 = std::numeric_limits<double>::quiet_NaN();
    if (e_transfer < 1.0) {
        xi1 = xi1_geometry;
        xi2 = xi2_geometry;
        while (xi2 <= xi1) {
            xi2 += kTwoPi;
        }
        xi2 += static_cast<double>(input.transfer_revolution) * kTwoPi;
    } else {
        if (input.transfer_revolution != 0) {
            return make_result(Problem1ResidualStatus::InvalidBranch);
        }
        xi1 = normalize_angle_minus_pi_pi(lambda1 - phi0_used);
        xi2 = normalize_angle_minus_pi_pi(input.encounter_global_angle - phi0_used);
        if (xi2 <= xi1) {
            return make_result(Problem1ResidualStatus::InvalidBranch);
        }
    }

    double deltaF_transfer = std::numeric_limits<double>::quiet_NaN();
    double deltaF_target = std::numeric_limits<double>::quiet_NaN();
    try {
        deltaF_transfer = orbit_F(e_transfer, xi2) - orbit_F(e_transfer, xi1);
        deltaF_target = orbit_F(launch.target_e, theta2_end) - orbit_F(launch.target_e, theta2_start);
    } catch (const std::domain_error&) {
        return make_result(Problem1ResidualStatus::InvalidBranch);
    }

    if (!is_finite(deltaF_transfer) || !is_finite(deltaF_target) ||
        !(deltaF_transfer > 0.0) || !(deltaF_target > 0.0)) {
        return make_result(Problem1ResidualStatus::InvalidTimeOfFlight);
    }

    const double transfer_time_scale_free = std::pow(p_transfer, 1.5) * deltaF_transfer;
    const double target_time_scale_free = std::pow(launch.target_p, 1.5) * deltaF_target;
    if (!is_finite(transfer_time_scale_free) || !is_finite(target_time_scale_free)) {
        return make_result(Problem1ResidualStatus::InvalidTimeOfFlight);
    }

    return Problem1ResidualResult{
        Problem1ResidualStatus::Success,
        transfer_time_scale_free - target_time_scale_free,
        r1,
        r2,
        xi1,
        xi2,
        theta2_start,
        theta2_end,
        e_raw,
        e_transfer,
        phi0_used,
        p_transfer,
        deltaF_transfer,
        deltaF_target,
        transfer_time_scale_free,
        target_time_scale_free,
    };
}

// 从求解输入和当前 (φ, k, q) 构造残差评估输入。
Problem1ResidualInput make_residual_input(
    const Problem1SolveInput& input,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    return Problem1ResidualInput{
        input.departure_planet,
        input.target_planet,
        input.launch_time_seconds_since_j2000,
        input.transfer_perihelion_angle,
        encounter_global_angle,
        transfer_revolution,
        target_revolution,
    };
}

// 残差计算成功且 residual 为有限值。
bool is_success_with_finite_residual(const Problem1ResidualResult& result) {
    return result.status == Problem1ResidualStatus::Success && is_finite(result.residual);
}

// 判断两个残差值是否异号（含零）；用于扫描阶段发现变号区间。
bool residual_sign_changed(double a, double b) {
    if (!is_finite(a) || !is_finite(b)) {
        return false;
    }
    if (a == 0.0 || b == 0.0) {
        return true;
    }
    return (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0);
}

// 粗扫样本；用于在相邻点之间估计残差导数并识别 fold/tangent 区间。
struct PhiScanSample {
    double phi = 0.0;
    Problem1ResidualResult result = make_result(Problem1ResidualStatus::InvalidBranch);
};

// fold 区间机制门控：二次极值相对端点深度比上限 ρ，及极值预测相对残差上限。
constexpr double kFoldPredictedExtremumRelativeDepthMax = 0.5;
constexpr double kFoldPredictedExtremumRelativeResidualMax = 1e-4;
constexpr int kFoldTernaryIterations = 48;

// 用相邻粗扫点做中心/单侧差分，估计 residual(phi) 在样本点处的导数。
double estimate_residual_derivative_at_sample(
    const std::vector<PhiScanSample>& samples,
    std::size_t index
) {
    if (samples.size() < 2 || index >= samples.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (index == 0) {
        return estimate_problem1_residual_derivative_at_left_endpoint(
            samples[0].phi,
            samples[0].result.residual,
            samples[1].phi,
            samples[1].result.residual);
    }
    if (index + 1 == samples.size()) {
        return estimate_problem1_residual_derivative_at_right_endpoint(
            samples[index - 1].phi,
            samples[index - 1].result.residual,
            samples[index].phi,
            samples[index].result.residual);
    }
    return estimate_problem1_residual_derivative_central(
        samples[index - 1].phi,
        samples[index - 1].result.residual,
        samples[index + 1].phi,
        samples[index + 1].result.residual);
}

// 计算残差的归一化尺度 max(|T_transfer|, |T_target|, 1)；
// 解决时间尺度很大时绝对残差不足以衡量精度的问题。
double compute_problem1_residual_scale(const Problem1ResidualResult& result) {
    if (result.status != Problem1ResidualStatus::Success ||
        !is_finite(result.transfer_time_scale_free) ||
        !is_finite(result.target_time_scale_free)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::max({std::abs(result.transfer_time_scale_free), std::abs(result.target_time_scale_free), 1.0});
}

// 相对残差 = |residual| / scale；用于过滤未真正收敛的伪候选。
double compute_problem1_relative_residual(const Problem1ResidualResult& result) {
    const double scale = compute_problem1_residual_scale(result);
    if (!is_finite(scale) || !(scale > 0.0) || !is_finite(result.residual)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::abs(result.residual) / scale;
}

// 从成功的残差结果构造候选解；不通过相对残差阈值则返回 nullopt。
std::optional<Problem1Candidate> make_candidate_from_residual_result(
    const Problem1SolveInput& input,
    double encounter_global_angle,
    const Problem1ResidualResult& residual_result,
    int transfer_revolution,
    int target_revolution,
    double root_bracket_width,
    int bisection_iterations,
    bool refined_by_bisection
) {
    if (!is_success_with_finite_residual(residual_result)) {
        return std::nullopt;
    }

    const double mu = planet_params::get_solar_system_physical_params().GM_sun;
    const double time_of_flight_seconds = residual_result.transfer_time_scale_free / std::sqrt(mu);
    if (!is_finite(time_of_flight_seconds) || !(time_of_flight_seconds > 0.0)) {
        return std::nullopt;
    }

    const double residual_scale = compute_problem1_residual_scale(residual_result);
    const double relative_residual = compute_problem1_relative_residual(residual_result);
    if (!is_finite(residual_scale) || !(residual_scale > 0.0) ||
        !is_finite(relative_residual) || relative_residual < 0.0 ||
        !is_finite(root_bracket_width) || root_bracket_width < 0.0 ||
        bisection_iterations < 0) {
        return std::nullopt;
    }
    if (relative_residual > input.max_candidate_relative_residual) {
        return std::nullopt;
    }

    return Problem1Candidate{
        residual_result,
        normalize_angle_0_2pi(encounter_global_angle),
        input.launch_time_seconds_since_j2000,
        time_of_flight_seconds,
        input.launch_time_seconds_since_j2000 + time_of_flight_seconds,
        residual_scale,
        relative_residual,
        root_bracket_width,
        bisection_iterations,
        refined_by_bisection,
        transfer_revolution,
        target_revolution,
    };
}

// 将区间细化结果转为候选解；统一供二分与三分精修复用。
std::optional<Problem1Candidate> make_candidate_from_phi_refinement(
    const Problem1SolveInput& input,
    const Problem1PhiRefinementResult& refinement,
    int transfer_revolution,
    int target_revolution,
    bool refined_by_bisection
) {
    if (!refinement.success) {
        return std::nullopt;
    }
    return make_candidate_from_residual_result(
        input,
        refinement.phi,
        refinement.residual_result,
        transfer_revolution,
        target_revolution,
        refinement.bracket_width,
        refinement.iterations_used,
        refined_by_bisection);
}

// 在已知变号区间 [left_phi, right_phi] 上用二分法细化遇合角根。
std::optional<Problem1PhiRefinementResult> bisect_problem1_residual_on_interval_with_launch_state(
    const Problem1SolveInput& input,
    const Problem1LaunchState& launch_state,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    int transfer_revolution,
    int target_revolution
) {
    if (!is_success_with_finite_residual(left_result) || !is_success_with_finite_residual(right_result)) {
        return std::nullopt;
    }
    if (!(right_phi > left_phi)) {
        return std::nullopt;
    }

    const double bracket_width = right_phi - left_phi;

    if (std::abs(left_result.residual) <= input.residual_tolerance) {
        return Problem1PhiRefinementResult{
            true,
            0,
            left_phi,
            left_result,
            bracket_width,
        };
    }
    if (std::abs(right_result.residual) <= input.residual_tolerance) {
        return Problem1PhiRefinementResult{
            true,
            0,
            right_phi,
            right_result,
            bracket_width,
        };
    }
    if (!residual_sign_changed(left_result.residual, right_result.residual)) {
        return std::nullopt;
    }

    double current_left_phi = left_phi;
    double current_right_phi = right_phi;
    Problem1ResidualResult current_left_result = left_result;
    Problem1ResidualResult current_right_result = right_result;
    double mid_phi = 0.5 * (current_left_phi + current_right_phi);
    Problem1ResidualResult mid_result = make_result(Problem1ResidualStatus::InvalidBranch);
    int iterations_used = 0;

    for (int iteration = 0; iteration < input.max_bisection_iterations; ++iteration) {
        iterations_used = iteration + 1;
        mid_phi = 0.5 * (current_left_phi + current_right_phi);
        mid_result = evaluate_problem1_residual_with_launch_state(
            make_residual_input(input, mid_phi, transfer_revolution, target_revolution),
            launch_state);
        if (!is_success_with_finite_residual(mid_result)) {
            return std::nullopt;
        }
        if (std::abs(mid_result.residual) <= input.residual_tolerance ||
            std::abs(current_right_phi - current_left_phi) <= input.phi_tolerance) {
            return Problem1PhiRefinementResult{
                true,
                iterations_used,
                mid_phi,
                mid_result,
                std::abs(current_right_phi - current_left_phi),
            };
        }

        if (residual_sign_changed(current_left_result.residual, mid_result.residual)) {
            current_right_phi = mid_phi;
            current_right_result = mid_result;
        } else if (residual_sign_changed(mid_result.residual, current_right_result.residual)) {
            current_left_phi = mid_phi;
            current_left_result = mid_result;
        } else {
            return std::nullopt;
        }
    }

    return Problem1PhiRefinementResult{
        true,
        iterations_used,
        mid_phi,
        mid_result,
        std::abs(current_right_phi - current_left_phi),
    };
}

std::optional<Problem1Candidate> refine_problem1_root_by_bisection(
    const Problem1SolveInput& input,
    const Problem1LaunchState& launch_state,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    int transfer_revolution,
    int target_revolution
) {
    const auto refined = bisect_problem1_residual_on_interval_with_launch_state(
        input,
        launch_state,
        left_phi,
        left_result,
        right_phi,
        right_result,
        transfer_revolution,
        target_revolution);
    if (!refined.has_value() || !refined->success) {
        return std::nullopt;
    }
    return make_candidate_from_phi_refinement(
        input,
        *refined,
        transfer_revolution,
        target_revolution,
        true);
}

// 在单峰假设下用三分法找 residual(phi) 的极值点。
std::optional<Problem1PhiRefinementResult> ternary_search_problem1_residual_extremum_on_interval_with_launch_state(
    const Problem1SolveInput& input,
    const Problem1LaunchState& launch_state,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    Problem1ResidualExtremumKind extremum_kind,
    int transfer_revolution,
    int target_revolution,
    int max_iterations
) {
    if (!is_success_with_finite_residual(left_result) || !is_success_with_finite_residual(right_result)) {
        return std::nullopt;
    }
    if (!(right_phi > left_phi) || max_iterations <= 0) {
        return std::nullopt;
    }

    const bool seek_minimum = extremum_kind == Problem1ResidualExtremumKind::Minimum;
    auto extremum_value = [&](const Problem1ResidualResult& result) {
        return seek_minimum ? result.residual : -result.residual;
    };

    double current_left_phi = left_phi;
    double current_right_phi = right_phi;
    int iterations_used = 0;

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        iterations_used = iteration + 1;
        const double span = current_right_phi - current_left_phi;
        if (span <= input.phi_tolerance) {
            break;
        }

        const double m1_phi = current_left_phi + span / 3.0;
        const double m2_phi = current_right_phi - span / 3.0;
        const Problem1ResidualResult m1_result = evaluate_problem1_residual_with_launch_state(
            make_residual_input(input, m1_phi, transfer_revolution, target_revolution),
            launch_state);
        const Problem1ResidualResult m2_result = evaluate_problem1_residual_with_launch_state(
            make_residual_input(input, m2_phi, transfer_revolution, target_revolution),
            launch_state);
        if (!is_success_with_finite_residual(m1_result) || !is_success_with_finite_residual(m2_result)) {
            return std::nullopt;
        }

        const double m1_value = extremum_value(m1_result);
        const double m2_value = extremum_value(m2_result);
        if (m1_value < m2_value) {
            current_right_phi = m2_phi;
        } else {
            current_left_phi = m1_phi;
        }
    }

    const double extremum_phi = 0.5 * (current_left_phi + current_right_phi);
    const Problem1ResidualResult extremum_result = evaluate_problem1_residual_with_launch_state(
        make_residual_input(input, extremum_phi, transfer_revolution, target_revolution),
        launch_state);
    if (!is_success_with_finite_residual(extremum_result)) {
        return std::nullopt;
    }

    return Problem1PhiRefinementResult{
        true,
        iterations_used,
        extremum_phi,
        extremum_result,
        right_phi - left_phi,
    };
}

// fold/tangent 区间：三分精修极值后收录切触根，并对极值两侧变号子区间二分。
std::vector<Problem1Candidate> refine_fold_interval_by_quadratic_extremum(
    const Problem1SolveInput& input,
    const Problem1LaunchState& launch_state,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    const Problem1QuadraticExtremumEstimate& quadratic_estimate,
    int transfer_revolution,
    int target_revolution
) {
    std::vector<Problem1Candidate> refined;

    const auto ternary_refined = ternary_search_problem1_residual_extremum_on_interval_with_launch_state(
        input,
        launch_state,
        left_phi,
        left_result,
        right_phi,
        right_result,
        quadratic_estimate.is_minimum ? Problem1ResidualExtremumKind::Minimum
                                      : Problem1ResidualExtremumKind::Maximum,
        transfer_revolution,
        target_revolution,
        kFoldTernaryIterations);
    if (!ternary_refined.has_value() || !ternary_refined->success) {
        return refined;
    }

    if (auto tangent_candidate = make_candidate_from_phi_refinement(
            input,
            *ternary_refined,
            transfer_revolution,
            target_revolution,
            true)) {
        refined.push_back(*tangent_candidate);
    }

    const Problem1ResidualResult& extremum_result = ternary_refined->residual_result;
    const double extremum_phi = ternary_refined->phi;

    if (residual_sign_changed(left_result.residual, extremum_result.residual)) {
        if (auto left_refined = refine_problem1_root_by_bisection(
                input,
                launch_state,
                left_phi,
                left_result,
                extremum_phi,
                extremum_result,
                transfer_revolution,
                target_revolution)) {
            refined.push_back(*left_refined);
        }
    }

    if (residual_sign_changed(extremum_result.residual, right_result.residual)) {
        if (auto right_refined = refine_problem1_root_by_bisection(
                input,
                launch_state,
                extremum_phi,
                extremum_result,
                right_phi,
                right_result,
                transfer_revolution,
                target_revolution)) {
            refined.push_back(*right_refined);
        }
    }

    return refined;
}

// 几何门控（导数异号 + 内点二次极值）与机制门控（ρ 或 |f_pred|/S）下的 fold 区间检测。
std::optional<Problem1QuadraticExtremumEstimate> detect_fold_interval_by_quadratic_extremum(
    double left_phi,
    const Problem1ResidualResult& left_result,
    double left_derivative,
    double right_phi,
    const Problem1ResidualResult& right_result,
    double right_derivative
) {
    if (!is_success_with_finite_residual(left_result) ||
        !is_success_with_finite_residual(right_result)) {
        return std::nullopt;
    }
    if (residual_sign_changed(left_result.residual, right_result.residual)) {
        return std::nullopt;
    }
    if (!is_finite(left_derivative) || !is_finite(right_derivative) ||
        left_derivative * right_derivative >= 0.0) {
        return std::nullopt;
    }

    const auto quadratic_estimate = estimate_problem1_residual_quadratic_extremum_on_interval(
        left_phi,
        left_result.residual,
        left_derivative,
        right_phi,
        right_result.residual,
        right_derivative);
    if (!quadratic_estimate.has_value()) {
        return std::nullopt;
    }

    const double endpoint_max = std::max(
        std::abs(left_result.residual),
        std::abs(right_result.residual));
    const double relative_depth =
        std::abs(quadratic_estimate->predicted_residual) /
        std::max(endpoint_max, kDefaultEpsilon);

    const double residual_scale = std::max(
        compute_problem1_residual_scale(left_result),
        compute_problem1_residual_scale(right_result));
    if (!is_finite(residual_scale) || !(residual_scale > 0.0)) {
        return std::nullopt;
    }
    const double relative_predicted =
        std::abs(quadratic_estimate->predicted_residual) / residual_scale;

    const bool passes_mechanism_gate =
        relative_depth < kFoldPredictedExtremumRelativeDepthMax ||
        relative_predicted < kFoldPredictedExtremumRelativeResidualMax;
    if (!passes_mechanism_gate) {
        return std::nullopt;
    }

    return quadratic_estimate;
}

// 去重添加候选：同 (k,q) 且遇合角接近时保留相对残差更小的。
void add_candidate_if_not_duplicate(std::vector<Problem1Candidate>* candidates, Problem1Candidate candidate) {
    for (Problem1Candidate& existing : *candidates) {
        if (existing.transfer_revolution != candidate.transfer_revolution ||
            existing.target_revolution != candidate.target_revolution) {
            continue;
        }
        const double angle_diff =
            normalize_angle_minus_pi_pi(existing.encounter_global_angle - candidate.encounter_global_angle);
        if (std::abs(angle_diff) < 1e-8) {
            if (candidate.relative_residual < existing.relative_residual ||
                (candidate.relative_residual == existing.relative_residual &&
                 (candidate.root_bracket_width < existing.root_bracket_width ||
                  (candidate.root_bracket_width == existing.root_bracket_width &&
                   std::abs(candidate.residual_result.residual) < std::abs(existing.residual_result.residual))))) {
                existing = candidate;
            }
            return;
        }
    }
    candidates->push_back(candidate);
}

// 处理一对相邻粗扫样本构成的区间：变号则二分；同号则二次极值门控后 fold 精修。
void process_phi_scan_interval(
    const Problem1SolveInput& input,
    const Problem1LaunchState& launch_state,
    const PhiScanSample& left_sample,
    const PhiScanSample& right_sample,
    double left_derivative,
    double right_derivative,
    int transfer_revolution,
    int target_revolution,
    std::vector<Problem1Candidate>* candidates
) {
    if (!is_success_with_finite_residual(left_sample.result) ||
        !is_success_with_finite_residual(right_sample.result)) {
        return;
    }

    const double left_phi = left_sample.phi;
    const double right_phi = right_sample.phi;
    const Problem1ResidualResult& left_result = left_sample.result;
    const Problem1ResidualResult& right_result = right_sample.result;

    if (std::abs(left_result.residual) <= input.residual_tolerance) {
        auto candidate = make_candidate_from_residual_result(
            input,
            left_phi,
            left_result,
            transfer_revolution,
            target_revolution,
            0.0,
            0,
            false);
        if (candidate.has_value()) {
            add_candidate_if_not_duplicate(candidates, *candidate);
        }
    }
    if (std::abs(right_result.residual) <= input.residual_tolerance) {
        auto candidate = make_candidate_from_residual_result(
            input,
            right_phi,
            right_result,
            transfer_revolution,
            target_revolution,
            0.0,
            0,
            false);
        if (candidate.has_value()) {
            add_candidate_if_not_duplicate(candidates, *candidate);
        }
    }

    if (residual_sign_changed(left_result.residual, right_result.residual)) {
        auto refined = refine_problem1_root_by_bisection(
            input,
            launch_state,
            left_phi,
            left_result,
            right_phi,
            right_result,
            transfer_revolution,
            target_revolution);
        if (refined.has_value()) {
            add_candidate_if_not_duplicate(candidates, *refined);
        }
        return;
    }

    const auto fold_estimate = detect_fold_interval_by_quadratic_extremum(
        left_phi,
        left_result,
        left_derivative,
        right_phi,
        right_result,
        right_derivative);
    if (!fold_estimate.has_value()) {
        return;
    }

    const auto fold_refined = refine_fold_interval_by_quadratic_extremum(
        input,
        launch_state,
        left_phi,
        left_result,
        right_phi,
        right_result,
        *fold_estimate,
        transfer_revolution,
        target_revolution);
    for (const Problem1Candidate& candidate : fold_refined) {
        add_candidate_if_not_duplicate(candidates, candidate);
    }
}

}  // namespace

Problem1ResidualResult evaluate_problem1_residual(const Problem1ResidualInput& input) {
    const Problem1LaunchState launch = compute_problem1_launch_state(
        input.departure_planet,
        input.target_planet,
        input.launch_time_seconds_since_j2000);
    return evaluate_problem1_residual_with_launch_state(input, launch);
}

// Problem 1 求解器：120 等分粗扫，变号区间二分；同号区间经二次极值门控后 fold 精修。
// 枚举所有 (k,q) 多圈分支，返回通过残差过滤的候选转移轨道列表。
std::vector<Problem1Candidate> solve_problem1(const Problem1SolveInput& input) {
    validate_problem1_solve_input(input);

    const Problem1LaunchState launch_state = compute_problem1_launch_state(
        input.departure_planet,
        input.target_planet,
        input.launch_time_seconds_since_j2000);

    std::vector<Problem1Candidate> candidates;
    for (int transfer_revolution = 0; transfer_revolution <= input.max_transfer_revolution; ++transfer_revolution) {
        for (int target_revolution = 0; target_revolution <= input.max_target_revolution; ++target_revolution) {
            std::vector<PhiScanSample> samples;
            samples.reserve(static_cast<std::size_t>(input.phi_scan_count));

            for (int i = 0; i < input.phi_scan_count; ++i) {
                const double phi = static_cast<double>(i) * kTwoPi / static_cast<double>(input.phi_scan_count);
                const Problem1ResidualResult result =
                    evaluate_problem1_residual_with_launch_state(
                        make_residual_input(input, phi, transfer_revolution, target_revolution),
                        launch_state);
                if (!is_success_with_finite_residual(result)) {
                    continue;
                }
                samples.push_back(PhiScanSample{phi, result});
            }

            if (samples.size() < 2) {
                continue;
            }

            for (std::size_t index = 1; index < samples.size(); ++index) {
                process_phi_scan_interval(
                    input,
                    launch_state,
                    samples[index - 1],
                    samples[index],
                    estimate_residual_derivative_at_sample(samples, index - 1),
                    estimate_residual_derivative_at_sample(samples, index),
                    transfer_revolution,
                    target_revolution,
                    &candidates);
            }

            const PhiScanSample& first_sample = samples.front();
            const PhiScanSample& last_sample = samples.back();
            PhiScanSample wrapped_first_sample = first_sample;
            wrapped_first_sample.phi = first_sample.phi + kTwoPi;
            process_phi_scan_interval(
                input,
                launch_state,
                last_sample,
                wrapped_first_sample,
                estimate_residual_derivative_at_sample(samples, samples.size() - 1),
                estimate_residual_derivative_at_sample(samples, 0),
                transfer_revolution,
                target_revolution,
                &candidates);
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const Problem1Candidate& a, const Problem1Candidate& b) {
        if (a.arrival_time_seconds_since_j2000 < b.arrival_time_seconds_since_j2000) {
            return true;
        }
        if (a.arrival_time_seconds_since_j2000 > b.arrival_time_seconds_since_j2000) {
            return false;
        }
        if (a.relative_residual < b.relative_residual) {
            return true;
        }
        if (a.relative_residual > b.relative_residual) {
            return false;
        }
        if (a.root_bracket_width < b.root_bracket_width) {
            return true;
        }
        if (a.root_bracket_width > b.root_bracket_width) {
            return false;
        }
        return std::abs(a.residual_result.residual) < std::abs(b.residual_result.residual);
    });

    return candidates;
}

double estimate_problem1_residual_derivative_central(
    double previous_phi,
    double previous_residual,
    double next_phi,
    double next_residual
) {
    const double delta_phi = next_phi - previous_phi;
    if (!spaceship_cpp::common::is_finite(previous_phi) ||
        !spaceship_cpp::common::is_finite(next_phi) ||
        !spaceship_cpp::common::is_finite(previous_residual) ||
        !spaceship_cpp::common::is_finite(next_residual) ||
        spaceship_cpp::common::near_zero(delta_phi, spaceship_cpp::common::kDefaultEpsilon)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (next_residual - previous_residual) / delta_phi;
}

double estimate_problem1_residual_derivative_at_left_endpoint(
    double left_phi,
    double left_residual,
    double right_phi,
    double right_residual
) {
    const double delta_phi = right_phi - left_phi;
    if (!spaceship_cpp::common::is_finite(left_phi) ||
        !spaceship_cpp::common::is_finite(right_phi) ||
        !spaceship_cpp::common::is_finite(left_residual) ||
        !spaceship_cpp::common::is_finite(right_residual) ||
        !(delta_phi > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (right_residual - left_residual) / delta_phi;
}

double estimate_problem1_residual_derivative_at_right_endpoint(
    double left_phi,
    double left_residual,
    double right_phi,
    double right_residual
) {
    return estimate_problem1_residual_derivative_at_left_endpoint(
        left_phi,
        left_residual,
        right_phi,
        right_residual);
}

std::optional<Problem1QuadraticExtremumEstimate>
estimate_problem1_residual_quadratic_extremum_on_interval(
    double left_phi,
    double left_residual,
    double left_derivative,
    double right_phi,
    double right_residual,
    double right_derivative
) {
    using spaceship_cpp::common::is_finite;
    using spaceship_cpp::common::kDefaultEpsilon;
    using spaceship_cpp::common::near_zero;

    const double interval_width = right_phi - left_phi;
    if (!is_finite(left_phi) || !is_finite(right_phi) || !is_finite(left_residual) ||
        !is_finite(right_residual) || !is_finite(left_derivative) ||
        !is_finite(right_derivative) || !(interval_width > 0.0)) {
        return std::nullopt;
    }

    const double derivative_span = right_derivative - left_derivative;
    if (near_zero(derivative_span, kDefaultEpsilon)) {
        return std::nullopt;
    }

    // 归一化 t in [0,1]：phi = left_phi + t*h；二次模型满足 f_hat(0)=f_L, f_hat'(0)=h*f'_L, f_hat'(1)=h*f'_R。
    // f_hat(t)=a t^2 + b t + c => t* = -f'_L / (f'_R - f'_L)。
    const double normalized_extremum = -left_derivative / derivative_span;
    if (!(normalized_extremum > 0.0 && normalized_extremum < 1.0)) {
        return std::nullopt;
    }

    const double quadratic_a = 0.5 * interval_width * derivative_span;
    const double quadratic_b = interval_width * left_derivative;
    const double quadratic_c = left_residual;
    const double predicted_residual =
        quadratic_a * normalized_extremum * normalized_extremum +
        quadratic_b * normalized_extremum +
        quadratic_c;

    if (!is_finite(predicted_residual)) {
        return std::nullopt;
    }

    return Problem1QuadraticExtremumEstimate{
        left_phi + normalized_extremum * interval_width,
        predicted_residual,
        quadratic_a > 0.0,
    };
}

std::optional<Problem1QuadraticExtremumEstimate> detect_problem1_fold_interval_by_quadratic_extremum(
    double left_phi,
    const Problem1ResidualResult& left_result,
    double left_derivative,
    double right_phi,
    const Problem1ResidualResult& right_result,
    double right_derivative
) {
    return detect_fold_interval_by_quadratic_extremum(
        left_phi,
        left_result,
        left_derivative,
        right_phi,
        right_result,
        right_derivative);
}

std::optional<Problem1PhiRefinementResult> bisect_problem1_residual_on_interval(
    const Problem1SolveInput& input,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    int transfer_revolution,
    int target_revolution
) {
    const Problem1LaunchState launch_state = compute_problem1_launch_state(
        input.departure_planet,
        input.target_planet,
        input.launch_time_seconds_since_j2000);
    return bisect_problem1_residual_on_interval_with_launch_state(
        input,
        launch_state,
        left_phi,
        left_result,
        right_phi,
        right_result,
        transfer_revolution,
        target_revolution);
}

std::optional<Problem1PhiRefinementResult> ternary_search_problem1_residual_extremum_on_interval(
    const Problem1SolveInput& input,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    Problem1ResidualExtremumKind extremum_kind,
    int transfer_revolution,
    int target_revolution,
    int max_iterations
) {
    const Problem1LaunchState launch_state = compute_problem1_launch_state(
        input.departure_planet,
        input.target_planet,
        input.launch_time_seconds_since_j2000);
    return ternary_search_problem1_residual_extremum_on_interval_with_launch_state(
        input,
        launch_state,
        left_phi,
        left_result,
        right_phi,
        right_result,
        extremum_kind,
        transfer_revolution,
        target_revolution,
        max_iterations);
}

}  // namespace spaceship_cpp::problem1
