/*
 * 文件作用：测试 Problem 1 残差评估。
 * 主要工作：检查有效输入、非法分支和轨道几何边界下的残差状态。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

// 中文说明：提供固定绝对容差的浮点比较，用于转移圆锥几何反推验证。
bool approx_equal(double a, double b, double eps = 1e-10) {
    return std::abs(a - b) <= eps;
}

// 中文说明：提供绝对/相对容差浮点比较，用于行星半径量级下的重建半径断言。
bool approx_equal_rel(double a, double b, double abs_tol, double rel_tol) {
    const double diff = std::abs(a - b);
    if (diff <= abs_tol) {
        return true;
    }
    return diff <= rel_tol * std::max(std::abs(a), std::abs(b));
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;

    {
        // 中文说明：验证由两点 (r,xi) 反推的转移离心率 e 与半通径 p 能重建原始半径。
        const double r1 = 1.0;
        const double r2 = 1.2;
        const double xi1 = 0.3;
        const double xi2 = 1.1;
        const double e = problem1::compute_transfer_e_from_two_points(r1, xi1, r2, xi2);
        const double p = problem1::compute_transfer_p_from_departure(r1, e, xi1);

        const double r1_rebuilt = p / (1.0 + e * std::cos(xi1));
        const double r2_rebuilt = p / (1.0 + e * std::cos(xi2));
        assert(approx_equal(r1_rebuilt, r1, 1e-10));
        assert(approx_equal(r2_rebuilt, r2, 1e-10));
    }

    {
        // 中文说明：验证出发/到达角相同（退化几何）时 compute_transfer_e_from_two_points 抛出 domain_error。
        bool threw = false;
        try {
            (void)problem1::compute_transfer_e_from_two_points(1.0, 0.0, 1.0, 0.0);
        } catch (const std::domain_error&) {
            threw = true;
        }
        assert(threw);
    }

    {
        // 中文说明：验证近日点角 phi0 未归一化到 [0,2pi) 时 e 取负值，但 p 与几何仍自洽（符号约定）。
        const double e_true = 0.4;
        const double p_true = 2.0;
        const double phi0_true = spaceship_cpp::common::kPi;
        const double raw_phi0 = phi0_true - spaceship_cpp::common::kPi;
        const double Phi1 = 0.2;
        const double Phi2 = 1.0;

        const double xi1_true = Phi1 - phi0_true;
        const double xi2_true = Phi2 - phi0_true;
        const double r1 = p_true / (1.0 + e_true * std::cos(xi1_true));
        const double r2 = p_true / (1.0 + e_true * std::cos(xi2_true));

        const double xi1_raw = Phi1 - raw_phi0;
        const double xi2_raw = Phi2 - raw_phi0;
        const double e_raw = problem1::compute_transfer_e_from_two_points(r1, xi1_raw, r2, xi2_raw);
        assert(e_raw < 0.0);
        assert(approx_equal(e_raw, -e_true, 1e-10));

        const double p_raw = problem1::compute_transfer_p_from_departure(r1, e_raw, xi1_raw);
        assert(approx_equal(p_raw, p_true, 1e-10));
    }

    {
        // 中文说明：验证 Earth→Mars 有效输入下残差成功时，转移几何、ΔF 与 residual=TOF_transfer-TOF_target 自洽。
        const problem1::Problem1ResidualInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            1.0,
            0,
            0,
        };
        const problem1::Problem1ResidualResult result = problem1::evaluate_problem1_residual(input);
        if (result.status == problem1::Problem1ResidualStatus::Success) {
            assert(std::isfinite(result.residual));
            assert(result.transfer_p > 0.0);
            assert(result.transfer_e >= 0.0);
            assert(std::isfinite(result.transfer_e_raw));
            assert(std::isfinite(result.transfer_perihelion_angle_used));
            assert(result.deltaF_transfer > 0.0);
            assert(result.deltaF_target > 0.0);
            assert(approx_equal(
                result.residual,
                result.transfer_time_scale_free - result.target_time_scale_free,
                1e-8));
            const double r1_rebuilt = result.transfer_p / (1.0 + result.transfer_e * std::cos(result.xi1));
            const double r2_rebuilt = result.transfer_p / (1.0 + result.transfer_e * std::cos(result.xi2));
            assert(approx_equal_rel(r1_rebuilt, result.r1, 1e-6, 1e-12));
            assert(approx_equal_rel(r2_rebuilt, result.r2, 1e-6, 1e-12));
        } else {
            assert(result.status != problem1::Problem1ResidualStatus::InvalidInput);
        }
    }

    {
        // 中文说明：验证目标全局角为 NaN 时 evaluate_problem1_residual 返回 InvalidInput 而非未定义行为。
        const problem1::Problem1ResidualInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            std::numeric_limits<double>::quiet_NaN(),
            0,
            0,
        };
        const problem1::Problem1ResidualResult result = problem1::evaluate_problem1_residual(input);
        assert(result.status == problem1::Problem1ResidualStatus::InvalidInput);
    }

    return 0;
}
