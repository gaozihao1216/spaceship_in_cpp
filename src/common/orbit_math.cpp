/*
 * 文件作用：实现轨道数学基础函数。
 * 主要工作：计算轨道半径、异常角转换和不同轨道类型下的时间积分。
 */
#include "spaceship_cpp/common/orbit_math.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace spaceship_cpp::common {

namespace {

// 返回 quiet_NaN，作为轨道公式失败时的统一哨兵值。
double nan_value() noexcept {
    return std::numeric_limits<double>::quiet_NaN();
}

// 校验容差参数是否为有限非负数；防止非法 tolerance 进入后续分支判断。
bool is_valid_tolerance(double tolerance) noexcept {
    return is_finite(tolerance) && tolerance >= 0.0;
}

// 校验偏心率 e 和异常角 xi 的基本有效性；统一 try_orbit_F 系列的前置检查。
bool validate_orbit_inputs(double e, double xi, double tolerance, std::string* message) {
    if (!is_finite(e)) {
        if (message != nullptr) {
            *message = "eccentricity must be finite";
        }
        return false;
    }
    if (!is_finite(xi)) {
        if (message != nullptr) {
            *message = "xi must be finite";
        }
        return false;
    }
    if (!is_valid_tolerance(tolerance)) {
        if (message != nullptr) {
            *message = "tolerance must be finite and non-negative";
        }
        return false;
    }
    if (e < 0.0) {
        if (message != nullptr) {
            *message = "eccentricity must be non-negative";
        }
        return false;
    }
    return true;
}

// 将 xi 归约到 [-π, π) 并记录整圈数 revolutions；
// 解决椭圆轨道 F(e,xi) 在大角度下需叠加整圈贡献的问题。
double reduce_angle_minus_pi_pi(double xi, long long* revolutions) noexcept {
    const double shifted_turns = std::floor((xi + kPi) / kTwoPi);
    if (revolutions != nullptr) {
        *revolutions = static_cast<long long>(shifted_turns);
    }
    return xi - shifted_turns * kTwoPi;
}

// 椭圆轨道 (e<1) 的时间积分 F(e,xi) 闭式解；
// 解决转移/目标轨道飞行时间计算的核心被积函数求值。
double elliptic_orbit_F(double e, double xi) {
    const double delta = 1.0 - e * e;
    const double delta_sqrt = std::sqrt(delta);
    const double delta_pow_3_2 = delta * delta_sqrt;

    long long revolutions = 0;
    const double reduced_xi = reduce_angle_minus_pi_pi(xi, &revolutions);
    const double sin_xi = std::sin(reduced_xi);
    const double cos_xi = std::cos(reduced_xi);
    const double denominator = 1.0 + e * cos_xi;
    const double E = std::atan2(delta_sqrt * sin_xi, e + cos_xi);

    const double local_value = E / delta_pow_3_2 - (e * sin_xi) / (delta * denominator);
    const double revolution_value = static_cast<double>(revolutions) * kTwoPi / delta_pow_3_2;
    return local_value + revolution_value;
}

// 抛物线轨道 (e≈1) 的时间积分 F(xi)；
// 处理 e=1 边界情形的专用公式。
double parabolic_orbit_F(double xi) {
    const double tangent = std::tan(0.5 * xi);
    return 0.5 * tangent + tangent * tangent * tangent / 6.0;
}

// 双曲线轨道 (e>1) 的时间积分 F(e,xi)；
// 解决超速转移等开放轨道段的飞行时间计算。
double hyperbolic_orbit_F(double e, double xi) {
    const double q = e * e - 1.0;
    const double sin_xi = std::sin(xi);
    const double cos_xi = std::cos(xi);
    const double denominator = 1.0 + e * cos_xi;
    const double z = std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * xi);
    const double H = 2.0 * std::atanh(z);
    return (e * sin_xi) / (q * denominator) - H / (q * std::sqrt(q));
}

}  // namespace

// 根据偏心率 e 分类轨道类型（椭圆/抛物/双曲）；
// 解决不同圆锥曲线需选用不同 F 公式的问题。
ConicType classify_conic(double e, double eccentricity_tolerance) {
    if (!is_finite(e)) {
        throw std::domain_error("eccentricity must be finite");
    }
    if (!is_valid_tolerance(eccentricity_tolerance)) {
        throw std::domain_error("eccentricity_tolerance must be finite and non-negative");
    }
    if (e < 0.0) {
        throw std::domain_error("eccentricity must be non-negative");
    }

    if (std::abs(e - 1.0) <= eccentricity_tolerance) {
        return ConicType::Parabolic;
    }
    if (e < 1.0) {
        return ConicType::Elliptic;
    }
    return ConicType::Hyperbolic;
}

// 判断 F(e,xi) 在当前 (e,xi) 是否有定义（避开 1+e·cos(xi)=0 奇点）；
// 解决积分前需排除物理不可达角度的问题。
bool is_orbit_F_defined(double e, double xi, double tolerance) {
    std::string message;
    if (!validate_orbit_inputs(e, xi, tolerance, &message)) {
        return false;
    }

    const ConicType conic_type = classify_conic(e, tolerance);
    if (conic_type == ConicType::Parabolic) {
        return std::abs(xi) < kPi - tolerance;
    }
    if (conic_type == ConicType::Hyperbolic) {
        const double xi_limit = std::acos(-1.0 / e);
        if (!(std::abs(xi) < xi_limit - tolerance)) {
            return false;
        }

        const double z = std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * xi);
        return std::abs(z) < 1.0 - tolerance;
    }
    return true;
}

// 轨道时间积分的被积函数 1/(1+e·cos(xi))²；
// 即 dF/dxi，用于牛顿迭代反解真近点角。
double orbit_F_integrand(double e, double xi) {
    if (!is_finite(e)) {
        throw std::domain_error("eccentricity must be finite");
    }
    if (!is_finite(xi)) {
        throw std::domain_error("xi must be finite");
    }
    if (e < 0.0) {
        throw std::domain_error("eccentricity must be non-negative");
    }

    const double denominator = 1.0 + e * std::cos(xi);
    if (!is_finite(denominator) || near_zero(denominator, kDefaultEpsilon)) {
        throw std::domain_error("orbit_F integrand is singular");
    }
    return 1.0 / (denominator * denominator);
}

// 安全计算 F(e,xi)：失败时返回 OrbitFResult{ok=false, message} 而非抛异常；
// 解决调用方需在无效分支上继续枚举而非中断的问题。
OrbitFResult try_orbit_F(double e, double xi, double tolerance) {
    std::string message;
    if (!validate_orbit_inputs(e, xi, tolerance, &message)) {
        return OrbitFResult{false, nan_value(), message};
    }

    const ConicType conic_type = classify_conic(e, tolerance);
    if (conic_type == ConicType::Parabolic) {
        if (std::abs(xi) >= kPi - tolerance) {
            return OrbitFResult{false, nan_value(), "parabolic orbit_F requires |xi| < pi"};
        }
        return OrbitFResult{true, parabolic_orbit_F(xi), ""};
    }

    if (conic_type == ConicType::Hyperbolic) {
        const double xi_limit = std::acos(-1.0 / e);
        if (std::abs(xi) >= xi_limit - tolerance) {
            return OrbitFResult{false, nan_value(), "hyperbolic orbit_F requires |xi| < acos(-1/e)"};
        }

        const double z = std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * xi);
        if (std::abs(z) >= 1.0 - tolerance) {
            return OrbitFResult{false, nan_value(), "hyperbolic orbit_F requires |z| < 1"};
        }
        return OrbitFResult{true, hyperbolic_orbit_F(e, xi), ""};
    }

    return OrbitFResult{true, elliptic_orbit_F(e, xi), ""};
}

// 直接计算 F(e,xi)，输入非法时抛 domain_error；
// 适合调用方已保证输入有效的热路径。
double orbit_F(double e, double xi, double tolerance) {
    const OrbitFResult result = try_orbit_F(e, xi, tolerance);
    if (!result.ok) {
        throw std::domain_error(result.message);
    }
    return result.value;
}

// 安全计算 ∂F/∂xi；用于牛顿法反解真近点角和敏感性分析。
OrbitFResult try_orbit_F_xi_derivative(double e, double xi, double tolerance) {
    const OrbitFResult domain_check = try_orbit_F(e, xi, tolerance);
    if (!domain_check.ok) {
        return OrbitFResult{false, std::numeric_limits<double>::quiet_NaN(), domain_check.message};
    }

    try {
        const double value = orbit_F_integrand(e, xi);
        if (!is_finite(value)) {
            return OrbitFResult{
                false,
                std::numeric_limits<double>::quiet_NaN(),
                "orbit_F_xi_derivative produced non-finite value",
            };
        }
        return OrbitFResult{true, value, ""};
    } catch (const std::exception& ex) {
        return OrbitFResult{false, std::numeric_limits<double>::quiet_NaN(), ex.what()};
    }
}

// 直接计算 ∂F/∂xi，失败时抛异常。
double orbit_F_xi_derivative(double e, double xi, double tolerance) {
    const OrbitFResult result = try_orbit_F_xi_derivative(e, xi, tolerance);
    if (!result.ok) {
        throw std::domain_error(result.message);
    }
    return result.value;
}

// ∂F/∂θ 与 ∂F/∂xi 等价（θ 即 xi），提供语义别名。
OrbitFResult try_orbit_F_theta_derivative(double e, double theta, double tolerance) {
    return try_orbit_F_xi_derivative(e, theta, tolerance);
}

double orbit_F_theta_derivative(double e, double theta, double tolerance) {
    return orbit_F_xi_derivative(e, theta, tolerance);
}

// 安全计算 ∂F/∂e；用于表格诊断中偏心率敏感性分析。
OrbitFResult try_orbit_F_e_derivative(double e, double theta, double tolerance) {
    const OrbitFResult domain_check = try_orbit_F(e, theta, tolerance);
    if (!domain_check.ok) {
        return OrbitFResult{false, std::numeric_limits<double>::quiet_NaN(), domain_check.message};
    }

    try {
        const ConicType conic_type = classify_conic(e, tolerance);
        if (conic_type == ConicType::Parabolic) {
            return OrbitFResult{
                false,
                std::numeric_limits<double>::quiet_NaN(),
                "orbit_F_e_derivative is unsupported for parabolic eccentricity",
            };
        }

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double denominator = 1.0 + e * cos_theta;
        if (!is_finite(denominator) || near_zero(denominator, kDefaultEpsilon)) {
            return OrbitFResult{
                false,
                std::numeric_limits<double>::quiet_NaN(),
                "orbit_F_e_derivative has singular denominator",
            };
        }

        const double F_value = orbit_F(e, theta, tolerance);
        double value = std::numeric_limits<double>::quiet_NaN();
        if (conic_type == ConicType::Elliptic) {
            const double one_minus_e2 = 1.0 - e * e;
            if (!(one_minus_e2 > tolerance)) {
                return OrbitFResult{
                    false,
                    std::numeric_limits<double>::quiet_NaN(),
                    "orbit_F_e_derivative is unstable near e = 1",
                };
            }
            value =
                3.0 * e * F_value / one_minus_e2 -
                sin_theta * (2.0 + e * cos_theta) /
                    (one_minus_e2 * denominator * denominator);
        } else {
            const double e2_minus_one = e * e - 1.0;
            if (!(e2_minus_one > tolerance)) {
                return OrbitFResult{
                    false,
                    std::numeric_limits<double>::quiet_NaN(),
                    "orbit_F_e_derivative is unstable near e = 1",
                };
            }
            value =
                sin_theta * (2.0 + e * cos_theta) /
                    (e2_minus_one * denominator * denominator) -
                3.0 * e * F_value / e2_minus_one;
        }

        if (!is_finite(value)) {
            return OrbitFResult{
                false,
                std::numeric_limits<double>::quiet_NaN(),
                "orbit_F_e_derivative produced non-finite value",
            };
        }
        return OrbitFResult{true, value, ""};
    } catch (const std::exception& ex) {
        return OrbitFResult{false, std::numeric_limits<double>::quiet_NaN(), ex.what()};
    }
}

// 直接计算 ∂F/∂e，失败时抛异常。
double orbit_F_e_derivative(double e, double theta, double tolerance) {
    const OrbitFResult result = try_orbit_F_e_derivative(e, theta, tolerance);
    if (!result.ok) {
        throw std::domain_error(result.message);
    }
    return result.value;
}

}  // namespace spaceship_cpp::common
