/*
 * 文件作用：声明项目通用数学常量和数值辅助函数。
 * 主要工作：提供角度归一化、有限性判断、近零判断和轨道异常角辅助计算。
 */
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>

namespace spaceship_cpp::common {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kHalfPi = 0.5 * kPi;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kRadToDeg = 180.0 / kPi;

constexpr double kDefaultEpsilon = 1e-12;

constexpr double square(double x) noexcept {
    return x * x;
}

inline bool is_finite(double x) noexcept {
    return std::isfinite(x);
}

inline bool near_zero(double x, double eps = kDefaultEpsilon) noexcept {
    return std::abs(x) <= eps;
}

inline bool nearly_equal(double a, double b, double eps = kDefaultEpsilon) noexcept {
    return std::abs(a - b) <= eps;
}

constexpr double deg_to_rad(double degrees) noexcept {
    return degrees * kDegToRad;
}

constexpr double rad_to_deg(double radians) noexcept {
    return radians * kRadToDeg;
}

double normalize_angle_0_2pi(double angle) noexcept;

double normalize_angle_minus_pi_pi(double angle) noexcept;

double positive_mod(double x, double period) noexcept;

double clamp(double x, double lo, double hi);

double safe_divide(double numerator, double denominator, double eps = kDefaultEpsilon);

}  // namespace spaceship_cpp::common
