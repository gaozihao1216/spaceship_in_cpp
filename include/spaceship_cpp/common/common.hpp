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
// 常用角度常量，统一避免各模块重复定义 pi 和角度转换比例。
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kHalfPi = 0.5 * kPi;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kRadToDeg = 180.0 / kPi;

constexpr double kDefaultEpsilon = 1e-12;

// 小型数值工具：解决公式里反复出现平方运算的问题。
constexpr double square(double x) noexcept {
    return x * x;
}

// 统一浮点有效性判断，避免 NaN/Inf 进入轨道公式后扩散。
inline bool is_finite(double x) noexcept {
    return std::isfinite(x);
}

// 判断数值是否可以按 0 处理，主要用于分母和几何退化检查。
inline bool near_zero(double x, double eps = kDefaultEpsilon) noexcept {
    return std::abs(x) <= eps;
}

// 绝对容差比较，用于测试和边界判断。
inline bool nearly_equal(double a, double b, double eps = kDefaultEpsilon) noexcept {
    return std::abs(a - b) <= eps;
}

// 角度单位转换，解决输入/输出中 degree 与 radian 的换算问题。
constexpr double deg_to_rad(double degrees) noexcept {
    return degrees * kDegToRad;
}

constexpr double rad_to_deg(double radians) noexcept {
    return radians * kRadToDeg;
}

// 将任意角归一化到 [0, 2pi)，适合存储全局方向角。
double normalize_angle_0_2pi(double angle) noexcept;

// 将角差归一化到 [-pi, pi)，适合比较两个方向的最短差值。
double normalize_angle_minus_pi_pi(double angle) noexcept;

// 正模运算，解决 C++ fmod 对负数返回负余数的问题。
double positive_mod(double x, double period) noexcept;

// 带输入校验的夹取函数，防止上界下界写反时静默出错。
double clamp(double x, double lo, double hi);

// 带分母阈值检查的除法，解决退化几何下直接除零的问题。
double safe_divide(double numerator, double denominator, double eps = kDefaultEpsilon);

}  // namespace spaceship_cpp::common
