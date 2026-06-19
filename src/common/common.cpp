/*
 * 文件作用：实现项目通用数学工具。
 * 主要工作：提供角度归一化、数值判断和异常角辅助函数的具体实现。
 */
#include "spaceship_cpp/common/common.hpp"

namespace spaceship_cpp::common {

// 将 x 映射到 [0, period) 区间；解决 C++ fmod 对负数返回负余数、
// 导致周期角/时间取模结果不一致的问题。
double positive_mod(double x, double period) noexcept {
    if (!is_finite(x) || !is_finite(period) || period <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double result = std::fmod(x, period);
    if (result < 0.0) {
        return result + period;
    }
    if (result >= period) {
        return 0.0;
    }
    return result;
}

// 将任意角归一化到 [0, 2π)；解决全局方向角跨多圈后数值过大、
// 不便存储和比较的问题。
double normalize_angle_0_2pi(double angle) noexcept {
    if (!is_finite(angle)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double normalized = positive_mod(angle, kTwoPi);
    if (!is_finite(normalized)) {
        return normalized;
    }

    return near_zero(normalized, kDefaultEpsilon) || nearly_equal(normalized, kTwoPi, kDefaultEpsilon)
        ? 0.0
        : normalized;
}

// 将角差归一化到 [-π, π)；解决比较两个方向角时，
// 应取最短旋转路径而非整圈差值的问题。
double normalize_angle_minus_pi_pi(double angle) noexcept {
    if (!is_finite(angle)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double shifted = normalize_angle_0_2pi(angle + kPi);
    if (!is_finite(shifted)) {
        return shifted;
    }

    const double normalized = shifted - kPi;
    return nearly_equal(normalized, kPi, kDefaultEpsilon) ? -kPi : normalized;
}

// 将 x 限制在 [lo, hi] 内；解决插值/搜索中数值越界，
// 并在上下界写反时主动抛错而非静默返回错误结果。
double clamp(double x, double lo, double hi) {
    if (!is_finite(x) || !is_finite(lo) || !is_finite(hi)) {
        throw std::invalid_argument("clamp requires finite inputs");
    }
    if (lo > hi) {
        throw std::invalid_argument("clamp requires lo <= hi");
    }
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

// 带分母阈值检查的除法；解决轨道几何退化（分母接近 0）时
// 直接除法产生 Inf/NaN 并污染后续计算的问题。
double safe_divide(double numerator, double denominator, double eps) {
    if (!is_finite(numerator) || !is_finite(denominator) || !is_finite(eps)) {
        throw std::invalid_argument("safe_divide requires finite inputs");
    }
    if (eps < 0.0) {
        throw std::invalid_argument("safe_divide requires non-negative epsilon");
    }
    if (std::abs(denominator) <= eps) {
        throw std::invalid_argument("safe_divide denominator too close to zero");
    }
    return numerator / denominator;
}

}  // namespace spaceship_cpp::common
