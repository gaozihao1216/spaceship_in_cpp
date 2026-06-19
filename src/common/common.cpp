/*
 * 文件作用：实现项目通用数学工具。
 * 主要工作：提供角度归一化、数值判断和异常角辅助函数的具体实现。
 */
#include "spaceship_cpp/common/common.hpp"

namespace spaceship_cpp::common {

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
