/*
 * 文件作用：基础构建和链接冒烟测试。
 * 主要工作：验证公共常量和基础工具能被测试程序正常包含与使用。
 */
#include "spaceship_cpp/common/common.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

int main() {
    namespace common = spaceship_cpp::common;

    static_assert(common::square(3.0) == 9.0);
    static_assert(common::deg_to_rad(180.0) == common::kPi);
    static_assert(common::rad_to_deg(common::kHalfPi) == 90.0);

    assert(common::is_finite(1.0));
    assert(!common::is_finite(std::numeric_limits<double>::infinity()));
    assert(common::near_zero(1e-13));
    assert(common::nearly_equal(1.0, 1.0 + 5e-13));

    assert(common::nearly_equal(common::positive_mod(5.5, 2.0), 1.5));
    assert(common::nearly_equal(common::positive_mod(-0.5, 2.0), 1.5));
    assert(std::isnan(common::positive_mod(1.0, 0.0)));

    assert(common::near_zero(common::normalize_angle_0_2pi(0.0)));
    assert(common::near_zero(common::normalize_angle_0_2pi(common::kTwoPi)));
    assert(common::nearly_equal(common::normalize_angle_0_2pi(-common::kHalfPi), 1.5 * common::kPi));

    assert(common::nearly_equal(common::normalize_angle_minus_pi_pi(common::kPi), -common::kPi));
    assert(common::nearly_equal(common::normalize_angle_minus_pi_pi(1.5 * common::kPi), -common::kHalfPi));
    assert(common::nearly_equal(common::normalize_angle_minus_pi_pi(-1.5 * common::kPi), common::kHalfPi));

    assert(common::clamp(-1.0, 0.0, 1.0) == 0.0);
    assert(common::clamp(0.5, 0.0, 1.0) == 0.5);
    assert(common::clamp(2.0, 0.0, 1.0) == 1.0);

    bool threw = false;
    try {
        (void)common::clamp(0.0, 2.0, 1.0);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    assert(common::safe_divide(6.0, 2.0) == 3.0);

    threw = false;
    try {
        (void)common::safe_divide(1.0, 1e-13);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    return 0;
}
