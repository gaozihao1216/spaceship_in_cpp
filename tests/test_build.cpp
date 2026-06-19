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

    // 中文说明：验证 compile-time 常量工具（square、deg/rad 转换）在编译期即可正确求值。
    static_assert(common::square(3.0) == 9.0);
    static_assert(common::deg_to_rad(180.0) == common::kPi);
    static_assert(common::rad_to_deg(common::kHalfPi) == 90.0);

    // 中文说明：验证 is_finite / near_zero / nearly_equal 能区分有限数、无穷大和近零/近相等浮点。
    assert(common::is_finite(1.0));
    assert(!common::is_finite(std::numeric_limits<double>::infinity()));
    assert(common::near_zero(1e-13));
    assert(common::nearly_equal(1.0, 1.0 + 5e-13));

    // 中文说明：验证 positive_mod 对正/负被除数给出 [0,mod) 余数，并在 mod=0 时返回 NaN 而非除零崩溃。
    assert(common::nearly_equal(common::positive_mod(5.5, 2.0), 1.5));
    assert(common::nearly_equal(common::positive_mod(-0.5, 2.0), 1.5));
    assert(std::isnan(common::positive_mod(1.0, 0.0)));

    // 中文说明：验证 normalize_angle_0_2pi 将任意角归一化到 [0,2pi)，含 0、2pi 与负角边界。
    assert(common::near_zero(common::normalize_angle_0_2pi(0.0)));
    assert(common::near_zero(common::normalize_angle_0_2pi(common::kTwoPi)));
    assert(common::nearly_equal(common::normalize_angle_0_2pi(-common::kHalfPi), 1.5 * common::kPi));

    // 中文说明：验证 normalize_angle_minus_pi_pi 将角差归一化到 [-pi,pi)，用于最短方向比较。
    assert(common::nearly_equal(common::normalize_angle_minus_pi_pi(common::kPi), -common::kPi));
    assert(common::nearly_equal(common::normalize_angle_minus_pi_pi(1.5 * common::kPi), -common::kHalfPi));
    assert(common::nearly_equal(common::normalize_angle_minus_pi_pi(-1.5 * common::kPi), common::kHalfPi));

    // 中文说明：验证 clamp 在区间内/外正确夹取，并在上下界颠倒时抛出 invalid_argument。
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

    // 中文说明：验证 safe_divide 在正常除法下返回正确商，并在分母过小时拒绝除零。
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
