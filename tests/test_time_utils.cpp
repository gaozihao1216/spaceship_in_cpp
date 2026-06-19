/*
 * 文件作用：测试时间转换工具。
 * 主要工作：验证天与秒之间的转换精度和符号处理。
 */
#include "spaceship_cpp/common/time_utils.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {

// 中文说明：提供固定绝对容差的浮点比较，用于 UTC 秒字段往返精度断言。
bool approx_equal(double a, double b, double eps = 1e-12) {
    return std::abs(a - b) <= eps;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;

    // 中文说明：验证 J2000 历元及前后各一天相对 J2000 的秒偏移（0、±86400）。
    assert(common::seconds_from_j2000({2000, 1, 1, 12, 0, 0.0}) == 0.0);
    assert(common::seconds_from_j2000({2000, 1, 2, 12, 0, 0.0}) == 86400.0);
    assert(common::seconds_from_j2000({1999, 12, 31, 12, 0, 0.0}) == -86400.0);

    {
        // 中文说明：验证 datetime_from_j2000_seconds(0) 精确还原 J2000 历元时刻。
        const common::DateTimeUtc dt = common::datetime_from_j2000_seconds(0.0);
        assert(dt.year == 2000 && dt.month == 1 && dt.day == 1);
        assert(dt.hour == 12 && dt.minute == 0);
        assert(approx_equal(dt.second, 0.0));
    }

    {
        // 中文说明：验证正向偏移一天后日期字段正确递增。
        const common::DateTimeUtc dt = common::datetime_from_j2000_seconds(86400.0);
        assert(dt.year == 2000 && dt.month == 1 && dt.day == 2);
        assert(dt.hour == 12 && dt.minute == 0);
        assert(approx_equal(dt.second, 0.0));
    }

    {
        // 中文说明：验证 datetime → seconds → datetime 往返在亚秒精度内保持一致。
        const common::DateTimeUtc original{2000, 1, 1, 12, 0, 1.25};
        const double seconds = common::seconds_from_j2000(original);
        const common::DateTimeUtc roundtrip = common::datetime_from_j2000_seconds(seconds);
        assert(roundtrip.year == original.year);
        assert(roundtrip.month == original.month);
        assert(roundtrip.day == original.day);
        assert(roundtrip.hour == original.hour);
        assert(roundtrip.minute == original.minute);
        assert(approx_equal(roundtrip.second, original.second, 1e-9));
    }

    {
        // 中文说明：验证历元前亚秒负偏移能正确映射到 11:59:59.x 并往返一致。
        const common::DateTimeUtc dt = common::datetime_from_j2000_seconds(-0.25);
        assert(dt.year == 2000 && dt.month == 1 && dt.day == 1);
        assert(dt.hour == 11 && dt.minute == 59);
        assert(approx_equal(dt.second, 59.75, 1e-9));
        assert(approx_equal(common::seconds_from_j2000(dt), -0.25, 1e-9));
    }

    {
        // 中文说明：验证跨日大负偏移（约 -12 小时）的日期/时分秒分解与往返一致性。
        const common::DateTimeUtc dt = common::datetime_from_j2000_seconds(-43200.25);
        assert(dt.year == 1999 && dt.month == 12 && dt.day == 31);
        assert(dt.hour == 23 && dt.minute == 59);
        assert(approx_equal(dt.second, 59.75, 1e-9));
        assert(approx_equal(common::seconds_from_j2000(dt), -43200.25, 1e-9));
    }

    // 中文说明：验证非法月份被拒绝；闰年规则下 2000-02-29 合法而 1900-02-29 非法。
    bool threw = false;
    try {
        (void)common::seconds_from_j2000({2000, 13, 1, 0, 0, 0.0});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    assert(common::is_valid_datetime_utc({2000, 2, 29, 0, 0, 0.0}));
    assert(!common::is_valid_datetime_utc({1900, 2, 29, 0, 0, 0.0}));

    return 0;
}
