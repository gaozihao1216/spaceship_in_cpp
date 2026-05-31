#include "spaceship_cpp/common/time_utils.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <cmath>
#include <stdexcept>

namespace spaceship_cpp::common {

namespace {

constexpr int kSecondsPerMinute = 60;
constexpr int kMinutesPerHour = 60;
constexpr int kHoursPerDay = 24;
constexpr int kSecondsPerHour = kSecondsPerMinute * kMinutesPerHour;
constexpr int kSecondsPerDay = kSecondsPerHour * kHoursPerDay;

constexpr long long days_from_civil(int year, unsigned month, unsigned day) noexcept {
    year -= month <= 2;
    const int era = (year >= 0 ? year : year - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(year - era * 400);
    const unsigned doy = (153 * (month + (month > 2 ? -3U : 9U)) + 2) / 5 + day - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return static_cast<long long>(era) * 146097LL + static_cast<long long>(doe) - 719468LL;
}

struct CivilDate {
    int year;
    unsigned month;
    unsigned day;
};

constexpr CivilDate civil_from_days(long long z) noexcept {
    z += 719468LL;
    const long long era = (z >= 0 ? z : z - 146096LL) / 146097LL;
    const unsigned doe = static_cast<unsigned>(z - era * 146097LL);
    const unsigned yoe = (doe - doe / 1460U + doe / 36524U - doe / 146096U) / 365U;
    const int year = static_cast<int>(yoe) + static_cast<int>(era) * 400;
    const unsigned doy = doe - (365U * yoe + yoe / 4U - yoe / 100U);
    const unsigned mp = (5U * doy + 2U) / 153U;
    const unsigned day = doy - (153U * mp + 2U) / 5U + 1U;
    const unsigned month = mp < 10U ? mp + 3U : mp - 9U;
    return CivilDate{year + (month <= 2U), month, day};
}

constexpr bool is_leap_year(int year) noexcept {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

constexpr int days_in_month(int year, int month) noexcept {
    switch (month) {
        case 1:
        case 3:
        case 5:
        case 7:
        case 8:
        case 10:
        case 12:
            return 31;
        case 4:
        case 6:
        case 9:
        case 11:
            return 30;
        case 2:
            return is_leap_year(year) ? 29 : 28;
        default:
            return 0;
    }
}

long long floor_divide(double value, double divisor) {
    return static_cast<long long>(std::floor(value / divisor));
}

}  // namespace

bool is_valid_datetime_utc(const DateTimeUtc& dt) {
    if (dt.month < 1 || dt.month > 12) {
        return false;
    }
    const int month_days = days_in_month(dt.year, dt.month);
    if (dt.day < 1 || dt.day > month_days) {
        return false;
    }
    if (dt.hour < 0 || dt.hour >= kHoursPerDay) {
        return false;
    }
    if (dt.minute < 0 || dt.minute >= kMinutesPerHour) {
        return false;
    }
    if (!is_finite(dt.second) || dt.second < 0.0 || dt.second >= 60.0) {
        return false;
    }
    return true;
}

double seconds_from_j2000(const DateTimeUtc& dt) {
    if (!is_valid_datetime_utc(dt)) {
        throw std::invalid_argument("invalid UTC datetime");
    }

    const long long epoch_days = days_from_civil(2000, 1, 1);
    const long long target_days = days_from_civil(dt.year, static_cast<unsigned>(dt.month), static_cast<unsigned>(dt.day));
    const long long day_delta = target_days - epoch_days;
    const double seconds_of_day =
        static_cast<double>(dt.hour * kSecondsPerHour + dt.minute * kSecondsPerMinute) + dt.second;
    return static_cast<double>(day_delta) * static_cast<double>(kSecondsPerDay) + seconds_of_day - 43200.0;
}

DateTimeUtc datetime_from_j2000_seconds(double seconds_since_j2000) {
    if (!is_finite(seconds_since_j2000)) {
        throw std::invalid_argument("seconds_since_j2000 must be finite");
    }

    const long long epoch_days = days_from_civil(2000, 1, 1);
    const double shifted_seconds = seconds_since_j2000 + 43200.0;
    const long long day_offset = floor_divide(shifted_seconds, static_cast<double>(kSecondsPerDay));
    double seconds_of_day = shifted_seconds - static_cast<double>(day_offset) * static_cast<double>(kSecondsPerDay);

    if (seconds_of_day < 0.0) {
        seconds_of_day += static_cast<double>(kSecondsPerDay);
    }
    if (seconds_of_day >= static_cast<double>(kSecondsPerDay) - 1e-12) {
        seconds_of_day = 0.0;
        const CivilDate next = civil_from_days(epoch_days + day_offset + 1);
        return DateTimeUtc{next.year, static_cast<int>(next.month), static_cast<int>(next.day), 0, 0, 0.0};
    }

    const CivilDate civil = civil_from_days(epoch_days + day_offset);
    const int hour = static_cast<int>(seconds_of_day / static_cast<double>(kSecondsPerHour));
    seconds_of_day -= static_cast<double>(hour * kSecondsPerHour);
    const int minute = static_cast<int>(seconds_of_day / static_cast<double>(kSecondsPerMinute));
    const double second = seconds_of_day - static_cast<double>(minute * kSecondsPerMinute);
    return DateTimeUtc{
        civil.year,
        static_cast<int>(civil.month),
        static_cast<int>(civil.day),
        hour,
        minute,
        second,
    };
}

}  // namespace spaceship_cpp::common
