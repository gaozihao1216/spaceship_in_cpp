#pragma once

namespace spaceship_cpp::common {

struct DateTimeUtc {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
};

constexpr const char* kJ2000EpochIso = "2000-01-01T12:00:00Z";
constexpr double kJ2000JulianDate = 2451545.0;

double seconds_from_j2000(const DateTimeUtc& dt);

DateTimeUtc datetime_from_j2000_seconds(double seconds_since_j2000);

bool is_valid_datetime_utc(const DateTimeUtc& dt);

}  // namespace spaceship_cpp::common
