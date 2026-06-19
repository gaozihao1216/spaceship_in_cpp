/*
 * 文件作用：Problem 1 求解热点剖析，验证时间残差 / orbit_F 是否占主导耗时。
 * 主要工作：对比完整求解、残差评估、几何段、行星查询、orbit_F 各段耗时占比。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
namespace common = spaceship_cpp::common;
namespace config = spaceship_cpp::config;
namespace planet_params = spaceship_cpp::planet_params;
namespace problem1 = spaceship_cpp::problem1;

struct TimingStats {
    double mean_us = 0.0;
    double total_ms = 0.0;
    int samples = 0;
};

template <typename Fn>
TimingStats bench(Fn&& fn, int warmup, int iterations) {
    for (int i = 0; i < warmup; ++i) {
        fn();
    }
    std::vector<double> samples_us;
    samples_us.reserve(static_cast<std::size_t>(iterations));
    for (int i = 0; i < iterations; ++i) {
        const auto t0 = Clock::now();
        fn();
        const auto t1 = Clock::now();
        samples_us.push_back(
            std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    TimingStats stats{};
    stats.samples = iterations;
    stats.mean_us =
        std::accumulate(samples_us.begin(), samples_us.end(), 0.0) /
        static_cast<double>(iterations);
    stats.total_ms = stats.mean_us * static_cast<double>(iterations) / 1000.0;
    return stats;
}

struct ResidualGeometrySample {
    planet_params::PlanetId departure{};
    planet_params::PlanetId target{};
    double launch_time = 0.0;
    double transfer_theta = 0.0;
    double phi = 0.0;
    int k = 0;
    int q = 0;

    double e_transfer = 0.0;
    double e_target = 0.0;
    double xi1 = 0.0;
    double xi2 = 0.0;
    double theta2_start = 0.0;
    double theta2_end = 0.0;
    bool valid = false;
};

// 复刻 evaluate_problem1_residual 的几何段（不含 orbit_F），用于分段计时。
ResidualGeometrySample extract_geometry_sample(
    const problem1::Problem1ResidualInput& input
) {
    ResidualGeometrySample sample{};
    sample.departure = input.departure_planet;
    sample.target = input.target_planet;
    sample.launch_time = input.launch_time_seconds_since_j2000;
    sample.transfer_theta = input.transfer_perihelion_angle;
    sample.phi = input.encounter_global_angle;
    sample.k = input.transfer_revolution;
    sample.q = input.target_revolution;

    const planet_params::PlanetParams& target =
        planet_params::get_planet_params(input.target_planet);
    const planet_params::PlanetState state1 =
        planet_params::planet_state_at_time(
            input.departure_planet, input.launch_time_seconds_since_j2000);
    const double r1 = state1.radius;
    const double lambda1 = state1.theta_global;
    const double theta2_start =
        planet_params::planet_true_anomaly_at_time(
            input.target_planet, input.launch_time_seconds_since_j2000);

    const double u2_base = common::normalize_angle_0_2pi(
        input.encounter_global_angle - target.orbit.theta_0);
    const double r2 = planet_params::planet_radius_at_true_anomaly(
        input.target_planet, u2_base);

    const double raw_phi0 = input.transfer_perihelion_angle;
    const double xi1_raw_base = common::normalize_angle_0_2pi(lambda1 - raw_phi0);
    const double xi2_raw_base =
        common::normalize_angle_0_2pi(input.encounter_global_angle - raw_phi0);

    double e_raw = 0.0;
    double p_transfer = 0.0;
    try {
        e_raw = problem1::compute_transfer_e_from_two_points(
            r1, xi1_raw_base, r2, xi2_raw_base);
        p_transfer = problem1::compute_transfer_p_from_departure(
            r1, e_raw, xi1_raw_base);
    } catch (const std::domain_error&) {
        return sample;
    }

    double e_transfer = e_raw;
    double phi0_used = raw_phi0;
    if (e_transfer < 0.0) {
        e_transfer = -e_transfer;
        phi0_used = raw_phi0 + common::kPi;
    }
    phi0_used = common::normalize_angle_0_2pi(phi0_used);

    double xi1_geometry = xi1_raw_base;
    double xi2_geometry = xi2_raw_base;
    if (e_raw < 0.0) {
        xi1_geometry = common::normalize_angle_0_2pi(lambda1 - phi0_used);
        xi2_geometry =
            common::normalize_angle_0_2pi(input.encounter_global_angle - phi0_used);
    }

    double theta2_end = u2_base;
    while (theta2_end <= theta2_start) {
        theta2_end += common::kTwoPi;
    }
    theta2_end += static_cast<double>(input.target_revolution) * common::kTwoPi;

    double xi1 = 0.0;
    double xi2 = 0.0;
    if (e_transfer < 1.0) {
        xi1 = xi1_geometry;
        xi2 = xi2_geometry;
        while (xi2 <= xi1) {
            xi2 += common::kTwoPi;
        }
        xi2 += static_cast<double>(input.transfer_revolution) * common::kTwoPi;
    } else {
        if (input.transfer_revolution != 0) {
            return sample;
        }
        xi1 = common::normalize_angle_minus_pi_pi(lambda1 - phi0_used);
        xi2 = common::normalize_angle_minus_pi_pi(
            input.encounter_global_angle - phi0_used);
        if (xi2 <= xi1) {
            return sample;
        }
    }

    sample.e_transfer = e_transfer;
    sample.e_target = target.orbit.e;
    sample.xi1 = xi1;
    sample.xi2 = xi2;
    sample.theta2_start = theta2_start;
    sample.theta2_end = theta2_end;
    sample.valid = true;
    return sample;
}

problem1::Problem1SolveInput make_input(
    planet_params::PlanetId dep,
    planet_params::PlanetId tgt,
    int max_k,
    int max_q,
    const config::Problem1SolveDefaults& defaults
) {
    auto input = config::make_problem1_solve_input(
        dep, tgt, 0.0, 0.5, defaults);
    input.max_transfer_revolution = max_k;
    input.max_target_revolution = max_q;
    return input;
}

void print_stats(const char* label, const TimingStats& stats, double pct = -1.0) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  " << label << ": mean=" << stats.mean_us << " us";
    if (pct >= 0.0) {
        std::cout << " (" << pct << "%)";
    }
    std::cout << ", n=" << stats.samples << '\n';
}

void profile_pair(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    int max_k,
    int max_q,
    const config::Problem1SolveDefaults& defaults
) {
    const auto input = make_input(departure, target, max_k, max_q, defaults);
    const auto solve_input = input;

    std::cout << "\n=== "
              << planet_params::planet_name(departure) << " -> "
              << planet_params::planet_name(target)
              << " (max_k=" << max_k << ", max_q=" << max_q
              << ", phi_scan=" << defaults.phi_scan_count << ") ===\n";

    // 找一个成功残差样本用于分段基准。
    ResidualGeometrySample geom{};
    problem1::Problem1ResidualInput residual_input{};
    for (int i = 0; i < defaults.phi_scan_count; ++i) {
        const double phi =
            static_cast<double>(i) * common::kTwoPi /
            static_cast<double>(defaults.phi_scan_count);
        residual_input = problem1::Problem1ResidualInput{
            departure,
            target,
            0.0,
            0.5,
            phi,
            0,
            0,
        };
        const auto result = problem1::evaluate_problem1_residual(residual_input);
        if (result.status == problem1::Problem1ResidualStatus::Success) {
            geom = extract_geometry_sample(residual_input);
            break;
        }
    }
    if (!geom.valid) {
        std::cout << "  warning: no successful residual sample found\n";
        return;
    }

    constexpr int kResidualIterations = 50000;
    constexpr int kOrbitFIterations = 200000;
    constexpr int kGeometryIterations = 50000;
    constexpr int kPlanetIterations = 50000;
    constexpr int kSolveIterations = 20;

    volatile double sink = 0.0;

    const TimingStats full_residual = bench(
        [&]() {
            const auto r = problem1::evaluate_problem1_residual(residual_input);
            sink += r.residual;
        },
        200,
        kResidualIterations);

    const TimingStats geometry_only = bench(
        [&]() {
            const auto g = extract_geometry_sample(residual_input);
            sink += g.xi2;
        },
        200,
        kGeometryIterations);

    const TimingStats planet_only = bench(
        [&]() {
            const auto s = planet_params::planet_state_at_time(
                departure, 0.0);
            const double t = planet_params::planet_true_anomaly_at_time(
                target, 0.0);
            const double r = planet_params::planet_radius_at_true_anomaly(
                target, t);
            sink += s.radius + r;
        },
        200,
        kPlanetIterations);

    const TimingStats inverse_orbit_f = bench(
        [&]() {
            const double t1 = planet_params::planet_true_anomaly_at_time(
                departure, 0.0);
            const double t2 = planet_params::planet_true_anomaly_at_time(
                target, 0.0);
            sink += t1 + t2;
        },
        200,
        kPlanetIterations);

    const TimingStats direct_orbit_f_four = bench(
        [&]() {
            const double d1 =
                common::orbit_F(geom.e_transfer, geom.xi2) -
                common::orbit_F(geom.e_transfer, geom.xi1);
            const double d2 =
                common::orbit_F(geom.e_target, geom.theta2_end) -
                common::orbit_F(geom.e_target, geom.theta2_start);
            sink += d1 + d2;
        },
        500,
        kOrbitFIterations);

    const TimingStats full_solve = bench(
        [&]() {
            const auto c = problem1::solve_problem1(solve_input);
            sink += static_cast<double>(c.size());
        },
        2,
        kSolveIterations);

    // 按单次残差评估等效归一化到同一基准。
    const double orbit_f_per_residual_us =
        direct_orbit_f_four.mean_us * static_cast<double>(kResidualIterations) /
        static_cast<double>(kOrbitFIterations);
    const double inverse_f_per_residual_us =
        inverse_orbit_f.mean_us * static_cast<double>(kResidualIterations) /
        static_cast<double>(kPlanetIterations);
    const double geometry_per_residual_us =
        geometry_only.mean_us * static_cast<double>(kResidualIterations) /
        static_cast<double>(kGeometryIterations);
    const double planet_per_residual_us =
        planet_only.mean_us * static_cast<double>(kResidualIterations) /
        static_cast<double>(kPlanetIterations);

    const double accounted_us =
        orbit_f_per_residual_us + geometry_per_residual_us + planet_per_residual_us;
    const double other_us = std::max(0.0, full_residual.mean_us - accounted_us);

    std::cout << "Segment breakdown (per evaluate_problem1_residual call):\n";
    print_stats(
        "full residual",
        full_residual,
        100.0);
    print_stats(
        "orbit_F delta in residual (4 direct calls)",
        TimingStats{orbit_f_per_residual_us, 0.0, 0},
        100.0 * orbit_f_per_residual_us / full_residual.mean_us);
    print_stats(
        "planet_true_anomaly inverse orbit_F (2 Newton solves)",
        TimingStats{inverse_f_per_residual_us, 0.0, 0},
        100.0 * inverse_f_per_residual_us / full_residual.mean_us);
    print_stats(
        "geometry (no orbit_F)",
        TimingStats{geometry_per_residual_us, 0.0, 0},
        100.0 * geometry_per_residual_us / full_residual.mean_us);
    print_stats(
        "planet state queries (3 calls)",
        TimingStats{planet_per_residual_us, 0.0, 0},
        100.0 * planet_per_residual_us / full_residual.mean_us);
    print_stats(
        "other/overhead (residual - above)",
        TimingStats{other_us, 0.0, 0},
        100.0 * other_us / full_residual.mean_us);

    const int branches =
        (max_k + 1) * (max_q + 1);
    const int scan_evals = branches * defaults.phi_scan_count;
    const double est_scan_ms =
        full_residual.mean_us * static_cast<double>(scan_evals) / 1000.0;

    const auto candidates = problem1::solve_problem1(solve_input);
    std::cout << "\nSolve summary:\n";
    std::cout << "  candidates found: " << candidates.size() << '\n';
    std::cout << "  min scan residual evals: " << scan_evals << '\n';
    std::cout << "  est. scan-only time: " << est_scan_ms << " ms\n";
    std::cout << "  measured full solve: mean="
              << full_solve.mean_us / 1000.0 << " ms"
              << " (n=" << full_solve.samples << ")\n";
    std::cout << "  scan fraction of solve (estimate): "
              << 100.0 * est_scan_ms / (full_solve.mean_us / 1000.0) << "%\n";

    std::cout << "\nKey sample orbit params for orbit_F:\n";
    std::cout << "  e_transfer=" << geom.e_transfer
              << ", xi1=" << geom.xi1 << ", xi2=" << geom.xi2 << '\n';
    std::cout << "  e_target=" << geom.e_target
              << ", theta_start=" << geom.theta2_start
              << ", theta_end=" << geom.theta2_end << '\n';
}

}  // namespace

int main() {
    const config::GlobalConfig& cfg = config::global_config();
    const auto& defaults = cfg.problem1_solve;

    std::cout << "=== Problem 1 profiling ===\n";
    std::cout << "Note: orbit_F is a closed-form time integral (not numerical quadrature),\n";
    std::cout << "      but each residual calls it 4 times (2 orbit segments x 2 endpoints).\n";

    profile_pair(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        1, 1, defaults);
    profile_pair(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        2, 2, defaults);
    profile_pair(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mercury,
        1, 1, defaults);

    return 0;
}
