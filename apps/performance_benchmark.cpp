/*
 * 文件作用：Problem 1/2 求解耗时与 BFS 分支扩展的性能评估程序。
 * 主要工作：
 *   - Problem 1：时间匹配转移求解（含 k/q 多圈分支）
 *   - Problem 2：在飞掠点联立弹弓方程 + 外向轨道几何，对 (theta', alpha) 求根
 *   - 模拟 BFS 扩展：首段纯 P1，中间段用 P2 联立求解弹弓后转移轨道
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
namespace common = spaceship_cpp::common;
namespace config = spaceship_cpp::config;
namespace planet_params = spaceship_cpp::planet_params;
namespace problem1 = spaceship_cpp::problem1;
namespace problem2 = spaceship_cpp::problem2;

struct TimingStats {
    int samples = 0;
    double mean_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
};

struct RevLimits {
    int max_k = 0;
    int max_q = 0;
};

struct ExpansionStats {
    planet_params::PlanetId departure{};
    double departure_time = 0.0;
    bool uses_problem2_joint_solve = false;
    int problem1_candidate_count = 0;
    int valid_edge_count = 0;
    double expand_wall_ms = 0.0;
};

constexpr double kLaunchTimeSeconds = 0.0;
constexpr double kFixedTransferTheta = 0.5;
constexpr int kTimedIterations = 10;
constexpr int kProblem2ThetaPrimeScanCount = 60;
constexpr int kProblem2AlphaScanCount = 60;
constexpr double kProblem2ResidualTolerance = 1e-4;

const std::array<planet_params::PlanetId, 4> kInnerPlanets{
    planet_params::PlanetId::Mercury,
    planet_params::PlanetId::Venus,
    planet_params::PlanetId::Earth,
    planet_params::PlanetId::Mars,
};

TimingStats compute_stats(const std::vector<double>& samples_ms) {
    TimingStats stats{};
    if (samples_ms.empty()) {
        return stats;
    }
    stats.samples = static_cast<int>(samples_ms.size());
    stats.mean_ms =
        std::accumulate(samples_ms.begin(), samples_ms.end(), 0.0) /
        static_cast<double>(samples_ms.size());
    stats.min_ms = *std::min_element(samples_ms.begin(), samples_ms.end());
    stats.max_ms = *std::max_element(samples_ms.begin(), samples_ms.end());
    return stats;
}

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

const char* planet_label(planet_params::PlanetId id) {
    return planet_params::planet_name(id);
}

problem1::Problem1SolveInput make_solve_input(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    double launch_time,
    double transfer_theta,
    const config::Problem1SolveDefaults& defaults,
    const RevLimits& rev
) {
    auto input = config::make_problem1_solve_input(
        departure,
        target,
        launch_time,
        transfer_theta,
        defaults);
    input.max_transfer_revolution = rev.max_k;
    input.max_target_revolution = rev.max_q;
    return input;
}

TimingStats benchmark_problem1_solve(
    const problem1::Problem1SolveInput& input,
    int iterations
) {
    for (int i = 0; i < 3; ++i) {
        (void)problem1::solve_problem1(input);
    }

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(iterations));
    for (int i = 0; i < iterations; ++i) {
        const auto start = Clock::now();
        (void)problem1::solve_problem1(input);
        const auto end = Clock::now();
        samples.push_back(elapsed_ms(start, end));
    }
    return compute_stats(samples);
}

struct Problem2SolveResult {
    int solution_count = 0;
    double wall_ms = 0.0;
};

// Problem 2 联立求解：在飞掠行星 J 已知入射轨道时，扫描 (theta', alpha)，
// 用弹弓残差 + 外向两点几何联立，找出满足 slingshot_residual≈0 的弹弓后转移轨道。
Problem2SolveResult solve_problem2_slingshot_transfer(
    planet_params::PlanetId flyby_planet,
    planet_params::PlanetId target_planet,
    double flyby_time,
    double incoming_e,
    double incoming_theta
) {
    const auto& params_j = planet_params::get_planet_params(flyby_planet);
    const auto& params_k = planet_params::get_planet_params(target_planet);
    const double phi = planet_params::planet_true_anomaly_at_time(flyby_planet, flyby_time);

    const auto start = Clock::now();
    int solution_count = 0;
    for (int i = 0; i < kProblem2ThetaPrimeScanCount; ++i) {
        const double theta_prime =
            static_cast<double>(i) * common::kTwoPi /
            static_cast<double>(kProblem2ThetaPrimeScanCount);
        for (int j = 0; j < kProblem2AlphaScanCount; ++j) {
            const double alpha =
                static_cast<double>(j) * common::kTwoPi /
                static_cast<double>(kProblem2AlphaScanCount);
            const auto result = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
                params_j.orbit.p,
                params_j.orbit.e,
                params_k.orbit.p,
                params_k.orbit.e,
                phi,
                alpha,
                incoming_e,
                incoming_theta,
                theta_prime);
            if (result.valid &&
                std::abs(result.slingshot_residual) <= kProblem2ResidualTolerance) {
                ++solution_count;
            }
        }
    }
    return Problem2SolveResult{
        solution_count,
        elapsed_ms(start, Clock::now()),
    };
}

TimingStats benchmark_problem2_joint_solve(
    planet_params::PlanetId flyby_planet,
    planet_params::PlanetId target_planet,
    double flyby_time,
    double incoming_e,
    double incoming_theta,
    int iterations
) {
    for (int i = 0; i < 2; ++i) {
        (void)solve_problem2_slingshot_transfer(
            flyby_planet, target_planet, flyby_time, incoming_e, incoming_theta);
    }

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(iterations));
    for (int i = 0; i < iterations; ++i) {
        const auto result = solve_problem2_slingshot_transfer(
            flyby_planet, target_planet, flyby_time, incoming_e, incoming_theta);
        samples.push_back(result.wall_ms);
    }
    return compute_stats(samples);
}

const problem1::Problem1Candidate* first_candidate(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    double launch_time,
    double transfer_theta,
    const config::Problem1SolveDefaults& defaults,
    const RevLimits& rev
) {
    const auto input = make_solve_input(
        departure, target, launch_time, transfer_theta, defaults, rev);
    const auto candidates = problem1::solve_problem1(input);
    if (candidates.empty()) {
        return nullptr;
    }
    return &candidates.front();
}

ExpansionStats simulate_bfs_expansion(
    planet_params::PlanetId departure,
    double departure_time,
    double incoming_e,
    double incoming_theta,
    bool has_incoming,
    double fixed_transfer_theta,
    const config::Problem1SolveDefaults& defaults,
    const RevLimits& rev
) {
    ExpansionStats stats{};
    stats.departure = departure;
    stats.departure_time = departure_time;
    stats.uses_problem2_joint_solve = has_incoming;

    const auto expand_start = Clock::now();
    for (planet_params::PlanetId target : kInnerPlanets) {
        if (target == departure) {
            continue;
        }

        if (!has_incoming) {
            // 首段：仅 Problem 1 时间匹配转移。
            const auto input = make_solve_input(
                departure,
                target,
                departure_time,
                fixed_transfer_theta,
                defaults,
                rev);
            const auto candidates = problem1::solve_problem1(input);
            stats.problem1_candidate_count += static_cast<int>(candidates.size());
            stats.valid_edge_count += static_cast<int>(candidates.size());
        } else {
            // 中间段：Problem 2 联立弹弓方程求弹弓后外向转移轨道。
            const auto p2 = solve_problem2_slingshot_transfer(
                departure,
                target,
                departure_time,
                incoming_e,
                incoming_theta);
            stats.problem1_candidate_count += 0;
            stats.valid_edge_count += p2.solution_count;
        }
    }
    stats.expand_wall_ms = elapsed_ms(expand_start, Clock::now());
    return stats;
}

void print_timing(const std::string& label, const TimingStats& stats) {
    std::cout << std::fixed << std::setprecision(3);
    std::cout
        << label
        << ": mean=" << stats.mean_ms << " ms"
        << ", min=" << stats.min_ms << " ms"
        << ", max=" << stats.max_ms << " ms"
        << ", n=" << stats.samples << '\n';
}

void print_expansion(const ExpansionStats& stats) {
    std::cout
        << "  expand from " << planet_label(stats.departure)
        << " @ t=" << stats.departure_time << " s"
        << ", mode=" << (stats.uses_problem2_joint_solve ? "P2_joint" : "P1_only")
        << ", p1_candidates=" << stats.problem1_candidate_count
        << ", valid_edges=" << stats.valid_edge_count
        << ", wall_time=" << stats.expand_wall_ms << " ms\n";
}

void run_benchmark_for_rev_limits(
    const config::Problem1SolveDefaults& defaults,
    const RevLimits& rev
) {
    std::cout << "\n######## max_k=" << rev.max_k
              << ", max_q=" << rev.max_q << " ########\n";

    std::cout << "--- Problem 1 full solve ---\n";
    const std::array<std::pair<planet_params::PlanetId, planet_params::PlanetId>, 3> kEarthPairs{
        std::pair{planet_params::PlanetId::Earth, planet_params::PlanetId::Mars},
        std::pair{planet_params::PlanetId::Earth, planet_params::PlanetId::Venus},
        std::pair{planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury},
    };

    std::vector<double> p1_means;
    std::vector<int> p1_counts;
    for (const auto& [departure, target] : kEarthPairs) {
        const auto input = make_solve_input(
            departure, target, kLaunchTimeSeconds, kFixedTransferTheta, defaults, rev);
        const auto candidates = problem1::solve_problem1(input);
        p1_counts.push_back(static_cast<int>(candidates.size()));
        const auto stats = benchmark_problem1_solve(input, kTimedIterations);
        p1_means.push_back(stats.mean_ms);
        std::cout << planet_label(departure) << " -> " << planet_label(target)
                  << ": candidates=" << candidates.size() << ", ";
        print_timing("solve", stats);
    }

    const double avg_p1_ms =
        std::accumulate(p1_means.begin(), p1_means.end(), 0.0) /
        static_cast<double>(p1_means.size());
    const double avg_p1_candidates =
        std::accumulate(p1_counts.begin(), p1_counts.end(), 0.0) /
        static_cast<double>(p1_counts.size());
    std::cout << "P1 average (Earth departures): " << avg_p1_ms
              << " ms, avg candidates/pair=" << avg_p1_candidates << '\n';

    const auto* earth_mars = first_candidate(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        kLaunchTimeSeconds,
        kFixedTransferTheta,
        defaults,
        rev);

    std::cout << "\n--- Problem 2 joint slingshot transfer solve ---\n";
    std::cout << "Scan grid: theta' x alpha = "
              << kProblem2ThetaPrimeScanCount << " x " << kProblem2AlphaScanCount
              << ", tolerance=" << kProblem2ResidualTolerance << '\n';

    if (earth_mars != nullptr) {
        print_timing(
            "P2 joint solve Mars->Venus (incoming from Earth->Mars)",
            benchmark_problem2_joint_solve(
                planet_params::PlanetId::Mars,
                planet_params::PlanetId::Venus,
                earth_mars->arrival_time_seconds_since_j2000,
                earth_mars->residual_result.transfer_e,
                earth_mars->residual_result.transfer_perihelion_angle_used,
                kTimedIterations));

        const auto p2_demo = solve_problem2_slingshot_transfer(
            planet_params::PlanetId::Mars,
            planet_params::PlanetId::Venus,
            earth_mars->arrival_time_seconds_since_j2000,
            earth_mars->residual_result.transfer_e,
            earth_mars->residual_result.transfer_perihelion_angle_used);
        std::cout << "  demo Mars->Venus P2 solutions found: "
                  << p2_demo.solution_count << '\n';
    } else {
        std::cout << "  skipped: no Earth->Mars P1 candidate for P2 demo\n";
    }

    std::cout << "\n--- BFS branch expansion (Earth -> Mercury goal) ---\n";
    std::cout << "Depth 0: P1 only. Depth>0: P2 joint slingshot+outgoing orbit solve.\n";

    std::vector<ExpansionStats> expansions;
    expansions.push_back(simulate_bfs_expansion(
        planet_params::PlanetId::Earth,
        kLaunchTimeSeconds,
        0.0,
        0.0,
        false,
        kFixedTransferTheta,
        defaults,
        rev));

    if (earth_mars != nullptr) {
        expansions.push_back(simulate_bfs_expansion(
            planet_params::PlanetId::Mars,
            earth_mars->arrival_time_seconds_since_j2000,
            earth_mars->residual_result.transfer_e,
            earth_mars->residual_result.transfer_perihelion_angle_used,
            true,
            kFixedTransferTheta,
            defaults,
            rev));
    }

    const auto* earth_venus = first_candidate(
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Venus,
        kLaunchTimeSeconds,
        kFixedTransferTheta,
        defaults,
        rev);
    if (earth_venus != nullptr) {
        expansions.push_back(simulate_bfs_expansion(
            planet_params::PlanetId::Venus,
            earth_venus->arrival_time_seconds_since_j2000,
            earth_venus->residual_result.transfer_e,
            earth_venus->residual_result.transfer_perihelion_angle_used,
            true,
            kFixedTransferTheta,
            defaults,
            rev));
    }

    double b_depth0 = 0.0;
    double b_depth_gt0 = 0.0;
    int depth0_n = 0;
    int depth_gt0_n = 0;
    double avg_expand_ms = 0.0;
    double avg_p1_expand_ms = 0.0;
    double avg_p2_expand_ms = 0.0;

    for (const ExpansionStats& e : expansions) {
        print_expansion(e);
        avg_expand_ms += e.expand_wall_ms;
        if (!e.uses_problem2_joint_solve) {
            b_depth0 += static_cast<double>(e.valid_edge_count);
            avg_p1_expand_ms += e.expand_wall_ms;
            ++depth0_n;
        } else {
            b_depth_gt0 += static_cast<double>(e.valid_edge_count);
            avg_p2_expand_ms += e.expand_wall_ms;
            ++depth_gt0_n;
        }
    }

    if (depth0_n > 0) {
        b_depth0 /= static_cast<double>(depth0_n);
        avg_p1_expand_ms /= static_cast<double>(depth0_n);
    }
    if (depth_gt0_n > 0) {
        b_depth_gt0 /= static_cast<double>(depth_gt0_n);
        avg_p2_expand_ms /= static_cast<double>(depth_gt0_n);
    }
    avg_expand_ms /= static_cast<double>(expansions.size());

    std::cout << "\nBranching summary:\n";
    std::cout << "  avg valid edges (depth=0, P1): " << b_depth0 << '\n';
    std::cout << "  avg valid edges (depth>0, P2 joint): " << b_depth_gt0 << '\n';
    std::cout << "  avg expand wall (depth=0): " << avg_p1_expand_ms << " ms\n";
    std::cout << "  avg expand wall (depth>0): " << avg_p2_expand_ms << " ms\n";

    const double t0 = avg_p1_expand_ms;
    const double t1 = avg_p2_expand_ms > 0.0 ? avg_p2_expand_ms : t0;
    const double b0 = b_depth0;
    const double b1 = b_depth_gt0 > 0.0 ? b_depth_gt0 : b0;

    std::cout << "\nBFS time estimate (Earth->Mercury, fixed theta):\n";
    for (int max_depth : {2, 3, 4}) {
        double nodes = 0.0;
        double time_ms = 0.0;
        double frontier = 1.0;
        for (int depth = 0; depth < max_depth; ++depth) {
            nodes += frontier;
            time_ms += frontier * ((depth == 0) ? t0 : t1);
            frontier *= (depth == 0) ? b0 : b1;
        }
        std::cout << "  max_depth=" << max_depth
                  << ": nodes~=" << static_cast<int>(nodes)
                  << ", time~=" << time_ms / 1000.0 << " s ("
                  << time_ms << " ms)\n";
    }
}

}  // namespace

int main() {
    const config::GlobalConfig& cfg = config::global_config();
    const auto& defaults = cfg.problem1_solve;

    std::cout << "=== Performance benchmark (corrected P1/P2 roles) ===\n";
    std::cout << "fixed launch_time=" << kLaunchTimeSeconds
              << " s, fixed transfer_theta=" << kFixedTransferTheta << " rad\n";
    std::cout << "problem1 phi_scan_count=" << defaults.phi_scan_count << '\n';
    std::cout << "inner planets: Mercury, Venus, Earth, Mars\n";
    std::cout << "Problem 2: joint slingshot solve via (theta', alpha) scan, "
                 "NOT post-filter on P1.\n";

    run_benchmark_for_rev_limits(defaults, RevLimits{1, 1});
    run_benchmark_for_rev_limits(defaults, RevLimits{2, 2});

    std::cout << "\nNote: full P2 root finder is not implemented; "
                 "benchmark uses grid scan as solve-cost proxy.\n";
    return 0;
}
