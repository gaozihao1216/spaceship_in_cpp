/*
 * 文件作用：用随机 (e, θ) 与内行星转移序列测量 P2 solve 耗时，并与“和善”单点 bench 对比。
 * 主要工作：在内四星 {Earth, Mercury, Venus, Mars} 上随机采样行星对与入射轨道，统计 scan/solve 分布。
 */
#include "spaceship_cpp/bfs/trajectory_search_config.hpp"
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_solve.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kDayInSeconds = 86400.0;
constexpr double kTenYearsSeconds = 10.0 * 365.25 * kDayInSeconds;
constexpr int kRandomP2Samples = 40;
constexpr int kRandomSequenceSamples = 24;
constexpr std::uint64_t kSeed = 20260619ULL;

using Clock = std::chrono::steady_clock;
namespace bfs = spaceship_cpp::bfs;
namespace config = spaceship_cpp::config;
namespace planet_params = spaceship_cpp::planet_params;
namespace problem1 = spaceship_cpp::problem1;
namespace problem2 = spaceship_cpp::problem2;
namespace common = spaceship_cpp::common;

using planet_params::PlanetId;

const std::array<PlanetId, 4> kInnerPlanets{
    PlanetId::Mercury,
    PlanetId::Venus,
    PlanetId::Earth,
    PlanetId::Mars,
};

const std::array<PlanetId, 3> kFlybyTargets{
    PlanetId::Mercury,
    PlanetId::Venus,
    PlanetId::Mars,
};

struct TimingSample {
    std::string category;
    std::string description;
    double scan_ms = 0.0;
    double solve_ms = 0.0;
    bool scan_ok = false;
    bool solve_ok = false;
    std::size_t solution_count = 0;
    bool has_profile = false;
    problem2::Problem2FlybyGSearchProfile profile{};
};

struct ProfileAggregate {
    problem2::Problem2FlybyGSearchProfile sum{};
    std::size_t count = 0;

    void add(const problem2::Problem2FlybyGSearchProfile& profile) {
        problem2::merge_problem2_flyby_G_search_profile(sum, profile);
        ++count;
    }

    double mean_ms(double (problem2::Problem2FlybyGSearchProfile::*field)) const {
        if (count == 0U) {
            return 0.0;
        }
        return sum.*field / static_cast<double>(count);
    }

    double mean_count(std::size_t problem2::Problem2FlybyGSearchProfile::*field) const {
        if (count == 0U) {
            return 0.0;
        }
        return static_cast<double>(sum.*field) / static_cast<double>(count);
    }
};

struct TimingStats {
    double min_ms = 0.0;
    double max_ms = 0.0;
    double mean_ms = 0.0;
    double median_ms = 0.0;
    std::size_t count = 0;
    std::size_t scan_ok_count = 0;
    std::size_t solve_ok_count = 0;
};

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

const char* planet_label(PlanetId id) {
    return planet_params::planet_name(id);
}

TimingStats compute_solve_stats(const std::vector<TimingSample>& samples) {
    TimingStats stats{};
    std::vector<double> solve_ms;
    solve_ms.reserve(samples.size());
    for (const TimingSample& sample : samples) {
        if (sample.scan_ok) {
            ++stats.scan_ok_count;
            solve_ms.push_back(sample.solve_ms);
        }
        if (sample.solve_ok) {
            ++stats.solve_ok_count;
        }
    }
    if (solve_ms.empty()) {
        return stats;
    }
    stats.count = solve_ms.size();
    std::sort(solve_ms.begin(), solve_ms.end());
    stats.min_ms = solve_ms.front();
    stats.max_ms = solve_ms.back();
    stats.median_ms = solve_ms[solve_ms.size() / 2U];
    stats.mean_ms =
        std::accumulate(solve_ms.begin(), solve_ms.end(), 0.0) /
        static_cast<double>(solve_ms.size());
    return stats;
}

void print_stats(const char* label, const TimingStats& stats) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << label << ": n=" << stats.count
              << " scan_ok=" << stats.scan_ok_count
              << " solve_ok=" << stats.solve_ok_count
              << " min=" << stats.min_ms
              << " median=" << stats.median_ms
              << " mean=" << stats.mean_ms
              << " max=" << stats.max_ms << " ms\n";
}

problem2::Problem2FlybySolveConfig make_p2_config(
    PlanetId flyby,
    PlanetId target,
    double flyby_time,
    const bfs::TrajectorySearchGlobalConfig& global,
    bool collect_profile
) {
    const config::GlobalConfig& defaults = config::global_config();
    problem2::Problem2FlybySolveConfig solve_config = config::make_problem2_flyby_solve_config(
        flyby,
        target,
        flyby_time,
        defaults.problem2_theta_prime_scan,
        defaults.problem2_route_a_newton,
        bfs::make_problem1_solve_defaults(global.problem1));
    solve_config.collect_g_search_profile = collect_profile;
    return solve_config;
}

void print_profile_breakdown(
    const char* label,
    const ProfileAggregate& aggregate,
    double mean_solve_ms
) {
    if (aggregate.count == 0U) {
        return;
    }
    const auto& sum = aggregate.sum;
    const double n = static_cast<double>(aggregate.count);
    const double mean_route_a = sum.route_a_ms / n;
    const double mean_case_c = sum.case_c_middle_ms / n;
    const double mean_incoming = sum.incoming_cache_ms / n;
    const double mean_enrich = sum.enrich_ms / n;
    const double mean_framework = mean_solve_ms - mean_route_a - mean_case_c - mean_incoming - mean_enrich;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << label << " (n=" << aggregate.count << ", mean solve_ms=" << mean_solve_ms << ")\n";
    std::cout << "  time_ms: route_a=" << mean_route_a
              << " (" << (100.0 * mean_route_a / std::max(mean_solve_ms, 1e-9)) << "%)"
              << ", case_c=" << mean_case_c
              << " (" << (100.0 * mean_case_c / std::max(mean_solve_ms, 1e-9)) << "%)"
              << ", framework=" << mean_framework
              << " (" << (100.0 * mean_framework / std::max(mean_solve_ms, 1e-9)) << "%)"
              << ", incoming_F=" << mean_incoming
              << ", enrich=" << mean_enrich << '\n';
    std::cout << std::setprecision(1);
    std::cout << "  counts/sample: interval_visits=" << (sum.interval_visits / n)
              << " equal=" << (sum.interval_equal / n)
              << " case_b=" << (sum.interval_case_b / n)
              << " case_c=" << (sum.interval_case_c / n)
              << " discarded=" << (sum.interval_discarded / n) << '\n';
    std::cout << "  calls/sample: route_a=" << (sum.route_a_calls / n)
              << " (iters=" << (sum.route_a_iterations / n) << ")"
              << ", case_c_middle=" << (sum.case_c_middle_calls / n)
              << ", case_b_probe=" << (sum.case_b_probe_calls / n)
              << ", g_newton=" << (sum.g_newton_calls / n)
              << " (iters=" << (sum.g_newton_iterations / n) << ")"
              << ", equal_branch_proc=" << (sum.branch_interval_process_calls / n) << '\n';
}

double percent_of(std::size_t part, std::size_t whole) {
    if (whole == 0U) {
        return 0.0;
    }
    return 100.0 * static_cast<double>(part) / static_cast<double>(whole);
}

void print_case_c_branch_analysis(const char* label, const ProfileAggregate& aggregate) {
    if (aggregate.count == 0U) {
        return;
    }
    const auto& p = aggregate.sum;
    const double n = static_cast<double>(aggregate.count);
    const double top_total =
        static_cast<double>(p.top_interval_equal + p.top_interval_case_b + p.top_interval_case_c) /
        n;

    std::cout << label << " Case C/B branch topology (aggregated over " << aggregate.count
              << " samples)\n";
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  top-level intervals/sample (depth=0): total=" << top_total
              << " equal=" << (p.top_interval_equal / n)
              << " case_b=" << (p.top_interval_case_b / n)
              << " case_c=" << (p.top_interval_case_c / n)
              << "  [expect ~63 per k-layer]\n";

    if (p.case_c_split_samples == 0U) {
        std::cout << "  Case C splits: none\n";
    } else {
        std::cout << std::setprecision(1);
        std::cout << "  Case C splits/sample: " << (static_cast<double>(p.case_c_split_samples) / n)
                  << ", mean endpoint gap="
                  << (static_cast<double>(p.case_c_endpoint_gap_sum) /
                      static_cast<double>(p.case_c_split_samples))
                  << '\n';
        std::cout << std::setprecision(1);
        std::cout << "  n_middle vs [n_min,n_max]: gt_max="
                  << percent_of(p.case_c_middle_gt_max_endpoints, p.case_c_split_samples) << "%"
                  << " in_range="
                  << percent_of(p.case_c_middle_in_endpoint_range, p.case_c_split_samples) << "%"
                  << " lt_min="
                  << percent_of(p.case_c_middle_lt_min_endpoints, p.case_c_split_samples) << "%\n";
        std::cout << "  after Case C split, left child: equal="
                  << percent_of(p.case_c_child_left_equal, p.case_c_split_samples) << "%"
                  << " case_b="
                  << percent_of(p.case_c_child_left_case_b, p.case_c_split_samples) << "%"
                  << " case_c="
                  << percent_of(p.case_c_child_left_case_c, p.case_c_split_samples) << "%\n";
        std::cout << "  after Case C split, right child: equal="
                  << percent_of(p.case_c_child_right_equal, p.case_c_split_samples) << "%"
                  << " case_b="
                  << percent_of(p.case_c_child_right_case_b, p.case_c_split_samples) << "%"
                  << " case_c="
                  << percent_of(p.case_c_child_right_case_c, p.case_c_split_samples) << "%\n";
    }

    if (p.case_b_split_samples == 0U) {
        std::cout << "  Case B splits: none\n";
    } else {
        std::cout << std::setprecision(1);
        std::cout << "  Case B splits/sample: " << (static_cast<double>(p.case_b_split_samples) / n)
                  << ", n_middle gt_max="
                  << percent_of(p.case_b_middle_gt_max_endpoints, p.case_b_split_samples) << "%"
                  << " in_range="
                  << percent_of(p.case_b_middle_in_endpoint_range, p.case_b_split_samples) << "%"
                  << " lt_min="
                  << percent_of(p.case_b_middle_lt_min_endpoints, p.case_b_split_samples) << "%\n";
        std::cout << "  after Case B split, child still Case C: left="
                  << percent_of(p.case_b_child_left_case_c, p.case_b_split_samples) << "%"
                  << " right="
                  << percent_of(p.case_b_child_right_case_c, p.case_b_split_samples) << "%\n";
    }
}

void compare_fast_vs_slow_profiles(const std::vector<TimingSample>& samples) {
    std::vector<const TimingSample*> with_profile;
    with_profile.reserve(samples.size());
    for (const TimingSample& sample : samples) {
        if (sample.scan_ok && sample.has_profile) {
            with_profile.push_back(&sample);
        }
    }
    if (with_profile.size() < 4U) {
        std::cout << "  (not enough profile samples for fast/slow split)\n";
        return;
    }

    std::sort(
        with_profile.begin(),
        with_profile.end(),
        [](const TimingSample* lhs, const TimingSample* rhs) {
            return lhs->solve_ms < rhs->solve_ms;
        });

    const std::size_t quartile = std::max<std::size_t>(1U, with_profile.size() / 4U);
    ProfileAggregate fast{};
    ProfileAggregate slow{};
    double fast_solve_sum = 0.0;
    double slow_solve_sum = 0.0;
    for (std::size_t i = 0; i < quartile; ++i) {
        fast.add(with_profile[i]->profile);
        fast_solve_sum += with_profile[i]->solve_ms;
    }
    for (std::size_t i = with_profile.size() - quartile; i < with_profile.size(); ++i) {
        slow.add(with_profile[i]->profile);
        slow_solve_sum += with_profile[i]->solve_ms;
    }

    const double fast_mean_solve = fast_solve_sum / static_cast<double>(quartile);
    const double slow_mean_solve = slow_solve_sum / static_cast<double>(quartile);

    std::cout << "\n--- G-search profile: fast vs slow quartile (random_p2) ---\n";
    print_profile_breakdown("FAST quartile", fast, fast_mean_solve);
    print_profile_breakdown("SLOW quartile", slow, slow_mean_solve);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  slow/fast ratios: solve_ms="
              << (slow_mean_solve / std::max(fast_mean_solve, 1e-9)) << "x"
              << ", route_a_ms="
              << (slow.mean_ms(&problem2::Problem2FlybyGSearchProfile::route_a_ms) /
                  std::max(fast.mean_ms(&problem2::Problem2FlybyGSearchProfile::route_a_ms), 1e-9))
              << "x"
              << ", case_c_ms="
              << (slow.mean_ms(&problem2::Problem2FlybyGSearchProfile::case_c_middle_ms) /
                  std::max(fast.mean_ms(&problem2::Problem2FlybyGSearchProfile::case_c_middle_ms), 1e-9))
              << "x"
              << ", route_a_calls="
              << (slow.mean_count(&problem2::Problem2FlybyGSearchProfile::route_a_calls) /
                  std::max(fast.mean_count(&problem2::Problem2FlybyGSearchProfile::route_a_calls), 1e-9))
              << "x"
              << ", case_c_calls="
              << (slow.mean_count(&problem2::Problem2FlybyGSearchProfile::case_c_middle_calls) /
                  std::max(fast.mean_count(&problem2::Problem2FlybyGSearchProfile::case_c_middle_calls), 1e-9))
              << "x"
              << ", interval_case_c="
              << (slow.mean_count(&problem2::Problem2FlybyGSearchProfile::interval_case_c) /
                  std::max(fast.mean_count(&problem2::Problem2FlybyGSearchProfile::interval_case_c), 1e-9))
              << "x\n";

    const TimingSample* slowest = with_profile.back();
    std::cout << "\n  slowest sample detail: " << slowest->description << '\n';
    print_profile_breakdown("  ", ProfileAggregate{.sum = slowest->profile, .count = 1U}, slowest->solve_ms);

    std::cout << "\n--- Case C middle branch count analysis ---\n";
    print_case_c_branch_analysis("FAST quartile", fast);
    print_case_c_branch_analysis("SLOW quartile", slow);
    ProfileAggregate slowest_agg{};
    slowest_agg.add(slowest->profile);
    print_case_c_branch_analysis("SLOWEST sample", slowest_agg);
}

TimingSample run_p2_cold(
    PlanetId flyby,
    PlanetId target,
    double flyby_time,
    double incoming_e,
    double incoming_theta,
    const bfs::TrajectorySearchGlobalConfig& global,
    std::string category,
    std::string description
) {
    TimingSample sample{};
    sample.category = std::move(category);
    sample.description = std::move(description);

    const problem2::Problem2FlybySolveConfig solve_config =
        make_p2_config(flyby, target, flyby_time, global, true);

    problem2::Problem2FlybySolveInput scan_input{};
    scan_input.flyby_planet = flyby;
    scan_input.target_planet = target;
    scan_input.flyby_time_seconds_since_j2000 = flyby_time;

    const Clock::time_point scan_start = Clock::now();
    const problem2::Problem2ThetaPrimeInitialScanResult scan =
        problem2::run_problem2_flyby_theta_prime_initial_scan(scan_input, solve_config);
    sample.scan_ms = elapsed_ms(scan_start, Clock::now());
    sample.scan_ok = scan.ok;
    if (!scan.ok) {
        return sample;
    }

    problem2::Problem2FlybySolveInput solve_input = scan_input;
    solve_input.incoming_eccentricity = incoming_e;
    solve_input.incoming_theta = incoming_theta;
    solve_input.incoming_theta_is_local = false;

    const Clock::time_point solve_start = Clock::now();
    const problem2::Problem2FlybySolveResult solved =
        problem2::solve_problem2_flyby_with_scan(solve_input, solve_config, scan);
    sample.solve_ms = elapsed_ms(solve_start, Clock::now());
    sample.solve_ok = solved.ok && !solved.solutions.empty();
    sample.solution_count = solved.solutions.size();
    if (solved.has_g_search_profile) {
        sample.has_profile = true;
        sample.profile = solved.g_search_profile;
    }
    return sample;
}

PlanetId pick_random_flyby_target(std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> dist(0U, kFlybyTargets.size() - 1U);
    return kFlybyTargets[dist(rng)];
}

PlanetId pick_random_inner(std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> dist(0U, kInnerPlanets.size() - 1U);
    return kInnerPlanets[dist(rng)];
}

double random_theta(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(0.0, common::kTwoPi);
    return dist(rng);
}

double random_eccentricity(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(0.05, 8.0);
    return dist(rng);
}

double random_flyby_time(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(0.0, kTenYearsSeconds);
    return dist(rng);
}

std::vector<TimingSample> run_easy_benchmark(
    const bfs::TrajectorySearchGlobalConfig& global
) {
    std::vector<TimingSample> samples;
    const auto p1_defaults = bfs::make_problem1_solve_defaults(global.problem1);
    const auto venus_leg0 = config::make_problem1_solve_input(
        PlanetId::Earth, PlanetId::Venus, 0.0, 1.0, p1_defaults);
    const auto venus_candidates = problem1::solve_problem1(venus_leg0);
    if (venus_candidates.empty()) {
        return samples;
    }
    const problem1::Problem1Candidate& candidate = venus_candidates.front();
    samples.push_back(run_p2_cold(
        PlanetId::Venus,
        PlanetId::Mercury,
        candidate.arrival_time_seconds_since_j2000,
        candidate.residual_result.transfer_e,
        candidate.residual_result.transfer_perihelion_angle_used,
        global,
        "easy_bench",
        "Venus->Mercury, leg0 first candidate (e,theta)"));
    return samples;
}

std::vector<TimingSample> run_random_p2_samples(
    const bfs::TrajectorySearchGlobalConfig& global,
    std::mt19937_64& rng
) {
    std::vector<TimingSample> samples;
    samples.reserve(static_cast<std::size_t>(kRandomP2Samples));
    for (int i = 0; i < kRandomP2Samples; ++i) {
        PlanetId flyby = pick_random_flyby_target(rng);
        PlanetId target = pick_random_flyby_target(rng);
        while (target == flyby) {
            target = pick_random_flyby_target(rng);
        }
        const double flyby_time = random_flyby_time(rng);
        const double incoming_e = random_eccentricity(rng);
        const double incoming_theta = random_theta(rng);
        std::ostringstream desc;
        desc << planet_label(flyby) << "->" << planet_label(target)
             << " t=" << (flyby_time / kDayInSeconds) << "d"
             << " e=" << incoming_e << " th=" << incoming_theta;
        samples.push_back(run_p2_cold(
            flyby,
            target,
            flyby_time,
            incoming_e,
            incoming_theta,
            global,
            "random_p2",
            desc.str()));
    }
    return samples;
}

std::vector<TimingSample> run_random_sequence_samples(
    const bfs::TrajectorySearchGlobalConfig& global,
    std::mt19937_64& rng
) {
    std::vector<TimingSample> samples;
    const auto p1_defaults = bfs::make_problem1_solve_defaults(global.problem1);
    std::uniform_int_distribution<int> leg_dist(1, 3);

    for (int seq_idx = 0; seq_idx < kRandomSequenceSamples; ++seq_idx) {
        const int extra_legs = leg_dist(rng);
        std::vector<PlanetId> sequence;
        sequence.push_back(PlanetId::Earth);
        for (int leg = 0; leg < extra_legs; ++leg) {
            sequence.push_back(pick_random_flyby_target(rng));
        }

        std::ostringstream seq_label;
        for (std::size_t i = 0; i < sequence.size(); ++i) {
            if (i > 0U) {
                seq_label << "->";
            }
            seq_label << planet_label(sequence[i]);
        }

        const double launch_time = 0.0;
        const double leg0_theta = random_theta(rng);
        const auto leg0_input = config::make_problem1_solve_input(
            PlanetId::Earth,
            sequence[1U],
            launch_time,
            leg0_theta,
            p1_defaults);
        const auto leg0_candidates = problem1::solve_problem1(leg0_input);

        double flyby_time = random_flyby_time(rng);
        if (!leg0_candidates.empty()) {
            flyby_time = leg0_candidates.front().arrival_time_seconds_since_j2000;
        }

        for (std::size_t leg = 1U; leg + 1U < sequence.size(); ++leg) {
            const PlanetId flyby = sequence[leg];
            const PlanetId target = sequence[leg + 1U];
            const double incoming_e = random_eccentricity(rng);
            const double incoming_theta = random_theta(rng);
            if (leg == 1U && !leg0_candidates.empty()) {
                flyby_time = leg0_candidates.front().arrival_time_seconds_since_j2000;
            } else if (leg > 1U) {
                flyby_time = random_flyby_time(rng);
            }
            std::ostringstream desc;
            desc << "seq[" << seq_idx << "] " << seq_label.str()
                 << " leg" << leg << " " << planet_label(flyby) << "->"
                 << planet_label(target)
                 << " rand(e,th)";
            samples.push_back(run_p2_cold(
                flyby,
                target,
                flyby_time,
                incoming_e,
                incoming_theta,
                global,
                "random_sequence",
                desc.str()));
        }
    }
    return samples;
}

void print_slowest(
    const std::vector<TimingSample>& samples,
    const char* category,
    std::size_t top_n
) {
    std::vector<TimingSample> filtered;
    for (const TimingSample& sample : samples) {
        if (sample.category == category && sample.scan_ok) {
            filtered.push_back(sample);
        }
    }
    std::sort(
        filtered.begin(),
        filtered.end(),
        [](const TimingSample& lhs, const TimingSample& rhs) {
            return lhs.solve_ms > rhs.solve_ms;
        });
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  slowest " << category << " (top " << top_n << "):\n";
    for (std::size_t i = 0; i < std::min(top_n, filtered.size()); ++i) {
        const TimingSample& s = filtered[i];
        std::cout << "    [" << i << "] solve_ms=" << s.solve_ms
                  << " scan_ok=" << s.scan_ok
                  << " solve_ok=" << s.solve_ok
                  << " n_sol=" << s.solution_count
                  << " | " << s.description << '\n';
    }
}

}  // namespace

int main() {
    auto global = bfs::default_trajectory_search_global_config();
    global.mission.launch_time_seconds_since_j2000 = 0.0;
    global.constraints.max_total_time_seconds = 40.0 * 365.25 * kDayInSeconds;
    global.constraints.max_launch_v_inf_mps = std::numeric_limits<double>::infinity();

    std::mt19937_64 rng(kSeed);

    std::cout << "=== P2 random (e,theta) + inner-planet sequence timing ===\n";
    std::cout << "seed=" << kSeed
              << " random_p2=" << kRandomP2Samples
              << " random_sequences=" << kRandomSequenceSamples << '\n';

    const std::vector<TimingSample> easy = run_easy_benchmark(global);
    const std::vector<TimingSample> random_p2 = run_random_p2_samples(global, rng);
    const std::vector<TimingSample> random_seq = run_random_sequence_samples(global, rng);

    std::cout << "\n--- solve_ms distribution ---\n";
    print_stats("easy_bench (Venus->Mercury leg0)", compute_solve_stats(easy));
    print_stats("random_p2 (random pair/t/e/theta)", compute_solve_stats(random_p2));
    print_stats("random_sequence (random seq + rand e,th)", compute_solve_stats(random_seq));

    const TimingStats easy_stats = compute_solve_stats(easy);
    const TimingStats random_stats = compute_solve_stats(random_p2);
    const TimingStats seq_stats = compute_solve_stats(random_seq);

    std::cout << "\n--- easy vs random comparison ---\n";
    std::cout << std::fixed << std::setprecision(2);
    if (easy_stats.count > 0U && random_stats.count > 0U) {
        std::cout << "random_p2 mean / easy mean = "
                  << (random_stats.mean_ms / easy_stats.mean_ms) << "x\n";
        std::cout << "random_p2 median / easy median = "
                  << (random_stats.median_ms / easy_stats.median_ms) << "x\n";
    }
    if (easy_stats.count > 0U && seq_stats.count > 0U) {
        std::cout << "random_sequence mean / easy mean = "
                  << (seq_stats.mean_ms / easy_stats.mean_ms) << "x\n";
    }
    if (random_stats.count > 0U) {
        std::cout << "random_p2 solve_ok rate = "
                  << (100.0 * static_cast<double>(random_stats.solve_ok_count) /
                      static_cast<double>(random_stats.scan_ok_count))
                  << "% (timing includes failed solves)\n";
    }

    std::cout << "\n--- slowest cases ---\n";
    print_slowest(random_p2, "random_p2", 5U);
    print_slowest(random_seq, "random_sequence", 5U);

    compare_fast_vs_slow_profiles(random_p2);

    ProfileAggregate all_random_p2{};
    for (const TimingSample& sample : random_p2) {
        if (sample.has_profile) {
            all_random_p2.add(sample.profile);
        }
    }
    std::cout << "\n--- Case C branch analysis (all random_p2) ---\n";
    print_case_c_branch_analysis("ALL", all_random_p2);

    std::cout << "\nConclusion: compare route_a_ms vs case_c_ms ratios in fast/slow quartiles;\n"
              << "the diverging bucket shows which G-search path drives tail latency.\n";

    return 0;
}
