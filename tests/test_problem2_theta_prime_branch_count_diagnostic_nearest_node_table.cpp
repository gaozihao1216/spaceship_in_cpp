#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr int kThetaSampleCount = 32;

struct ResidualBranch {
    bool problem1_valid = false;
    bool slingshot_valid = false;
    std::string problem1_invalid_reason;
    std::string slingshot_invalid_reason;
    double theta_prime = 0.0;
    double alpha = 0.0;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double problem1_residual_seconds = 0.0;
    double slingshot_residual = 0.0;
};

struct SampleSummary {
    double theta_prime = 0.0;
    std::vector<ResidualBranch> branches;
    std::map<int, int> raw_count_by_k;
    std::map<int, int> slingshot_valid_count_by_k;
};

struct TopologyStats {
    int stable_interval_count = 0;
    int topology_change_interval_count = 0;
    int sign_change_pair_count = 0;
    int residual_invalid_pair_count = 0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

std::vector<ResidualBranch> make_residual_branches(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double phi,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem2 = spaceship_cpp::problem2;
    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);
    std::vector<ResidualBranch> out;
    out.reserve(branches.size());
    for (const auto& branch : branches) {
        ResidualBranch rb{};
        rb.theta_prime = theta_prime;
        rb.alpha = branch.target_arrival_true_anomaly;
        rb.transfer_revolution = branch.transfer_revolution;
        rb.target_revolution = branch.target_revolution;
        rb.time_of_flight_seconds = branch.time_of_flight_seconds;
        rb.target_time_seconds = branch.target_time_seconds;
        rb.problem1_residual_seconds = branch.residual_seconds;
        if (!branch.valid) {
            rb.problem1_invalid_reason =
                branch.invalid_reason.empty() ? "problem1_branch_invalid" : branch.invalid_reason;
            out.push_back(rb);
            continue;
        }
        rb.problem1_valid = true;
        const auto residual = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            phi,
            branch.target_arrival_true_anomaly,
            incoming_e,
            incoming_theta,
            theta_prime);
        if (!residual.valid) {
            rb.slingshot_invalid_reason =
                residual.invalid_reason.empty() ? "slingshot_invalid" : residual.invalid_reason;
            out.push_back(rb);
            continue;
        }
        rb.slingshot_valid = true;
        rb.slingshot_residual = residual.slingshot_residual;
        out.push_back(rb);
    }
    std::sort(out.begin(), out.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.time_of_flight_seconds) <
               std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.time_of_flight_seconds);
    });
    return out;
}

SampleSummary build_direct_summary(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    auto branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet, phi, beta, theta_A, 1, 1);
    SampleSummary summary{};
    summary.theta_prime = theta_prime;
    summary.branches = make_residual_branches(
        departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime, branches);
    for (const auto& branch : summary.branches) {
        if (branch.problem1_valid) {
            summary.raw_count_by_k[branch.transfer_revolution] += 1;
        }
        if (branch.slingshot_valid) {
            summary.slingshot_valid_count_by_k[branch.transfer_revolution] += 1;
        }
    }
    return summary;
}

SampleSummary build_table_summary(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    bool* fallback_used
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    (void)beta;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto query = problem1::query_problem1_from_2deg_nearest_node(
        loader, departure_planet, target_planet, phi, beta, theta_A, 1, 1, options);
    *fallback_used = query.used_direct_solve_fallback;
    SampleSummary summary{};
    summary.theta_prime = theta_prime;
    summary.branches = make_residual_branches(
        departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime, query.branches);
    for (const auto& branch : summary.branches) {
        if (branch.problem1_valid) {
            summary.raw_count_by_k[branch.transfer_revolution] += 1;
        }
        if (branch.slingshot_valid) {
            summary.slingshot_valid_count_by_k[branch.transfer_revolution] += 1;
        }
    }
    return summary;
}

TopologyStats compute_topology_stats(const std::vector<SampleSummary>& samples) {
    TopologyStats stats{};
    for (std::size_t i = 0; i + 1 < samples.size(); ++i) {
        const auto& left = samples[i];
        const auto& right = samples[i + 1];
        if (left.raw_count_by_k == right.raw_count_by_k) {
            stats.stable_interval_count += 1;
        } else {
            stats.topology_change_interval_count += 1;
        }
        std::map<int, std::vector<ResidualBranch>> left_by_k;
        std::map<int, std::vector<ResidualBranch>> right_by_k;
        for (const auto& branch : left.branches) {
            if (branch.problem1_valid) {
                left_by_k[branch.transfer_revolution].push_back(branch);
            }
        }
        for (const auto& branch : right.branches) {
            if (branch.problem1_valid) {
                right_by_k[branch.transfer_revolution].push_back(branch);
            }
        }
        for (auto& [k, branches] : left_by_k) {
            const auto other_it = right_by_k.find(k);
            if (other_it == right_by_k.end()) {
                continue;
            }
            auto& other = other_it->second;
            std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
            });
            std::sort(other.begin(), other.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
            });
            const int pair_count = std::min<int>(branches.size(), other.size());
            for (int rank = 0; rank < pair_count; ++rank) {
                const auto& a = branches[static_cast<std::size_t>(rank)];
                const auto& b = other[static_cast<std::size_t>(rank)];
                if (!a.slingshot_valid || !b.slingshot_valid) {
                    stats.residual_invalid_pair_count += 1;
                    continue;
                }
                if ((a.slingshot_residual < 0.0 && b.slingshot_residual > 0.0) ||
                    (a.slingshot_residual > 0.0 && b.slingshot_residual < 0.0)) {
                    stats.sign_change_pair_count += 1;
                }
            }
        }
    }
    return stats;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_branch_count_nearest_node_table_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    problem1::Problem1NearestNodeQueryOptions options{};
    options.residual_tolerance_seconds = 1e-2;
    options.max_newton_iterations = 80;
    options.fallback_direct_solve = true;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double incoming_e = 0.3;
    const double incoming_theta = 0.4;

    std::vector<SampleSummary> direct_samples;
    std::vector<SampleSummary> table_samples;
    int table_fallback_count = 0;
    for (int i = 0; i < kThetaSampleCount; ++i) {
        const double theta_prime = kTwoPi * static_cast<double>(i) / static_cast<double>(kThetaSampleCount);
        direct_samples.push_back(build_direct_summary(
            departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta, theta_prime));
        bool fallback_used = false;
        table_samples.push_back(build_table_summary(
            loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
            theta_prime, &fallback_used));
        if (fallback_used) {
            table_fallback_count += 1;
        }
    }

    const TopologyStats direct_stats = compute_topology_stats(direct_samples);
    const TopologyStats table_stats = compute_topology_stats(table_samples);
    const bool ok = table_fallback_count < kThetaSampleCount &&
                    table_stats.stable_interval_count >= 0 &&
                    table_stats.topology_change_interval_count >= 0;

    std::cout << "Problem2ThetaPrimeBranchCountNearestNodeTableSummary\n";
    std::cout << "theta_sample_count=" << kThetaSampleCount << '\n';
    std::cout << "direct_stable_interval_count=" << direct_stats.stable_interval_count << '\n';
    std::cout << "table_stable_interval_count=" << table_stats.stable_interval_count << '\n';
    std::cout << "direct_topology_change_interval_count=" << direct_stats.topology_change_interval_count << '\n';
    std::cout << "table_topology_change_interval_count=" << table_stats.topology_change_interval_count << '\n';
    std::cout << "direct_sign_change_pair_count=" << direct_stats.sign_change_pair_count << '\n';
    std::cout << "table_sign_change_pair_count=" << table_stats.sign_change_pair_count << '\n';
    std::cout << "direct_residual_invalid_pair_count=" << direct_stats.residual_invalid_pair_count << '\n';
    std::cout << "table_residual_invalid_pair_count=" << table_stats.residual_invalid_pair_count << '\n';
    std::cout << "table_fallback_count=" << table_fallback_count << '\n';
    std::cout << "diagnostic_ok=" << (ok ? 1 : 0) << '\n';
    return ok ? 0 : 1;
}
