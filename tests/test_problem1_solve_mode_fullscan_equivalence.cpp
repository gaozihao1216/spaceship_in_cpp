#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

struct PairCase {
    spaceship_cpp::planet_params::PlanetId from;
    spaceship_cpp::planet_params::PlanetId to;
};

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    const int value = std::atoi(raw);
    return value > 0 ? value : fallback;
}

unsigned env_uint(const char* name, unsigned fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    const unsigned long value = std::strtoul(raw, nullptr, 10);
    return value > 0 ? static_cast<unsigned>(value) : fallback;
}

double angle_diff(double lhs, double rhs) {
    double diff = std::fmod(lhs - rhs + spaceship_cpp::common::kPi, spaceship_cpp::common::kTwoPi);
    if (diff < 0.0) diff += spaceship_cpp::common::kTwoPi;
    return std::abs(diff - spaceship_cpp::common::kPi);
}

bool branches_equal(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& lhs,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& rhs,
    std::string* reason
) {
    if (lhs.size() != rhs.size()) {
        *reason = "branch_count_mismatch";
        return false;
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const auto& a = lhs[i];
        const auto& b = rhs[i];
        if (a.valid != b.valid ||
            a.transfer_revolution != b.transfer_revolution ||
            a.target_revolution != b.target_revolution) {
            *reason = "branch_identity_mismatch";
            return false;
        }
        if (angle_diff(a.target_arrival_true_anomaly, b.target_arrival_true_anomaly) > 1e-12 ||
            angle_diff(a.encounter_global_angle, b.encounter_global_angle) > 1e-12 ||
            std::abs(a.target_time_seconds - b.target_time_seconds) > 1e-9 ||
            std::abs(a.residual_seconds - b.residual_seconds) > 1e-12) {
            *reason = "branch_value_mismatch";
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const int sample_count = env_int("ROOT_TABLE_5DEG_TEST_SAMPLE_COUNT", 20);
    const unsigned seed = env_uint("ROOT_TABLE_5DEG_TEST_SEED", 12345);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> angle_dist(0.0, common::kTwoPi);

    const std::vector<PairCase> pairs{
        {planet::PlanetId::Earth, planet::PlanetId::Mars},
        {planet::PlanetId::Mars, planet::PlanetId::Earth},
        {planet::PlanetId::Venus, planet::PlanetId::Mercury},
    };

    int query_count = 0;
    int mismatch_count = 0;
    int total_branch_count = 0;

    for (const auto& pair : pairs) {
        int pair_mismatch_count = 0;
        int pair_branch_count = 0;
        for (int i = 0; i < sample_count; ++i) {
            const double nu_A = angle_dist(rng);
            const double nu_B = angle_dist(rng);
            const double theta_A = angle_dist(rng);

            const auto old_branches = problem1::solve_problem1_from_departure_anomalies(
                pair.from, pair.to, nu_A, nu_B, theta_A, 1, 1);

            problem1::Problem1SolveOptions options{};
            options.mode = problem1::Problem1SolveMode::FullScan2880;
            options.max_transfer_revolution = 1;
            options.max_target_revolution = 1;
            const auto mode_result = problem1::solve_problem1_from_departure_anomalies_with_mode(
                pair.from, pair.to, nu_A, nu_B, theta_A, options);

            std::string reason;
            if (!branches_equal(old_branches, mode_result.branches, &reason) ||
                mode_result.mode_profile.mode_requested != problem1::Problem1SolveMode::FullScan2880 ||
                mode_result.mode_profile.mode_used != problem1::Problem1SolveMode::FullScan2880 ||
                mode_result.mode_profile.adaptive_attempted ||
                mode_result.mode_profile.fallback_used) {
                mismatch_count += 1;
                pair_mismatch_count += 1;
                std::cout << "Problem1SolveModeFullScanEquivalenceMismatch\n";
                std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
                std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
                std::cout << "sample_index=" << i << '\n';
                std::cout << "reason=" << reason << '\n';
                std::cout << "nu_A=" << nu_A << '\n';
                std::cout << "nu_B=" << nu_B << '\n';
                std::cout << "theta_A=" << theta_A << '\n';
                std::cout << "old_branch_count=" << old_branches.size() << '\n';
                std::cout << "new_branch_count=" << mode_result.branches.size() << '\n';
            }

            query_count += 1;
            pair_branch_count += static_cast<int>(old_branches.size());
            total_branch_count += static_cast<int>(old_branches.size());
        }

        std::cout << "Problem1SolveModeFullScanEquivalencePairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "query_count=" << sample_count << '\n';
        std::cout << "branch_count=" << pair_branch_count << '\n';
        std::cout << "mismatch_count=" << pair_mismatch_count << '\n';
    }

    std::cout << "Problem1SolveModeFullScanEquivalenceSummary\n";
    std::cout << "query_count=" << query_count << '\n';
    std::cout << "total_branch_count=" << total_branch_count << '\n';
    std::cout << "mismatch_count=" << mismatch_count << '\n';
    std::cout << "equivalence_ok=" << (mismatch_count == 0 ? 1 : 0) << '\n';
    return mismatch_count == 0 ? 0 : 1;
}
