#include "problem1_5deg_test_utils.hpp"

#include <iomanip>
#include <iostream>

namespace {

problem1_5deg_test::QueryPoint perturb_query(
    const problem1_5deg_test::QueryPoint& query,
    int axis,
    double delta
) {
    auto out = query;
    if (axis == 0) {
        out.nu_A = spaceship_cpp::common::normalize_angle_0_2pi(out.nu_A + delta);
    } else if (axis == 1) {
        out.nu_B = spaceship_cpp::common::normalize_angle_0_2pi(out.nu_B + delta);
    } else {
        out.theta_A = spaceship_cpp::common::normalize_angle_0_2pi(out.theta_A + delta);
    }
    return out;
}

bool is_boundary_like(
    const problem1_5deg_test::QueryPoint& query,
    int base_branch_count
) {
    constexpr double degree = spaceship_cpp::common::kPi / 180.0;
    const double perturbations[] = {0.1 * degree, 0.2 * degree, 0.5 * degree, 1.0 * degree, 2.5 * degree};
    for (double h : perturbations) {
        for (int axis = 0; axis < 3; ++axis) {
            for (double sign : {-1.0, 1.0}) {
                const auto shifted = perturb_query(query, axis, sign * h);
                const auto branches = problem1_5deg_test::solve_direct(shifted);
                if (static_cast<int>(branches.size()) != base_branch_count) {
                    return true;
                }
            }
        }
    }
    return false;
}

}  // namespace

int main() {
    namespace util = problem1_5deg_test;
    namespace planet = spaceship_cpp::planet_params;

    std::cout << std::setprecision(12) << std::scientific;
    const int sample_count = util::env_int("ROOT_TABLE_5DEG_TEST_SAMPLE_COUNT", 100);
    const unsigned seed = util::env_uint("ROOT_TABLE_5DEG_TEST_SEED", 12345);
    auto pairs = util::selected_pairs_from_env();
    std::mt19937_64 rng(seed);

    long long boundary_like_case_count = 0;
    long long non_boundary_case_count = 0;
    long long boundary_direct_branch_count = 0;
    long long non_boundary_direct_branch_count = 0;
    long long boundary_missing_count = 0;
    long long non_boundary_missing_count = 0;
    long long boundary_matched_count = 0;
    long long non_boundary_matched_count = 0;

    for (const auto& pair : pairs) {
        long long pair_boundary_cases = 0;
        long long pair_non_boundary_cases = 0;
        long long pair_boundary_missing = 0;
        long long pair_non_boundary_missing = 0;
        long long pair_boundary_direct = 0;
        long long pair_non_boundary_direct = 0;
        const auto queries = util::make_random_queries_for_pair(pair, sample_count, &rng);
        for (const auto& query : queries) {
            auto direct = util::solve_direct(query);
            const bool boundary_like = is_boundary_like(query, static_cast<int>(direct.size()));
            const auto simulated = util::simulate_5deg_nearest_refine(query);
            const auto match = util::compare_branches(direct, simulated.refined_branches);
            if (boundary_like) {
                boundary_like_case_count += 1;
                pair_boundary_cases += 1;
                boundary_direct_branch_count += match.direct_branch_count;
                boundary_missing_count += match.missing_count;
                boundary_matched_count += match.matched_branch_count;
                pair_boundary_direct += match.direct_branch_count;
                pair_boundary_missing += match.missing_count;
            } else {
                non_boundary_case_count += 1;
                pair_non_boundary_cases += 1;
                non_boundary_direct_branch_count += match.direct_branch_count;
                non_boundary_missing_count += match.missing_count;
                non_boundary_matched_count += match.matched_branch_count;
                pair_non_boundary_direct += match.direct_branch_count;
                pair_non_boundary_missing += match.missing_count;
            }
        }

        std::cout << "Problem1FiveDegBranchBoundaryStressPairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "boundary_like_case_count=" << pair_boundary_cases << '\n';
        std::cout << "non_boundary_case_count=" << pair_non_boundary_cases << '\n';
        std::cout << "boundary_case_missing_rate_5deg=" << static_cast<double>(pair_boundary_missing) / static_cast<double>(std::max<long long>(1, pair_boundary_direct)) << '\n';
        std::cout << "non_boundary_case_missing_rate_5deg=" << static_cast<double>(pair_non_boundary_missing) / static_cast<double>(std::max<long long>(1, pair_non_boundary_direct)) << '\n';
    }

    std::cout << "Problem1FiveDegBranchBoundaryStressGlobalSummary\n";
    std::cout << "boundary_like_case_count=" << boundary_like_case_count << '\n';
    std::cout << "non_boundary_case_count=" << non_boundary_case_count << '\n';
    std::cout << "boundary_direct_branch_count=" << boundary_direct_branch_count << '\n';
    std::cout << "non_boundary_direct_branch_count=" << non_boundary_direct_branch_count << '\n';
    std::cout << "boundary_matched_count=" << boundary_matched_count << '\n';
    std::cout << "non_boundary_matched_count=" << non_boundary_matched_count << '\n';
    std::cout << "boundary_missing_count=" << boundary_missing_count << '\n';
    std::cout << "non_boundary_missing_count=" << non_boundary_missing_count << '\n';
    std::cout << "boundary_case_missing_rate_5deg=" << static_cast<double>(boundary_missing_count) / static_cast<double>(std::max<long long>(1, boundary_direct_branch_count)) << '\n';
    std::cout << "non_boundary_case_missing_rate_5deg=" << static_cast<double>(non_boundary_missing_count) / static_cast<double>(std::max<long long>(1, non_boundary_direct_branch_count)) << '\n';
    std::cout << "problem1_5deg_branch_boundary_stress_ok=1\n";
    return 0;
}
