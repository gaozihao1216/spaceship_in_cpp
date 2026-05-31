#include "problem1_5deg_test_utils.hpp"

#include <iomanip>
#include <iostream>

int main() {
    namespace util = problem1_5deg_test;
    namespace planet = spaceship_cpp::planet_params;

    std::cout << std::setprecision(12) << std::scientific;
    const int sample_count = util::env_int("ROOT_TABLE_5DEG_TEST_SAMPLE_COUNT", 100);
    const unsigned seed = util::env_uint("ROOT_TABLE_5DEG_TEST_SEED", 12345);
    auto pairs = util::selected_pairs_from_env();
    std::mt19937_64 rng(seed);

    long long total_query_count = 0;
    long long total_direct_branch_count = 0;
    long long total_grid_branch_count = 0;
    long long total_matched_branch_count = 0;
    long long total_missing_count = 0;
    long long total_new_count = 0;
    long long branch_count_change_count = 0;
    long long branch_order_mismatch_count = 0;
    double theta_error_sum = 0.0;
    double alpha_error_sum = 0.0;
    double arrival_error_sum = 0.0;
    double theta_error_max = 0.0;
    double alpha_error_max = 0.0;
    double arrival_error_max = 0.0;

    for (const auto& pair : pairs) {
        long long pair_queries = 0;
        long long pair_branch_count_change = 0;
        long long pair_direct = 0;
        long long pair_grid = 0;
        long long pair_missing = 0;
        const auto queries = util::make_random_queries_for_pair(pair, sample_count, &rng);
        for (const auto& query : queries) {
            auto direct = util::solve_direct(query);
            const auto grid_query = util::nearest_5deg_query(query);
            auto grid = util::solve_direct(grid_query);
            const auto match = util::compare_branches(direct, grid);

            total_query_count += 1;
            pair_queries += 1;
            total_direct_branch_count += match.direct_branch_count;
            total_grid_branch_count += match.candidate_branch_count;
            total_matched_branch_count += match.matched_branch_count;
            total_missing_count += match.missing_count;
            total_new_count += match.new_count;
            branch_order_mismatch_count += match.branch_order_mismatch_count;
            theta_error_sum += match.theta_prime_error_sum;
            alpha_error_sum += match.alpha_error_sum;
            arrival_error_sum += match.arrival_time_error_sum;
            theta_error_max = std::max(theta_error_max, match.theta_prime_error_max);
            alpha_error_max = std::max(alpha_error_max, match.alpha_error_max);
            arrival_error_max = std::max(arrival_error_max, match.arrival_time_error_max);

            pair_direct += match.direct_branch_count;
            pair_grid += match.candidate_branch_count;
            pair_missing += match.missing_count;
            if (match.direct_branch_count != match.candidate_branch_count) {
                branch_count_change_count += 1;
                pair_branch_count_change += 1;
            }
        }

        std::cout << "Problem1FiveDegNearestSensitivityPairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "query_count=" << pair_queries << '\n';
        std::cout << "direct_branch_count=" << pair_direct << '\n';
        std::cout << "grid_branch_count=" << pair_grid << '\n';
        std::cout << "branch_count_change_count=" << pair_branch_count_change << '\n';
        std::cout << "missing_branch_count=" << pair_missing << '\n';
        std::cout << "branch_count_change_rate=" << static_cast<double>(pair_branch_count_change) / static_cast<double>(std::max<long long>(1, pair_queries)) << '\n';
        std::cout << "missing_rate=" << static_cast<double>(pair_missing) / static_cast<double>(std::max<long long>(1, pair_direct)) << '\n';
    }

    const double matched_den = static_cast<double>(std::max<long long>(1, total_matched_branch_count));
    std::cout << "Problem1FiveDegNearestSensitivityGlobalSummary\n";
    std::cout << "total_query_count=" << total_query_count << '\n';
    std::cout << "total_direct_branch_count=" << total_direct_branch_count << '\n';
    std::cout << "total_grid_branch_count=" << total_grid_branch_count << '\n';
    std::cout << "matched_branch_count=" << total_matched_branch_count << '\n';
    std::cout << "missing_count=" << total_missing_count << '\n';
    std::cout << "new_count=" << total_new_count << '\n';
    std::cout << "branch_count_change_count=" << branch_count_change_count << '\n';
    std::cout << "branch_count_change_rate=" << static_cast<double>(branch_count_change_count) / static_cast<double>(std::max<long long>(1, total_query_count)) << '\n';
    std::cout << "missing_rate=" << static_cast<double>(total_missing_count) / static_cast<double>(std::max<long long>(1, total_direct_branch_count)) << '\n';
    std::cout << "new_branch_rate=" << static_cast<double>(total_new_count) / static_cast<double>(std::max<long long>(1, total_grid_branch_count)) << '\n';
    std::cout << "branch_order_mismatch_count=" << branch_order_mismatch_count << '\n';
    std::cout << "mean_abs_theta_prime_error=" << theta_error_sum / matched_den << '\n';
    std::cout << "max_abs_theta_prime_error=" << theta_error_max << '\n';
    std::cout << "mean_abs_alpha_error=" << alpha_error_sum / matched_den << '\n';
    std::cout << "max_abs_alpha_error=" << alpha_error_max << '\n';
    std::cout << "mean_abs_arrival_time_error=" << arrival_error_sum / matched_den << '\n';
    std::cout << "max_abs_arrival_time_error=" << arrival_error_max << '\n';
    std::cout << "problem1_5deg_nearest_sensitivity_ok=1\n";
    return 0;
}
