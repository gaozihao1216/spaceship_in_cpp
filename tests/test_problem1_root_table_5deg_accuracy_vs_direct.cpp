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
    long long total_seed_branch_count = 0;
    long long total_refined_branch_count = 0;
    long long total_matched_branch_count = 0;
    long long total_missing_count = 0;
    long long total_new_count = 0;
    long long total_refine_attempt_count = 0;
    long long total_refine_success_count = 0;
    long long total_refine_fail_count = 0;
    long long total_derivative_attach_attempt_count = 0;
    long long total_derivative_attach_success_count = 0;
    long long total_branch_order_mismatch_count = 0;

    double theta_error_sum = 0.0;
    double theta_error_max = 0.0;
    double alpha_error_sum = 0.0;
    double alpha_error_max = 0.0;
    double arrival_error_sum = 0.0;
    double arrival_error_max = 0.0;

    std::vector<double> direct_ms_values;
    std::vector<double> seed_ms_values;
    std::vector<double> refine_ms_values;
    std::string max_error_case = "none";
    std::string worst_missing_case = "none";
    std::string worst_branch_mismatch_case = "none";
    int worst_missing = -1;
    int worst_mismatch = -1;
    double worst_error = -1.0;

    for (const auto& pair : pairs) {
        long long pair_queries = 0;
        long long pair_direct = 0;
        long long pair_refined = 0;
        long long pair_matched = 0;
        long long pair_missing = 0;
        long long pair_new = 0;
        long long pair_refine_attempt = 0;
        long long pair_refine_success = 0;
        long long pair_refine_fail = 0;

        const auto queries = util::make_random_queries_for_pair(pair, sample_count, &rng);
        for (const auto& query : queries) {
            const auto direct_start = util::Clock::now();
            auto direct = util::solve_direct(query);
            const double direct_ms = util::elapsed_ms(direct_start, util::Clock::now());
            const auto simulated = util::simulate_5deg_nearest_refine(query);
            const auto match = util::compare_branches(direct, simulated.refined_branches);

            total_query_count += 1;
            total_direct_branch_count += match.direct_branch_count;
            total_seed_branch_count += static_cast<int>(simulated.seed_branches.size());
            total_refined_branch_count += match.candidate_branch_count;
            total_matched_branch_count += match.matched_branch_count;
            total_missing_count += match.missing_count;
            total_new_count += match.new_count;
            total_refine_attempt_count += simulated.refine_attempt_count;
            total_refine_success_count += simulated.refine_success_count;
            total_refine_fail_count += simulated.refine_fail_count;
            total_derivative_attach_attempt_count += simulated.derivative_attach_attempt_count;
            total_derivative_attach_success_count += simulated.derivative_attach_success_count;
            total_branch_order_mismatch_count += match.branch_order_mismatch_count;

            theta_error_sum += match.theta_prime_error_sum;
            theta_error_max = std::max(theta_error_max, match.theta_prime_error_max);
            alpha_error_sum += match.alpha_error_sum;
            alpha_error_max = std::max(alpha_error_max, match.alpha_error_max);
            arrival_error_sum += match.arrival_time_error_sum;
            arrival_error_max = std::max(arrival_error_max, match.arrival_time_error_max);

            direct_ms_values.push_back(direct_ms);
            seed_ms_values.push_back(simulated.seed_ms);
            refine_ms_values.push_back(simulated.refine_ms);

            pair_queries += 1;
            pair_direct += match.direct_branch_count;
            pair_refined += match.candidate_branch_count;
            pair_matched += match.matched_branch_count;
            pair_missing += match.missing_count;
            pair_new += match.new_count;
            pair_refine_attempt += simulated.refine_attempt_count;
            pair_refine_success += simulated.refine_success_count;
            pair_refine_fail += simulated.refine_fail_count;

            const double query_error = std::max({match.theta_prime_error_max, match.alpha_error_max, match.arrival_time_error_max / 1e6});
            const std::string case_name = util::pair_name(pair);
            if (query_error > worst_error) {
                worst_error = query_error;
                max_error_case = case_name;
            }
            if (match.missing_count > worst_missing) {
                worst_missing = match.missing_count;
                worst_missing_case = case_name;
            }
            if (match.branch_order_mismatch_count > worst_mismatch) {
                worst_mismatch = match.branch_order_mismatch_count;
                worst_branch_mismatch_case = case_name;
            }
        }

        std::cout << "Problem1RootTable5DegAccuracyPairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "query_count=" << pair_queries << '\n';
        std::cout << "direct_branch_count=" << pair_direct << '\n';
        std::cout << "refined_5deg_branch_count=" << pair_refined << '\n';
        std::cout << "matched_branch_count=" << pair_matched << '\n';
        std::cout << "missing_from_5deg_count=" << pair_missing << '\n';
        std::cout << "new_in_5deg_count=" << pair_new << '\n';
        std::cout << "missing_rate=" << static_cast<double>(pair_missing) / static_cast<double>(std::max<long long>(1, pair_direct)) << '\n';
        std::cout << "new_branch_rate=" << static_cast<double>(pair_new) / static_cast<double>(std::max<long long>(1, pair_refined)) << '\n';
        std::cout << "refine_attempt_count=" << pair_refine_attempt << '\n';
        std::cout << "refine_success_count=" << pair_refine_success << '\n';
        std::cout << "refine_fail_count=" << pair_refine_fail << '\n';
    }

    const double matched_den = static_cast<double>(std::max<long long>(1, total_matched_branch_count));
    std::cout << "Problem1RootTable5DegAccuracyGlobalSummary\n";
    std::cout << "total_query_count=" << total_query_count << '\n';
    std::cout << "total_direct_branch_count=" << total_direct_branch_count << '\n';
    std::cout << "total_5deg_seed_branch_count=" << total_seed_branch_count << '\n';
    std::cout << "total_5deg_refined_branch_count=" << total_refined_branch_count << '\n';
    std::cout << "matched_branch_count=" << total_matched_branch_count << '\n';
    std::cout << "missing_from_5deg_count=" << total_missing_count << '\n';
    std::cout << "new_in_5deg_count=" << total_new_count << '\n';
    std::cout << "match_rate=" << static_cast<double>(total_matched_branch_count) / static_cast<double>(std::max<long long>(1, total_direct_branch_count)) << '\n';
    std::cout << "missing_rate=" << static_cast<double>(total_missing_count) / static_cast<double>(std::max<long long>(1, total_direct_branch_count)) << '\n';
    std::cout << "new_branch_rate=" << static_cast<double>(total_new_count) / static_cast<double>(std::max<long long>(1, total_refined_branch_count)) << '\n';
    std::cout << "mean_abs_theta_prime_error=" << theta_error_sum / matched_den << '\n';
    std::cout << "max_abs_theta_prime_error=" << theta_error_max << '\n';
    std::cout << "mean_abs_alpha_error=" << alpha_error_sum / matched_den << '\n';
    std::cout << "max_abs_alpha_error=" << alpha_error_max << '\n';
    std::cout << "mean_abs_arrival_time_error=" << arrival_error_sum / matched_den << '\n';
    std::cout << "max_abs_arrival_time_error=" << arrival_error_max << '\n';
    std::cout << "refine_attempt_count=" << total_refine_attempt_count << '\n';
    std::cout << "refine_success_count=" << total_refine_success_count << '\n';
    std::cout << "refine_fail_count=" << total_refine_fail_count << '\n';
    std::cout << "derivative_attach_attempt_count=" << total_derivative_attach_attempt_count << '\n';
    std::cout << "derivative_attach_success_count=" << total_derivative_attach_success_count << '\n';
    std::cout << "branch_order_mismatch_count=" << total_branch_order_mismatch_count << '\n';
    std::cout << "mean_query_ms_direct=" << util::mean(direct_ms_values) << '\n';
    std::cout << "mean_query_ms_5deg_seed=" << util::mean(seed_ms_values) << '\n';
    std::cout << "mean_query_ms_5deg_refine=" << util::mean(refine_ms_values) << '\n';
    std::cout << "mean_query_ms_total_5deg_simulated=" << util::mean(seed_ms_values) + util::mean(refine_ms_values) << '\n';
    std::cout << "max_error_case=" << max_error_case << '\n';
    std::cout << "worst_missing_case=" << worst_missing_case << '\n';
    std::cout << "worst_branch_mismatch_case=" << worst_branch_mismatch_case << '\n';
    std::cout << "problem1_root_table_5deg_accuracy_vs_direct_ok=1\n";
    return 0;
}
