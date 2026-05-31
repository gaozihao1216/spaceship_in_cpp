#include "problem1_5deg_test_utils.hpp"

#include <iomanip>
#include <iostream>

namespace {

void print_distribution(const std::string& name, const std::vector<double>& values) {
    std::cout << name << "_mean_ms=" << problem1_5deg_test::mean(values) << '\n';
    std::cout << name << "_median_ms=" << problem1_5deg_test::percentile(values, 0.50) << '\n';
    std::cout << name << "_p90_ms=" << problem1_5deg_test::percentile(values, 0.90) << '\n';
    std::cout << name << "_p99_ms=" << problem1_5deg_test::percentile(values, 0.99) << '\n';
    std::cout << name << "_max_ms=" << problem1_5deg_test::percentile(values, 1.0) << '\n';
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

    std::vector<double> direct_ms;
    std::vector<double> seed_ms;
    std::vector<double> refine_ms;
    std::vector<double> derivative_attach_ms;
    std::vector<double> dedup_ms;
    std::vector<double> matching_ms;
    std::vector<double> simulated_total_ms;
    std::vector<double> estimated_real_5deg_without_io_ms;

    for (const auto& pair : pairs) {
        std::vector<double> pair_direct_ms;
        std::vector<double> pair_simulated_ms;
        const auto queries = util::make_random_queries_for_pair(pair, sample_count, &rng);
        for (const auto& query : queries) {
            const auto direct_start = util::Clock::now();
            auto direct = util::solve_direct(query);
            const double direct_elapsed = util::elapsed_ms(direct_start, util::Clock::now());

            const auto simulated = util::simulate_5deg_nearest_refine(query);
            const auto match_start = util::Clock::now();
            (void)util::compare_branches(direct, simulated.refined_branches);
            const double matching_elapsed = util::elapsed_ms(match_start, util::Clock::now());
            const double simulated_total = simulated.seed_ms + simulated.refine_ms +
                                           simulated.derivative_attach_ms + simulated.dedup_ms + matching_elapsed;
            const double estimated_real = simulated_total - simulated.seed_ms;

            direct_ms.push_back(direct_elapsed);
            seed_ms.push_back(simulated.seed_ms);
            refine_ms.push_back(simulated.refine_ms);
            derivative_attach_ms.push_back(simulated.derivative_attach_ms);
            dedup_ms.push_back(simulated.dedup_ms);
            matching_ms.push_back(matching_elapsed);
            simulated_total_ms.push_back(simulated_total);
            estimated_real_5deg_without_io_ms.push_back(estimated_real);

            pair_direct_ms.push_back(direct_elapsed);
            pair_simulated_ms.push_back(simulated_total);
        }

        std::cout << "Problem1FiveDegSpeedBreakdownPairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "query_count=" << pair_direct_ms.size() << '\n';
        std::cout << "mean_direct_solve_ms=" << util::mean(pair_direct_ms) << '\n';
        std::cout << "mean_simulated_5deg_total_ms=" << util::mean(pair_simulated_ms) << '\n';
    }

    std::cout << "Problem1FiveDegSpeedBreakdownGlobalSummary\n";
    std::cout << "query_count=" << direct_ms.size() << '\n';
    print_distribution("direct_solve", direct_ms);
    print_distribution("nearest_5deg_node_solve", seed_ms);
    print_distribution("refine", refine_ms);
    print_distribution("derivative_attach", derivative_attach_ms);
    print_distribution("dedup", dedup_ms);
    print_distribution("matching", matching_ms);
    print_distribution("simulated_5deg_total", simulated_total_ms);
    print_distribution("estimated_real_5deg_query_without_io", estimated_real_5deg_without_io_ms);
    std::cout << "estimated_table_read_time_ms_assumed=0\n";
    std::cout << "problem1_5deg_speed_breakdown_ok=1\n";
    return 0;
}
