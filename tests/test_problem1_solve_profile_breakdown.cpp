#include "problem1_5deg_test_utils.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>

namespace {

struct QueryRecord {
    problem1_5deg_test::QueryPoint query;
    spaceship_cpp::problem1::Problem1SolveWithDiagnosticResult result;
    double solve_ms = 0.0;
    int branch_count = 0;
    int valid_branch_count = 0;
    int kq_combo_count = 0;
};

struct BinStats {
    int query_count = 0;
    double solve_ms_sum = 0.0;
};

std::string bin_name(int branch_count) {
    if (branch_count == 0) return "branch_count_0";
    if (branch_count <= 3) return "branch_count_1_to_3";
    if (branch_count <= 8) return "branch_count_4_to_8";
    return "branch_count_gt_8";
}

int valid_branch_count(const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches) {
    int count = 0;
    for (const auto& branch : branches) {
        if (branch.valid) count += 1;
    }
    return count;
}

int kq_combo_count(const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches) {
    std::set<std::pair<int, int>> combos;
    for (const auto& branch : branches) {
        if (!branch.valid) continue;
        combos.insert({branch.transfer_revolution, branch.target_revolution});
    }
    return static_cast<int>(combos.size());
}

void print_distribution(const std::string& prefix, const std::vector<double>& values) {
    std::cout << prefix << "_mean=" << problem1_5deg_test::mean(values) << '\n';
    std::cout << prefix << "_median=" << problem1_5deg_test::percentile(values, 0.50) << '\n';
    std::cout << prefix << "_p90=" << problem1_5deg_test::percentile(values, 0.90) << '\n';
    std::cout << prefix << "_p99=" << problem1_5deg_test::percentile(values, 0.99) << '\n';
    std::cout << prefix << "_max=" << problem1_5deg_test::percentile(values, 1.0) << '\n';
}

void print_hot_case(
    const problem1_5deg_test::PairCase& pair,
    int rank,
    const QueryRecord& record
) {
    namespace planet = spaceship_cpp::planet_params;
    std::cout << "Problem1SolveProfileHotCase\n";
    std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
    std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
    std::cout << "rank=" << rank << '\n';
    std::cout << "nu_A=" << record.query.nu_A << '\n';
    std::cout << "nu_B=" << record.query.nu_B << '\n';
    std::cout << "theta_A=" << record.query.theta_A << '\n';
    std::cout << "solve_ms=" << record.solve_ms << '\n';
    std::cout << "branch_count=" << record.branch_count << '\n';
    std::cout << "valid_branch_count=" << record.valid_branch_count << '\n';
    std::cout << "kq_combo_count=" << record.kq_combo_count << '\n';
    std::cout << "alpha_scan_samples=" << record.result.diagnostic.alpha_scan_samples << '\n';
    std::cout << "residual_evaluations=" << record.result.diagnostic.residual_evaluations << '\n';
    std::cout << "bisection_refinements=" << record.result.diagnostic.bisection_refinements << '\n';
    for (std::size_t i = 0; i < record.result.branches.size(); ++i) {
        const auto& branch = record.result.branches[i];
        std::cout << "Problem1SolveProfileHotCaseBranch\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "rank=" << rank << '\n';
        std::cout << "branch_index=" << i << '\n';
        std::cout << "transfer_revolution=" << branch.transfer_revolution << '\n';
        std::cout << "target_revolution=" << branch.target_revolution << '\n';
        std::cout << "theta_prime=" << branch.target_arrival_true_anomaly << '\n';
        std::cout << "alpha=" << branch.encounter_global_angle << '\n';
        std::cout << "target_time_seconds=" << branch.target_time_seconds << '\n';
        std::cout << "residual_seconds=" << branch.residual_seconds << '\n';
    }
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet = spaceship_cpp::planet_params;
    namespace util = problem1_5deg_test;

    std::cout << std::setprecision(12) << std::scientific;
    const int sample_count = util::env_int("ROOT_TABLE_5DEG_TEST_SAMPLE_COUNT", 100);
    const unsigned seed = util::env_uint("ROOT_TABLE_5DEG_TEST_SEED", 12345);
    const auto pairs = util::selected_pairs_from_env();
    std::mt19937_64 rng(seed);

    std::map<std::string, BinStats> global_bins;
    long long global_query_count = 0;
    double global_solve_ms_sum = 0.0;
    long long global_branch_count_sum = 0;
    long long global_alpha_scan_samples_sum = 0;
    long long global_residual_evaluations_sum = 0;
    long long global_bisection_refinements_sum = 0;
    double global_root_scanning_ms_sum = 0.0;
    double global_residual_eval_ms_sum = 0.0;
    double global_sorting_ms_sum = 0.0;

    for (const auto& pair : pairs) {
        std::vector<QueryRecord> records;
        std::vector<double> solve_ms_values;
        std::vector<double> branch_count_values;
        std::vector<double> root_scan_ms_values;
        std::vector<double> residual_eval_ms_values;
        std::vector<double> sorting_ms_values;
        std::map<std::string, BinStats> pair_bins;

        const auto queries = util::make_random_queries_for_pair(pair, sample_count, &rng);
        records.reserve(queries.size());
        for (const auto& query : queries) {
            const auto start = util::Clock::now();
            auto solved = problem1::solve_problem1_from_departure_anomalies_diagnostic(
                query.from,
                query.to,
                query.nu_A,
                query.nu_B,
                query.theta_A,
                1,
                1);
            const double solve_ms = util::elapsed_ms(start, util::Clock::now());

            QueryRecord record{};
            record.query = query;
            record.result = std::move(solved);
            record.solve_ms = solve_ms;
            record.branch_count = static_cast<int>(record.result.branches.size());
            record.valid_branch_count = valid_branch_count(record.result.branches);
            record.kq_combo_count = kq_combo_count(record.result.branches);
            records.push_back(record);

            solve_ms_values.push_back(solve_ms);
            branch_count_values.push_back(static_cast<double>(record.branch_count));
            root_scan_ms_values.push_back(record.result.diagnostic.root_scanning_seconds * 1000.0);
            residual_eval_ms_values.push_back(record.result.diagnostic.residual_evaluation_seconds * 1000.0);
            sorting_ms_values.push_back(record.result.diagnostic.sorting_conversion_seconds * 1000.0);

            const std::string bin = bin_name(record.branch_count);
            pair_bins[bin].query_count += 1;
            pair_bins[bin].solve_ms_sum += solve_ms;
            global_bins[bin].query_count += 1;
            global_bins[bin].solve_ms_sum += solve_ms;

            global_query_count += 1;
            global_solve_ms_sum += solve_ms;
            global_branch_count_sum += record.branch_count;
            global_alpha_scan_samples_sum += record.result.diagnostic.alpha_scan_samples;
            global_residual_evaluations_sum += record.result.diagnostic.residual_evaluations;
            global_bisection_refinements_sum += record.result.diagnostic.bisection_refinements;
            global_root_scanning_ms_sum += record.result.diagnostic.root_scanning_seconds * 1000.0;
            global_residual_eval_ms_sum += record.result.diagnostic.residual_evaluation_seconds * 1000.0;
            global_sorting_ms_sum += record.result.diagnostic.sorting_conversion_seconds * 1000.0;
        }

        const auto slowest_it = std::max_element(records.begin(), records.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.solve_ms < rhs.solve_ms;
        });

        std::cout << "Problem1SolveProfilePairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "query_count=" << records.size() << '\n';
        std::cout << "mean_solve_ms=" << util::mean(solve_ms_values) << '\n';
        std::cout << "median_solve_ms=" << util::percentile(solve_ms_values, 0.50) << '\n';
        std::cout << "p90_solve_ms=" << util::percentile(solve_ms_values, 0.90) << '\n';
        std::cout << "p99_solve_ms=" << util::percentile(solve_ms_values, 0.99) << '\n';
        std::cout << "max_solve_ms=" << util::percentile(solve_ms_values, 1.0) << '\n';
        std::cout << "mean_branch_count=" << util::mean(branch_count_values) << '\n';
        std::cout << "p90_branch_count=" << util::percentile(branch_count_values, 0.90) << '\n';
        std::cout << "p99_branch_count=" << util::percentile(branch_count_values, 0.99) << '\n';
        std::cout << "max_branch_count=" << util::percentile(branch_count_values, 1.0) << '\n';
        std::cout << "mean_root_scanning_ms=" << util::mean(root_scan_ms_values) << '\n';
        std::cout << "mean_residual_evaluation_ms=" << util::mean(residual_eval_ms_values) << '\n';
        std::cout << "mean_sorting_conversion_ms=" << util::mean(sorting_ms_values) << '\n';
        if (slowest_it != records.end()) {
            std::cout << "slowest_query_nu_A=" << slowest_it->query.nu_A << '\n';
            std::cout << "slowest_query_nu_B=" << slowest_it->query.nu_B << '\n';
            std::cout << "slowest_query_theta_A=" << slowest_it->query.theta_A << '\n';
            std::cout << "slowest_query_branch_count=" << slowest_it->branch_count << '\n';
        }

        for (const auto& [bin, stats] : pair_bins) {
            std::cout << "Problem1SolveProfileBranchCountBinPairSummary\n";
            std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
            std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
            std::cout << "branch_count_bin=" << bin << '\n';
            std::cout << "query_count=" << stats.query_count << '\n';
            std::cout << "mean_solve_ms=" << stats.solve_ms_sum / static_cast<double>(std::max(1, stats.query_count)) << '\n';
        }

        std::sort(records.begin(), records.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.solve_ms > rhs.solve_ms;
        });
        const int hot_count = std::min<int>(5, records.size());
        for (int i = 0; i < hot_count; ++i) {
            print_hot_case(pair, i + 1, records[static_cast<std::size_t>(i)]);
        }
    }

    std::cout << "Problem1SolveProfileBranchCountBinGlobalSummary\n";
    for (const auto& [bin, stats] : global_bins) {
        std::cout << "branch_count_bin=" << bin << '\n';
        std::cout << "query_count=" << stats.query_count << '\n';
        std::cout << "mean_solve_ms=" << stats.solve_ms_sum / static_cast<double>(std::max(1, stats.query_count)) << '\n';
    }

    std::cout << "Problem1SolveProfileGlobalSummary\n";
    std::cout << "query_count=" << global_query_count << '\n';
    std::cout << "mean_solve_ms=" << global_solve_ms_sum / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_branch_count=" << static_cast<double>(global_branch_count_sum) / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_alpha_scan_samples=" << static_cast<double>(global_alpha_scan_samples_sum) / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_residual_evaluations=" << static_cast<double>(global_residual_evaluations_sum) / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_bisection_refinements=" << static_cast<double>(global_bisection_refinements_sum) / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_root_scanning_ms=" << global_root_scanning_ms_sum / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_residual_evaluation_ms=" << global_residual_eval_ms_sum / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "mean_sorting_conversion_ms=" << global_sorting_ms_sum / static_cast<double>(std::max<long long>(1, global_query_count)) << '\n';
    std::cout << "problem1_solve_profile_breakdown_ok=1\n";
    return 0;
}
