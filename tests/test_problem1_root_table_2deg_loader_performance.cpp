#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

std::vector<long long> deterministic_indices(long long total_nodes, int count) {
    std::vector<long long> indices;
    indices.reserve(static_cast<std::size_t>(count));
    std::uint64_t state = 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < count; ++i) {
        state = state * 2862933555777941757ULL + 3037000493ULL;
        indices.push_back(static_cast<long long>(state % static_cast<std::uint64_t>(total_nodes)));
    }
    return indices;
}

double percentile(std::vector<double> values, double fraction) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const std::size_t index = std::min<std::size_t>(
        values.size() - 1,
        static_cast<std::size_t>(fraction * static_cast<double>(values.size() - 1)));
    return values[index];
}

void run_round(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const std::vector<long long>& indices,
    const std::string& round
) {
    using Clock = std::chrono::steady_clock;
    std::vector<double> query_ms;
    query_ms.reserve(indices.size());
    long long branch_count_sum = 0;

    const auto total_start = Clock::now();
    for (const long long linear_index : indices) {
        const auto query_start = Clock::now();
        const auto node = loader.load_node_by_linear_index(linear_index);
        const auto query_end = Clock::now();
        assert(node.valid);
        query_ms.push_back(std::chrono::duration<double, std::milli>(query_end - query_start).count());
        branch_count_sum += static_cast<long long>(node.branches.size());
    }
    const auto total_end = Clock::now();

    const double total_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    const double mean_ms = query_ms.empty()
        ? 0.0
        : std::accumulate(query_ms.begin(), query_ms.end(), 0.0) / static_cast<double>(query_ms.size());
    const double max_ms = query_ms.empty() ? 0.0 : *std::max_element(query_ms.begin(), query_ms.end());

    std::cout << "Problem1RootTable2DegLoaderPerformanceSummary\n";
    std::cout << "round=" << round << '\n';
    std::cout << "query_count=" << indices.size() << '\n';
    std::cout << "total_ms=" << total_ms << '\n';
    std::cout << "mean_ms=" << mean_ms << '\n';
    std::cout << "max_ms=" << max_ms << '\n';
    std::cout << "p50_ms=" << percentile(query_ms, 0.50) << '\n';
    std::cout << "p90_ms=" << percentile(query_ms, 0.90) << '\n';
    std::cout << "p99_ms=" << percentile(query_ms, 0.99) << '\n';
    std::cout << "branch_count_sum=" << branch_count_sum << '\n';
    std::cout << "loader_performance_ok=1\n";
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem1_root_table_2deg_loader_performance_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto indices = deterministic_indices(loader.total_nodes(), 1000);
    run_round(loader, indices, "cold");
    run_round(loader, indices, "warm");
    return 0;
}
