/*
 * 文件作用：测试 Problem 1 表格诊断工具。
 * 主要工作：验证诊断采样、统计汇总和直接求解对照流程可以正常运行。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_diagnostics.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const problem1::Problem1TableConfig config{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        0.0,
        spaceship_cpp::common::kPi,
        2,
        0.0,
        spaceship_cpp::common::kPi,
        2,
        0.0,
        spaceship_cpp::common::kHalfPi,
        3,
        1,
        1,
    };

    const problem1::Problem1Table table = problem1::build_problem1_table(config);
    const problem1::Problem1TableSummary summary = problem1::summarize_problem1_table(table);
    const problem1::Problem1TableBranchDiagnostics branch_diagnostics =
        problem1::summarize_problem1_table_branches(table);

    assert(summary.cell_count == static_cast<int>(table.cells().size()));
    assert(summary.valid_geometry_cell_count >= 0);
    assert(summary.valid_geometry_cell_count <= summary.cell_count);
    assert(summary.total_branch_count >= 0);
    assert(summary.valid_branch_count >= 0);
    assert(summary.valid_branch_count <= summary.total_branch_count);

    if (summary.valid_geometry_cell_count > 0) {
        assert(std::isfinite(summary.min_departure_radius));
        assert(std::isfinite(summary.max_departure_radius));
        assert(std::isfinite(summary.min_transfer_e));
        assert(std::isfinite(summary.max_transfer_e));
    }

    assert(branch_diagnostics.cell_count == summary.cell_count);
    assert(branch_diagnostics.valid_cell_count + branch_diagnostics.invalid_cell_count == summary.cell_count);
    assert(!branch_diagnostics.conic_type_counts.empty());
    assert(!branch_diagnostics.branch_count_distribution.empty());
    assert(branch_diagnostics.min_transfer_revolution == 0);
    assert(branch_diagnostics.max_transfer_revolution == 1);
    assert(branch_diagnostics.min_target_revolution == 0);
    assert(branch_diagnostics.max_target_revolution == 1);
    assert(!branch_diagnostics.branch_count_by_transfer_revolution.empty());
    assert(!branch_diagnostics.branch_count_by_target_revolution.empty());

    std::filesystem::create_directories("results");
    const std::filesystem::path output_path = "results/test_problem1_table_diagnostics.csv";
    const std::filesystem::path summary_output_path = "results/test_problem1_table_diagnostics_summary.csv";
    problem1::write_problem1_table_csv(table, output_path.string());
    problem1::write_problem1_table_branch_summary_csv(table, summary_output_path.string());
    assert(std::filesystem::exists(output_path));
    assert(std::filesystem::file_size(output_path) > 0);
    assert(std::filesystem::exists(summary_output_path));
    assert(std::filesystem::file_size(summary_output_path) > 0);

    return 0;
}
