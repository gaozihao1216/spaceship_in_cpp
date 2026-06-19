/*
 * 文件作用：声明 Problem 1 诊断辅助接口。
 * 主要工作：汇总直接求解和表格采样结果，生成用于分析的统计数据。
 */
#pragma once

#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <map>
#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

struct Problem1SolveSummary {
    // 直接求解候选数量和飞行时间范围。
    int candidate_count;
    double min_time_of_flight_seconds;
    double max_time_of_flight_seconds;
    // 候选残差范围，用于快速判断本次 solve 是否可信。
    double min_relative_residual;
    double max_relative_residual;
};

// 汇总候选解，解决 diagnostics app 不直接遍历细节字段的问题。
Problem1SolveSummary summarize_problem1_candidates(
    const std::vector<Problem1Candidate>& candidates
);

// 将候选解导出 CSV，便于用外部工具画图或做误差分析。
void write_problem1_candidates_csv(
    const std::vector<Problem1Candidate>& candidates,
    const std::string& output_path
);

struct Problem1TableSummary {
    // 表格单元数量和有效几何数量。
    int cell_count;
    int valid_geometry_cell_count;
    int total_branch_count;
    int valid_branch_count;

    // 端点半径范围，用来检查表格覆盖的行星轨道区间是否合理。
    double min_departure_radius;
    double max_departure_radius;
    double min_target_radius;
    double max_target_radius;

    // 有效分支飞行时间范围。
    double min_time_of_flight_seconds;
    double max_time_of_flight_seconds;
    double min_time_of_flight_days;
    double max_time_of_flight_days;

    // 转移轨道参数范围，帮助发现异常几何。
    double min_transfer_e;
    double max_transfer_e;
    double min_transfer_p;
    double max_transfer_p;
};

struct Problem1TableBranchDiagnostics {
    // 分支层面的有效/无效统计。
    int cell_count;
    int valid_cell_count;
    int invalid_cell_count;
    int min_transfer_revolution;
    int max_transfer_revolution;
    int min_target_revolution;
    int max_target_revolution;
    // 按失败原因、圆锥类型、分支数量等维度聚合，定位表格问题。
    std::map<std::string, int> invalid_reason_counts;
    std::map<std::string, int> conic_type_counts;
    std::map<int, int> branch_count_distribution;
    std::map<int, int> branch_count_by_transfer_revolution;
    std::map<int, int> branch_count_by_target_revolution;
    std::map<std::string, int> valid_branch_count_by_pair;
    std::map<std::string, int> invalid_branch_count_by_pair;
    std::map<std::string, int> invalid_reason_counts_by_pair;
};

// 汇总整个表格的几何和时间范围。
Problem1TableSummary summarize_problem1_table(const Problem1Table& table);

// 汇总表格分支拓扑和失败原因分布。
Problem1TableBranchDiagnostics summarize_problem1_table_branches(const Problem1Table& table);

// 将表格单元导出 CSV，解决诊断结果需要离线查看的问题。
void write_problem1_table_csv(
    const Problem1Table& table,
    const std::string& output_path
);

void write_problem1_table_branch_summary_csv(
    const Problem1Table& table,
    const std::string& output_path
);

}  // namespace spaceship_cpp::problem1
