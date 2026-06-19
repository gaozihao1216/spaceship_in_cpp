/*
 * 文件作用：实现 Problem 1 诊断统计工具。
 * 主要工作：对直接求解和表格分支进行采样、汇总和一致性检查。
 */
#include "spaceship_cpp/problem1/problem1_diagnostics.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>

namespace spaceship_cpp::problem1 {

namespace {

constexpr double kSecondsPerDay = 86400.0;

double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

double tof_days(double time_of_flight_seconds) {
    return time_of_flight_seconds / kSecondsPerDay;
}

const char* conic_type_label(Problem1TransferConicType conic_type) {
    switch (conic_type) {
        case Problem1TransferConicType::Invalid:
            return "invalid";
        case Problem1TransferConicType::Elliptic:
            return "elliptic";
        case Problem1TransferConicType::Parabolic:
            return "parabolic";
        case Problem1TransferConicType::Hyperbolic:
            return "hyperbolic";
    }
    return "unknown";
}

const char* non_empty_label(const std::string& value, const char* fallback) {
    // 中文注释：CSV 里把空字符串归并到稳定标签，便于后续统计和筛选。
    return value.empty() ? fallback : value.c_str();
}

std::string branch_pair_label(int transfer_revolution, int target_revolution) {
    // 中文注释：统一 (k,q) 标签格式，便于 summary CSV 和测试代码共享同一语义。
    return "k=" + std::to_string(transfer_revolution) + ";q=" + std::to_string(target_revolution);
}

void write_nan_branch_columns(std::ofstream& output) {
    output
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan,"
        << "nan";
}

}  // namespace

// 汇总候选解的飞行时间范围和相对残差范围。
Problem1SolveSummary summarize_problem1_candidates(const std::vector<Problem1Candidate>& candidates) {
    if (candidates.empty()) {
        const double nan = nan_value();
        return Problem1SolveSummary{0, nan, nan, nan, nan};
    }

    Problem1SolveSummary summary{
        static_cast<int>(candidates.size()),
        candidates.front().time_of_flight_seconds,
        candidates.front().time_of_flight_seconds,
        candidates.front().relative_residual,
        candidates.front().relative_residual,
    };

    for (const Problem1Candidate& candidate : candidates) {
        summary.min_time_of_flight_seconds =
            std::min(summary.min_time_of_flight_seconds, candidate.time_of_flight_seconds);
        summary.max_time_of_flight_seconds =
            std::max(summary.max_time_of_flight_seconds, candidate.time_of_flight_seconds);
        summary.min_relative_residual =
            std::min(summary.min_relative_residual, candidate.relative_residual);
        summary.max_relative_residual =
            std::max(summary.max_relative_residual, candidate.relative_residual);
    }

    return summary;
}

// 将候选解列表导出为 CSV，便于外部工具画图或误差分析。
void write_problem1_candidates_csv(
    const std::vector<Problem1Candidate>& candidates,
    const std::string& output_path
) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("failed to open output CSV file");
    }

    output
        << "index,"
        << "transfer_revolution,"
        << "target_revolution,"
        << "encounter_global_angle,"
        << "launch_time_seconds_since_j2000,"
        << "time_of_flight_seconds,"
        << "arrival_time_seconds_since_j2000,"
        << "transfer_e_raw,"
        << "transfer_e,"
        << "transfer_perihelion_angle_used,"
        << "transfer_p,"
        << "relative_residual,"
        << "residual_scale,"
        << "root_bracket_width,"
        << "bisection_iterations,"
        << "refined_by_bisection,"
        << "residual,"
        << "deltaF_transfer,"
        << "deltaF_target\n";

    for (std::size_t index = 0; index < candidates.size(); ++index) {
        const Problem1Candidate& candidate = candidates[index];
        output
            << index << ','
            << candidate.transfer_revolution << ','
            << candidate.target_revolution << ','
            << candidate.encounter_global_angle << ','
            << candidate.launch_time_seconds_since_j2000 << ','
            << candidate.time_of_flight_seconds << ','
            << candidate.arrival_time_seconds_since_j2000 << ','
            << candidate.residual_result.transfer_e_raw << ','
            << candidate.residual_result.transfer_e << ','
            << candidate.residual_result.transfer_perihelion_angle_used << ','
            << candidate.residual_result.transfer_p << ','
            << candidate.relative_residual << ','
            << candidate.residual_scale << ','
            << candidate.root_bracket_width << ','
            << candidate.bisection_iterations << ','
            << (candidate.refined_by_bisection ? 1 : 0) << ','
            << candidate.residual_result.residual << ','
            << candidate.residual_result.deltaF_transfer << ','
            << candidate.residual_result.deltaF_target << '\n';
    }
}

// 汇总表格的几何有效性、分支数量和飞行时间/偏心率范围。
Problem1TableSummary summarize_problem1_table(const Problem1Table& table) {
    const std::vector<Problem1TableCell>& cells = table.cells();
    const double nan = nan_value();

    Problem1TableSummary summary{};
    summary.cell_count = static_cast<int>(cells.size());
    summary.valid_geometry_cell_count = 0;
    summary.total_branch_count = 0;
    summary.valid_branch_count = 0;
    summary.min_departure_radius = nan;
    summary.max_departure_radius = nan;
    summary.min_target_radius = nan;
    summary.max_target_radius = nan;
    summary.min_time_of_flight_seconds = nan;
    summary.max_time_of_flight_seconds = nan;
    summary.min_time_of_flight_days = nan;
    summary.max_time_of_flight_days = nan;
    summary.min_transfer_e = nan;
    summary.max_transfer_e = nan;
    summary.min_transfer_p = nan;
    summary.max_transfer_p = nan;

    bool has_geometry_stats = false;
    bool has_tof_stats = false;
    for (const Problem1TableCell& cell : cells) {
        summary.total_branch_count += static_cast<int>(cell.time_of_flight_branches.size());
        if (!cell.valid) {
            continue;
        }
        ++summary.valid_geometry_cell_count;
        if (!has_geometry_stats) {
            summary.min_departure_radius = cell.departure_radius;
            summary.max_departure_radius = cell.departure_radius;
            summary.min_target_radius = cell.target_radius;
            summary.max_target_radius = cell.target_radius;
            summary.min_transfer_e = cell.transfer_e;
            summary.max_transfer_e = cell.transfer_e;
            summary.min_transfer_p = cell.transfer_p;
            summary.max_transfer_p = cell.transfer_p;
            has_geometry_stats = true;
        } else {
            summary.min_departure_radius = std::min(summary.min_departure_radius, cell.departure_radius);
            summary.max_departure_radius = std::max(summary.max_departure_radius, cell.departure_radius);
            summary.min_target_radius = std::min(summary.min_target_radius, cell.target_radius);
            summary.max_target_radius = std::max(summary.max_target_radius, cell.target_radius);
            summary.min_transfer_e = std::min(summary.min_transfer_e, cell.transfer_e);
            summary.max_transfer_e = std::max(summary.max_transfer_e, cell.transfer_e);
            summary.min_transfer_p = std::min(summary.min_transfer_p, cell.transfer_p);
            summary.max_transfer_p = std::max(summary.max_transfer_p, cell.transfer_p);
        }

        for (const Problem1TimeOfFlightBranch& branch : cell.time_of_flight_branches) {
            if (!branch.valid) {
                continue;
            }
            ++summary.valid_branch_count;
            const double days = tof_days(branch.time_of_flight_seconds);
            if (!has_tof_stats) {
                summary.min_time_of_flight_seconds = branch.time_of_flight_seconds;
                summary.max_time_of_flight_seconds = branch.time_of_flight_seconds;
                summary.min_time_of_flight_days = days;
                summary.max_time_of_flight_days = days;
                has_tof_stats = true;
            } else {
                summary.min_time_of_flight_seconds =
                    std::min(summary.min_time_of_flight_seconds, branch.time_of_flight_seconds);
                summary.max_time_of_flight_seconds =
                    std::max(summary.max_time_of_flight_seconds, branch.time_of_flight_seconds);
                summary.min_time_of_flight_days = std::min(summary.min_time_of_flight_days, days);
                summary.max_time_of_flight_days = std::max(summary.max_time_of_flight_days, days);
            }
        }
    }

    return summary;
}

// 按 (k,q) 对、圆锥类型、失败原因等维度聚合分支拓扑统计。
Problem1TableBranchDiagnostics summarize_problem1_table_branches(const Problem1Table& table) {
    Problem1TableBranchDiagnostics diagnostics{};
    diagnostics.cell_count = static_cast<int>(table.cells().size());
    diagnostics.valid_cell_count = 0;
    diagnostics.invalid_cell_count = 0;
    diagnostics.min_transfer_revolution = std::numeric_limits<int>::max();
    diagnostics.max_transfer_revolution = std::numeric_limits<int>::min();
    diagnostics.min_target_revolution = std::numeric_limits<int>::max();
    diagnostics.max_target_revolution = std::numeric_limits<int>::min();

    for (const Problem1TableCell& cell : table.cells()) {
        if (cell.valid) {
            ++diagnostics.valid_cell_count;
        } else {
            ++diagnostics.invalid_cell_count;
            const std::string key = cell.invalid_reason.empty() ? "none" : cell.invalid_reason;
            ++diagnostics.invalid_reason_counts[key];
        }

        ++diagnostics.conic_type_counts[conic_type_label(cell.conic_type)];
        ++diagnostics.branch_count_distribution[static_cast<int>(cell.time_of_flight_branches.size())];
        for (const Problem1TimeOfFlightBranch& branch : cell.time_of_flight_branches) {
            diagnostics.min_transfer_revolution =
                std::min(diagnostics.min_transfer_revolution, branch.transfer_revolution);
            diagnostics.max_transfer_revolution =
                std::max(diagnostics.max_transfer_revolution, branch.transfer_revolution);
            diagnostics.min_target_revolution =
                std::min(diagnostics.min_target_revolution, branch.target_revolution);
            diagnostics.max_target_revolution =
                std::max(diagnostics.max_target_revolution, branch.target_revolution);
            ++diagnostics.branch_count_by_transfer_revolution[branch.transfer_revolution];
            ++diagnostics.branch_count_by_target_revolution[branch.target_revolution];
            const std::string pair_key =
                branch_pair_label(branch.transfer_revolution, branch.target_revolution);
            if (branch.valid) {
                ++diagnostics.valid_branch_count_by_pair[pair_key];
            } else {
                ++diagnostics.invalid_branch_count_by_pair[pair_key];
                const std::string invalid_key =
                    pair_key + "|reason=" + (branch.invalid_reason.empty() ? "none" : branch.invalid_reason);
                ++diagnostics.invalid_reason_counts_by_pair[invalid_key];
            }
        }
    }

    if (diagnostics.branch_count_by_transfer_revolution.empty()) {
        diagnostics.min_transfer_revolution = 0;
        diagnostics.max_transfer_revolution = -1;
    }
    if (diagnostics.branch_count_by_target_revolution.empty()) {
        diagnostics.min_target_revolution = 0;
        diagnostics.max_target_revolution = -1;
    }

    return diagnostics;
}

// 将表格所有单元及其分支导出为扁平 CSV（每行一个分支）。
void write_problem1_table_csv(const Problem1Table& table, const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("failed to open output CSV file");
    }

    output
        << "schema_version,"
        << "departure_planet,"
        << "target_planet,"
        << "departure_true_anomaly_index,"
        << "target_true_anomaly_index,"
        << "transfer_theta_departure_index,"
        << "nu_A_input,"
        << "nu_B_input,"
        << "theta_A_input,"
        << "nu_A,"
        << "nu_B,"
        << "theta_A,"
        << "theta_B,"
        << "r_A,"
        << "r_B,"
        << "lambda_A,"
        << "lambda_B,"
        << "delta_lambda,"
        << "transfer_perihelion_angle_global_raw,"
        << "transfer_perihelion_angle_global_used,"
        << "transfer_e_raw,"
        << "transfer_e,"
        << "transfer_p,"
        << "transfer_a,"
        << "conic_type,"
        << "normalized_negative_e,"
        << "cell_valid,"
        << "cell_invalid_reason,"
        << "branch_index,"
        << "branch_count,"
        << "branch_valid,"
        << "transfer_revolution,"
        << "target_revolution,"
        << "theta_B_branch,"
        << "target_theta_start,"
        << "target_theta_end_branch,"
        << "deltaF_transfer,"
        << "deltaF_target,"
        << "time_of_flight_scale_free,"
        << "target_time_of_flight_scale_free,"
        << "time_of_flight_seconds,"
        << "target_time_of_flight_seconds,"
        << "residual_scale_free,"
        << "residual_seconds,"
        << "time_of_flight_days,"
        << "branch_invalid_reason\n";

    const Problem1TableMetadata& metadata = table.metadata();
    const char* departure_name = planet_params::planet_name(table.config().departure_planet);
    const char* target_name = planet_params::planet_name(table.config().target_planet);

    for (const Problem1TableCell& cell : table.cells()) {
        if (cell.time_of_flight_branches.empty()) {
            output
                << metadata.schema_version << ','
                << departure_name << ','
                << target_name << ','
                << cell.departure_true_anomaly_index << ','
                << cell.target_true_anomaly_index << ','
                << cell.transfer_theta_departure_index << ','
                << cell.departure_true_anomaly_input << ','
                << cell.target_true_anomaly_input << ','
                << cell.transfer_theta_departure_input << ','
                << cell.departure_true_anomaly << ','
                << cell.target_true_anomaly << ','
                << cell.transfer_theta_departure << ','
                << cell.transfer_theta_arrival << ','
                << cell.departure_radius << ','
                << cell.target_radius << ','
                << cell.departure_global_angle << ','
                << cell.target_global_angle << ','
                << cell.delta_global_angle << ','
                << cell.transfer_perihelion_angle_global_raw << ','
                << cell.transfer_perihelion_angle_global_used << ','
                << cell.transfer_e_raw << ','
                << cell.transfer_e << ','
                << cell.transfer_p << ','
                << cell.transfer_a << ','
                << conic_type_label(cell.conic_type) << ','
                << (cell.normalized_negative_e ? 1 : 0) << ','
                << (cell.valid ? 1 : 0) << ','
                << cell.invalid_reason << ','
                << -1 << ','
                << 0 << ',';
            write_nan_branch_columns(output);
            output << '\n';
            continue;
        }

        for (std::size_t branch_index = 0; branch_index < cell.time_of_flight_branches.size(); ++branch_index) {
            const Problem1TimeOfFlightBranch& branch = cell.time_of_flight_branches[branch_index];
            output
                << metadata.schema_version << ','
                << departure_name << ','
                << target_name << ','
                << cell.departure_true_anomaly_index << ','
                << cell.target_true_anomaly_index << ','
                << cell.transfer_theta_departure_index << ','
                << cell.departure_true_anomaly_input << ','
                << cell.target_true_anomaly_input << ','
                << cell.transfer_theta_departure_input << ','
                << cell.departure_true_anomaly << ','
                << cell.target_true_anomaly << ','
                << cell.transfer_theta_departure << ','
                << cell.transfer_theta_arrival << ','
                << cell.departure_radius << ','
                << cell.target_radius << ','
                << cell.departure_global_angle << ','
                << cell.target_global_angle << ','
                << cell.delta_global_angle << ','
                << cell.transfer_perihelion_angle_global_raw << ','
                << cell.transfer_perihelion_angle_global_used << ','
                << cell.transfer_e_raw << ','
                << cell.transfer_e << ','
                << cell.transfer_p << ','
                << cell.transfer_a << ','
                << conic_type_label(cell.conic_type) << ','
                << (cell.normalized_negative_e ? 1 : 0) << ','
                << (cell.valid ? 1 : 0) << ','
                << cell.invalid_reason << ','
                << branch_index << ','
                << cell.time_of_flight_branches.size() << ','
                << (branch.valid ? 1 : 0) << ','
                << branch.transfer_revolution << ','
                << branch.target_revolution << ','
                << branch.theta_arrival_branch << ','
                << branch.target_true_anomaly_start << ','
                << branch.target_true_anomaly_end_branch << ','
                << branch.deltaF_transfer << ','
                << branch.deltaF_target << ','
                << branch.time_of_flight_scale_free << ','
                << branch.target_time_of_flight_scale_free << ','
                << branch.time_of_flight_seconds << ','
                << branch.target_time_of_flight_seconds << ','
                << branch.residual_scale_free << ','
                << branch.residual_seconds << ','
                << tof_days(branch.time_of_flight_seconds) << ','
                << branch.invalid_reason << '\n';
        }
    }
}

// 将分支拓扑聚合统计导出为 summary CSV，便于快速判断表格质量。
void write_problem1_table_branch_summary_csv(
    const Problem1Table& table,
    const std::string& output_path
) {
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("failed to open output summary CSV file");
    }

    const Problem1TableBranchDiagnostics diagnostics = summarize_problem1_table_branches(table);
    const Problem1TableMetadata& metadata = table.metadata();

    output
        << "row_type,"
        << "schema_version,"
        << "departure_planet,"
        << "target_planet,"
        << "key,"
        << "count\n";

    // 中文注释：先写整体计数，方便脚本不扫描整表也能快速判断表质量。
    output
        << "overall,"
        << metadata.schema_version << ','
        << metadata.departure_planet_name << ','
        << metadata.target_planet_name << ','
        << "cell_count,"
        << diagnostics.cell_count << '\n';
    output
        << "overall,"
        << metadata.schema_version << ','
        << metadata.departure_planet_name << ','
        << metadata.target_planet_name << ','
        << "valid_cell_count,"
        << diagnostics.valid_cell_count << '\n';
    output
        << "overall,"
        << metadata.schema_version << ','
        << metadata.departure_planet_name << ','
        << metadata.target_planet_name << ','
        << "invalid_cell_count,"
        << diagnostics.invalid_cell_count << '\n';
    output
        << "overall,"
        << metadata.schema_version << ','
        << metadata.departure_planet_name << ','
        << metadata.target_planet_name << ','
        << "transfer_revolution_range,"
        << diagnostics.min_transfer_revolution << ':' << diagnostics.max_transfer_revolution << '\n';
    output
        << "overall,"
        << metadata.schema_version << ','
        << metadata.departure_planet_name << ','
        << metadata.target_planet_name << ','
        << "target_revolution_range,"
        << diagnostics.min_target_revolution << ':' << diagnostics.max_target_revolution << '\n';

    // 中文注释：invalid reason、conic type、branch 数量分布分别单独成段，便于后续聚合统计。
    for (const auto& [reason, count] : diagnostics.invalid_reason_counts) {
        output
            << "invalid_reason,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << non_empty_label(reason, "none") << ','
            << count << '\n';
    }

    for (const auto& [conic_type, count] : diagnostics.conic_type_counts) {
        output
            << "conic_type,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << non_empty_label(conic_type, "unknown") << ','
            << count << '\n';
    }

    for (const auto& [branch_count, count] : diagnostics.branch_count_distribution) {
        output
            << "branch_count_distribution,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << branch_count << ','
            << count << '\n';
    }

    for (const auto& [transfer_revolution, count] : diagnostics.branch_count_by_transfer_revolution) {
        output
            << "branch_count_by_transfer_revolution,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << transfer_revolution << ','
            << count << '\n';
    }

    for (const auto& [target_revolution, count] : diagnostics.branch_count_by_target_revolution) {
        output
            << "branch_count_by_target_revolution,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << target_revolution << ','
            << count << '\n';
    }

    for (const auto& [pair_key, count] : diagnostics.valid_branch_count_by_pair) {
        output
            << "valid_branch_count_by_pair,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << pair_key << ','
            << count << '\n';
    }

    for (const auto& [pair_key, count] : diagnostics.invalid_branch_count_by_pair) {
        output
            << "invalid_branch_count_by_pair,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << pair_key << ','
            << count << '\n';
    }

    for (const auto& [pair_reason_key, count] : diagnostics.invalid_reason_counts_by_pair) {
        output
            << "invalid_reason_by_pair,"
            << metadata.schema_version << ','
            << metadata.departure_planet_name << ','
            << metadata.target_planet_name << ','
            << pair_reason_key << ','
            << count << '\n';
    }
}

}  // namespace spaceship_cpp::problem1
