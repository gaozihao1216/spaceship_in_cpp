/*
 * 文件作用：运行 Problem 1 endpoint transfer-time table 的命令行诊断程序。
 * 主要工作：构造表格配置、采样表格单元，并输出几何和飞行时间分支统计。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_diagnostics.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <filesystem>
#include <iostream>
#include <string>

namespace {

void run_table_diagnostics(
    spaceship_cpp::planet_params::PlanetId departure,
    spaceship_cpp::planet_params::PlanetId target,
    const spaceship_cpp::config::Problem1TableDefaults& defaults
) {
    namespace config = spaceship_cpp::config;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const problem1::Problem1TableConfig table_config =
        config::make_problem1_table_config(departure, target, defaults);

    const problem1::Problem1Table table = problem1::build_problem1_table(table_config);
    const problem1::Problem1TableSummary summary = problem1::summarize_problem1_table(table);

    const std::string departure_name = planet_params::planet_name(departure);
    const std::string target_name = planet_params::planet_name(target);
    const std::string output_path =
        "results/problem1_table_" + departure_name + "_to_" + target_name + ".csv";
    const std::string summary_output_path =
        "results/problem1_table_" + departure_name + "_to_" + target_name + "_summary.csv";

    problem1::write_problem1_table_csv(table, output_path);
    problem1::write_problem1_table_branch_summary_csv(table, summary_output_path);

    std::cout
        << departure_name << " -> " << target_name << '\n'
        << "  cell_count: " << summary.cell_count << '\n'
        << "  valid_geometry_cell_count: " << summary.valid_geometry_cell_count << '\n'
        << "  total_branch_count: " << summary.total_branch_count << '\n'
        << "  valid_branch_count: " << summary.valid_branch_count << '\n'
        << "  min_time_of_flight_days: " << summary.min_time_of_flight_days << '\n'
        << "  max_time_of_flight_days: " << summary.max_time_of_flight_days << '\n'
        << "  min_transfer_e: " << summary.min_transfer_e << '\n'
        << "  max_transfer_e: " << summary.max_transfer_e << '\n'
        << "  min_departure_radius: " << summary.min_departure_radius << '\n'
        << "  max_departure_radius: " << summary.max_departure_radius << '\n'
        << "  csv: " << output_path << '\n'
        << "  summary_csv: " << summary_output_path << '\n';
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;

    std::filesystem::create_directories("results");
    const config::GlobalConfig& cfg = config::global_config();

    config::Problem1TableDefaults table_defaults = cfg.problem1_table_smoke;
    table_defaults.departure_true_anomaly_count = cfg.problem1_diagnostics.table_departure_true_anomaly_count;
    table_defaults.target_true_anomaly_count = cfg.problem1_diagnostics.table_target_true_anomaly_count;
    table_defaults.transfer_theta_departure_count =
        cfg.problem1_diagnostics.table_transfer_theta_departure_count;

    run_table_diagnostics(
        spaceship_cpp::planet_params::PlanetId::Earth,
        spaceship_cpp::planet_params::PlanetId::Mars,
        table_defaults);
    run_table_diagnostics(
        spaceship_cpp::planet_params::PlanetId::Earth,
        spaceship_cpp::planet_params::PlanetId::Venus,
        table_defaults);
    run_table_diagnostics(
        spaceship_cpp::planet_params::PlanetId::Earth,
        spaceship_cpp::planet_params::PlanetId::Mercury,
        table_defaults);

    return 0;
}
