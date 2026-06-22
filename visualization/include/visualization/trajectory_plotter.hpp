/*
 * 文件作用：轨迹可视化模块公共接口（占位）。
 * 主要工作：未来将 TrajectorySearchOutput 转为绘图数据或导出文件；当前未实现。
 */
#pragma once

#include "spaceship_cpp/bfs/trajectory_solution.hpp"

#include <string>

namespace spaceship_cpp::visualization {

struct PlotExportResult {
    bool ok = false;
    std::string error_message;
    std::string output_path;
};

// 占位：将搜索结果导出为可视化中间格式（如 JSON）。尚未实现。
PlotExportResult export_trajectory_for_plot(
    const bfs::TrajectorySearchOutput& trajectory,
    const std::string& output_path
);

}  // namespace spaceship_cpp::visualization
