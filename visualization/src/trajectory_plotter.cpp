/*
 * 文件作用：轨迹可视化占位实现。
 */
#include "visualization/trajectory_plotter.hpp"

namespace spaceship_cpp::visualization {

PlotExportResult export_trajectory_for_plot(
    const bfs::TrajectorySearchOutput& /*trajectory*/,
    const std::string& /*output_path*/
) {
    return PlotExportResult{
        .ok = false,
        .error_message = "visualization_not_implemented",
        .output_path = {},
    };
}

}  // namespace spaceship_cpp::visualization
