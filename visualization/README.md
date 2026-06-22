# 轨迹可视化（框架）

本目录用于 **Earth → 目标行星** 转移轨道的图像展示，与 `spaceship_cpp::bfs::search_best_trajectory` 的输出对接。

## 当前状态

**仅占位框架，尚未实现绘图。**

## 计划数据流

```
TrajectorySearchInput (t_launch, Earth, Mercury)
        │
        ▼
  search_best_trajectory()   ← src/bfs/bfs.cpp
        │
        ▼
  TrajectorySearchOutput
    · visit_sequence
    · legs[]: (e, p, theta, t_dep, t_arr, …)
        │
        ▼
  visualization/             ← 本目录（待实现）
    · 读取 legs 绘制日心转移椭圆
    · 叠加行星轨道与到达时刻位置
```

## 目录规划

| 路径 | 用途 |
|------|------|
| `include/visualization/trajectory_plotter.hpp` | C++ 导出接口占位（可选 JSON/CSV 中转） |
| `src/trajectory_plotter.cpp` | 占位实现 |
| `python/plot_trajectory.py` | Python + matplotlib 绘图脚本占位 |

## 轨道参数说明（供可视化使用）

每段 `TransferLegDescriptor` 中：

- `eccentricity` (e)、`semi_latus_rectum_au` (p)、`perihelion_angle_global_rad` (θ)：定义该段 **出射转移椭圆**（日心惯性系）
- `departure_time_seconds_since_j2000` / `arrival_time_seconds_since_j2000`：在该时间段内绘制弧段
- `visit_sequence`：行星访问顺序，用于标注 flyby 节点

详细任务说明见 [`docs/zh/earth_to_mercury.md`](../docs/zh/earth_to_mercury.md)。

## 后续实现建议

1. 在 `trajectory_search` 或单独工具中将 `TrajectorySearchOutput` 序列化为 JSON
2. `python/plot_trajectory.py` 读取 JSON，用 `(e,p,θ)` 生成椭圆点列并绘图
3. 可选：将 `visualization` 编为独立 CMake target，依赖 `spaceship_cpp`
