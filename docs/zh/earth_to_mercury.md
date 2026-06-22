# Earth → Mercury 轨迹搜索

本项目核心任务：**给定地球发射时刻，搜索到达水星（或经内行星飞掠）的最佳转移路径**，并输出可供可视化的轨道参数。

## 问题定义

| 输入 | 说明 |
|------|------|
| `launch_time_seconds_since_j2000` | 从 Earth 出发的发射时刻（J2000 秒） |
| `departure_planet` | 出发行星，当前固定为 **Earth** |
| `destination_planet` | 目标行星，默认 **Mercury** |

| 输出 | 说明 |
|------|------|
| `visit_sequence` | 完整行星访问顺序，如 `Earth → Mercury` 或 `Earth → Venus → Mercury` |
| `legs[]` | 每段转移的 `(e, p, theta)` 及时间 |
| `score` | `v∞,launch + v∞,arrival`（m/s，越小越好） |

## 公共 API（BFS 入口）

**头文件**：`include/spaceship_cpp/bfs/bfs.hpp`  
**实现**：`src/bfs/bfs.cpp`

```cpp
#include "spaceship_cpp/bfs/bfs.hpp"

spaceship_cpp::bfs::TrajectorySearchInput input{};
input.launch_time_seconds_since_j2000 = 0.0;
input.departure_planet = PlanetId::Earth;
input.destination_planet = PlanetId::Mercury;

const auto result = spaceship_cpp::bfs::search_best_trajectory(input);
if (result.found_solution) {
    // result.visit_sequence, result.legs[i].eccentricity / semi_latus_rectum_au / ...
}
```

`search_best_trajectory` 内部执行 **Steps 1–4** 管线，并**只保留 `visit_sequence` 终点为目标行星**的解。

### 与分步 API 的关系

| 层级 | 位置 | 用途 |
|------|------|------|
| **任务 API** | `bfs.hpp` / `bfs.cpp` | 应用与可视化对接：起终点 → 最佳轨道描述 |
| **管线 Steps 1–4** | `free_path_bfs.cpp`, `step3_*.cpp`, `step4_*.cpp` | 内部实现；测试与调试可单独调用 |
| **Problem 1/2** | `problem1/`, `problem2/` | 单段转移与飞掠求解预言机 |

## 单段轨道参数（可视化用）

`TransferLegDescriptor` 中每段转移提供：

```
e   = eccentricity
p   = semi_latus_rectum_au   (AU)
θ   = perihelion_angle_global_rad   (rad，日心惯性系近日点幅角)
k,q = 转移/目标多圈分支编号
t_dep, t_arr, tof = 出发/到达时刻与飞行时间（秒，J2000）
```

在日心平面极坐标下，该段出射椭圆可写为：

\[
r(\nu) = \frac{p}{1 + e\cos(\nu - \theta)}
\]

其中 \(\nu\) 为真近点角（与 `theta` 同一全局参考方向）。

## 命令行运行

```bash
cmake --build build -j4 --target trajectory_search
./build/trajectory_search
```

`apps/trajectory_search.cpp` 默认搜索 **Earth → Mercury**，打印 `TrajectorySearchOutput`。

## 可视化（占位）

目录 [`visualization/`](../../visualization/README.md) 预留绘图框架，尚未实现。计划流程：

1. `search_best_trajectory` → `TrajectorySearchOutput`
2. 导出 JSON（待实现 `visualization/export`）
3. `visualization/python/plot_trajectory.py` 读取并绘图

## 默认搜索配置

来自 `config::global_config().trajectory_search`：

- leg0 θ 粗扫 60 点，细扫 180 点
- 每首段目标 3 个 θ 种子
- `max_search_legs = 6`（最多 6 段转移，允许多行星飞掠路径）
- Problem 1：`k = q = 1`

可在调用 `search_best_trajectory(input, custom_config)` 时覆盖。

## 相关文档

- [BFS 搜索逻辑（分步实现）](bfs_search.md)
- [代码结构](code_structure.md)
