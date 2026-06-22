# spaceship_in_cpp

C++ 重写版行星际轨迹搜索工程。当前核心任务：**给定地球发射时刻，搜索到达水星的最佳转移路径**（直达或经内行星飞掠）。

## 模块

| 模块 | 说明 |
|------|------|
| `problem1` | 两点时间匹配转移求解 |
| `problem2` | 飞掠弹弓 + 外向转移联立求解 |
| `bfs` | **轨迹搜索入口**（`search_best_trajectory`） |
| `visualization` | 轨道绘图（占位框架） |
| `planet_params` / `trajectory` | 行星状态与轨道物理 |

## 快速开始

```bash
cmake -S . -B build
cmake --build build -j4 --target trajectory_search
./build/trajectory_search
```

## 文档

- 任务与 API：[`docs/zh/earth_to_mercury.md`](docs/zh/earth_to_mercury.md)
- 文档索引：[`docs/zh/README.md`](docs/zh/README.md)
- 可视化占位：[`visualization/README.md`](visualization/README.md)

## 公共 API 示例

```cpp
#include "spaceship_cpp/bfs/bfs.hpp"

spaceship_cpp::bfs::TrajectorySearchInput input{};
input.launch_time_seconds_since_j2000 = 0.0;
input.destination_planet = PlanetId::Mercury;

const auto result = spaceship_cpp::bfs::search_best_trajectory(input);
// result.visit_sequence, result.legs[i].eccentricity / semi_latus_rectum_au / ...
```

## 测试

```bash
cmake --build build -j4
ctest --test-dir build --output-on-failure
```
