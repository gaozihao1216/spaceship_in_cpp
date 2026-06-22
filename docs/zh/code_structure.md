# 代码结构

本文档描述在移除根表（root-table）预计算路径之后，当前源码的目录布局。

## 顶层目录

- `CMakeLists.txt`：构建 `spaceship_cpp` 静态库、`trajectory_search` 主程序与测试套件。
- `README.md`：项目简介与构建命令。
- `apps/trajectory_search.cpp`：Earth → Mercury 命令行入口（调用 `search_best_trajectory`）。
- `visualization/`：轨迹可视化占位框架（尚未接入构建）。
- `scripts/build.sh`：项目构建脚本。
- `scripts/run_tests.sh`：测试运行脚本。

## 公共头文件（`include/`）

- `include/spaceship_cpp/common/common.hpp`：共享常量、角度归一化、数值辅助函数、轨道异常角工具。
- `include/spaceship_cpp/common/orbit_math.hpp`：Problem 1 和轨迹模块共用的轨道数学工具。
- `include/spaceship_cpp/common/time_utils.hpp`：时间单位转换工具。
- `include/spaceship_cpp/config/global_config.hpp`：全局默认配置访问接口。
- `include/spaceship_cpp/core_types/core_types.hpp`：共享核心类型定义（当前为占位）。
- `include/spaceship_cpp/planet_params/planet_params.hpp`：行星编号、物理参数、指定时刻行星状态查询。
- `include/spaceship_cpp/problem1/orbit_time_integral.hpp`：Problem 1 使用的轨道时间积分辅助函数（转发自 common）。
- `include/spaceship_cpp/problem1/problem1.hpp`：Problem 1 直接残差评估与二分法求解器 API。
- `include/spaceship_cpp/problem1/problem1_diagnostics.hpp`：Problem 1 诊断结果类型与汇总/导出工具。
- `include/spaceship_cpp/problem1/problem1_table.hpp`：端点几何与转移时间表格 API（实验性保留）。
- `include/spaceship_cpp/problem2/problem2.hpp`：Problem 2 模块的总入口头文件。
- `include/spaceship_cpp/problem2/problem2_slingshot.hpp`：Problem 2 飞掠弹弓残差方程。
- `include/spaceship_cpp/trajectory/flyby_physics.hpp`：物理飞掠可行性计算。
- `include/spaceship_cpp/trajectory/orbit_velocity.hpp`：轨道速度与行星相对速度计算。
- `include/spaceship_cpp/bfs/bfs.hpp`：**BFS 公共入口**（`search_best_trajectory`）。
- `include/spaceship_cpp/bfs/trajectory_solution.hpp`：任务输入/输出与 `TransferLegDescriptor`。
- `include/spaceship_cpp/bfs/trajectory_search_config.hpp`：Steps 1–4 管线配置。
- `include/spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp`：全局近日点角与 Problem 2 局部角之间的转换。
- `include/spaceship_cpp/bfs/trajectory_search_state.hpp`：轨迹搜索的基础状态与边数据结构。

## 源码布局（`src/`）

- `src/common/common.cpp`：通用数学工具实现（角度归一化、clamp、safe_divide 等）。
- `src/common/orbit_math.cpp`：轨道时间积分 F(e,ξ) 及导数实现。
- `src/common/time_utils.cpp`：UTC 与 J2000 秒偏移互转实现。
- `src/config/global_config.cpp`：全局默认配置实现。
- `src/core_types/core_types.cpp`：核心类型占位实现。
- `src/planet_params/planet_params.cpp`：八大行星参数表与按时间求行星状态。
- `src/problem1/problem1.cpp`：Problem 1 残差评估与扫描+二分求解器实现。
- `src/problem1/diagnostics/problem1_diagnostics.cpp`：Problem 1 诊断统计与 CSV 导出。
- `src/problem1/tables/problem1_table.cpp`：端点几何与转移时间表格实现。
- `src/problem2/problem2.cpp`：Problem 2 模块占位实现。
- `src/problem2/flyby/problem2_slingshot.cpp`：Problem 2 弹弓残差方程实现。
- `src/trajectory/flyby_physics.cpp`：飞掠物理可行性判断实现。
- `src/trajectory/orbit_velocity.cpp`：轨道速度计算实现。
- `src/bfs/bfs.cpp`：**BFS 公共 API 实现**（封装 Steps 1–4，按目标行星筛选最优解）。
- `src/bfs/free_path_bfs.cpp`：Step 2 自由路径 BFS。
- `src/bfs/step3_top_k_sequences.cpp`：Step 3 Top-K 序列筛选。
- `src/bfs/step4_fixed_sequence_fine_search.cpp`：Step 4 细 θ + 固定序列精搜。
- `src/bfs/fixed_sequence_bfs.cpp`：固定行星序列 BFS。
- `src/bfs/leg0_theta_feasibility.cpp`：Step 1 leg0 θ 可行性扫描。
- `src/bfs/trajectory_search_config.cpp`：管线默认配置工厂。
- `src/bfs/adapters/problem2_angle_frame_adapter.cpp`：角度坐标系适配器实现。

## 测试（`tests/`）

测试按所验证的模块/主题分子目录存放（详见 `tests/README.md`）：

- `tests/CMakeLists.txt`：通过子目录注册全部测试目标。
- `tests/common/`：构建冒烟测试、时间工具、轨道时间积分（`orbit_math`）。
- `tests/config/`：全局默认配置。
- `tests/planet_params/`：行星参数表与指定时刻状态查询。
- `tests/problem1/`：残差评估、直接求解器、端点表格、表格理论、表格诊断。
- `tests/problem2/`：弹弓残差方程。
- `tests/trajectory/`：轨道速度辅助函数。

## 已移除的根表路径

根表生成器、2 度表格加载器、最近节点查询层，以及基于表格的 BFS/Problem2 搜索模块已从活跃代码树和 CMake 构建目标中移除。

当前 Problem 1 工作应使用**直接求解器**或**保留的端点转移时间表格（endpoint transfer-time table）**。
