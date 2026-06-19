# 代码结构

本文档描述在移除根表（root-table）预计算路径之后，当前源码的目录布局。

## 顶层目录

- `CMakeLists.txt`：构建 `spaceship_cpp` 静态库、两个诊断程序（apps）和测试套件。
- `README.md`：项目简介与构建命令。
- `apps/problem1_solve_diagnostics.cpp`：Problem 1 直接求解器的命令行诊断程序。
- `apps/problem1_table_diagnostics.cpp`：Problem 1 端点表格（endpoint table）的命令行诊断程序。
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
- `include/spaceship_cpp/bfs/bfs.hpp`：BFS 模块占位入口头文件。
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
- `src/bfs/bfs.cpp`：BFS 模块占位实现。
- `src/bfs/adapters/problem2_angle_frame_adapter.cpp`：角度坐标系适配器实现。

## 测试（`tests/`）

- `tests/CMakeLists.txt`：注册与当前非根表代码路径匹配的测试目标。
- `tests/test_build.cpp`：基础构建/链接冒烟测试。
- `tests/test_global_config.cpp`：全局配置检查。
- `tests/test_orbit_time_integral.cpp`：轨道时间积分检查。
- `tests/test_orbit_velocity_helpers.cpp`：轨道速度辅助函数检查。
- `tests/test_planet_params.cpp`：行星参数检查。
- `tests/test_planet_position.cpp`：行星指定时刻状态检查。
- `tests/test_problem1_residual.cpp`：Problem 1 残差函数检查。
- `tests/test_problem1_solve.cpp`：Problem 1 直接求解器检查。
- `tests/test_problem1_table.cpp`：端点表格检查。
- `tests/test_problem1_table_diagnostics.cpp`：Problem 1 表格诊断检查。
- `tests/test_problem1_table_theory.cpp`：Problem 1 表格理论性质检查。
- `tests/test_problem2_slingshot_equations.cpp`：Problem 2 弹弓方程检查。
- `tests/test_time_utils.cpp`：时间工具检查。

## 已移除的根表路径

根表生成器、2 度表格加载器、最近节点查询层，以及基于表格的 BFS/Problem2 搜索模块已从活跃代码树和 CMake 构建目标中移除。

当前 Problem 1 工作应使用**直接求解器**或**保留的端点转移时间表格（endpoint transfer-time table）**。
