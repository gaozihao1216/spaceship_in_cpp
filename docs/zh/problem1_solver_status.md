# Problem 1 求解器现状

## 当前求解器边界

活跃实现位于：

- `include/spaceship_cpp/problem1/problem1.hpp`
- `src/problem1/problem1.cpp`

这是 Problem 1 的**直接残差求解器**。

## 当前直接求解流程

1. 在 [0, 2π) 上扫描遇合角（encounter angle）候选点
2. 对每个 (k, q) 多圈分支检测残差符号变化区间
3. 对接受的变号区间用二分法细化根
4. 返回按飞行时间排序的、物理上有效的转移候选解

## 当前角色

- Problem 1 的**正确性基准**（ground truth）
- 测试和诊断程序使用的转移时间候选来源
- 移除根表预计算路径后的**替代方向**

## 与表格模块的关系

保留的 `Problem1Table` 模块是一个较小的**端点几何 + 转移时间表格**实验。

它**不负责**根的预计算（precomputation），仅用于：

- 离散采样 (ν_A, ν_B, θ_A) 三维空间
- 枚举每个几何点下的 (k, q) 飞行时间分支
- 为后续插值/查询实验提供数据基础

## 核心数学问题（一句话）

给定出发行星、目标行星、发射时刻和转移轨道近日点方向，求使**转移轨道飞行时间 = 目标行星轨道飞行时间**的遇合角 φ——即对残差函数求根。

## 关键输入参数

| 参数 | 含义 |
|------|------|
| `departure_planet` / `target_planet` | 出发/目标行星 |
| `launch_time_seconds_since_j2000` | 发射时刻（J2000 秒偏移） |
| `transfer_perihelion_angle` | 转移轨道近日点全局角 |
| `max_transfer_revolution` (k) | 转移轨道允许额外绕日圈数 |
| `max_target_revolution` (q) | 目标行星允许额外绕日圈数 |
| `phi_scan_count` | 遇合角初始扫描点数 |
| `phi_tolerance` | 二分收敛角度容差 |
| `max_candidate_relative_residual` | 候选解相对残差过滤阈值 |

## 诊断程序

- `apps/problem1_solve_diagnostics.cpp`：批量运行直接求解，输出 CSV
- `apps/problem1_table_diagnostics.cpp`：构建端点表格，输出统计和 CSV
