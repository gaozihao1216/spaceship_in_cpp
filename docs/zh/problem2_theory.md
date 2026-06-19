# Problem 2 理论框架

## 目标

Problem 2 应被定义为一个**建立在 Problem 1 转移预言机（oracle）之上的图搜索问题**，然后再编写生产级搜索代码。

换句话说：Problem 1 负责回答「从 A 到 B 能不能转移、花多长时间」；Problem 2 负责在行星网络上找一条满足约束的最优/可行路径。

## 数学对象

### 行星集合

设

`P = {Earth, Mars, ...}`

为候选行星或遭遇体（encounter body）的集合。

### 遭遇状态（Encounter State）

一个遭遇状态至少应包含：

- 当前所在行星
- 当前时刻
- 累计 Δv（速度增量）
- 路径历史
- 分支元数据（k, q 等）

根据精度需求，还可能需要：

- 入射转移几何（incoming transfer geometry）
- 出射转移几何（outgoing transfer geometry）
- 上一条转移轨道参数 (e, p, θ)
- 紧凑的入射速度或分支类别状态

### 转移边（Transfer Edge）

Problem 1 提供从行星 A 到行星 B 的一条转移边。

每条边应包含：

- 出发行星 A
- 目标行星 B
- 出发异常角或出发时刻
- 到达异常角或到达时刻
- 转移飞行时间（TOF）
- 转移多圈数 k
- 目标多圈数 q
- 遭遇角 α（encounter alpha）
- 转移轨道参数 (e, p, θ)
- 残差与有效性诊断信息
- 估计的 Δv 贡献

### Problem 1 求解器作为预言机

Problem 2 **不应在内部重新求解全部轨道几何**。

它应调用 Problem 1 预言机，策略为：

1. 尝试 Route B 缓存的 Hessian 查询
2. 若结果有效且分支完整，使用 Route B 结果
3. 否则回退到 Route A
4. 若 Route A 也无效，则该边不存在

> 注意：Problem 2 不应依赖在线 Hessian Route B（当前代码中 Route B 尚未实现）。

## 抽象搜索图

### 节点（Node）

图节点或标签包含：

- 当前行星编号
- 当前时刻或时刻桶（time bucket）
- 累计代价
- 路径序列
- 分支元数据

### 有向边（Directed Edge）

一条有向边是由 Problem 1 预言机生成的、从当前行星到下一颗行星的可行转移。

### 代价（Cost）

候选代价模型包括：

- Δv（速度增量）
- 时间惩罚
- 任务目标惩罚
- 风险或无效性惩罚

### 可行性过滤（Feasibility）

可行性过滤器可能包括：

- 时间窗口约束
- 最大总飞行时间
- 最大多圈数 (k, q)
- 最大遭遇次数
- 最小/最大单次转移时间
- 每颗行星的遭遇约束

## 搜索算法选择

### BFS（广度优先搜索）

若边代价近似均匀，BFS 可以接受。

### Dijkstra

若代价按 Δv 或时间加权，图是有权图，Dijkstra 是自然基线。

### A*

若有目标行星且存在有意义的下界启发式，应考虑 A*。

### 动态规划 / 标签设定（Label Setting）

若状态包含连续时间，搜索可能需要离散化或标签设定公式。

### 多标签最短路径（Multi-Label Shortest Path）

若同一「行星-时刻」状态可持有多个 Pareto 最优标签，则 Problem 2 不是单标签最短路径问题，而是**多标签最短路径问题**。

这就是为什么 Problem 2 不自动等于「只是 BFS」。状态可能依赖于：

- 时间
- Δv
- 入射几何
- 分支标识
- 路径约束

## 支配剪枝（Dominance Pruning）

第一版支配规则草案：

标签 `L1` 支配 `L2`，当且仅当：

- 同一当前行星
- 相近的当前时刻桶
- `L1.accumulated_delta_v <= L2.accumulated_delta_v`
- `L1.elapsed_time <= L2.elapsed_time`
- `L1.remaining_flexibility` 不劣于 L2
- 且至少一个不等式严格成立

但该规则**仅当入射几何不实质影响下一次转移代价时**才安全。

若入射几何改变下游可行性或 Δv，则不能只按「行星 + 时刻」合并标签——状态必须保留入射分支类别或速度状态信息。

## Problem 1 预言机策略

未来 `query_transfer(A, B, departure_state)` 的策略：

1. 尝试 Route B 缓存 Hessian 查询
2. 若有效且分支完整 → 使用 Route B 结果
3. 否则 → 回退 Route A
4. 若 Route A 无效 → 无边

这使 Problem 2 依赖一个稳定的预言机边界，而不是把 Route B 在线 Hessian 工作嵌入搜索本身。

## 最小实现路线图

### Version 0

- 固定行星列表
- 固定最大搜索深度
- 固定离散发射时刻网格
- 每条边使用 Problem 1 Route A
- 目标仅是找到可行路径

### Version 1

- 加入 Route B 缓存查询
- 回退 Route A
- 加入代价（总时间或近似 Δv）

### Version 2

- 多标签 Dijkstra
- 支配剪枝
- 时间窗口约束

### Version 3

- 更真实的 Δv 模型
- 入射/出射速度匹配
- 引力辅助或遭遇连续性约束

## 开放问题

1. Problem 2 的真正优化目标是什么？
   - 最小化 Δv？
   - 最小化时间？
   - 最大化可达行星数？
   - 仅找任意可行路径？

2. 中间遭遇点是否允许脉冲修正（impulsive correction）？

3. 每次遭遇时速度是否必须连续，还是几何可达性就足够？

4. 是否允许等待（在同一行星轨道上等相位）？

5. 是否允许重复访问同一颗行星？

6. `max_transfer_revolution` 和 `max_target_revolution` 是否有全局约束？

7. 输出是单一最优解，还是 Pareto 前沿？

## 与当前代码的对应关系

| 理论概念 | 当前代码位置 |
|----------|-------------|
| Problem 1 预言机 | `problem1::solve_problem1` / `evaluate_problem1_residual` |
| 转移边数据结构 | `bfs::TrajectorySearchEdge` |
| 搜索状态 | `bfs::TrajectorySearchState` |
| 飞掠物理过滤 | `trajectory::evaluate_flyby_physical_feasibility` |
| 弹弓残差（Problem 2 几何） | `problem2::evaluate_problem2_slingshot_residual` |
| 角度坐标适配 | `bfs::global_periapsis_angle_to_problem2_local` |
