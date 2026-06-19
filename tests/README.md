# 测试目录说明

测试按**验证目标**分子目录，与 `src/` 模块划分对应。

```
tests/
├── CMakeLists.txt          # 公共 add_spaceship_test 宏 + 子目录入口
├── common/                 # 通用工具与构建冒烟
├── config/                 # 全局配置
├── planet_params/          # 行星参数与状态
├── problem1/               # Problem 1 求解、表格与诊断
├── problem2/               # Problem 2 弹弓方程
└── trajectory/             # 轨道速度与飞掠相关工具
```

## 各目录测试目的

| 目录 | 测试目标 | 文件 |
|------|----------|------|
| `common/` | 链接冒烟、时间工具、轨道时间积分 F(e,ξ) | `test_build`, `test_time_utils`, `test_orbit_time_integral` |
| `config/` | 全局默认配置与工厂函数 | `test_global_config` |
| `planet_params/` | 行星常数表、按时刻求位置/真近点角 | `test_planet_params`, `test_planet_position` |
| `problem1/` | 残差、直接求解器、端点表格、理论约束、诊断导出 | `test_problem1_*`（5 个） |
| `problem2/` | 飞掠弹弓不变量与残差方程 | `test_problem2_slingshot_equations` |
| `trajectory/` | 日心轨道速度与相对行星速度 | `test_orbit_velocity_helpers` |

运行全部测试：

```bash
bash scripts/run_tests.sh
```
