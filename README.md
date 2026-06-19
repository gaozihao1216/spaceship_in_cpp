# spaceship_in_cpp

这是与原 Python 项目 spaceship 同级的独立 C++ 重写工程。

当前模块划分：

- problem1：Problem 1 直接求解、诊断和 endpoint transfer-time table
- problem2：Problem 2 弹弓残差方程
- core_types：核心类型
- planet_params：行星参数和行星状态
- common：共用函数
- trajectory：轨道速度和飞掠物理工具
- bfs：BFS 相关基础状态、角度适配器与搜索逻辑文档（`docs/zh/bfs_search.md`）

更详细的当前文件结构见 `docs/zh/code_structure.md`（英文：`docs/en/code_structure.md`）。完整文档索引见 `docs/README.md`。

构建方式：

```bash
bash scripts/build.sh
bash scripts/run_tests.sh
```
