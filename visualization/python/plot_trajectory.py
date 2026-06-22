#!/usr/bin/env python3
"""
轨迹可视化占位脚本。

用法（待实现）:
    python plot_trajectory.py trajectory.json

输入 JSON 应由 TrajectorySearchOutput 序列化得到，包含 visit_sequence 与各 leg 的 (e, p, theta)。
"""

import sys


def main() -> int:
    print("visualization not implemented yet")
    print("see visualization/README.md and docs/zh/earth_to_mercury.md")
    if len(sys.argv) > 1:
        print(f"would plot: {sys.argv[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
