# Earth → Mercury Trajectory Search

This project’s primary mission: **given an Earth launch time, find the best transfer path to Mercury** (direct or via inner-planet flybys) and return orbit parameters for visualization.

## Public API

- **Header**: `include/spaceship_cpp/bfs/bfs.hpp`
- **Implementation**: `src/bfs/bfs.cpp`
- **Function**: `search_best_trajectory(TrajectorySearchInput, TrajectorySearchGlobalConfig)`

Input: launch time, departure planet (Earth), destination planet (Mercury).  
Output: `TrajectorySearchOutput` with `visit_sequence`, per-leg `(e, p, theta)`, times, and score.

See the Chinese doc for full detail: [`docs/zh/earth_to_mercury.md`](../zh/earth_to_mercury.md).

## Visualization (stub)

Directory [`visualization/`](../../visualization/README.md) — placeholder only; plotting not implemented yet.

## Run

```bash
cmake --build build -j4 --target trajectory_search
./build/trajectory_search
```
