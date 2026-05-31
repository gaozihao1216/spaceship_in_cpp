#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

int main() {
    std::ifstream input("../spaceship/doc/code_theory_framework.md");
    if (!input) {
        // 中文注释：不同测试运行目录下相对路径可能不同，这里允许做有限的本地回退查找。
        input.clear();
        input.open("spaceship/doc/code_theory_framework.md");
    }
    if (!input) {
        input.clear();
        input.open("../../spaceship/doc/code_theory_framework.md");
    }
    assert(input);

    std::ostringstream buffer;
    buffer << input.rdbuf();
    const std::string text = buffer.str();

    // 中文注释：文档必须明确说明 Problem1Table 不是纯 3D 连续表，而是 3 个连续角轴加 2 个离散 (k,q) 分支轴。
    assert(text.find("Table_{A \\to B}(\\nu_A,\\nu_B,\\theta_A; k, q)") != std::string::npos);
    assert(text.find("three continuous angle axes and two discrete transfer-revolution branch axes") !=
        std::string::npos);
    assert(text.find("launch-phase-specific") != std::string::npos);
    assert(text.find("representative only") != std::string::npos);
    assert(text.find("q does not participate in interpolation admissibility") != std::string::npos);
    assert(text.find("q is evaluated online with real `t_depart` or `nu_B_depart`") != std::string::npos);
    assert(text.find("geometry + `k`") != std::string::npos);
    assert(text.find("admissibility checks transfer-side k validity only") != std::string::npos);
    assert(text.find("representative q branch validity must not affect interpolation admissibility") !=
        std::string::npos);
    assert(text.find("不同 `q`") != std::string::npos);
    assert(text.find("pure transfer branch view") != std::string::npos);
    assert(text.find("`1x1x1`") != std::string::npos);
    assert(text.find("Future Interpolator Stub") != std::string::npos);
    assert(text.find("transfer time(`k`)") != std::string::npos);
    assert(text.find("`q` 在线计算") != std::string::npos);
    assert(text.find("sorted solution list") != std::string::npos);
    assert(text.find("`k,q` 不是表格维度") != std::string::npos);
    assert(text.find("experimental / transitional helper") != std::string::npos);
    assert(text.find("Problem1RootTable Draft") != std::string::npos);
    assert(text.find("solutions_sorted_by_time_of_flight") != std::string::npos);
    assert(text.find("solve_problem1_from_departure_anomalies") != std::string::npos);
    assert(text.find("time_of_flight_seconds") != std::string::npos);
    assert(text.find("fake launch-time inversion") != std::string::npos);
    assert(text.find("encounter_global_angle") != std::string::npos);
    assert(text.find("target_arrival_true_anomaly") != std::string::npos);
    assert(text.find("First Root-Table Query Layer") != std::string::npos);
    assert(text.find("Route A：nearest-node linear Taylor seed + Newton refinement") != std::string::npos);
    assert(text.find("Route B：nearest-node quadratic raw approximation") != std::string::npos);
    assert(text.find("nearest-node baseline") != std::string::npos);
    assert(text.find("TangentFiniteDifference") != std::string::npos);
    assert(text.find("NewtonRefinedFiniteDifference") != std::string::npos);
    assert(text.find("Tangent finite-difference Hessian") != std::string::npos);
    assert(text.find("alpha_{tangent}") != std::string::npos);
    assert(text.find("admissibility gate") != std::string::npos);
    assert(text.find("admissible_for_fast_approximation") != std::string::npos);
    assert(text.find("fallback 到 Route A") != std::string::npos);
    assert(text.find("faster raw approximation candidate") != std::string::npos);
    assert(text.find("对照 / 验证工具") != std::string::npos);
    assert(text.find("8 节点 Hessian fitting") != std::string::npos);
    assert(text.find("完整 branch-safe trilinear interpolation") != std::string::npos);
    assert(text.find("query_problem1_root_table_fast_with_fallback") != std::string::npos);
    assert(text.find("Problem1RootTable Resolution Sweep") != std::string::npos);
    assert(text.find("8 x 8 x 12") != std::string::npos || text.find("8x8x12") != std::string::npos);
    assert(text.find("smoke-test / functionality-test") != std::string::npos);
    assert(text.find("RUN_EXPENSIVE_ROOT_TABLE_SWEEP=1") != std::string::npos);
    assert(text.find("multi-node seed pool") != std::string::npos);
    assert(text.find("Virtual 2-Degree Nearest-Node Experiment") != std::string::npos);
    assert(text.find("180^3 = 5,832,000") != std::string::npos);
    assert(text.find("avg node solve time * 180^3") != std::string::npos);
    assert(text.find("nearest-node feasibility probe") != std::string::npos);
    assert(text.find("same `(k,q)` seed") != std::string::npos);
    assert(text.find("Newton basin") != std::string::npos);
    assert(text.find("audit Newton refinement") != std::string::npos);
    assert(text.find("profile `solve_problem1_from_departure_anomalies(...)`") != std::string::npos);
    assert(text.find("Residual Derivative Fallback and Newton Audit") != std::string::npos);
    assert(text.find("negative-e normalization branch") != std::string::npos);
    assert(text.find("finite-difference residual derivative fallback") != std::string::npos);
    assert(text.find("Analytic Derivative Canonicalization for Negative E") != std::string::npos);
    assert(text.find("e = |e_raw|") != std::string::npos);
    assert(text.find("theta_B_star") != std::string::npos);

    std::ifstream problem1_status("../spaceship/doc/problem1_solver_status.md");
    if (!problem1_status) {
        problem1_status.clear();
        problem1_status.open("spaceship/doc/problem1_solver_status.md");
    }
    if (!problem1_status) {
        problem1_status.clear();
        problem1_status.open("../../spaceship/doc/problem1_solver_status.md");
    }
    assert(problem1_status);
    std::ostringstream problem1_status_buffer;
    problem1_status_buffer << problem1_status.rdbuf();
    const std::string problem1_status_text = problem1_status_buffer.str();
    assert(problem1_status_text.find("Route A") != std::string::npos);
    assert(problem1_status_text.find("Route B") != std::string::npos);
    assert(problem1_status_text.find("ProjectedTangentFiniteDifference") != std::string::npos);
    assert(problem1_status_text.find("precomputed Hessian") != std::string::npos);
    assert(problem1_status_text.find("Route A fallback") != std::string::npos);
    assert(problem1_status_text.find("residual_tolerance_seconds = 1e-2") != std::string::npos);
    assert(problem1_status_text.find("online Hessian is not the target path") != std::string::npos);

    std::ifstream problem2_theory("../spaceship/doc/problem2_theory.md");
    if (!problem2_theory) {
        problem2_theory.clear();
        problem2_theory.open("spaceship/doc/problem2_theory.md");
    }
    if (!problem2_theory) {
        problem2_theory.clear();
        problem2_theory.open("../../spaceship/doc/problem2_theory.md");
    }
    assert(problem2_theory);
    std::ostringstream problem2_theory_buffer;
    problem2_theory_buffer << problem2_theory.rdbuf();
    const std::string problem2_text = problem2_theory_buffer.str();
    assert(problem2_text.find("Problem 1 oracle") != std::string::npos);
    assert(problem2_text.find("BFS") != std::string::npos);
    assert(problem2_text.find("Dijkstra") != std::string::npos);
    assert(problem2_text.find("A*") != std::string::npos);
    assert(problem2_text.find("multi-label shortest path") != std::string::npos);
    assert(problem2_text.find("dominance") != std::string::npos);
    assert(problem2_text.find("Route B cached Hessian query") != std::string::npos);
    assert(problem2_text.find("fallback Route A") != std::string::npos);
    assert(problem2_text.find("Open Questions") != std::string::npos);

    return 0;
}
