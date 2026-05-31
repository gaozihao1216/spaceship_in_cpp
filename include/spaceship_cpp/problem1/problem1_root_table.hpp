#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

struct Problem1RootTableCell;

struct Problem1SolutionBranch {
    bool valid = false;
    // 中文注释：当前 root-solution-table 与 solve_problem1 保持一致，主角变量保存的是 encounter 全局角。
    double encounter_global_angle = 0.0;
    // 中文注释：目标行星到达 true anomaly 由 encounter_global_angle 减去目标轨道 theta_0 得到。
    double target_arrival_true_anomaly = 0.0;
    int transfer_revolution = 0;
    int target_revolution = 0;
    // 中文注释：root table 的核心时间语义是“自 departure 起经过的飞行时间”，而不是绝对到达时刻。
    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double residual_seconds = 0.0;
    // 中文注释：只有从带 absolute launch_time 的 solve_problem1(...) candidate 转换时，
    // 这个绝对 J2000 到达时刻才有意义；纯 root-table cell 不应依赖它。
    double arrival_time_seconds_since_j2000 = 0.0;

    double transfer_e = 0.0;
    double transfer_p = 0.0;
    double transfer_a = 0.0;
    double theta_B = 0.0;
    bool derivatives_available = false;
    double d_encounter_global_angle_d_nu_A = 0.0;
    double d_encounter_global_angle_d_nu_B = 0.0;
    double d_encounter_global_angle_d_theta_A = 0.0;
    double root_bracket_width = std::numeric_limits<double>::quiet_NaN();
    int bisection_iterations = -1;
    std::string adaptive_source_reason;
    double adaptive_ternary_phi_star = std::numeric_limits<double>::quiet_NaN();
    double adaptive_ternary_residual = std::numeric_limits<double>::quiet_NaN();
    double adaptive_ternary_residual_sq = std::numeric_limits<double>::quiet_NaN();
    double adaptive_local_interval_left_phi = std::numeric_limits<double>::quiet_NaN();
    double adaptive_local_interval_right_phi = std::numeric_limits<double>::quiet_NaN();
    int adaptive_local_subdivision_index = -1;
    std::string invalid_reason;
};

struct Problem1RootResidualResult {
    Problem1ResidualStatus status = Problem1ResidualStatus::InvalidInput;
    bool valid = false;
    double residual_scale_free = 0.0;
    double residual_seconds = 0.0;
    double transfer_time_scale_free = 0.0;
    double target_time_scale_free = 0.0;
    double transfer_time_seconds = 0.0;
    double target_time_seconds = 0.0;
    double encounter_global_angle = 0.0;
    double target_arrival_true_anomaly = 0.0;
    double transfer_e_raw = 0.0;
    double transfer_e = 0.0;
    double transfer_p = 0.0;
    double transfer_a = 0.0;
    double theta_B = 0.0;
    std::string invalid_reason;
};

struct Problem1RootResidualDerivatives {
    bool valid = false;
    double R_alpha = 0.0;
    double R_nu_A = 0.0;
    double R_nu_B = 0.0;
    double R_theta_A = 0.0;
    double d_alpha_d_nu_A = 0.0;
    double d_alpha_d_nu_B = 0.0;
    double d_alpha_d_theta_A = 0.0;
    double F_alpha = 0.0;
    std::string invalid_reason;
};

enum class Problem1RootDerivativeMode {
    AnalyticOnly,
    FiniteDifferenceOnly,
    AnalyticWithFiniteDifferenceFallback,
};

struct Problem1RootLinearPrediction {
    bool valid = false;
    std::string method;
    std::string invalid_reason;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double predicted_encounter_global_angle = 0.0;
};

struct Problem1RootQSheetSelectionResult {
    int selected_q = 0;
    bool selection_failed = false;
    bool q_changed = false;
    double selected_continuity_error = std::numeric_limits<double>::quiet_NaN();
    double source_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
};

struct Problem1RootNewtonTraceStep {
    int iteration = 0;
    double alpha = 0.0;
    double residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double residual_scale_free = std::numeric_limits<double>::quiet_NaN();
    double R_alpha = std::numeric_limits<double>::quiet_NaN();
    double delta_alpha = std::numeric_limits<double>::quiet_NaN();
    bool derivative_valid = false;
    std::string derivative_source;
    bool alpha_normalized = false;
    bool step_clamped = false;
    bool residual_increased = false;
    std::string reason;
};

struct Problem1RootNewtonDiagnostic {
    bool valid = false;
    bool converged = false;
    bool derivative_failed = false;
    bool derivative_attach_failed_after_convergence = false;
    bool step_clamped = false;
    bool residual_increased = false;
    bool likely_wrong_root = false;
    int fallback_used_count = 0;
    int finite_difference_success_count = 0;
    int iterations = 0;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double initial_alpha = 0.0;
    double final_alpha = 0.0;
    double initial_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double final_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    std::vector<Problem1RootNewtonTraceStep> trace;
    std::string invalid_reason;
};

struct Problem1RootRefinementResult {
    bool valid = false;
    Problem1SolutionBranch branch;
    Problem1RootNewtonDiagnostic diagnostic;
};

struct Problem1RootNearestNode {
    bool valid = false;
    int i = 0;
    int j = 0;
    int k = 0;
    double node_nu_A = 0.0;
    double node_nu_B = 0.0;
    double node_theta_A = 0.0;
    const Problem1RootTableCell* cell = nullptr;
    std::string invalid_reason;
};

struct Problem1SolveDiagnostic {
    bool valid = false;
    int enumerated_branch_pairs = 0;
    long long alpha_scan_samples = 0;
    long long residual_evaluations = 0;
    long long bisection_refinements = 0;
    double residual_evaluation_seconds = 0.0;
    double root_scanning_seconds = 0.0;
    double sorting_conversion_seconds = 0.0;
    int final_branch_count = 0;
    std::string invalid_reason;
};

struct Problem1SolveWithDiagnosticResult {
    std::vector<Problem1SolutionBranch> branches;
    Problem1SolveDiagnostic diagnostic;
};

struct Problem1BranchCompareOptions {
    double angle_tolerance = 1e-8;
    double target_anomaly_tolerance = 1e-8;
    double time_tolerance_seconds = 1e-4;
    double relative_time_tolerance = 1e-10;
    int max_detail_count = 5;
};

struct Problem1BranchCompareDetail {
    double nu_A_depart = std::numeric_limits<double>::quiet_NaN();
    double nu_B_depart = std::numeric_limits<double>::quiet_NaN();
    double theta_A = std::numeric_limits<double>::quiet_NaN();
    int transfer_revolution = 0;
    int target_revolution = 0;
    double encounter_global_angle = std::numeric_limits<double>::quiet_NaN();
    double target_arrival_true_anomaly = std::numeric_limits<double>::quiet_NaN();
    double time_of_flight_seconds = std::numeric_limits<double>::quiet_NaN();
    double target_time_seconds = std::numeric_limits<double>::quiet_NaN();
    double residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double relative_residual = std::numeric_limits<double>::quiet_NaN();
};

struct Problem1BranchCompareReport {
    int matched_count = 0;
    int missing_count = 0;
    int extra_count = 0;
    double max_angle_error = 0.0;
    double max_target_anomaly_error = 0.0;
    double max_time_error = 0.0;
    double max_residual_error = 0.0;
    std::vector<Problem1BranchCompareDetail> missing_branches;
    std::vector<Problem1BranchCompareDetail> extra_branches;
};

struct Problem1AdaptiveExtraBranchDiagnostic {
    Problem1BranchCompareDetail branch;
    double root_bracket_width = std::numeric_limits<double>::quiet_NaN();
    int bisection_iterations = -1;
    std::string adaptive_source_reason;
    double ternary_phi_star = std::numeric_limits<double>::quiet_NaN();
    double ternary_residual = std::numeric_limits<double>::quiet_NaN();
    double ternary_residual_sq = std::numeric_limits<double>::quiet_NaN();
    double local_interval_left_phi = std::numeric_limits<double>::quiet_NaN();
    double local_interval_right_phi = std::numeric_limits<double>::quiet_NaN();
    int local_subdivision_index = -1;
    double sweep_min_abs_residual = std::numeric_limits<double>::quiet_NaN();
    double sweep_phi_at_min_abs_residual = std::numeric_limits<double>::quiet_NaN();
    bool sweep_sign_change_exists = false;
    int sweep_valid_sample_count = 0;
    int sweep_invalid_sample_count = 0;
    double sweep_residual_at_adaptive_phi = std::numeric_limits<double>::quiet_NaN();
    double sweep_left_residual = std::numeric_limits<double>::quiet_NaN();
    double sweep_right_residual = std::numeric_limits<double>::quiet_NaN();
    int nearest_baseline_transfer_revolution = -1;
    int nearest_baseline_target_revolution = -1;
    double nearest_baseline_angle_diff = std::numeric_limits<double>::quiet_NaN();
    double nearest_baseline_target_anomaly_diff = std::numeric_limits<double>::quiet_NaN();
    double nearest_baseline_time_diff = std::numeric_limits<double>::quiet_NaN();
    double nearest_baseline_residual_diff = std::numeric_limits<double>::quiet_NaN();
    std::string nearest_baseline_match_failure_reason;
    bool high_res_oracle_enabled = false;
    int high_res_scan_count = 0;
    bool high_res_oracle_matched = false;
    int high_res_oracle_branch_count = 0;
};

enum class Problem1SolveMode {
    FullScan2880,
    AdaptiveScanWithFallback,
    AdaptiveScanDebugCompare,
};

struct Problem1AdaptiveScanOptions {
    int coarse_sample_count = 360;
    int fine_subdivision_count = 32;
    int interval_expand_steps = 1;
    int max_candidate_interval_count = 64;
    int max_residual_evaluation_budget = 1200;
};

struct Problem1SolveOptions {
    Problem1SolveMode mode = Problem1SolveMode::FullScan2880;
    int max_transfer_revolution = 1;
    int max_target_revolution = 1;
    Problem1AdaptiveScanOptions adaptive;
    bool fallback_to_full_scan = true;
};

struct Problem1SolveModeProfile {
    Problem1SolveMode mode_requested = Problem1SolveMode::FullScan2880;
    Problem1SolveMode mode_used = Problem1SolveMode::FullScan2880;
    bool adaptive_attempted = false;
    bool fallback_used = false;
    bool debug_compare_enabled = false;
    int debug_compare_matched_count = 0;
    int debug_compare_missing_from_adaptive = 0;
    int debug_compare_new_in_adaptive = 0;
    double debug_compare_max_alpha_diff = 0.0;
    double debug_compare_max_theta_prime_diff = 0.0;
    double debug_compare_max_target_time_diff = 0.0;
    double debug_compare_max_residual_diff = 0.0;
    int residual_evaluation_count = 0;
    int branch_count = 0;
    double total_ms = 0.0;
    long long adaptive_coarse_scan_count = 0;
    long long adaptive_interval_total = 0;
    long long adaptive_sign_change_interval_count = 0;
    long long adaptive_near_zero_interval_count = 0;
    long long adaptive_local_min_interval_count = 0;
    long long adaptive_local_max_interval_count = 0;
    long long adaptive_valid_boundary_interval_count = 0;
    long long adaptive_wrap_interval_count = 0;
    long long adaptive_rapid_change_interval_count = 0;
    long long adaptive_residual_eval_count = 0;
    long long adaptive_refined_interval_count = 0;
    long long adaptive_bisection_interval_count = 0;
    long long adaptive_ternary_interval_count = 0;
    long long adaptive_local_fine_scan_interval_count = 0;
    long long adaptive_candidate_count_before_dedup = 0;
    long long adaptive_candidate_count_after_dedup = 0;
    long long adaptive_ternary_accept_count = 0;
    long long adaptive_ternary_reject_count = 0;
    long long adaptive_local_fine_scan_root_count = 0;
    bool adaptive_fallback_to_fullscan = false;
    double adaptive_coarse_scan_seconds = 0.0;
    double adaptive_interval_collection_seconds = 0.0;
    double adaptive_interval_refine_seconds = 0.0;
    double adaptive_local_fine_scan_seconds = 0.0;
    double adaptive_ternary_seconds = 0.0;
    double adaptive_bisection_seconds = 0.0;
    double adaptive_candidate_dedup_seconds = 0.0;
    double adaptive_sorting_seconds = 0.0;
    std::vector<Problem1BranchCompareDetail> debug_compare_missing_branches;
    std::vector<Problem1BranchCompareDetail> debug_compare_extra_branches;
    std::vector<Problem1AdaptiveExtraBranchDiagnostic> debug_compare_extra_diagnostics;
};

struct Problem1SolveWithModeResult {
    std::vector<Problem1SolutionBranch> branches;
    Problem1SolveDiagnostic diagnostic;
    Problem1SolveModeProfile mode_profile;
};

// 中文注释：当前 residual gate 仍是 Route B 的实验参数；Projected 模式每个 stencil 只做一次 projection，
// 如果 query 端 residual_after_correction_too_large，仍应 fallback 到 Route A。
inline constexpr double kProblem1RootExperimentalTangentResidualTolerance = 1.0;

struct Problem1RootApproximationDiagnostics {
    std::string hessian_method;
    bool hessian_valid = false;
    double tangent_residual_max_scale_free = std::numeric_limits<double>::quiet_NaN();
    double tangent_residual_max_seconds = std::numeric_limits<double>::quiet_NaN();
    double hessian_step = std::numeric_limits<double>::quiet_NaN();
    double raw_residual_scale_free = std::numeric_limits<double>::quiet_NaN();
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double alpha_linear = std::numeric_limits<double>::quiet_NaN();
    double alpha_quadratic = std::numeric_limits<double>::quiet_NaN();
    int source_target_revolution = 0;
    int selected_target_revolution = 0;
    bool q_sheet_selection_changed = false;
    bool q_sheet_selection_failed = false;
    double selected_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
    double source_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
    bool admissible_for_fast_approximation = false;
    std::string admissibility_reason;
};

struct Problem1RootApproximationResult {
    bool valid = false;
    std::string method;
    std::string invalid_reason;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double predicted_encounter_global_angle = 0.0;
    double target_arrival_true_anomaly = 0.0;
    double residual_scale_free = 0.0;
    double residual_seconds = 0.0;
    double transfer_time_seconds = 0.0;
    double target_time_seconds = 0.0;
    double transfer_e = 0.0;
    double transfer_p = 0.0;
    double transfer_a = 0.0;
    double theta_B = 0.0;
    Problem1RootApproximationDiagnostics diagnostics;
};

enum class Problem1RootHessianMethod {
    // 中文注释：legacy/diagnostic/ablation mode，保留旧的 tangent residual pre-gate 行为。
    TangentFiniteDifference,
    // 中文注释：Route B 默认 Hessian mode；不调用 Newton，只做一次 stencil-level residual projection。
    ProjectedTangentFiniteDifference,
    NewtonRefinedFiniteDifference,
};

struct Problem1RootHessian {
    bool valid = false;
    double H_nu_A_nu_A = 0.0;
    double H_nu_B_nu_B = 0.0;
    double H_theta_A_theta_A = 0.0;
    double H_nu_A_nu_B = 0.0;
    double H_nu_A_theta_A = 0.0;
    double H_nu_B_theta_A = 0.0;
    double tangent_residual_max_scale_free = std::numeric_limits<double>::quiet_NaN();
    double tangent_residual_max_seconds = std::numeric_limits<double>::quiet_NaN();
    std::string method;
    std::string invalid_reason;
};

struct Problem1RootQuadraticPrediction {
    bool valid = false;
    std::string method;
    std::string invalid_reason;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double predicted_encounter_global_angle = 0.0;
    std::string hessian_method;
    bool hessian_valid = false;
    double tangent_residual_max_scale_free = std::numeric_limits<double>::quiet_NaN();
    double tangent_residual_max_seconds = std::numeric_limits<double>::quiet_NaN();
    double hessian_step = std::numeric_limits<double>::quiet_NaN();
};

struct Problem1RootTableQueryResult {
    bool valid = false;
    std::string method;
    std::string invalid_reason;
    std::vector<Problem1SolutionBranch> branches;
};

struct Problem1RootCellAdmissibilityResult {
    bool admissible = false;
    std::string reason;
    int corner_count = 0;
    std::map<int, int> reference_root_count_by_k;
    std::map<int, int> min_root_count_by_k;
    std::map<int, int> max_root_count_by_k;
};

struct Problem1RouteBSafeQueryResult {
    bool valid = false;
    bool fallback_required = false;
    bool branch_count_complete = false;
    std::string method;
    std::string reason;
    Problem1RootCellAdmissibilityResult cell_admissibility;
    std::map<int, int> expected_count_by_k;
    std::map<int, int> candidate_count_by_k;
    std::map<int, int> missing_count_by_k;
    std::map<int, int> extra_count_by_k;
    std::map<int, bool> incomplete_by_k;
    std::vector<Problem1RootApproximationResult> approximations;
};

struct Problem1RootTableCell {
    double nu_A_depart = 0.0;
    double nu_B_depart = 0.0;
    double theta_A = 0.0;
    bool solved = false;
    // 中文注释：root table 的自然排序键是 time_of_flight_seconds。
    std::vector<Problem1SolutionBranch> solutions_sorted_by_time_of_flight;
    std::string invalid_reason;
};

struct Problem1RootTableConfig {
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;

    double nu_A_start;
    double nu_A_step;
    int nu_A_count;

    double nu_B_depart_start;
    double nu_B_depart_step;
    int nu_B_depart_count;

    double theta_A_start;
    double theta_A_step;
    int theta_A_count;

    int max_transfer_revolution;
    int max_target_revolution;

    std::string schema_version = "problem1_root_table_draft_v0";
};

class Problem1RootTable {
public:
    explicit Problem1RootTable(Problem1RootTableConfig config);

    const Problem1RootTableConfig& config() const;
    std::size_t flat_index(int nu_A_index, int nu_B_depart_index, int theta_A_index) const;
    const Problem1RootTableCell& at(int nu_A_index, int nu_B_depart_index, int theta_A_index) const;
    const std::vector<Problem1RootTableCell>& cells() const;
    std::vector<Problem1RootTableCell>& mutable_cells();

private:
    Problem1RootTableConfig config_;
    std::vector<Problem1RootTableCell> cells_;
};

Problem1SolutionBranch problem1_solution_branch_from_candidate(
    planet_params::PlanetId target_planet,
    const Problem1Candidate& candidate
);

std::vector<Problem1SolutionBranch> convert_problem1_candidates_to_solution_branches(
    planet_params::PlanetId target_planet,
    const std::vector<Problem1Candidate>& candidates
);

// 中文注释：未来完整 root-solution-table 需要这个底层接口，直接以 departure anomalies 为输入；
// 它不应依赖 absolute launch_time，也不应通过 synodic period / relative phase 反推假时间。
std::vector<Problem1SolutionBranch> solve_problem1_from_departure_anomalies(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int max_transfer_revolution,
    int max_target_revolution
);

Problem1SolveWithDiagnosticResult solve_problem1_from_departure_anomalies_diagnostic(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int max_transfer_revolution,
    int max_target_revolution
);

Problem1SolveWithModeResult solve_problem1_from_departure_anomalies_with_mode(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1SolveOptions& options
);

Problem1BranchCompareReport compare_problem1_solution_branches_by_identity(
    const std::vector<Problem1SolutionBranch>& candidate_branches,
    const std::vector<Problem1SolutionBranch>& baseline_branches,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1BranchCompareOptions& options = Problem1BranchCompareOptions{}
);

Problem1RootResidualResult evaluate_problem1_root_residual(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
);

double problem1_residual_seconds_to_scale_free(double residual_seconds);

// 中文注释：q 是 target-time sheet label，不是 branch identity。这里通过 target_time_seconds
// 对 source branch 的连续性来选择 q，而不是按 residual 最小值选 q。
Problem1RootQSheetSelectionResult select_q_by_target_time_sheet_continuity(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int transfer_revolution,
    double alpha_linear,
    const Problem1SolutionBranch& source_branch,
    int max_target_revolution
);

Problem1RootResidualDerivatives evaluate_problem1_root_residual_derivatives(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
);

Problem1RootResidualDerivatives estimate_problem1_root_residual_derivatives_finite_difference(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution,
    double step
);

Problem1RootResidualDerivatives evaluate_problem1_root_residual_derivatives_with_mode(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution,
    Problem1RootDerivativeMode mode,
    double finite_difference_step
);

Problem1SolutionBranch attach_problem1_root_derivatives(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1SolutionBranch& branch
);

Problem1SolutionBranch attach_problem1_root_derivatives_with_mode(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1SolutionBranch& branch,
    Problem1RootDerivativeMode mode,
    double finite_difference_step
);

Problem1SolutionBranch refine_problem1_root_branch_newton(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    int max_iterations,
    // 中文注释：residual_tolerance 的单位是 scale-free residual；如果调用侧想按 seconds 指定阈值，
    // 必须先用 problem1_residual_seconds_to_scale_free(...) 转换。alpha_step_tolerance 的单位是 radians。
    double residual_tolerance,
    double alpha_step_tolerance
);

// 中文注释：推荐外部调用侧使用 seconds 版本，避免把秒制残差容差误传给 scale-free Newton 接口。
Problem1SolutionBranch refine_problem1_root_branch_newton_seconds(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance
);

Problem1RootRefinementResult refine_problem1_root_branch_newton_diagnostic(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    int max_iterations,
    // 中文注释：residual_tolerance 的单位是 scale-free residual；如果调用侧想按 seconds 指定阈值，
    // 必须先用 problem1_residual_seconds_to_scale_free(...) 转换。alpha_step_tolerance 的单位是 radians。
    double residual_tolerance,
    double alpha_step_tolerance,
    double max_alpha_step = std::numeric_limits<double>::infinity(),
    bool enable_backtracking = false,
    int max_backtracking_steps = 0,
    Problem1RootDerivativeMode derivative_mode = Problem1RootDerivativeMode::AnalyticOnly,
    double finite_difference_step = 1e-6
);

// 中文注释：推荐外部调用侧使用 seconds 版本，避免把秒制残差容差误传给 scale-free Newton 接口。
Problem1RootRefinementResult refine_problem1_root_branch_newton_diagnostic_seconds(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance,
    double max_alpha_step = std::numeric_limits<double>::infinity(),
    bool enable_backtracking = false,
    int max_backtracking_steps = 0,
    Problem1RootDerivativeMode derivative_mode = Problem1RootDerivativeMode::AnalyticOnly,
    double finite_difference_step = 1e-6
);

// 中文注释：residual-first Newton wrapper。每轮先检查 public residual_seconds 是否达标；
// 达标则立即成功返回。只有残差未达标时才计算导数和 Newton step。
Problem1RootRefinementResult refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance,
    Problem1RootDerivativeMode derivative_mode = Problem1RootDerivativeMode::AnalyticOnly,
    double finite_difference_step = 1e-6
);

Problem1SolutionBranch refine_problem1_root_branch_newton_residual_first_seconds(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance,
    Problem1RootDerivativeMode derivative_mode = Problem1RootDerivativeMode::AnalyticOnly,
    double finite_difference_step = 1e-6
);

Problem1RootLinearPrediction predict_problem1_root_branch_linear_from_node(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
);

Problem1SolutionBranch refine_problem1_root_branch_from_linear_prediction(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_iterations,
    double residual_tolerance,
    double alpha_step_tolerance
);

// 中文注释：推荐外部调用侧使用 seconds 版本，避免把秒制残差容差误传给 scale-free Newton 接口。
Problem1SolutionBranch refine_problem1_root_branch_from_linear_prediction_seconds(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance
);

Problem1RootNearestNode find_nearest_problem1_root_table_node(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
);

Problem1RootCellAdmissibilityResult evaluate_problem1_root_cell_admissibility(
    const Problem1RootTable& table,
    int nu_A_index,
    int nu_B_depart_index,
    int theta_A_index,
    int max_transfer_revolution,
    int max_target_revolution
);

Problem1RootApproximationResult evaluate_problem1_root_linear_approximation_from_node(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
);

Problem1RootApproximationResult evaluate_problem1_root_linear_route_b_approximation_from_node(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_target_revolution
);

Problem1RootHessian estimate_problem1_root_hessian_finite_difference(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double step
);

Problem1RootHessian estimate_problem1_root_hessian_tangent_finite_difference(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double step,
    double tangent_residual_tolerance
);

Problem1RootHessian estimate_problem1_root_hessian_projected_tangent_finite_difference(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double step,
    double residual_tolerance_seconds,
    double projection_derivative_min_abs
);

Problem1RootQuadraticPrediction predict_problem1_root_branch_quadratic_from_node(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method = Problem1RootHessianMethod::ProjectedTangentFiniteDifference,
    double tangent_residual_tolerance = kProblem1RootExperimentalTangentResidualTolerance
);

Problem1RootApproximationResult evaluate_problem1_root_quadratic_approximation_from_node(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method = Problem1RootHessianMethod::ProjectedTangentFiniteDifference,
    double tangent_residual_tolerance = kProblem1RootExperimentalTangentResidualTolerance,
    int max_target_revolution = -1
);

Problem1RootApproximationResult evaluate_problem1_root_quadratic_route_b_approximation_from_node(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_target_revolution,
    double hessian_step,
    Problem1RootHessianMethod hessian_method = Problem1RootHessianMethod::ProjectedTangentFiniteDifference,
    double tangent_residual_tolerance = kProblem1RootExperimentalTangentResidualTolerance
);

Problem1RootTableQueryResult query_problem1_root_table_linear_newton(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_iterations,
    double residual_tolerance,
    double alpha_step_tolerance
);

// 中文注释：推荐外部调用侧使用 seconds 版本，避免把秒制残差容差误传给 scale-free Newton 接口。
Problem1RootTableQueryResult query_problem1_root_table_linear_newton_seconds(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance
);

std::vector<Problem1RootApproximationResult> query_problem1_root_table_quadratic_raw(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method = Problem1RootHessianMethod::ProjectedTangentFiniteDifference,
    double tangent_residual_tolerance = kProblem1RootExperimentalTangentResidualTolerance
);

Problem1RouteBSafeQueryResult query_problem1_root_table_route_b_safe(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method = Problem1RootHessianMethod::ProjectedTangentFiniteDifference,
    double tangent_residual_tolerance = kProblem1RootExperimentalTangentResidualTolerance
);

Problem1RouteBSafeQueryResult query_problem1_root_table_route_b_linear_safe(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
);

Problem1RouteBSafeQueryResult query_problem1_root_table_route_b_quadratic_safe(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method = Problem1RootHessianMethod::ProjectedTangentFiniteDifference,
    double tangent_residual_tolerance = kProblem1RootExperimentalTangentResidualTolerance
);

Problem1RootTable build_problem1_root_table(const Problem1RootTableConfig& config);

}  // namespace spaceship_cpp::problem1
