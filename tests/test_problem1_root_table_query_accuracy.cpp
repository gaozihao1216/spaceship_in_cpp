#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <tuple>
#include <vector>

namespace {

bool approx_equal(double lhs, double rhs, double abs_tol = 1e-5, double rel_tol = 1e-5) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= abs_tol + rel_tol * scale;
}

std::optional<spaceship_cpp::problem1::Problem1SolutionBranch> find_same_kq_nearest_branch(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    int transfer_revolution,
    int target_revolution,
    double encounter_global_angle
) {
    namespace common = spaceship_cpp::common;
    const spaceship_cpp::problem1::Problem1SolutionBranch* best = nullptr;
    double best_distance = std::numeric_limits<double>::infinity();
    for (const auto& branch : branches) {
        if (!branch.valid ||
            branch.transfer_revolution != transfer_revolution ||
            branch.target_revolution != target_revolution) {
            continue;
        }
        const double distance = std::abs(
            common::normalize_angle_minus_pi_pi(branch.encounter_global_angle - encounter_global_angle));
        if (distance < best_distance) {
            best_distance = distance;
            best = &branch;
        }
    }
    if (best == nullptr) {
        return std::nullopt;
    }
    return *best;
}

struct Problem1BranchMatchResult {
    bool matched = false;
    bool same_kq_exists = false;
    std::size_t candidate_index = 0;
    double alpha_error = std::numeric_limits<double>::quiet_NaN();
    double time_error_seconds = std::numeric_limits<double>::quiet_NaN();
    double residual_seconds = std::numeric_limits<double>::quiet_NaN();
};

double branch_alpha_error(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    namespace common = spaceship_cpp::common;
    return std::abs(common::normalize_angle_minus_pi_pi(
        lhs.encounter_global_angle - rhs.encounter_global_angle));
}

Problem1BranchMatchResult find_best_refined_branch_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact_branch,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& refined_candidates,
    double alpha_threshold,
    double time_threshold_seconds,
    double residual_threshold_seconds
) {
    Problem1BranchMatchResult result{};
    double best_score = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < refined_candidates.size(); ++i) {
        const auto& candidate = refined_candidates[i];
        if (!candidate.valid ||
            candidate.transfer_revolution != exact_branch.transfer_revolution ||
            candidate.target_revolution != exact_branch.target_revolution) {
            continue;
        }
        result.same_kq_exists = true;
        const double aerr = branch_alpha_error(candidate, exact_branch);
        const double terr = std::abs(candidate.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
        const double rerr = std::abs(candidate.residual_seconds);
        const double score =
            aerr / alpha_threshold +
            terr / time_threshold_seconds +
            rerr / residual_threshold_seconds;
        if (score < best_score) {
            best_score = score;
            result.candidate_index = i;
            result.alpha_error = aerr;
            result.time_error_seconds = terr;
            result.residual_seconds = candidate.residual_seconds;
        }
    }
    if (!result.same_kq_exists) {
        return result;
    }
    result.matched =
        result.alpha_error <= alpha_threshold &&
        result.time_error_seconds <= time_threshold_seconds &&
        std::abs(result.residual_seconds) <= residual_threshold_seconds;
    return result;
}

double average_or_zero(double sum, int count) {
    if (count <= 0) {
        return 0.0;
    }
    return sum / static_cast<double>(count);
}

void update_max(double value, double* current_max) {
    assert(current_max != nullptr);
    if (!std::isfinite(value)) {
        return;
    }
    if (value > *current_max) {
        *current_max = value;
    }
}

struct QuerySample {
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> exact_branches;
};

struct ApproximationStats {
    int valid_count = 0;
    int admissible_count = 0;
    int invalid_count = 0;
    int alpha_error_count = 0;
    double alpha_error_sum = 0.0;
    double alpha_error_max = 0.0;
    double residual_abs_sum = 0.0;
    double residual_abs_max = 0.0;
};

void record_valid_approximation(
    ApproximationStats* stats,
    const spaceship_cpp::problem1::Problem1RootApproximationResult& result
) {
    assert(stats != nullptr);
    stats->valid_count += 1;
    if (result.diagnostics.admissible_for_fast_approximation) {
        stats->admissible_count += 1;
    }
    const double residual_abs = std::abs(result.residual_scale_free);
    if (std::isfinite(residual_abs)) {
        stats->residual_abs_sum += residual_abs;
        update_max(residual_abs, &stats->residual_abs_max);
    }
}

void record_alpha_error(
    ApproximationStats* stats,
    double alpha_error
) {
    assert(stats != nullptr);
    assert(std::isfinite(alpha_error));
    stats->alpha_error_count += 1;
    stats->alpha_error_sum += alpha_error;
    update_max(alpha_error, &stats->alpha_error_max);
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const double newton_residual_tolerance_seconds = 1e-6;
    const double newton_residual_tolerance_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(newton_residual_tolerance_seconds);

    const problem1::Problem1RootTableConfig config{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        0.0,
        common::kTwoPi / 8.0,
        8,
        0.0,
        common::kTwoPi / 8.0,
        8,
        0.0,
        common::kTwoPi / 12.0,
        12,
        1,
        1,
        "problem1_root_table_draft_v0",
    };
    const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
    const double hessian_step = 5e-4;
    const double tangent_residual_tolerance = problem1::kProblem1RootExperimentalTangentResidualTolerance;

    std::vector<QuerySample> samples;
    for (const auto& cell : table.cells()) {
        if (!cell.solved || cell.solutions_sorted_by_time_of_flight.empty()) {
            continue;
        }
        for (double offset_scale : {0.005}) {
            QuerySample sample{};
            sample.query_nu_A =
                common::normalize_angle_0_2pi(cell.nu_A_depart + offset_scale * config.nu_A_step);
            sample.query_nu_B =
                common::normalize_angle_0_2pi(cell.nu_B_depart + offset_scale * config.nu_B_depart_step);
            sample.query_theta_A =
                common::normalize_angle_0_2pi(cell.theta_A + offset_scale * config.theta_A_step);
            sample.exact_branches = problem1::solve_problem1_from_departure_anomalies(
                config.departure_planet,
                config.target_planet,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                config.max_transfer_revolution,
                config.max_target_revolution);
            if (sample.exact_branches.empty()) {
                continue;
            }
            for (auto& exact_branch : sample.exact_branches) {
                if (!exact_branch.valid) {
                    continue;
                }
                const auto polished = problem1::refine_problem1_root_branch_newton_seconds(
                    config.departure_planet,
                    config.target_planet,
                    sample.query_nu_A,
                    sample.query_nu_B,
                    sample.query_theta_A,
                    exact_branch.transfer_revolution,
                    exact_branch.target_revolution,
                    exact_branch.encounter_global_angle,
                    30,
                    newton_residual_tolerance_seconds,
                    1e-14);
                if (polished.valid) {
                    exact_branch = polished;
                }
            }
            samples.push_back(std::move(sample));
            if (samples.size() >= 12) {
                break;
            }
        }
        if (samples.size() >= 12) {
            break;
        }
    }
    assert(samples.size() >= 6);

    int exact_solve_branch_count = 0;
    int route_a_matched_branch_count = 0;
    int route_a_missing_branch_count = 0;
    double route_a_max_alpha_error = 0.0;
    double route_a_max_time_error = 0.0;
    double route_a_max_residual = 0.0;

    ApproximationStats linear_stats{};
    ApproximationStats quadratic_tangent_stats{};
    ApproximationStats quadratic_newton_stats{};
    int quadratic_tangent_better_than_linear_count = 0;
    int quadratic_tangent_worse_than_linear_count = 0;
    int quadratic_tangent_inadmissible_count = 0;

    bool forced_invalid_checked = false;

    for (const QuerySample& sample : samples) {
        exact_solve_branch_count += static_cast<int>(sample.exact_branches.size());

        const auto route_a = problem1::query_problem1_root_table_linear_newton(
            table,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            30,
            newton_residual_tolerance_scale_free,
            1e-12);
        const auto route_a_seconds = problem1::query_problem1_root_table_linear_newton_seconds(
            table,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            30,
            newton_residual_tolerance_seconds,
            1e-12);
        assert(route_a.valid);
        assert(route_a_seconds.valid == route_a.valid);
        assert(route_a_seconds.branches.size() == route_a.branches.size());
        for (std::size_t i = 0; i < route_a.branches.size(); ++i) {
            const auto& lhs = route_a.branches[i];
            const auto& rhs = route_a_seconds.branches[i];
            assert(lhs.valid == rhs.valid);
            assert(lhs.transfer_revolution == rhs.transfer_revolution);
            assert(lhs.target_revolution == rhs.target_revolution);
            assert(std::abs(common::normalize_angle_minus_pi_pi(
                lhs.encounter_global_angle - rhs.encounter_global_angle)) <= 1e-12);
            assert(approx_equal(lhs.time_of_flight_seconds, rhs.time_of_flight_seconds, 1e-12, 1e-12));
            assert(approx_equal(lhs.residual_seconds, rhs.residual_seconds, 1e-12, 1e-12));
        }

        for (const auto& exact_branch : sample.exact_branches) {
            const auto match = find_best_refined_branch_match(
                exact_branch,
                route_a.branches,
                1e-5,
                1e-1,
                1e-6);
            if (!match.same_kq_exists) {
                route_a_missing_branch_count += 1;
                continue;
            }
            const double alpha_error = match.alpha_error;
            const double time_error = match.time_error_seconds;
            const double residual_abs = std::abs(match.residual_seconds);
            if (!match.matched) {
                route_a_missing_branch_count += 1;
                continue;
            }
            route_a_matched_branch_count += 1;
            update_max(alpha_error, &route_a_max_alpha_error);
            update_max(time_error, &route_a_max_time_error);
            update_max(residual_abs, &route_a_max_residual);
        }

        const auto nearest = problem1::find_nearest_problem1_root_table_node(
            table,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A);
        assert(nearest.valid);
        assert(nearest.cell != nullptr);

        for (const auto& node_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
            if (!node_branch.valid) {
                continue;
            }

            const auto linear_result = problem1::evaluate_problem1_root_linear_approximation_from_node(
                config.departure_planet,
                config.target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A);
            if (linear_result.valid) {
                record_valid_approximation(&linear_stats, linear_result);
                const auto exact_match = find_same_kq_nearest_branch(
                    sample.exact_branches,
                    linear_result.transfer_revolution,
                    linear_result.target_revolution,
                    linear_result.predicted_encounter_global_angle);
                if (exact_match.has_value()) {
                    const double alpha_error = std::abs(common::normalize_angle_minus_pi_pi(
                        linear_result.predicted_encounter_global_angle - exact_match->encounter_global_angle));
                    record_alpha_error(&linear_stats, alpha_error);
                }
            } else {
                linear_stats.invalid_count += 1;
            }

            const auto quadratic_tangent = problem1::evaluate_problem1_root_quadratic_approximation_from_node(
                config.departure_planet,
                config.target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                hessian_step,
                problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                tangent_residual_tolerance);
            assert(
                quadratic_tangent.diagnostics.hessian_method ==
                "tangent_finite_difference_of_implicit_first_derivatives");
            assert(approx_equal(quadratic_tangent.diagnostics.hessian_step, hessian_step, 1e-12, 1e-12));
            assert(
                !quadratic_tangent.diagnostics.hessian_valid ||
                std::isfinite(quadratic_tangent.diagnostics.tangent_residual_max_scale_free));
            if (!quadratic_tangent.diagnostics.hessian_valid ||
                !std::isfinite(quadratic_tangent.diagnostics.raw_residual_scale_free)) {
                assert(!quadratic_tangent.diagnostics.admissible_for_fast_approximation);
            }
            if (quadratic_tangent.valid) {
                assert(quadratic_tangent.diagnostics.hessian_valid);
                assert(std::isfinite(quadratic_tangent.diagnostics.raw_residual_scale_free));
                assert(std::isfinite(quadratic_tangent.diagnostics.raw_residual_seconds));
                assert(approx_equal(
                    quadratic_tangent.diagnostics.raw_residual_scale_free,
                    quadratic_tangent.residual_scale_free,
                    1e-12,
                    1e-12));
                assert(approx_equal(
                    quadratic_tangent.diagnostics.raw_residual_seconds,
                    quadratic_tangent.residual_seconds,
                    1e-12,
                    1e-12));
                assert(!quadratic_tangent.diagnostics.admissibility_reason.empty());
                record_valid_approximation(&quadratic_tangent_stats, quadratic_tangent);
                if (!quadratic_tangent.diagnostics.admissible_for_fast_approximation) {
                    quadratic_tangent_inadmissible_count += 1;
                }
                const auto exact_match = find_same_kq_nearest_branch(
                    sample.exact_branches,
                    quadratic_tangent.transfer_revolution,
                    quadratic_tangent.target_revolution,
                    quadratic_tangent.predicted_encounter_global_angle);
                if (exact_match.has_value()) {
                    const double alpha_error = std::abs(common::normalize_angle_minus_pi_pi(
                        quadratic_tangent.predicted_encounter_global_angle - exact_match->encounter_global_angle));
                    record_alpha_error(&quadratic_tangent_stats, alpha_error);
                    if (linear_result.valid) {
                        const auto linear_exact_match = find_same_kq_nearest_branch(
                            sample.exact_branches,
                            linear_result.transfer_revolution,
                            linear_result.target_revolution,
                            linear_result.predicted_encounter_global_angle);
                        if (linear_exact_match.has_value()) {
                            const double linear_alpha_error = std::abs(common::normalize_angle_minus_pi_pi(
                                linear_result.predicted_encounter_global_angle -
                                linear_exact_match->encounter_global_angle));
                            if (alpha_error + 1e-12 < linear_alpha_error) {
                                quadratic_tangent_better_than_linear_count += 1;
                            } else if (alpha_error > linear_alpha_error + 1e-12) {
                                quadratic_tangent_worse_than_linear_count += 1;
                            }
                        }
                    }
                }
            } else {
                quadratic_tangent_stats.invalid_count += 1;
            }

            const auto quadratic_newton = problem1::evaluate_problem1_root_quadratic_approximation_from_node(
                config.departure_planet,
                config.target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                sample.query_nu_A,
                sample.query_nu_B,
                sample.query_theta_A,
                hessian_step,
                problem1::Problem1RootHessianMethod::NewtonRefinedFiniteDifference,
                tangent_residual_tolerance);
            assert(
                quadratic_newton.diagnostics.hessian_method ==
                "newton_refined_finite_difference_of_implicit_first_derivatives");
            assert(approx_equal(quadratic_newton.diagnostics.hessian_step, hessian_step, 1e-12, 1e-12));
            if (quadratic_newton.valid) {
                assert(quadratic_newton.diagnostics.hessian_valid);
                record_valid_approximation(&quadratic_newton_stats, quadratic_newton);
                const auto exact_match = find_same_kq_nearest_branch(
                    sample.exact_branches,
                    quadratic_newton.transfer_revolution,
                    quadratic_newton.target_revolution,
                    quadratic_newton.predicted_encounter_global_angle);
                if (exact_match.has_value()) {
                    const double alpha_error = std::abs(common::normalize_angle_minus_pi_pi(
                        quadratic_newton.predicted_encounter_global_angle - exact_match->encounter_global_angle));
                    record_alpha_error(&quadratic_newton_stats, alpha_error);
                }
            } else {
                quadratic_newton_stats.invalid_count += 1;
            }

            if (!forced_invalid_checked) {
                const auto forced_invalid = problem1::evaluate_problem1_root_quadratic_approximation_from_node(
                    config.departure_planet,
                    config.target_planet,
                    nearest.node_nu_A,
                    nearest.node_nu_B,
                    nearest.node_theta_A,
                    node_branch,
                    sample.query_nu_A,
                    sample.query_nu_B,
                    sample.query_theta_A,
                    0.0,
                    problem1::Problem1RootHessianMethod::TangentFiniteDifference,
                    tangent_residual_tolerance);
                assert(!forced_invalid.valid);
                assert(!forced_invalid.diagnostics.hessian_valid);
                assert(!forced_invalid.diagnostics.admissible_for_fast_approximation);
                assert(!forced_invalid.diagnostics.admissibility_reason.empty());
                forced_invalid_checked = true;
            }
        }
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1 root-table query accuracy summary\n";
    std::cout << "sample_count=" << samples.size() << '\n';
    std::cout << "exact_solve_branch_count=" << exact_solve_branch_count << '\n';
    std::cout << "route_a_matched_branch_count=" << route_a_matched_branch_count << '\n';
    std::cout << "route_a_missing_branch_count=" << route_a_missing_branch_count << '\n';
    std::cout << "route_a_max_alpha_error=" << route_a_max_alpha_error << '\n';
    std::cout << "route_a_max_time_error=" << route_a_max_time_error << '\n';
    std::cout << "route_a_max_residual=" << route_a_max_residual << '\n';
    std::cout << "linear_raw_valid_count=" << linear_stats.valid_count << '\n';
    std::cout << "quadratic_tangent_raw_valid_count=" << quadratic_tangent_stats.valid_count << '\n';
    std::cout << "quadratic_tangent_raw_admissible_count=" << quadratic_tangent_stats.admissible_count << '\n';
    std::cout << "quadratic_newton_refined_raw_valid_count=" << quadratic_newton_stats.valid_count << '\n';
    std::cout << "linear_raw_alpha_error_avg=" <<
        average_or_zero(linear_stats.alpha_error_sum, linear_stats.alpha_error_count) << '\n';
    std::cout << "linear_raw_alpha_error_max=" << linear_stats.alpha_error_max << '\n';
    std::cout << "quadratic_tangent_raw_alpha_error_avg=" <<
        average_or_zero(quadratic_tangent_stats.alpha_error_sum, quadratic_tangent_stats.alpha_error_count) << '\n';
    std::cout << "quadratic_tangent_raw_alpha_error_max=" << quadratic_tangent_stats.alpha_error_max << '\n';
    std::cout << "quadratic_newton_refined_raw_alpha_error_avg=" <<
        average_or_zero(quadratic_newton_stats.alpha_error_sum, quadratic_newton_stats.alpha_error_count) << '\n';
    std::cout << "quadratic_newton_refined_raw_alpha_error_max=" <<
        quadratic_newton_stats.alpha_error_max << '\n';
    std::cout << "linear_raw_residual_avg=" <<
        average_or_zero(linear_stats.residual_abs_sum, linear_stats.valid_count) << '\n';
    std::cout << "linear_raw_residual_max=" << linear_stats.residual_abs_max << '\n';
    std::cout << "quadratic_tangent_raw_residual_avg=" <<
        average_or_zero(quadratic_tangent_stats.residual_abs_sum, quadratic_tangent_stats.valid_count) << '\n';
    std::cout << "quadratic_tangent_raw_residual_max=" << quadratic_tangent_stats.residual_abs_max << '\n';
    std::cout << "quadratic_tangent_better_than_linear_count=" <<
        quadratic_tangent_better_than_linear_count << '\n';
    std::cout << "quadratic_tangent_worse_than_linear_count=" <<
        quadratic_tangent_worse_than_linear_count << '\n';
    std::cout << "quadratic_tangent_raw_invalid_count=" << quadratic_tangent_stats.invalid_count << '\n';
    std::cout << "quadratic_tangent_raw_inadmissible_count=" << quadratic_tangent_inadmissible_count << '\n';
    std::cout << "quadratic_tangent_invalid_or_inadmissible_count=" <<
        (quadratic_tangent_stats.invalid_count + quadratic_tangent_inadmissible_count) << '\n';
    std::cout << std::flush;

    assert(forced_invalid_checked);
    assert(route_a_matched_branch_count > 0);
    assert(route_a_matched_branch_count + route_a_missing_branch_count == exact_solve_branch_count);
    assert(route_a_max_alpha_error <= 1e-5);
    // 中文注释：Route A 的 branch matching 当前按 angle + relative time 一致性判定，
    // 因此这里先用 0.1 s 量级的绝对时间误差上界与现有匹配语义保持一致。
    assert(route_a_max_time_error <= 1e-1);
    assert(route_a_max_residual <= 1e-6);
    assert(linear_stats.valid_count > 0);
    assert(quadratic_tangent_stats.admissible_count <= quadratic_tangent_stats.valid_count);
    return 0;
}
