#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kStrictAlphaThreshold = 1e-5;
constexpr double kStrictTimeThresholdSeconds = 1e-1;
constexpr double kStrictResidualThresholdSeconds = 1e-6;
constexpr double kMediumAlphaThreshold = 1e-3;
constexpr double kMediumTimeThresholdSeconds = 1e3;
constexpr double kLooseAlphaThreshold = 1e-2;
constexpr double kLooseTimeThresholdSeconds = 1e5;

struct Problem1RootVirtualGridNode {
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double nu_A_node = 0.0;
    double nu_B_node = 0.0;
    double theta_A_node = 0.0;
    double dnu_A = 0.0;
    double dnu_B = 0.0;
    double dtheta_A = 0.0;
};

struct NewtonConfig {
    std::string name;
    int max_iterations = 30;
    double residual_tolerance_seconds = 1e-6;
    double alpha_step_tolerance = 1e-12;
    bool enable_backtracking = false;
    int max_backtracking_steps = 0;
};

double wrapped_alpha_distance(double a1, double a2) {
    return std::abs(normalize_angle_minus_pi_pi(a1 - a2));
}

bool is_strict_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <=
            kStrictAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <=
            kStrictTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kStrictResidualThresholdSeconds;
}

bool is_medium_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <=
            kMediumAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <=
            kMediumTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= 1e-4;
}

bool is_loose_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <=
            kLooseAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <=
            kLooseTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= 1e-2;
}

Problem1RootVirtualGridNode make_virtual_node(
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int nu_A_index,
    int nu_B_index,
    int theta_A_index
) {
    Problem1RootVirtualGridNode node{};
    node.nu_A_index = nu_A_index;
    node.nu_B_index = nu_B_index;
    node.theta_A_index = theta_A_index;
    node.nu_A_node = normalize_angle_0_2pi(static_cast<double>(nu_A_index) * kVirtualGridStepRadians);
    node.nu_B_node = normalize_angle_0_2pi(static_cast<double>(nu_B_index) * kVirtualGridStepRadians);
    node.theta_A_node = normalize_angle_0_2pi(static_cast<double>(theta_A_index) * kVirtualGridStepRadians);
    node.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

void print_trace_excerpt(const spaceship_cpp::problem1::Problem1RootRefinementResult& refined) {
    const auto& trace = refined.diagnostic.trace;
    const std::size_t n = trace.size();
    if (n == 0) {
        std::cout << "  trace=<empty>\n";
        return;
    }
    std::cout << "  trace_first_last_steps\n";
    for (std::size_t i = 0; i < n; ++i) {
        const bool emit = i < 5 || i + 5 >= n;
        if (!emit) {
            continue;
        }
        const auto& step = trace[i];
        std::cout << "    iteration=" << step.iteration
                  << " alpha=" << step.alpha
                  << " residual_seconds=" << step.residual_seconds
                  << " R_alpha=" << step.R_alpha
                  << " delta_alpha=" << step.delta_alpha
                  << " derivative_valid=" << step.derivative_valid
                  << " residual_increased=" << step.residual_increased
                  << " reason=" << step.reason << '\n';
    }
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const planet_params::PlanetId departure_planet = planet_params::PlanetId::Earth;
    const planet_params::PlanetId target_planet = planet_params::PlanetId::Mars;
    const double query_nu_A = 4.555309e+00;
    const double query_nu_B = 5.986479e+00;
    const double query_theta_A = 2.949606e+00;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const int mismatch_k = 1;
    const Problem1RootVirtualGridNode node =
        make_virtual_node(query_nu_A, query_nu_B, query_theta_A, 131, 172, 85);

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "virtual_2deg_newton_rescue_debug\n";
    std::cout << "group=mid_cell\n";
    std::cout << "query_nu_A=" << query_nu_A << '\n';
    std::cout << "query_nu_B=" << query_nu_B << '\n';
    std::cout << "query_theta_A=" << query_theta_A << '\n';
    std::cout << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
    std::cout << "node_offset_degrees=("
              << (node.dnu_A * 180.0 / kPi) << ","
              << (node.dnu_B * 180.0 / kPi) << ","
              << (node.dtheta_A * 180.0 / kPi) << ")\n";
    std::cout << "mismatch_k=" << mismatch_k << '\n';

    auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        max_transfer_revolution,
        max_target_revolution);
    for (auto& branch : exact_branches) {
        if (!branch.valid) {
            continue;
        }
        const auto polished = problem1::refine_problem1_root_branch_newton_seconds(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            branch.transfer_revolution,
            branch.target_revolution,
            branch.encounter_global_angle,
            30,
            1e-6,
            1e-14);
        if (polished.valid) {
            branch = polished;
        }
    }

    auto node_branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        node.nu_A_node,
        node.nu_B_node,
        node.theta_A_node,
        max_transfer_revolution,
        max_target_revolution);

    std::vector<problem1::Problem1SolutionBranch> exact_k_branches;
    std::vector<std::pair<int, problem1::Problem1SolutionBranch>> node_k_branches;
    for (const auto& branch : exact_branches) {
        if (branch.valid && branch.transfer_revolution == mismatch_k) {
            exact_k_branches.push_back(branch);
        }
    }
    for (std::size_t i = 0; i < node_branches.size(); ++i) {
        if (node_branches[i].valid && node_branches[i].transfer_revolution == mismatch_k) {
            node_k_branches.push_back({static_cast<int>(i), node_branches[i]});
        }
    }
    std::sort(exact_k_branches.begin(), exact_k_branches.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
    });
    std::sort(node_k_branches.begin(), node_k_branches.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.second.time_of_flight_seconds < rhs.second.time_of_flight_seconds;
    });

    std::cout << "exact_branches_in_k\n";
    for (std::size_t i = 0; i < exact_k_branches.size(); ++i) {
        const auto& b = exact_k_branches[i];
        std::cout << "  exact_sorted_index=" << i
                  << " k=" << b.transfer_revolution
                  << " q=" << b.target_revolution
                  << " time_of_flight_seconds=" << b.time_of_flight_seconds
                  << " alpha=" << b.encounter_global_angle
                  << " residual_seconds=" << b.residual_seconds << '\n';
    }
    std::cout << "nearest_node_branches_in_k\n";
    for (const auto& [idx, b] : node_k_branches) {
        std::cout << "  node_original_index=" << idx
                  << " k=" << b.transfer_revolution
                  << " q=" << b.target_revolution
                  << " time_of_flight_seconds=" << b.time_of_flight_seconds
                  << " alpha=" << b.encounter_global_angle
                  << " residual_seconds=" << b.residual_seconds << '\n';
    }

    const std::vector<int> source_indices{2, 3};
    const std::vector<NewtonConfig> configs{
        {"current", 30, 1e-6, 1e-12, false, 0},
        {"backtracking", 30, 1e-6, 1e-12, true, 8},
        {"more_iterations", 80, 1e-6, 1e-14, true, 12},
        {"relaxed_residual", 80, 1e-5, 1e-14, true, 12},
    };

    for (int source_index : source_indices) {
        const auto& source = node_branches.at(static_cast<std::size_t>(source_index));
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            node.nu_A_node,
            node.nu_B_node,
            node.theta_A_node,
            source,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        const double raw_alpha = normalize_angle_0_2pi(
            attached.encounter_global_angle +
            attached.d_encounter_global_angle_d_nu_A * node.dnu_A +
            attached.d_encounter_global_angle_d_nu_B * node.dnu_B +
            attached.d_encounter_global_angle_d_theta_A * node.dtheta_A);
        const auto raw_residual = problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            raw_alpha,
            source.transfer_revolution,
            source.target_revolution);

        const problem1::Problem1SolutionBranch* reference_exact = nullptr;
        double best_reference_gap = std::numeric_limits<double>::infinity();
        for (const auto& exact : exact_k_branches) {
            const double gap = std::abs(exact.time_of_flight_seconds - source.time_of_flight_seconds);
            if (gap < best_reference_gap) {
                best_reference_gap = gap;
                reference_exact = &exact;
            }
        }
        assert(reference_exact != nullptr);

        std::cout << "source_branch\n";
        std::cout << "source_branch_index=" << source_index << '\n';
        std::cout << "source_time=" << source.time_of_flight_seconds << '\n';
        std::cout << "source_alpha=" << source.encounter_global_angle << '\n';
        std::cout << "raw_predicted_alpha=" << raw_alpha << '\n';
        std::cout << "raw_residual_seconds=" << raw_residual.residual_seconds << '\n';
        std::cout << "reference_exact_time_by_source=" << reference_exact->time_of_flight_seconds << '\n';
        std::cout << "reference_exact_alpha_by_source=" << reference_exact->encounter_global_angle << '\n';

        for (const auto& config : configs) {
            const auto refined = problem1::refine_problem1_root_branch_newton_diagnostic_seconds(
                departure_planet,
                target_planet,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                source.transfer_revolution,
                source.target_revolution,
                raw_alpha,
                config.max_iterations,
                config.residual_tolerance_seconds,
                config.alpha_step_tolerance,
                std::numeric_limits<double>::infinity(),
                config.enable_backtracking,
                config.max_backtracking_steps,
                problem1::Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
            const problem1::Problem1SolutionBranch* nearest_exact = reference_exact;
            if (refined.valid) {
                double best_gap = std::numeric_limits<double>::infinity();
                for (const auto& exact : exact_k_branches) {
                    const double gap = std::abs(exact.time_of_flight_seconds - refined.branch.time_of_flight_seconds);
                    if (gap < best_gap) {
                        best_gap = gap;
                        nearest_exact = &exact;
                    }
                }
            }
            const double time_error = refined.valid
                ? std::abs(refined.branch.time_of_flight_seconds - nearest_exact->time_of_flight_seconds)
                : std::numeric_limits<double>::quiet_NaN();
            const double alpha_error = refined.valid
                ? wrapped_alpha_distance(refined.branch.encounter_global_angle, nearest_exact->encounter_global_angle)
                : std::numeric_limits<double>::quiet_NaN();
            const double final_residual_abs = std::abs(refined.diagnostic.final_residual_seconds);
            const bool final_residual_within_requested_tolerance =
                std::isfinite(final_residual_abs) &&
                final_residual_abs <= config.residual_tolerance_seconds;
            const bool strict = refined.valid && is_strict_match(refined.branch, *nearest_exact);
            const bool medium = refined.valid && is_medium_match(refined.branch, *nearest_exact);
            const bool loose = refined.valid && is_loose_match(refined.branch, *nearest_exact);

            std::cout << "config=" << config.name << '\n';
            std::cout << "  source_branch_index=" << source_index << '\n';
            std::cout << "  source_time=" << source.time_of_flight_seconds << '\n';
            std::cout << "  source_alpha=" << source.encounter_global_angle << '\n';
            std::cout << "  raw_predicted_alpha=" << raw_alpha << '\n';
            std::cout << "  raw_residual_seconds=" << raw_residual.residual_seconds << '\n';
            std::cout << "  refined_valid=" << refined.valid << '\n';
            std::cout << "  diagnostic_converged=" << refined.diagnostic.converged << '\n';
            std::cout << "  invalid_reason=" << refined.diagnostic.invalid_reason << '\n';
            std::cout << "  iterations=" << refined.diagnostic.iterations << '\n';
            std::cout << "  initial_residual_seconds=" << refined.diagnostic.initial_residual_seconds << '\n';
            std::cout << "  final_residual_seconds=" << refined.diagnostic.final_residual_seconds << '\n';
            std::cout << "  final_residual_abs_seconds=" << final_residual_abs << '\n';
            std::cout << "  final_residual_within_requested_tolerance="
                      << final_residual_within_requested_tolerance << '\n';
            std::cout << "  final_alpha="
                      << (refined.valid ? refined.branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN()) << '\n';
            std::cout << "  final_time_of_flight_seconds="
                      << (refined.valid ? refined.branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN()) << '\n';
            std::cout << "  nearest_exact_root_by_time=" << nearest_exact->time_of_flight_seconds << '\n';
            std::cout << "  time_error_to_nearest_exact=" << time_error << '\n';
            std::cout << "  alpha_wrapped_error_to_nearest_exact=" << alpha_error << '\n';
            std::cout << "  strict_match=" << strict << '\n';
            std::cout << "  medium_match=" << medium << '\n';
            std::cout << "  loose_match=" << loose << '\n';
            if (!refined.valid && final_residual_within_requested_tolerance) {
                std::cout << "  unexpected_invalid_after_residual_reached=1\n";
            } else if (!refined.valid) {
                std::cout << "  invalid_explanation=residual_above_requested_tolerance\n";
            }
            print_trace_excerpt(refined);
        }
    }

    return 0;
}
