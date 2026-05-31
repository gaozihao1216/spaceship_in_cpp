#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepDegrees = 2.0;
constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kNearNodeOffsetDegrees = 0.1;
constexpr double kMidCellOffsetDegrees = 1.0;
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

struct QuerySample {
    std::string group_name;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

struct CandidateDebug {
    spaceship_cpp::problem1::Problem1SolutionBranch node_branch;
    bool differentiated = false;
    spaceship_cpp::problem1::Problem1SolutionBranch differentiated_branch;
    bool raw_valid = false;
    double raw_alpha = std::numeric_limits<double>::quiet_NaN();
    double raw_alpha_error = std::numeric_limits<double>::quiet_NaN();
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double raw_transfer_time_seconds = std::numeric_limits<double>::quiet_NaN();
    double raw_time_error_seconds = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1RootRefinementResult refined;
};

double radians_to_degrees(double radians) {
    return radians * 180.0 / kPi;
}

double alpha_error(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return std::abs(normalize_angle_minus_pi_pi(lhs.encounter_global_angle - rhs.encounter_global_angle));
}

bool is_strict_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return alpha_error(candidate, exact) <= kStrictAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= kStrictTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kStrictResidualThresholdSeconds;
}

bool is_medium_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return alpha_error(candidate, exact) <= kMediumAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= kMediumTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= 1e-4;
}

bool is_loose_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return alpha_error(candidate, exact) <= kLooseAlphaThreshold &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <= kLooseTimeThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= 1e-2;
}

Problem1RootVirtualGridNode find_nearest_virtual_root_table_node(
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double step_radians = kVirtualGridStepRadians
) {
    const int axis_count = static_cast<int>(std::llround(kTwoPi / step_radians));
    const auto nearest_index = [&](double angle) {
        long long index = std::llround(normalize_angle_0_2pi(angle) / step_radians);
        index %= axis_count;
        if (index < 0) {
            index += axis_count;
        }
        return static_cast<int>(index);
    };

    Problem1RootVirtualGridNode node{};
    node.nu_A_index = nearest_index(query_nu_A);
    node.nu_B_index = nearest_index(query_nu_B);
    node.theta_A_index = nearest_index(query_theta_A);
    node.nu_A_node = normalize_angle_0_2pi(static_cast<double>(node.nu_A_index) * step_radians);
    node.nu_B_node = normalize_angle_0_2pi(static_cast<double>(node.nu_B_index) * step_radians);
    node.theta_A_node = normalize_angle_0_2pi(static_cast<double>(node.theta_A_index) * step_radians);
    node.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

std::vector<QuerySample> build_samples_per_group(int count_per_group) {
    std::vector<QuerySample> samples;
    const double near_offset = kNearNodeOffsetDegrees * kPi / 180.0;
    const double mid_offset = kMidCellOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < count_per_group; ++i) {
        const double base_nu_A = normalize_angle_0_2pi(0.37 + static_cast<double>(i) * 2.3999632297);
        const double base_nu_B = normalize_angle_0_2pi(1.11 + static_cast<double>(i) * 1.7548776662);
        const double base_theta_A = normalize_angle_0_2pi(0.23 + static_cast<double>(i) * 0.9182736455);
        const auto near_node = find_nearest_virtual_root_table_node(base_nu_A, base_nu_B, base_theta_A);
        samples.push_back({"near_node",
                           normalize_angle_0_2pi(near_node.nu_A_node + near_offset),
                           normalize_angle_0_2pi(near_node.nu_B_node + near_offset),
                           normalize_angle_0_2pi(near_node.theta_A_node + near_offset)});

        const double base2_nu_A = normalize_angle_0_2pi(0.91 + static_cast<double>(i) * 1.8123456789);
        const double base2_nu_B = normalize_angle_0_2pi(1.41 + static_cast<double>(i) * 2.2718281828);
        const double base2_theta_A = normalize_angle_0_2pi(0.63 + static_cast<double>(i) * 1.1447298860);
        const auto mid_node = find_nearest_virtual_root_table_node(base2_nu_A, base2_nu_B, base2_theta_A);
        samples.push_back({"mid_cell",
                           normalize_angle_0_2pi(mid_node.nu_A_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.nu_B_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.theta_A_node + mid_offset)});
    }

    const double earth_period =
        spaceship_cpp::planet_params::planet_orbital_period(spaceship_cpp::planet_params::PlanetId::Earth);
    const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};
    for (int i = 0; i < count_per_group; ++i) {
        const double launch_fraction = std::fmod(0.17 + 0.31 * static_cast<double>(i), 1.0);
        const double launch_time = launch_fraction * earth_period;
        const auto departure_state = spaceship_cpp::planet_params::planet_state_at_time(
            spaceship_cpp::planet_params::PlanetId::Earth, launch_time);
        const auto target_state = spaceship_cpp::planet_params::planet_state_at_time(
            spaceship_cpp::planet_params::PlanetId::Mars, launch_time);
        const double transfer_perihelion_angle =
            transfer_perihelion_angles[static_cast<std::size_t>(i) % transfer_perihelion_angles.size()];
        samples.push_back({"physical_launch",
                           departure_state.varphi,
                           target_state.varphi,
                           normalize_angle_0_2pi(departure_state.theta_global - transfer_perihelion_angle)});
    }
    return samples;
}

std::map<std::pair<int, int>, int> count_kq_distribution(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::map<std::pair<int, int>, int> counts;
    for (const auto& branch : branches) {
        if (!branch.valid) {
            continue;
        }
        counts[{branch.transfer_revolution, branch.target_revolution}] += 1;
    }
    return counts;
}

int count_multi_root_kq(const std::map<std::pair<int, int>, int>& counts) {
    int total = 0;
    for (const auto& [kq, count] : counts) {
        if (count > 1) {
            total += 1;
        }
    }
    return total;
}

void print_kq_distribution(
    const std::string& label,
    const std::map<std::pair<int, int>, int>& counts
) {
    std::cout << label << "_kq_distribution\n";
    for (const auto& [kq, count] : counts) {
        std::cout << "  (k=" << kq.first << ",q=" << kq.second << ")=" << count << '\n';
    }
}

struct ExactBranchDebugClassification {
    std::string primary;
    bool multi_root_same_kq_ambiguous = false;
};

void print_trace_excerpt(const spaceship_cpp::problem1::Problem1RootRefinementResult& refined) {
    const auto& trace = refined.diagnostic.trace;
    if (trace.empty()) {
        std::cout << "    trace: <empty>\n";
        return;
    }
    std::cout << "    trace_excerpt\n";
    for (std::size_t i = 0; i < trace.size(); ++i) {
        const bool print =
            i < 5 || i + 3 >= trace.size();
        if (!print) {
            continue;
        }
        const auto& step = trace[i];
        std::cout << "      iteration=" << step.iteration
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
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const double newton_residual_tolerance_seconds = 1e-6;
    const auto samples = build_samples_per_group(12);

    std::map<std::string, bool> emitted_group;
    bool printed_trace_for_first_failed_exact_branch = false;
    int emitted_sample_count = 0;

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Virtual 2-degree failure microscope\n";
    std::cout << "newton_residual_tolerance_seconds=" << newton_residual_tolerance_seconds << '\n';

    for (const QuerySample& sample : samples) {
        if (emitted_group[sample.group_name]) {
            continue;
        }

        const auto node = find_nearest_virtual_root_table_node(
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A);
        const auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            max_transfer_revolution,
            max_target_revolution);
        const auto node_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            node.nu_A_node,
            node.nu_B_node,
            node.theta_A_node,
            max_transfer_revolution,
            max_target_revolution);

        bool sample_has_failure = false;
        for (const auto& exact_branch : exact_branches) {
            if (!exact_branch.valid) {
                continue;
            }
            std::vector<const problem1::Problem1SolutionBranch*> same_kq_node_candidates;
            for (const auto& node_branch : node_branches) {
                if (!node_branch.valid) {
                    continue;
                }
                if (node_branch.transfer_revolution == exact_branch.transfer_revolution &&
                    node_branch.target_revolution == exact_branch.target_revolution) {
                    same_kq_node_candidates.push_back(&node_branch);
                }
            }
            if (same_kq_node_candidates.empty()) {
                sample_has_failure = true;
                break;
            }

            bool any_strict = false;
            for (const auto* node_branch_ptr : same_kq_node_candidates) {
                const auto differentiated = problem1::attach_problem1_root_derivatives_with_mode(
                    departure_planet,
                    target_planet,
                    node.nu_A_node,
                    node.nu_B_node,
                    node.theta_A_node,
                    *node_branch_ptr,
                    problem1::Problem1RootDerivativeMode::AnalyticOnly,
                    1e-6);
                if (!differentiated.valid || !differentiated.derivatives_available) {
                    continue;
                }
                const double alpha_linear = normalize_angle_0_2pi(
                    differentiated.encounter_global_angle +
                    differentiated.d_encounter_global_angle_d_nu_A * node.dnu_A +
                    differentiated.d_encounter_global_angle_d_nu_B * node.dnu_B +
                    differentiated.d_encounter_global_angle_d_theta_A * node.dtheta_A);
                const auto refined = problem1::refine_problem1_root_branch_newton_diagnostic_seconds(
                    departure_planet,
                    target_planet,
                    sample.query_nu_A,
                    sample.query_nu_B,
                    sample.query_theta_A,
                    node_branch_ptr->transfer_revolution,
                    node_branch_ptr->target_revolution,
                    alpha_linear,
                    30,
                    newton_residual_tolerance_seconds,
                    1e-12,
                    std::numeric_limits<double>::infinity(),
                    true,
                    8,
                    problem1::Problem1RootDerivativeMode::AnalyticOnly,
                    1e-6);
                if (refined.valid && is_strict_match(refined.branch, exact_branch)) {
                    any_strict = true;
                    break;
                }
            }
            if (!any_strict) {
                sample_has_failure = true;
                break;
            }
        }

        if (!sample_has_failure) {
            continue;
        }

        emitted_group[sample.group_name] = true;
        emitted_sample_count += 1;

        std::cout << "\n=== debug_sample_begin ===\n";
        std::cout << "group=" << sample.group_name << '\n';
        std::cout << "query_nu_A=" << sample.query_nu_A << '\n';
        std::cout << "query_nu_B=" << sample.query_nu_B << '\n';
        std::cout << "query_theta_A=" << sample.query_theta_A << '\n';
        std::cout << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
        std::cout << "node_nu_A=" << node.nu_A_node << '\n';
        std::cout << "node_nu_B=" << node.nu_B_node << '\n';
        std::cout << "node_theta_A=" << node.theta_A_node << '\n';
        std::cout << "dnu_A_degrees=" << radians_to_degrees(node.dnu_A) << '\n';
        std::cout << "dnu_B_degrees=" << radians_to_degrees(node.dnu_B) << '\n';
        std::cout << "dtheta_A_degrees=" << radians_to_degrees(node.dtheta_A) << '\n';
        std::cout << "exact_branch_count=" << exact_branches.size() << '\n';
        std::cout << "node_branch_count=" << node_branches.size() << '\n';

        const auto exact_kq = count_kq_distribution(exact_branches);
        const auto node_kq = count_kq_distribution(node_branches);
        print_kq_distribution("exact", exact_kq);
        print_kq_distribution("node", node_kq);

        int exact_kq_missing_in_node_count = 0;
        for (const auto& [kq, count] : exact_kq) {
            if (!node_kq.contains(kq)) {
                exact_kq_missing_in_node_count += 1;
            }
        }
        int node_kq_not_in_exact_count = 0;
        for (const auto& [kq, count] : node_kq) {
            if (!exact_kq.contains(kq)) {
                node_kq_not_in_exact_count += 1;
            }
        }
        std::cout << "exact_kq_missing_in_node_count=" << exact_kq_missing_in_node_count << '\n';
        std::cout << "node_kq_not_in_exact_count=" << node_kq_not_in_exact_count << '\n';
        std::cout << "multi_root_kq_exact_count=" << count_multi_root_kq(exact_kq) << '\n';
        std::cout << "multi_root_kq_node_count=" << count_multi_root_kq(node_kq) << '\n';

        for (std::size_t exact_index = 0; exact_index < exact_branches.size(); ++exact_index) {
            const auto& exact_branch = exact_branches[exact_index];
            if (!exact_branch.valid) {
                continue;
            }
            std::cout << "\n  exact_branch_index=" << exact_index
                      << " k=" << exact_branch.transfer_revolution
                      << " q=" << exact_branch.target_revolution
                      << " alpha=" << exact_branch.encounter_global_angle
                      << " tof=" << exact_branch.time_of_flight_seconds
                      << " residual_seconds=" << exact_branch.residual_seconds << '\n';

            std::vector<CandidateDebug> candidates;
            for (const auto& node_branch : node_branches) {
                if (!node_branch.valid ||
                    node_branch.transfer_revolution != exact_branch.transfer_revolution ||
                    node_branch.target_revolution != exact_branch.target_revolution) {
                    continue;
                }
                CandidateDebug candidate{};
                candidate.node_branch = node_branch;
                candidate.differentiated_branch = problem1::attach_problem1_root_derivatives_with_mode(
                    departure_planet,
                    target_planet,
                    node.nu_A_node,
                    node.nu_B_node,
                    node.theta_A_node,
                    node_branch,
                    problem1::Problem1RootDerivativeMode::AnalyticOnly,
                    1e-6);
                candidate.differentiated = candidate.differentiated_branch.valid &&
                    candidate.differentiated_branch.derivatives_available;

                if (candidate.differentiated) {
                    candidate.raw_alpha = normalize_angle_0_2pi(
                        candidate.differentiated_branch.encounter_global_angle +
                        candidate.differentiated_branch.d_encounter_global_angle_d_nu_A * node.dnu_A +
                        candidate.differentiated_branch.d_encounter_global_angle_d_nu_B * node.dnu_B +
                        candidate.differentiated_branch.d_encounter_global_angle_d_theta_A * node.dtheta_A);
                    const auto residual = problem1::evaluate_problem1_root_residual(
                        departure_planet,
                        target_planet,
                        sample.query_nu_A,
                        sample.query_nu_B,
                        sample.query_theta_A,
                        candidate.raw_alpha,
                        node_branch.transfer_revolution,
                        node_branch.target_revolution);
                    if (residual.valid) {
                        candidate.raw_valid = true;
                        candidate.raw_alpha_error = std::abs(normalize_angle_minus_pi_pi(
                            candidate.raw_alpha - exact_branch.encounter_global_angle));
                        candidate.raw_residual_seconds = residual.residual_seconds;
                        candidate.raw_transfer_time_seconds = residual.transfer_time_seconds;
                        candidate.raw_time_error_seconds =
                            std::abs(residual.transfer_time_seconds - exact_branch.time_of_flight_seconds);
                    }
                    candidate.refined = problem1::refine_problem1_root_branch_newton_diagnostic_seconds(
                        departure_planet,
                        target_planet,
                        sample.query_nu_A,
                        sample.query_nu_B,
                        sample.query_theta_A,
                        node_branch.transfer_revolution,
                        node_branch.target_revolution,
                        candidate.raw_alpha,
                        30,
                        newton_residual_tolerance_seconds,
                        1e-12,
                        std::numeric_limits<double>::infinity(),
                        true,
                        8,
                        problem1::Problem1RootDerivativeMode::AnalyticOnly,
                        1e-6);
                }
                candidates.push_back(candidate);
            }

            ExactBranchDebugClassification classification{};
            const auto exact_kq_it = exact_kq.find(
                {exact_branch.transfer_revolution, exact_branch.target_revolution});
            classification.multi_root_same_kq_ambiguous =
                (exact_kq_it != exact_kq.end() && exact_kq_it->second > 1) ||
                static_cast<int>(candidates.size()) > 1;

            if (candidates.empty()) {
                classification.primary = "no_same_kq_seed";
                std::cout << "    classification=" << classification.primary
                          << " multi_root_same_kq_ambiguous=" << classification.multi_root_same_kq_ambiguous << '\n';
                continue;
            }

            bool any_strict = false;
            bool any_medium = false;
            bool any_loose = false;
            bool any_invalid = false;
            bool any_raw_far = false;

            const double angular_distance_norm_degrees = std::sqrt(
                std::pow(radians_to_degrees(node.dnu_A), 2.0) +
                std::pow(radians_to_degrees(node.dnu_B), 2.0) +
                std::pow(radians_to_degrees(node.dtheta_A), 2.0));

            for (std::size_t candidate_index = 0; candidate_index < candidates.size(); ++candidate_index) {
                const auto& candidate = candidates[candidate_index];
                std::cout << "    candidate_index=" << candidate_index
                          << " node_alpha=" << candidate.node_branch.encounter_global_angle
                          << " node_tof=" << candidate.node_branch.time_of_flight_seconds
                          << " raw_alpha=" << candidate.raw_alpha
                          << " raw_alpha_error=" << candidate.raw_alpha_error
                          << " raw_residual_seconds=" << candidate.raw_residual_seconds
                          << " raw_transfer_time_seconds=" << candidate.raw_transfer_time_seconds
                          << " raw_time_error_seconds=" << candidate.raw_time_error_seconds
                          << " node_to_query_angular_distance_norm_degrees=" << angular_distance_norm_degrees << '\n';

                if (!candidate.raw_valid ||
                    candidate.raw_alpha_error > 0.1 ||
                    std::abs(candidate.raw_residual_seconds) > 1e5) {
                    any_raw_far = true;
                }

                const auto& refined = candidate.refined;
                bool strict = false;
                bool medium = false;
                bool loose = false;
                double final_alpha_error = std::numeric_limits<double>::quiet_NaN();
                double final_time_error = std::numeric_limits<double>::quiet_NaN();
                if (refined.valid) {
                    final_alpha_error = alpha_error(refined.branch, exact_branch);
                    final_time_error =
                        std::abs(refined.branch.time_of_flight_seconds - exact_branch.time_of_flight_seconds);
                    strict = is_strict_match(refined.branch, exact_branch);
                    medium = is_medium_match(refined.branch, exact_branch);
                    loose = is_loose_match(refined.branch, exact_branch);
                    any_strict = any_strict || strict;
                    any_medium = any_medium || medium;
                    any_loose = any_loose || loose;
                } else {
                    any_invalid = true;
                }

                std::cout << "      refined_valid=" << refined.valid
                          << " converged=" << refined.diagnostic.converged
                          << " invalid_reason=" << refined.diagnostic.invalid_reason
                          << " iterations=" << refined.diagnostic.iterations
                          << " initial_residual_seconds=" << refined.diagnostic.initial_residual_seconds
                          << " final_residual_seconds=" << refined.diagnostic.final_residual_seconds
                          << " final_alpha=" << refined.diagnostic.final_alpha
                          << " final_alpha_error=" << final_alpha_error
                          << " final_time_error=" << final_time_error
                          << " strict_match=" << strict
                          << " medium_match=" << medium
                          << " loose_match=" << loose << '\n';
            }

            if (any_strict) {
                classification.primary = "newton_converged_correct_root";
            } else if (any_medium || any_loose) {
                classification.primary = "matching_threshold_too_strict";
            } else if (any_invalid) {
                classification.primary = "newton_invalid";
            } else if (any_raw_far) {
                classification.primary = "same_kq_seed_exists_but_raw_far";
            } else {
                classification.primary = "newton_converged_but_wrong_root";
            }

            std::cout << "    classification=" << classification.primary
                      << " multi_root_same_kq_ambiguous=" << classification.multi_root_same_kq_ambiguous << '\n';

            if (!printed_trace_for_first_failed_exact_branch &&
                classification.primary != "newton_converged_correct_root") {
                printed_trace_for_first_failed_exact_branch = true;
                std::size_t best_index = 0;
                double best_score = std::numeric_limits<double>::infinity();
                for (std::size_t candidate_index = 0; candidate_index < candidates.size(); ++candidate_index) {
                    const auto& candidate = candidates[candidate_index];
                    double score = candidate.raw_alpha_error;
                    if (!std::isfinite(score)) {
                        score = std::numeric_limits<double>::infinity();
                    }
                    if (score < best_score) {
                        best_score = score;
                        best_index = candidate_index;
                    }
                }
                std::cout << "    detailed_trace_for_first_failed_exact_branch candidate_index=" << best_index << '\n';
                print_trace_excerpt(candidates[best_index].refined);
            }
        }

        std::cout << "=== debug_sample_end ===\n";
        if (emitted_sample_count >= 3) {
            break;
        }
    }

    assert(emitted_sample_count > 0);
    return 0;
}
