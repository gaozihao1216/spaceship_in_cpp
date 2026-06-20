/*
 * 文件作用：实现 Problem 2 出射段 θ' 初扫与 branch 导数估计。
 * 主要工作：离散扫描 Problem 1 多根并在配对 branch 上用差分估计 dφ/dθ'、de/dθ'。
 */
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include "spaceship_cpp/bfs/problem2_angle_frame_adapter.hpp"
#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spaceship_cpp::problem2 {
namespace {

using spaceship_cpp::bfs::problem2_local_periapsis_angle_to_global;
using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::problem1::Problem1Candidate;
using spaceship_cpp::problem1::Problem1ResidualStatus;

double quiet_nan() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool same_revolution_branch(
    const Problem2OutgoingBranchSolution& left,
    const Problem2OutgoingBranchSolution& right
) {
    return left.transfer_revolution == right.transfer_revolution &&
        left.target_revolution == right.target_revolution;
}

double phi_gap_between(const Problem2OutgoingBranchSolution& left, const Problem2OutgoingBranchSolution& right) {
    return std::abs(normalize_angle_minus_pi_pi(right.encounter_global_angle - left.encounter_global_angle));
}

Problem2OutgoingBranchSolution make_branch_solution_from_candidate(const Problem1Candidate& candidate) {
    Problem2OutgoingBranchSolution solution{};
    solution.transfer_revolution = candidate.transfer_revolution;
    solution.target_revolution = candidate.target_revolution;
    solution.encounter_global_angle = candidate.encounter_global_angle;
    solution.outgoing_eccentricity = candidate.residual_result.transfer_e;
    solution.outgoing_semi_latus_rectum = candidate.residual_result.transfer_p;
    solution.relative_problem1_residual = candidate.relative_residual;
    solution.dphi_dtheta_prime = quiet_nan();
    solution.de_dtheta_prime = quiet_nan();
    return solution;
}

void sort_solutions_by_phi(std::vector<Problem2OutgoingBranchSolution>& solutions) {
    std::sort(
        solutions.begin(),
        solutions.end(),
        [](const Problem2OutgoingBranchSolution& left, const Problem2OutgoingBranchSolution& right) {
            if (left.transfer_revolution != right.transfer_revolution) {
                return left.transfer_revolution < right.transfer_revolution;
            }
            if (left.target_revolution != right.target_revolution) {
                return left.target_revolution < right.target_revolution;
            }
            return left.encounter_global_angle < right.encounter_global_angle;
        });
}

void bucket_solutions_by_transfer_revolution(
    const std::vector<Problem2OutgoingBranchSolution>& solutions,
    int max_transfer_revolution,
    std::vector<std::vector<Problem2OutgoingBranchSolution>>& solutions_by_k
) {
    solutions_by_k.clear();
    if (max_transfer_revolution < 0) {
        return;
    }
    solutions_by_k.resize(static_cast<std::size_t>(max_transfer_revolution) + 1U);
    for (const auto& solution : solutions) {
        if (solution.transfer_revolution < 0 || solution.transfer_revolution > max_transfer_revolution) {
            continue;
        }
        solutions_by_k[static_cast<std::size_t>(solution.transfer_revolution)].push_back(solution);
    }
    for (auto& layer : solutions_by_k) {
        sort_solutions_by_phi(layer);
    }
}

void merge_flat_solutions_from_k_layers(Problem2ThetaPrimeNodeSnapshot& snapshot) {
    snapshot.solutions.clear();
    for (const auto& layer : snapshot.solutions_by_k) {
        snapshot.solutions.insert(snapshot.solutions.end(), layer.begin(), layer.end());
    }
    sort_solutions_by_phi(snapshot.solutions);
}

bool same_outgoing_branch_identity(
    const Problem2OutgoingBranchSolution& left,
    const Problem2OutgoingBranchSolution& right
) {
    return left.transfer_revolution == right.transfer_revolution &&
        left.target_revolution == right.target_revolution &&
        left.encounter_global_angle == right.encounter_global_angle;
}

void copy_solution_derivatives(
    const Problem2OutgoingBranchSolution& source,
    Problem2OutgoingBranchSolution& destination
) {
    destination.dphi_dtheta_prime = source.dphi_dtheta_prime;
    destination.de_dtheta_prime = source.de_dtheta_prime;
    destination.has_dphi_dtheta_prime = source.has_dphi_dtheta_prime;
    destination.has_de_dtheta_prime = source.has_de_dtheta_prime;
}

}  // namespace

std::vector<double> build_uniform_theta_prime_local_grid(int theta_prime_count) {
    std::vector<double> grid;
    if (theta_prime_count <= 0) {
        return grid;
    }
    grid.reserve(static_cast<std::size_t>(theta_prime_count));
    const double step = kTwoPi / static_cast<double>(theta_prime_count);
    for (int index = 0; index < theta_prime_count; ++index) {
        grid.push_back(-spaceship_cpp::common::kPi + static_cast<double>(index) * step);
    }
    return grid;
}

double estimate_theta_prime_derivative_central(
    double left_theta_prime,
    double left_value,
    double right_theta_prime,
    double right_value
) {
    const double delta_theta = right_theta_prime - left_theta_prime;
    if (!is_finite(delta_theta) || delta_theta == 0.0) {
        return quiet_nan();
    }
    return (right_value - left_value) / delta_theta;
}

double estimate_theta_prime_derivative_forward(
    double left_theta_prime,
    double left_value,
    double right_theta_prime,
    double right_value
) {
    return estimate_theta_prime_derivative_central(
        left_theta_prime,
        left_value,
        right_theta_prime,
        right_value);
}

double estimate_theta_prime_derivative_backward(
    double left_theta_prime,
    double left_value,
    double right_theta_prime,
    double right_value
) {
    return estimate_theta_prime_derivative_central(
        left_theta_prime,
        left_value,
        right_theta_prime,
        right_value);
}

std::optional<std::size_t> find_best_matching_outgoing_branch_index(
    const Problem2OutgoingBranchSolution& reference,
    const std::vector<Problem2OutgoingBranchSolution>& candidates,
    double max_phi_gap
) {
    if (!(max_phi_gap > 0.0)) {
        return std::nullopt;
    }

    std::optional<std::size_t> best_index;
    double best_gap = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        const auto& candidate = candidates[index];
        if (!same_revolution_branch(reference, candidate)) {
            continue;
        }
        const double gap = phi_gap_between(reference, candidate);
        if (!(gap <= max_phi_gap) || gap >= best_gap) {
            continue;
        }
        best_gap = gap;
        best_index = index;
    }
    return best_index;
}

std::vector<Problem2OutgoingBranchPair> pair_outgoing_branch_solutions_by_phi(
    const std::vector<Problem2OutgoingBranchSolution>& left_solutions,
    const std::vector<Problem2OutgoingBranchSolution>& right_solutions,
    double max_phi_gap
) {
    std::vector<Problem2OutgoingBranchPair> pairs;
    if (!(max_phi_gap > 0.0)) {
        return pairs;
    }

    std::unordered_set<std::size_t> used_right_indices;
    for (std::size_t left_index = 0; left_index < left_solutions.size(); ++left_index) {
        std::optional<std::size_t> best_right_index;
        double best_gap = std::numeric_limits<double>::infinity();
        for (std::size_t right_index = 0; right_index < right_solutions.size(); ++right_index) {
            if (used_right_indices.count(right_index) != 0) {
                continue;
            }
            const auto& left = left_solutions[left_index];
            const auto& right = right_solutions[right_index];
            if (!same_revolution_branch(left, right)) {
                continue;
            }
            const double gap = phi_gap_between(left, right);
            if (!(gap <= max_phi_gap) || gap >= best_gap) {
                continue;
            }
            best_gap = gap;
            best_right_index = right_index;
        }
        if (!best_right_index.has_value()) {
            continue;
        }
        used_right_indices.insert(*best_right_index);
        pairs.push_back(
            Problem2OutgoingBranchPair{
                .left_index = left_index,
                .right_index = *best_right_index,
                .phi_gap = best_gap,
            });
    }
    return pairs;
}

Problem2ThetaPrimeNodeSnapshot evaluate_problem2_theta_prime_node(
    const Problem2ThetaPrimeScanConfig& config,
    double theta_prime_local,
    const problem1::Problem1LaunchState& launch_state
) {
    Problem2ThetaPrimeNodeSnapshot snapshot{};
    snapshot.theta_prime_local = theta_prime_local;
    snapshot.theta_prime_global =
        problem2_local_periapsis_angle_to_global(config.flyby_planet, theta_prime_local);

    problem1::Problem1SolveInput solve_input = config.problem1_solve;
    solve_input.departure_planet = config.flyby_planet;
    solve_input.target_planet = config.target_planet;
    solve_input.launch_time_seconds_since_j2000 = config.flyby_time_seconds_since_j2000;
    solve_input.transfer_perihelion_angle = snapshot.theta_prime_global;

    const std::vector<Problem1Candidate> candidates =
        problem1::solve_problem1_with_launch_state(solve_input, launch_state);
    snapshot.solutions.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (candidate.residual_result.status != Problem1ResidualStatus::Success ||
            !is_finite(candidate.encounter_global_angle) ||
            !is_finite(candidate.residual_result.transfer_e) ||
            !is_finite(candidate.residual_result.transfer_p)) {
            continue;
        }
        snapshot.solutions.push_back(make_branch_solution_from_candidate(candidate));
    }
    sort_solutions_by_phi(snapshot.solutions);
    bucket_solutions_by_transfer_revolution(
        snapshot.solutions,
        config.problem1_solve.max_transfer_revolution,
        snapshot.solutions_by_k);
    return snapshot;
}

Problem2ThetaPrimeNodeSnapshot evaluate_problem2_theta_prime_node(
    const Problem2ThetaPrimeScanConfig& config,
    double theta_prime_local
) {
    const problem1::Problem1LaunchState launch_state = problem1::compute_problem1_launch_state(
        config.flyby_planet,
        config.target_planet,
        config.flyby_time_seconds_since_j2000);
    return evaluate_problem2_theta_prime_node(config, theta_prime_local, launch_state);
}

Problem2ThetaPrimeInitialScanResult run_problem2_theta_prime_initial_scan(
    const Problem2ThetaPrimeScanConfig& config
) {
    Problem2ThetaPrimeInitialScanResult result{};
    if (config.theta_prime_count < 3) {
        result.error_message = "theta_prime_count_must_be_at_least_3";
        return result;
    }
    if (!is_finite(config.flyby_time_seconds_since_j2000)) {
        result.error_message = "non_finite_flyby_time";
        return result;
    }
    if (!(config.branch_phi_pairing_max_gap > 0.0)) {
        result.error_message = "branch_phi_pairing_max_gap_must_be_positive";
        return result;
    }

    const std::vector<double> grid = build_uniform_theta_prime_local_grid(config.theta_prime_count);
    const problem1::Problem1LaunchState launch_state = problem1::compute_problem1_launch_state(
        config.flyby_planet,
        config.target_planet,
        config.flyby_time_seconds_since_j2000);
    result.max_transfer_revolution = config.problem1_solve.max_transfer_revolution;
    result.nodes.resize(grid.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int grid_index = 0; grid_index < static_cast<int>(grid.size()); ++grid_index) {
        result.nodes[static_cast<std::size_t>(grid_index)] =
            evaluate_problem2_theta_prime_node(config, grid[static_cast<std::size_t>(grid_index)], launch_state);
    }

    attach_problem2_theta_prime_solution_derivatives_by_k(
        result.nodes,
        result.max_transfer_revolution,
        config.branch_phi_pairing_max_gap);
    result.ok = true;
    return result;
}

void attach_problem2_theta_prime_solution_derivatives(
    std::vector<Problem2ThetaPrimeNodeSnapshot>& nodes,
    double branch_phi_pairing_max_gap
) {
    if (nodes.size() < 2 || !(branch_phi_pairing_max_gap > 0.0)) {
        return;
    }

    for (auto& snapshot : nodes) {
        for (auto& solution : snapshot.solutions) {
            solution.dphi_dtheta_prime = quiet_nan();
            solution.de_dtheta_prime = quiet_nan();
            solution.has_dphi_dtheta_prime = false;
            solution.has_de_dtheta_prime = false;
        }
    }

    const std::size_t node_count = nodes.size();
    for (std::size_t node_index = 0; node_index < node_count; ++node_index) {
        auto& center_node = nodes[node_index];
        for (auto& center_solution : center_node.solutions) {
            const bool has_left_neighbor = node_index > 0;
            const bool has_right_neighbor = node_index + 1 < node_count;

            std::optional<std::size_t> left_match_index;
            std::optional<std::size_t> right_match_index;
            if (has_left_neighbor) {
                left_match_index = find_best_matching_outgoing_branch_index(
                    center_solution,
                    nodes[node_index - 1].solutions,
                    branch_phi_pairing_max_gap);
            }
            if (has_right_neighbor) {
                right_match_index = find_best_matching_outgoing_branch_index(
                    center_solution,
                    nodes[node_index + 1].solutions,
                    branch_phi_pairing_max_gap);
            }

            if (has_left_neighbor && has_right_neighbor &&
                left_match_index.has_value() && right_match_index.has_value()) {
                const auto& left_solution = nodes[node_index - 1].solutions[*left_match_index];
                const auto& right_solution = nodes[node_index + 1].solutions[*right_match_index];
                const double left_theta = nodes[node_index - 1].theta_prime_local;
                const double right_theta = nodes[node_index + 1].theta_prime_local;

                const double dphi = estimate_theta_prime_derivative_central(
                    left_theta,
                    left_solution.encounter_global_angle,
                    right_theta,
                    right_solution.encounter_global_angle);
                const double de = estimate_theta_prime_derivative_central(
                    left_theta,
                    left_solution.outgoing_eccentricity,
                    right_theta,
                    right_solution.outgoing_eccentricity);
                if (is_finite(dphi)) {
                    center_solution.dphi_dtheta_prime = dphi;
                    center_solution.has_dphi_dtheta_prime = true;
                }
                if (is_finite(de)) {
                    center_solution.de_dtheta_prime = de;
                    center_solution.has_de_dtheta_prime = true;
                }
                continue;
            }

            if (node_index == 0 && has_right_neighbor && right_match_index.has_value()) {
                const auto& right_solution = nodes[1].solutions[*right_match_index];
                const double dphi = estimate_theta_prime_derivative_forward(
                    center_node.theta_prime_local,
                    center_solution.encounter_global_angle,
                    nodes[1].theta_prime_local,
                    right_solution.encounter_global_angle);
                const double de = estimate_theta_prime_derivative_forward(
                    center_node.theta_prime_local,
                    center_solution.outgoing_eccentricity,
                    nodes[1].theta_prime_local,
                    right_solution.outgoing_eccentricity);
                if (is_finite(dphi)) {
                    center_solution.dphi_dtheta_prime = dphi;
                    center_solution.has_dphi_dtheta_prime = true;
                }
                if (is_finite(de)) {
                    center_solution.de_dtheta_prime = de;
                    center_solution.has_de_dtheta_prime = true;
                }
                continue;
            }

            if (node_index + 1 == node_count && has_left_neighbor && left_match_index.has_value()) {
                const auto& left_solution = nodes[node_index - 1].solutions[*left_match_index];
                const double dphi = estimate_theta_prime_derivative_backward(
                    nodes[node_index - 1].theta_prime_local,
                    left_solution.encounter_global_angle,
                    center_node.theta_prime_local,
                    center_solution.encounter_global_angle);
                const double de = estimate_theta_prime_derivative_backward(
                    nodes[node_index - 1].theta_prime_local,
                    left_solution.outgoing_eccentricity,
                    center_node.theta_prime_local,
                    center_solution.outgoing_eccentricity);
                if (is_finite(dphi)) {
                    center_solution.dphi_dtheta_prime = dphi;
                    center_solution.has_dphi_dtheta_prime = true;
                }
                if (is_finite(de)) {
                    center_solution.de_dtheta_prime = de;
                    center_solution.has_de_dtheta_prime = true;
                }
            }
        }
    }
}

void attach_problem2_theta_prime_solution_derivatives_by_k(
    std::vector<Problem2ThetaPrimeNodeSnapshot>& nodes,
    int max_transfer_revolution,
    double branch_phi_pairing_max_gap
) {
    if (nodes.empty() || max_transfer_revolution < 0 || !(branch_phi_pairing_max_gap > 0.0)) {
        return;
    }

    for (int transfer_revolution = 0; transfer_revolution <= max_transfer_revolution; ++transfer_revolution) {
        std::vector<Problem2ThetaPrimeNodeSnapshot> layer_nodes;
        layer_nodes.reserve(nodes.size());
        for (const auto& node : nodes) {
            Problem2ThetaPrimeNodeSnapshot layer_node{};
            layer_node.theta_prime_local = node.theta_prime_local;
            layer_node.theta_prime_global = node.theta_prime_global;
            if (transfer_revolution < static_cast<int>(node.solutions_by_k.size())) {
                layer_node.solutions = node.solutions_by_k[static_cast<std::size_t>(transfer_revolution)];
            }
            layer_nodes.push_back(std::move(layer_node));
        }

        attach_problem2_theta_prime_solution_derivatives(layer_nodes, branch_phi_pairing_max_gap);

        for (std::size_t node_index = 0; node_index < nodes.size(); ++node_index) {
            auto& node = nodes[node_index];
            if (transfer_revolution >= static_cast<int>(node.solutions_by_k.size())) {
                continue;
            }
            auto& layer_solutions = node.solutions_by_k[static_cast<std::size_t>(transfer_revolution)];
            const auto& derived_layer = layer_nodes[node_index].solutions;
            for (auto& solution : layer_solutions) {
                for (const auto& derived_solution : derived_layer) {
                    if (same_outgoing_branch_identity(solution, derived_solution)) {
                        copy_solution_derivatives(derived_solution, solution);
                        break;
                    }
                }
            }
        }
    }

    for (auto& node : nodes) {
        merge_flat_solutions_from_k_layers(node);
    }
}

}  // namespace spaceship_cpp::problem2
