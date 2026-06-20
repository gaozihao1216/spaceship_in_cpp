/*
 * 文件作用：测试 Problem 2 θ' 初扫与 dφ/dθ'、de/dθ' 导数估计。
 * 主要工作：用细步长数值参考验证粗网格差分在 branch 配对后的误差界。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

namespace {

using spaceship_cpp::config::Problem2ThetaPrimeScanDefaults;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;
using spaceship_cpp::problem2::Problem2ThetaPrimeScanConfig;

constexpr double kDayInSeconds = 86400.0;

std::optional<Problem2OutgoingBranchSolution> find_solution_with_derivatives(
    const spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult& scan,
    std::size_t& out_node_index,
    std::size_t& out_solution_index
) {
    for (std::size_t node_index = 1; node_index + 1 < scan.nodes.size(); ++node_index) {
        const auto& node = scan.nodes[node_index];
        for (std::size_t solution_index = 0; solution_index < node.solutions.size(); ++solution_index) {
            const auto& solution = node.solutions[solution_index];
            if (solution.has_dphi_dtheta_prime && solution.has_de_dtheta_prime) {
                out_node_index = node_index;
                out_solution_index = solution_index;
                return solution;
            }
        }
    }
    return std::nullopt;
}

std::optional<Problem2OutgoingBranchSolution> matched_solution_at_theta_prime(
    const Problem2ThetaPrimeScanConfig& config,
    double theta_prime_local,
    const Problem2OutgoingBranchSolution& reference,
    double max_phi_gap
) {
    const auto snapshot =
        spaceship_cpp::problem2::evaluate_problem2_theta_prime_node(config, theta_prime_local);
    const auto match_index = spaceship_cpp::problem2::find_best_matching_outgoing_branch_index(
        reference,
        snapshot.solutions,
        max_phi_gap);
    if (!match_index.has_value()) {
        return std::nullopt;
    }
    return snapshot.solutions[*match_index];
}

double relative_error(double estimated, double reference) {
    return std::abs(estimated - reference) / std::max(std::abs(reference), 1.0);
}

void verify_branch_derivatives_at_theta_prime(
    const Problem2ThetaPrimeScanConfig& config,
    double theta_prime_local,
    const Problem2OutgoingBranchSolution& center_solution,
    double coarse_step,
    double reference_step
) {
    const auto minus_solution = matched_solution_at_theta_prime(
        config,
        theta_prime_local - reference_step,
        center_solution,
        config.branch_phi_pairing_max_gap);
    const auto plus_solution = matched_solution_at_theta_prime(
        config,
        theta_prime_local + reference_step,
        center_solution,
        config.branch_phi_pairing_max_gap);
    const auto minus_coarse = matched_solution_at_theta_prime(
        config,
        theta_prime_local - coarse_step,
        center_solution,
        config.branch_phi_pairing_max_gap);
    const auto plus_coarse = matched_solution_at_theta_prime(
        config,
        theta_prime_local + coarse_step,
        center_solution,
        config.branch_phi_pairing_max_gap);
    assert(minus_solution.has_value());
    assert(plus_solution.has_value());
    assert(minus_coarse.has_value());
    assert(plus_coarse.has_value());

    const double phi_span = std::max(
        std::abs(plus_coarse->encounter_global_angle - minus_coarse->encounter_global_angle),
        std::abs(plus_coarse->outgoing_eccentricity - minus_coarse->outgoing_eccentricity));
    assert(phi_span > 1e-14);

    const double reference_dphi = spaceship_cpp::problem2::estimate_theta_prime_derivative_central(
        theta_prime_local - reference_step,
        minus_solution->encounter_global_angle,
        theta_prime_local + reference_step,
        plus_solution->encounter_global_angle);
    const double reference_de = spaceship_cpp::problem2::estimate_theta_prime_derivative_central(
        theta_prime_local - reference_step,
        minus_solution->outgoing_eccentricity,
        theta_prime_local + reference_step,
        plus_solution->outgoing_eccentricity);
    assert(std::isfinite(reference_dphi));
    assert(std::isfinite(reference_de));

    const double coarse_dphi = spaceship_cpp::problem2::estimate_theta_prime_derivative_central(
        theta_prime_local - coarse_step,
        minus_coarse->encounter_global_angle,
        theta_prime_local + coarse_step,
        plus_coarse->encounter_global_angle);
    const double coarse_de = spaceship_cpp::problem2::estimate_theta_prime_derivative_central(
        theta_prime_local - coarse_step,
        minus_coarse->outgoing_eccentricity,
        theta_prime_local + coarse_step,
        plus_coarse->outgoing_eccentricity);
    assert(std::isfinite(coarse_dphi));
    assert(std::isfinite(coarse_de));

    const double half_step = 0.5 * coarse_step;
    const auto minus_half = matched_solution_at_theta_prime(
        config,
        theta_prime_local - half_step,
        center_solution,
        config.branch_phi_pairing_max_gap);
    const auto plus_half = matched_solution_at_theta_prime(
        config,
        theta_prime_local + half_step,
        center_solution,
        config.branch_phi_pairing_max_gap);
    assert(minus_half.has_value());
    assert(plus_half.has_value());

    const double half_step_dphi = spaceship_cpp::problem2::estimate_theta_prime_derivative_central(
        theta_prime_local - half_step,
        minus_half->encounter_global_angle,
        theta_prime_local + half_step,
        plus_half->encounter_global_angle);
    const double half_step_de = spaceship_cpp::problem2::estimate_theta_prime_derivative_central(
        theta_prime_local - half_step,
        minus_half->outgoing_eccentricity,
        theta_prime_local + half_step,
        plus_half->outgoing_eccentricity);

    const double richardson_phi_bound =
        std::abs(coarse_dphi - half_step_dphi) / 3.0 + 1e-9;
    const double richardson_e_bound =
        std::abs(coarse_de - half_step_de) / 3.0 + 1e-9;
    double allowed_phi_error = std::max(
        richardson_phi_bound * 4.0,
        0.01 * std::max(std::abs(reference_dphi), 1.0));
    double allowed_e_error = std::max(
        richardson_e_bound * 4.0,
        0.01 * std::max(std::abs(reference_de), 1.0));

    const double center_phi = center_solution.encounter_global_angle;
    const double second_dphi = (
        plus_solution->encounter_global_angle -
        2.0 * center_phi +
        minus_solution->encounter_global_angle) /
        (reference_step * reference_step);
    const double second_de = (
        plus_solution->outgoing_eccentricity -
        2.0 * center_solution.outgoing_eccentricity +
        minus_solution->outgoing_eccentricity) /
        (reference_step * reference_step);
    if (std::isfinite(second_dphi)) {
        allowed_phi_error = std::max(
            allowed_phi_error,
            0.5 * std::abs(second_dphi) * coarse_step * coarse_step + 1e-9);
    }
    if (std::isfinite(second_de)) {
        allowed_e_error = std::max(
            allowed_e_error,
            0.5 * std::abs(second_de) * coarse_step * coarse_step + 1e-9);
    }

    assert(std::abs(center_solution.dphi_dtheta_prime - reference_dphi) <= allowed_phi_error);
    assert(std::abs(center_solution.de_dtheta_prime - reference_de) <= allowed_e_error);
    assert(relative_error(center_solution.dphi_dtheta_prime, coarse_dphi) <= 1e-12);
    assert(relative_error(center_solution.de_dtheta_prime, coarse_de) <= 1e-12);
}

Problem2ThetaPrimeScanConfig make_test_scan_config(
    PlanetId flyby_planet,
    PlanetId target_planet,
    double flyby_time_seconds,
    int theta_prime_count
) {
    const auto& defaults = spaceship_cpp::config::global_config();
    const Problem2ThetaPrimeScanDefaults scan_defaults{
        theta_prime_count,
        defaults.problem2_theta_prime_scan.branch_phi_pairing_max_gap,
    };
    return spaceship_cpp::config::make_problem2_theta_prime_scan_config(
        flyby_planet,
        target_planet,
        flyby_time_seconds,
        scan_defaults,
        defaults.problem1_solve);
}

bool try_verify_derivatives_for_scenario(
    PlanetId flyby_planet,
    PlanetId target_planet,
    double flyby_time_seconds
) {
    const auto config = make_test_scan_config(
        flyby_planet,
        target_planet,
        flyby_time_seconds,
        32);
    const auto scan = spaceship_cpp::problem2::run_problem2_theta_prime_initial_scan(config);
    if (!scan.ok || scan.nodes.size() < 3U) {
        return false;
    }

    std::size_t node_index = 0;
    std::size_t solution_index = 0;
    const auto center_solution = find_solution_with_derivatives(scan, node_index, solution_index);
    if (!center_solution.has_value()) {
        return false;
    }

    const double coarse_step =
        scan.nodes[node_index + 1].theta_prime_local - scan.nodes[node_index].theta_prime_local;
    const double reference_step = coarse_step / 40.0;
    verify_branch_derivatives_at_theta_prime(
        config,
        scan.nodes[node_index].theta_prime_local,
        *center_solution,
        coarse_step,
        reference_step);
    return true;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace problem2 = spaceship_cpp::problem2;

    {
        const auto grid = problem2::build_uniform_theta_prime_local_grid(64);
        assert(grid.size() == 64U);
        assert(common::nearly_equal(grid.front(), -common::kPi, 1e-12));
        const double expected_step = common::kTwoPi / 64.0;
        for (std::size_t index = 1; index < grid.size(); ++index) {
            assert(common::nearly_equal(grid[index] - grid[index - 1], expected_step, 1e-12));
        }
    }

    {
        const auto config = make_test_scan_config(PlanetId::Earth, PlanetId::Mars, 0.0, 64);
        const auto scan = problem2::run_problem2_theta_prime_initial_scan(config);
        assert(scan.ok);
        assert(scan.nodes.size() == 64U);

        std::size_t total_solutions = 0;
        for (const auto& node : scan.nodes) {
            for (const auto& solution : node.solutions) {
                assert(std::isfinite(solution.encounter_global_angle));
                assert(std::isfinite(solution.outgoing_eccentricity));
                assert(std::isfinite(solution.outgoing_semi_latus_rectum));
                ++total_solutions;
            }
        }
        assert(total_solutions > 0U);
    }

    {
        bool verified = false;
        const std::pair<PlanetId, PlanetId> planet_pairs[] = {
            {PlanetId::Earth, PlanetId::Mars},
            {PlanetId::Mars, PlanetId::Earth},
            {PlanetId::Venus, PlanetId::Earth},
        };
        const double flyby_times[] = {0.0, 100.0 * kDayInSeconds};
        for (const auto& [flyby_planet, target_planet] : planet_pairs) {
            for (const double flyby_time : flyby_times) {
                if (try_verify_derivatives_for_scenario(flyby_planet, target_planet, flyby_time)) {
                    verified = true;
                    std::cout << "derivative_verified_for_flyby_scenario\n";
                    break;
                }
            }
            if (verified) {
                break;
            }
        }
        assert(verified);
    }

    {
        const auto config = make_test_scan_config(PlanetId::Earth, PlanetId::Mars, 0.0, 32);
        auto scan = problem2::run_problem2_theta_prime_initial_scan(config);
        assert(scan.ok);

        std::size_t derivative_count = 0;
        for (const auto& node : scan.nodes) {
            for (const auto& solution : node.solutions) {
                if (solution.has_dphi_dtheta_prime && solution.has_de_dtheta_prime) {
                    assert(std::isfinite(solution.dphi_dtheta_prime));
                    assert(std::isfinite(solution.de_dtheta_prime));
                    ++derivative_count;
                }
            }
        }
        assert(derivative_count > 0U);
    }

    std::cout << "problem2_theta_prime_derivative_ok\n";
    return 0;
}
