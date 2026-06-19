/*
 * 文件作用：实现 Problem 1 的直接残差评估和求解器。
 * 主要工作：构造转移轨道残差，扫描遇合角根区间，并用二分法细化候选解。
 */
#include "spaceship_cpp/problem1/problem1.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace spaceship_cpp::problem1 {

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::near_zero;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::common::orbit_F;

Problem1ResidualResult make_result(Problem1ResidualStatus status) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return Problem1ResidualResult{
        status,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
        nan,
    };
}

bool positive_orbit_denominator(double e, double angle) {
    const double denominator = 1.0 + e * std::cos(angle);
    return is_finite(denominator) && denominator > kDefaultEpsilon;
}

void validate_problem1_solve_input(const Problem1SolveInput& input) {
    if (!is_finite(input.launch_time_seconds_since_j2000) ||
        !is_finite(input.transfer_perihelion_angle) ||
        input.max_transfer_revolution < 0 ||
        input.max_target_revolution < 0 ||
        input.phi_scan_count < 3 ||
        !is_finite(input.phi_tolerance) ||
        !(input.phi_tolerance > 0.0) ||
        !is_finite(input.residual_tolerance) ||
        input.residual_tolerance < 0.0 ||
        input.max_bisection_iterations <= 0 ||
        !is_finite(input.max_candidate_relative_residual) ||
        !(input.max_candidate_relative_residual > 0.0)) {
        throw std::invalid_argument("invalid Problem1SolveInput");
    }
}

Problem1ResidualInput make_residual_input(
    const Problem1SolveInput& input,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    return Problem1ResidualInput{
        input.departure_planet,
        input.target_planet,
        input.launch_time_seconds_since_j2000,
        input.transfer_perihelion_angle,
        encounter_global_angle,
        transfer_revolution,
        target_revolution,
    };
}

bool is_success_with_finite_residual(const Problem1ResidualResult& result) {
    return result.status == Problem1ResidualStatus::Success && is_finite(result.residual);
}

bool residual_sign_changed(double a, double b) {
    if (!is_finite(a) || !is_finite(b)) {
        return false;
    }
    if (a == 0.0 || b == 0.0) {
        return true;
    }
    return (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0);
}

double compute_problem1_residual_scale(const Problem1ResidualResult& result) {
    if (result.status != Problem1ResidualStatus::Success ||
        !is_finite(result.transfer_time_scale_free) ||
        !is_finite(result.target_time_scale_free)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::max({std::abs(result.transfer_time_scale_free), std::abs(result.target_time_scale_free), 1.0});
}

double compute_problem1_relative_residual(const Problem1ResidualResult& result) {
    const double scale = compute_problem1_residual_scale(result);
    if (!is_finite(scale) || !(scale > 0.0) || !is_finite(result.residual)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::abs(result.residual) / scale;
}

std::optional<Problem1Candidate> make_candidate_from_residual_result(
    const Problem1SolveInput& input,
    double encounter_global_angle,
    const Problem1ResidualResult& residual_result,
    int transfer_revolution,
    int target_revolution,
    double root_bracket_width,
    int bisection_iterations,
    bool refined_by_bisection
) {
    if (!is_success_with_finite_residual(residual_result)) {
        return std::nullopt;
    }

    const double mu = planet_params::get_solar_system_physical_params().GM_sun;
    const double time_of_flight_seconds = residual_result.transfer_time_scale_free / std::sqrt(mu);
    if (!is_finite(time_of_flight_seconds) || !(time_of_flight_seconds > 0.0)) {
        return std::nullopt;
    }

    const double residual_scale = compute_problem1_residual_scale(residual_result);
    const double relative_residual = compute_problem1_relative_residual(residual_result);
    if (!is_finite(residual_scale) || !(residual_scale > 0.0) ||
        !is_finite(relative_residual) || relative_residual < 0.0 ||
        !is_finite(root_bracket_width) || root_bracket_width < 0.0 ||
        bisection_iterations < 0) {
        return std::nullopt;
    }
    if (relative_residual > input.max_candidate_relative_residual) {
        return std::nullopt;
    }

    return Problem1Candidate{
        residual_result,
        normalize_angle_0_2pi(encounter_global_angle),
        input.launch_time_seconds_since_j2000,
        time_of_flight_seconds,
        input.launch_time_seconds_since_j2000 + time_of_flight_seconds,
        residual_scale,
        relative_residual,
        root_bracket_width,
        bisection_iterations,
        refined_by_bisection,
        transfer_revolution,
        target_revolution,
    };
}

std::optional<Problem1Candidate> refine_problem1_root_by_bisection(
    const Problem1SolveInput& input,
    double left_phi,
    const Problem1ResidualResult& left_result,
    double right_phi,
    const Problem1ResidualResult& right_result,
    int transfer_revolution,
    int target_revolution
) {
    if (!is_success_with_finite_residual(left_result) || !is_success_with_finite_residual(right_result)) {
        return std::nullopt;
    }

    if (std::abs(left_result.residual) <= input.residual_tolerance) {
        return make_candidate_from_residual_result(
            input,
            left_phi,
            left_result,
            transfer_revolution,
            target_revolution,
            std::abs(right_phi - left_phi),
            0,
            true);
    }
    if (std::abs(right_result.residual) <= input.residual_tolerance) {
        return make_candidate_from_residual_result(
            input,
            right_phi,
            right_result,
            transfer_revolution,
            target_revolution,
            std::abs(right_phi - left_phi),
            0,
            true);
    }
    if (!residual_sign_changed(left_result.residual, right_result.residual)) {
        return std::nullopt;
    }

    double current_left_phi = left_phi;
    double current_right_phi = right_phi;
    Problem1ResidualResult current_left_result = left_result;
    Problem1ResidualResult current_right_result = right_result;
    double mid_phi = 0.5 * (current_left_phi + current_right_phi);
    Problem1ResidualResult mid_result = make_result(Problem1ResidualStatus::InvalidBranch);
    int iterations_used = 0;

    for (int iteration = 0; iteration < input.max_bisection_iterations; ++iteration) {
        iterations_used = iteration + 1;
        mid_phi = 0.5 * (current_left_phi + current_right_phi);
        mid_result = evaluate_problem1_residual(
            make_residual_input(input, mid_phi, transfer_revolution, target_revolution));
        if (!is_success_with_finite_residual(mid_result)) {
            return std::nullopt;
        }
        if (std::abs(mid_result.residual) <= input.residual_tolerance ||
            std::abs(current_right_phi - current_left_phi) <= input.phi_tolerance) {
            return make_candidate_from_residual_result(
                input,
                mid_phi,
                mid_result,
                transfer_revolution,
                target_revolution,
                std::abs(current_right_phi - current_left_phi),
                iterations_used,
                true);
        }

        if (residual_sign_changed(current_left_result.residual, mid_result.residual)) {
            current_right_phi = mid_phi;
            current_right_result = mid_result;
        } else if (residual_sign_changed(mid_result.residual, current_right_result.residual)) {
            current_left_phi = mid_phi;
            current_left_result = mid_result;
        } else {
            return std::nullopt;
        }
    }

    return make_candidate_from_residual_result(
        input,
        mid_phi,
        mid_result,
        transfer_revolution,
        target_revolution,
        std::abs(current_right_phi - current_left_phi),
        iterations_used,
        true);
}

void add_candidate_if_not_duplicate(std::vector<Problem1Candidate>* candidates, Problem1Candidate candidate) {
    for (Problem1Candidate& existing : *candidates) {
        if (existing.transfer_revolution != candidate.transfer_revolution ||
            existing.target_revolution != candidate.target_revolution) {
            continue;
        }
        const double angle_diff =
            normalize_angle_minus_pi_pi(existing.encounter_global_angle - candidate.encounter_global_angle);
        if (std::abs(angle_diff) < 1e-8) {
            if (candidate.relative_residual < existing.relative_residual ||
                (candidate.relative_residual == existing.relative_residual &&
                 (candidate.root_bracket_width < existing.root_bracket_width ||
                  (candidate.root_bracket_width == existing.root_bracket_width &&
                   std::abs(candidate.residual_result.residual) < std::abs(existing.residual_result.residual))))) {
                existing = candidate;
            }
            return;
        }
    }
    candidates->push_back(candidate);
}

}  // namespace

double compute_transfer_e_from_two_points(double r1, double xi1, double r2, double xi2) {
    if (!is_finite(r1) || !is_finite(r2) || !is_finite(xi1) || !is_finite(xi2)) {
        throw std::domain_error("compute_transfer_e_from_two_points requires finite inputs");
    }
    if (!(r1 > 0.0) || !(r2 > 0.0)) {
        throw std::domain_error("compute_transfer_e_from_two_points requires positive radii");
    }

    const double denominator = r1 * std::cos(xi1) - r2 * std::cos(xi2);
    if (near_zero(denominator, kDefaultEpsilon)) {
        throw std::domain_error("compute_transfer_e_from_two_points has singular denominator");
    }

    return (r2 - r1) / denominator;
}

double compute_transfer_p_from_departure(double r1, double e_transfer, double xi1) {
    if (!is_finite(r1) || !is_finite(e_transfer) || !is_finite(xi1)) {
        throw std::domain_error("compute_transfer_p_from_departure requires finite inputs");
    }
    if (!(r1 > 0.0)) {
        throw std::domain_error("compute_transfer_p_from_departure requires positive radius");
    }

    const double factor = 1.0 + e_transfer * std::cos(xi1);
    if (!(is_finite(factor) && factor > kDefaultEpsilon)) {
        throw std::domain_error("compute_transfer_p_from_departure requires positive orbit factor");
    }
    return r1 * factor;
}

Problem1ResidualResult evaluate_problem1_residual(const Problem1ResidualInput& input) {
    if (!is_finite(input.launch_time_seconds_since_j2000) ||
        !is_finite(input.transfer_perihelion_angle) ||
        !is_finite(input.encounter_global_angle) ||
        input.transfer_revolution < 0 ||
        input.target_revolution < 0) {
        return make_result(Problem1ResidualStatus::InvalidInput);
    }

    const planet_params::PlanetParams& target = planet_params::get_planet_params(input.target_planet);

    const planet_params::PlanetState state1 =
        planet_params::planet_state_at_time(input.departure_planet, input.launch_time_seconds_since_j2000);
    const double r1 = state1.radius;
    const double lambda1 = state1.theta_global;
    const double theta2_start =
        planet_params::planet_true_anomaly_at_time(input.target_planet, input.launch_time_seconds_since_j2000);

    const double u2_base = normalize_angle_0_2pi(input.encounter_global_angle - target.orbit.theta_0);
    const double r2 = planet_params::planet_radius_at_true_anomaly(input.target_planet, u2_base);

    const double raw_phi0 = input.transfer_perihelion_angle;
    const double xi1_raw_base = normalize_angle_0_2pi(lambda1 - raw_phi0);
    const double xi2_raw_base = normalize_angle_0_2pi(input.encounter_global_angle - raw_phi0);

    double e_raw = std::numeric_limits<double>::quiet_NaN();
    double p_transfer = std::numeric_limits<double>::quiet_NaN();
    try {
        e_raw = compute_transfer_e_from_two_points(r1, xi1_raw_base, r2, xi2_raw_base);
        p_transfer = compute_transfer_p_from_departure(r1, e_raw, xi1_raw_base);
    } catch (const std::domain_error&) {
        return make_result(Problem1ResidualStatus::SingularGeometry);
    }

    double e_transfer = e_raw;
    double phi0_used = raw_phi0;
    if (e_transfer < 0.0) {
        e_transfer = -e_transfer;
        phi0_used = raw_phi0 + spaceship_cpp::common::kPi;
    }
    phi0_used = normalize_angle_0_2pi(phi0_used);

    double xi1_geometry = xi1_raw_base;
    double xi2_geometry = xi2_raw_base;
    if (e_raw < 0.0) {
        xi1_geometry = normalize_angle_0_2pi(lambda1 - phi0_used);
        xi2_geometry = normalize_angle_0_2pi(input.encounter_global_angle - phi0_used);
    }

    if (!is_finite(e_transfer) || !(e_transfer >= 0.0) ||
        !is_finite(p_transfer) || !(p_transfer > 0.0) ||
        !positive_orbit_denominator(e_transfer, xi1_geometry) ||
        !positive_orbit_denominator(e_transfer, xi2_geometry)) {
        return make_result(Problem1ResidualStatus::InvalidTransferOrbit);
    }

    double theta2_end = u2_base;
    while (theta2_end <= theta2_start) {
        theta2_end += kTwoPi;
    }
    theta2_end += static_cast<double>(input.target_revolution) * kTwoPi;

    double xi1 = std::numeric_limits<double>::quiet_NaN();
    double xi2 = std::numeric_limits<double>::quiet_NaN();
    if (e_transfer < 1.0) {
        xi1 = xi1_geometry;
        xi2 = xi2_geometry;
        while (xi2 <= xi1) {
            xi2 += kTwoPi;
        }
        xi2 += static_cast<double>(input.transfer_revolution) * kTwoPi;
    } else {
        if (input.transfer_revolution != 0) {
            return make_result(Problem1ResidualStatus::InvalidBranch);
        }
        xi1 = normalize_angle_minus_pi_pi(lambda1 - phi0_used);
        xi2 = normalize_angle_minus_pi_pi(input.encounter_global_angle - phi0_used);
        if (xi2 <= xi1) {
            return make_result(Problem1ResidualStatus::InvalidBranch);
        }
    }

    double deltaF_transfer = std::numeric_limits<double>::quiet_NaN();
    double deltaF_target = std::numeric_limits<double>::quiet_NaN();
    try {
        deltaF_transfer = orbit_F(e_transfer, xi2) - orbit_F(e_transfer, xi1);
        deltaF_target = orbit_F(target.orbit.e, theta2_end) - orbit_F(target.orbit.e, theta2_start);
    } catch (const std::domain_error&) {
        return make_result(Problem1ResidualStatus::InvalidBranch);
    }

    if (!is_finite(deltaF_transfer) || !is_finite(deltaF_target) ||
        !(deltaF_transfer > 0.0) || !(deltaF_target > 0.0)) {
        return make_result(Problem1ResidualStatus::InvalidTimeOfFlight);
    }

    const double transfer_time_scale_free = std::pow(p_transfer, 1.5) * deltaF_transfer;
    const double target_time_scale_free = std::pow(target.orbit.p, 1.5) * deltaF_target;
    if (!is_finite(transfer_time_scale_free) || !is_finite(target_time_scale_free)) {
        return make_result(Problem1ResidualStatus::InvalidTimeOfFlight);
    }

    return Problem1ResidualResult{
        Problem1ResidualStatus::Success,
        transfer_time_scale_free - target_time_scale_free,
        r1,
        r2,
        xi1,
        xi2,
        theta2_start,
        theta2_end,
        e_raw,
        e_transfer,
        phi0_used,
        p_transfer,
        deltaF_transfer,
        deltaF_target,
        transfer_time_scale_free,
        target_time_scale_free,
    };
}

std::vector<Problem1Candidate> solve_problem1(const Problem1SolveInput& input) {
    validate_problem1_solve_input(input);

    std::vector<Problem1Candidate> candidates;
    for (int transfer_revolution = 0; transfer_revolution <= input.max_transfer_revolution; ++transfer_revolution) {
        for (int target_revolution = 0; target_revolution <= input.max_target_revolution; ++target_revolution) {
            bool has_previous_valid = false;
            bool has_first_valid = false;
            double previous_phi = 0.0;
            Problem1ResidualResult previous_result = make_result(Problem1ResidualStatus::InvalidBranch);
            double first_phi = 0.0;
            Problem1ResidualResult first_result = make_result(Problem1ResidualStatus::InvalidBranch);
            double last_valid_phi = 0.0;
            Problem1ResidualResult last_valid_result = make_result(Problem1ResidualStatus::InvalidBranch);

            for (int i = 0; i < input.phi_scan_count; ++i) {
                const double phi = static_cast<double>(i) * kTwoPi / static_cast<double>(input.phi_scan_count);
                const Problem1ResidualResult result =
                    evaluate_problem1_residual(make_residual_input(input, phi, transfer_revolution, target_revolution));

                if (!is_success_with_finite_residual(result)) {
                    has_previous_valid = false;
                    continue;
                }

                if (!has_first_valid) {
                    has_first_valid = true;
                    first_phi = phi;
                    first_result = result;
                }

                if (std::abs(result.residual) <= input.residual_tolerance) {
                    auto candidate = make_candidate_from_residual_result(
                        input,
                        phi,
                        result,
                        transfer_revolution,
                        target_revolution,
                        0.0,
                        0,
                        false);
                    if (candidate.has_value()) {
                        add_candidate_if_not_duplicate(&candidates, *candidate);
                    }
                }

                if (has_previous_valid && residual_sign_changed(previous_result.residual, result.residual)) {
                    auto refined = refine_problem1_root_by_bisection(
                        input,
                        previous_phi,
                        previous_result,
                        phi,
                        result,
                        transfer_revolution,
                        target_revolution);
                    if (refined.has_value()) {
                        add_candidate_if_not_duplicate(&candidates, *refined);
                    }
                }

                previous_phi = phi;
                previous_result = result;
                last_valid_phi = phi;
                last_valid_result = result;
                has_previous_valid = true;
            }

            if (has_first_valid && has_previous_valid &&
                residual_sign_changed(last_valid_result.residual, first_result.residual)) {
                auto refined = refine_problem1_root_by_bisection(
                    input,
                    last_valid_phi,
                    last_valid_result,
                    first_phi + kTwoPi,
                    first_result,
                    transfer_revolution,
                    target_revolution);
                if (refined.has_value()) {
                    add_candidate_if_not_duplicate(&candidates, *refined);
                }
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const Problem1Candidate& a, const Problem1Candidate& b) {
        if (a.arrival_time_seconds_since_j2000 < b.arrival_time_seconds_since_j2000) {
            return true;
        }
        if (a.arrival_time_seconds_since_j2000 > b.arrival_time_seconds_since_j2000) {
            return false;
        }
        if (a.relative_residual < b.relative_residual) {
            return true;
        }
        if (a.relative_residual > b.relative_residual) {
            return false;
        }
        if (a.root_bracket_width < b.root_bracket_width) {
            return true;
        }
        if (a.root_bracket_width > b.root_bracket_width) {
            return false;
        }
        return std::abs(a.residual_result.residual) < std::abs(b.residual_result.residual);
    });

    return candidates;
}

}  // namespace spaceship_cpp::problem1
