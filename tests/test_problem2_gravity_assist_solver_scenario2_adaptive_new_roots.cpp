#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem2/problem2_gravity_assist_solver.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kEncounterTimeFactor = 0.25;
constexpr double kIncomingE = 0.4;
constexpr double kIncomingTheta = 1.0;
constexpr double kGeometryDenominatorEpsilon = 1e-12;
constexpr double kBoundaryRoundoffAbsToleranceM = 1e-3;
constexpr double kBoundaryRoundoffRelTolerance = 1e-12;

struct GeometryCheck {
    bool valid = false;
    std::string invalid_reason;
    double r0 = std::numeric_limits<double>::quiet_NaN();
    double r1 = std::numeric_limits<double>::quiet_NaN();
    double denominator = std::numeric_limits<double>::quiet_NaN();
    double e_prime = std::numeric_limits<double>::quiet_NaN();
    double p_prime = std::numeric_limits<double>::quiet_NaN();
    double encounter_factor = std::numeric_limits<double>::quiet_NaN();
    double target_factor = std::numeric_limits<double>::quiet_NaN();
    double encounter_radius_error = std::numeric_limits<double>::quiet_NaN();
    double target_radius_error = std::numeric_limits<double>::quiet_NaN();
};

struct DirectValidation {
    bool matched = false;
    bool validated = false;
    int matched_k = 0;
    int matched_q = 0;
    double alpha_diff = std::numeric_limits<double>::quiet_NaN();
    double time_diff_seconds = std::numeric_limits<double>::quiet_NaN();
    double direct_problem1_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double direct_slingshot_residual = std::numeric_limits<double>::quiet_NaN();
    std::string direct_residual_source;
};

struct RootSetComparison {
    int matched_root_count = 0;
    int new_root_in_b_count = 0;
    int missing_root_from_a_count = 0;
    double max_theta_diff_matched = 0.0;
    double max_alpha_diff_matched = 0.0;
    std::vector<int> unmatched_b_indices;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

std::string residual_source_name(spaceship_cpp::problem2::Problem2ResidualSource source) {
    if (source == spaceship_cpp::problem2::Problem2ResidualSource::Strict) {
        return "Strict";
    }
    return "BoundaryAmbiguousRoundoff";
}

double wrapped_angle_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

bool same_root(
    const spaceship_cpp::problem2::Problem2GravityAssistSolution& lhs,
    const spaceship_cpp::problem2::Problem2GravityAssistSolution& rhs,
    double* theta_diff = nullptr,
    double* alpha_diff = nullptr
) {
    const double local_theta_diff = std::abs(lhs.theta_prime - rhs.theta_prime);
    const double local_alpha_diff = wrapped_angle_distance(lhs.alpha, rhs.alpha);
    if (theta_diff != nullptr) {
        *theta_diff = local_theta_diff;
    }
    if (alpha_diff != nullptr) {
        *alpha_diff = local_alpha_diff;
    }
    return lhs.transfer_revolution == rhs.transfer_revolution &&
           lhs.target_revolution == rhs.target_revolution &&
           local_theta_diff <= 1e-6 &&
           local_alpha_diff <= 1e-6 &&
           std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds) <= 1e3;
}

RootSetComparison compare_root_sets(
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& roots_a,
    const std::vector<spaceship_cpp::problem2::Problem2GravityAssistSolution>& roots_b
) {
    RootSetComparison comparison{};
    std::vector<bool> used_b(roots_b.size(), false);
    std::vector<bool> matched_a(roots_a.size(), false);

    for (std::size_t a = 0; a < roots_a.size(); ++a) {
        for (std::size_t b = 0; b < roots_b.size(); ++b) {
            double theta_diff = 0.0;
            double alpha_diff = 0.0;
            if (!used_b[b] && same_root(roots_a[a], roots_b[b], &theta_diff, &alpha_diff)) {
                used_b[b] = true;
                matched_a[a] = true;
                comparison.matched_root_count += 1;
                comparison.max_theta_diff_matched = std::max(comparison.max_theta_diff_matched, theta_diff);
                comparison.max_alpha_diff_matched = std::max(comparison.max_alpha_diff_matched, alpha_diff);
                break;
            }
        }
    }

    for (std::size_t b = 0; b < roots_b.size(); ++b) {
        if (!used_b[b]) {
            comparison.unmatched_b_indices.push_back(static_cast<int>(b));
        }
    }
    comparison.new_root_in_b_count = static_cast<int>(comparison.unmatched_b_indices.size());
    comparison.missing_root_from_a_count =
        static_cast<int>(std::count(matched_a.begin(), matched_a.end(), false));
    return comparison;
}

GeometryCheck compute_geometry(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
) {
    namespace problem2 = spaceship_cpp::problem2;
    GeometryCheck g{};
    g.r0 = problem2::problem2_orbit_radius(R_J, e_J, phi);
    if (!is_finite(g.r0) || !(g.r0 > 0.0)) {
        g.invalid_reason = "invalid_encounter_radius";
        return g;
    }
    g.r1 = problem2::problem2_orbit_radius(R_K, e_K, alpha);
    if (!is_finite(g.r1) || !(g.r1 > 0.0)) {
        g.invalid_reason = "invalid_target_radius";
        return g;
    }
    const double cos_encounter = std::cos(phi - theta_prime);
    const double cos_target = std::cos(alpha - theta_prime);
    g.denominator = g.r0 * cos_encounter - g.r1 * cos_target;
    if (!is_finite(g.denominator) || std::abs(g.denominator) <= kGeometryDenominatorEpsilon) {
        g.invalid_reason = "geometry_denominator_too_small";
        return g;
    }
    g.e_prime = (g.r1 - g.r0) / g.denominator;
    if (!is_finite(g.e_prime)) {
        g.invalid_reason = "non_finite_outgoing_eccentricity";
        return g;
    }
    g.encounter_factor = 1.0 + g.e_prime * cos_encounter;
    g.p_prime = g.r0 * g.encounter_factor;
    if (!is_finite(g.p_prime) || !(g.p_prime > 0.0)) {
        g.invalid_reason = "non_positive_outgoing_semi_latus_rectum";
        return g;
    }
    if (!is_finite(g.encounter_factor) || !(g.encounter_factor > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }
    g.target_factor = 1.0 + g.e_prime * cos_target;
    if (!is_finite(g.target_factor) || !(g.target_factor > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }
    const double encounter_radius = g.p_prime / g.encounter_factor;
    g.encounter_radius_error = encounter_radius - g.r0;
    if (!is_finite(encounter_radius) || std::abs(g.encounter_radius_error) > 1e-10) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }
    const double target_radius = g.p_prime / g.target_factor;
    g.target_radius_error = target_radius - g.r1;
    if (!is_finite(target_radius) || std::abs(g.target_radius_error) > 1e-10) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }
    g.valid = true;
    return g;
}

bool is_boundary_ambiguous(const std::string& strict_reason, const GeometryCheck& g) {
    if (strict_reason != "outgoing_orbit_does_not_pass_target" &&
        strict_reason != "outgoing_orbit_does_not_pass_encounter") {
        return false;
    }
    if (!is_finite(g.denominator) || !(std::abs(g.denominator) > kGeometryDenominatorEpsilon) ||
        !is_finite(g.p_prime) || !(g.p_prime > 0.0) ||
        !is_finite(g.encounter_factor) || !(g.encounter_factor > 0.0) ||
        !is_finite(g.target_factor) || !(g.target_factor > 0.0)) {
        return false;
    }
    const auto roundoff_ok = [](double abs_error, double radius) {
        return (is_finite(abs_error) && abs_error <= kBoundaryRoundoffAbsToleranceM) ||
               (is_finite(abs_error) && is_finite(radius) && radius > 0.0 &&
                abs_error / radius <= kBoundaryRoundoffRelTolerance);
    };
    return roundoff_ok(std::abs(g.encounter_radius_error), g.r0) &&
           roundoff_ok(std::abs(g.target_radius_error), g.r1);
}

DirectValidation validate_direct(
    const spaceship_cpp::problem2::Problem2GravityAssistSolution& root,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - root.theta_prime);

    const auto direct_branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet, phi, beta, theta_A, 1, 1);

    DirectValidation validation{};
    double best_score = std::numeric_limits<double>::infinity();
    problem1::Problem1SolutionBranch best{};
    for (const auto& branch : direct_branches) {
        if (!branch.valid || branch.transfer_revolution != root.transfer_revolution) {
            continue;
        }
        const double alpha_diff = wrapped_angle_distance(branch.target_arrival_true_anomaly, root.alpha);
        const double time_diff = std::abs(branch.time_of_flight_seconds - root.time_of_flight_seconds);
        const double score = alpha_diff + time_diff / 1.0e8;
        if (score < best_score) {
            best_score = score;
            best = branch;
            validation.matched = true;
        }
    }

    if (!validation.matched) {
        return validation;
    }

    validation.matched_k = best.transfer_revolution;
    validation.matched_q = best.target_revolution;
    validation.alpha_diff = wrapped_angle_distance(best.target_arrival_true_anomaly, root.alpha);
    validation.time_diff_seconds = std::abs(best.time_of_flight_seconds - root.time_of_flight_seconds);
    validation.direct_problem1_residual_seconds = best.residual_seconds;

    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);
    const auto strict = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
        departure_params.orbit.p,
        departure_params.orbit.e,
        target_params.orbit.p,
        target_params.orbit.e,
        phi,
        best.target_arrival_true_anomaly,
        kIncomingE,
        kIncomingTheta,
        root.theta_prime);
    if (strict.valid) {
        validation.direct_residual_source = "Strict";
        validation.direct_slingshot_residual = strict.slingshot_residual;
    } else {
        const GeometryCheck geometry = compute_geometry(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            phi,
            best.target_arrival_true_anomaly,
            root.theta_prime);
        const std::string strict_reason = strict.invalid_reason.empty() ? "slingshot_invalid" : strict.invalid_reason;
        if (is_boundary_ambiguous(strict_reason, geometry)) {
            const auto relaxed = problem2::evaluate_problem2_slingshot_residual(
                phi,
                departure_params.orbit.e,
                kIncomingE,
                kIncomingTheta,
                geometry.e_prime,
                root.theta_prime);
            if (relaxed.valid) {
                validation.direct_residual_source = "BoundaryAmbiguousRoundoff";
                validation.direct_slingshot_residual = relaxed.residual;
            }
        }
    }

    validation.validated = validation.matched &&
                           validation.alpha_diff <= 1e-7 &&
                           validation.time_diff_seconds <= 1e3 &&
                           is_finite(validation.direct_slingshot_residual) &&
                           std::abs(validation.direct_slingshot_residual) <= 1e-7;
    return validation;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    namespace problem2 = spaceship_cpp::problem2;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_scenario2_adaptive_new_roots_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = kEncounterTimeFactor * planet_params::planet_orbital_period(departure_planet);

    problem2::Problem2GravityAssistSolverOptions baseline_options{};
    baseline_options.theta_sample_count = 64;
    baseline_options.topology_adaptive_enabled = false;
    baseline_options.allow_roundoff_boundary_relaxed_residual = true;

    auto adaptive64_options = baseline_options;
    adaptive64_options.topology_adaptive_enabled = true;
    adaptive64_options.topology_max_depth = 10;
    adaptive64_options.topology_epsilon = 1e-5;

    auto adaptive128_options = adaptive64_options;
    adaptive128_options.theta_sample_count = 128;

    const auto baseline64 = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, kIncomingE, kIncomingTheta, baseline_options);
    const auto adaptive64 = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, kIncomingE, kIncomingTheta, adaptive64_options);
    const auto adaptive128 = problem2::solve_problem2_gravity_assist_with_table(
        loader, departure_planet, target_planet, encounter_time, kIncomingE, kIncomingTheta, adaptive128_options);

    assert(baseline64.ok);
    assert(adaptive64.ok);
    assert(adaptive128.ok);

    const RootSetComparison baseline_vs_adaptive64 =
        compare_root_sets(baseline64.solutions, adaptive64.solutions);
    const RootSetComparison adaptive64_vs_128 =
        compare_root_sets(adaptive64.solutions, adaptive128.solutions);

    std::cout << "Problem2Scenario2AdaptiveNewRootComparison\n";
    std::cout << "baseline_solution_count=" << baseline64.solutions.size() << '\n';
    std::cout << "adaptive_solution_count=" << adaptive64.solutions.size() << '\n';
    std::cout << "matched_root_count=" << baseline_vs_adaptive64.matched_root_count << '\n';
    std::cout << "new_adaptive_root_count=" << baseline_vs_adaptive64.new_root_in_b_count << '\n';
    std::cout << "missing_baseline_root_count=" << baseline_vs_adaptive64.missing_root_from_a_count << '\n';

    assert(baseline64.solutions.size() == 12);
    assert(adaptive64.solutions.size() == 16);
    assert(baseline_vs_adaptive64.matched_root_count == 12);
    assert(baseline_vs_adaptive64.new_root_in_b_count == 4);
    assert(baseline_vs_adaptive64.missing_root_from_a_count == 0);

    int validation_success_count = 0;
    int validation_failure_count = 0;
    double max_abs_new_root_direct_residual = 0.0;
    int output_root_index = 0;
    for (const int adaptive_index : baseline_vs_adaptive64.unmatched_b_indices) {
        const auto& root = adaptive64.solutions[static_cast<std::size_t>(adaptive_index)];
        assert(root.origin_was_topology_change);
        assert(std::abs(root.slingshot_residual) <= 1e-8);

        std::cout << "Problem2Scenario2AdaptiveNewRoot\n";
        std::cout << "root_index=" << output_root_index << '\n';
        std::cout << "theta_prime=" << root.theta_prime << '\n';
        std::cout << "alpha=" << root.alpha << '\n';
        std::cout << "k=" << root.transfer_revolution << '\n';
        std::cout << "q=" << root.target_revolution << '\n';
        std::cout << "time_of_flight_seconds=" << root.time_of_flight_seconds << '\n';
        std::cout << "target_time_seconds=" << root.target_time_seconds << '\n';
        std::cout << "outgoing_eccentricity=" << root.outgoing_eccentricity << '\n';
        std::cout << "outgoing_semi_latus_rectum=" << root.outgoing_semi_latus_rectum << '\n';
        std::cout << "slingshot_residual=" << root.slingshot_residual << '\n';
        std::cout << "problem1_residual_seconds=" << root.problem1_residual_seconds << '\n';
        std::cout << "residual_source=" << residual_source_name(root.residual_source) << '\n';
        std::cout << "boundary_ambiguous=" << (root.boundary_ambiguous ? 1 : 0) << '\n';
        std::cout << "origin_was_topology_change=" << (root.origin_was_topology_change ? 1 : 0) << '\n';
        std::cout << "bisection_iterations=" << root.bisection_iterations << '\n';
        std::cout << "final_theta_width=" << root.final_theta_width << '\n';

        const DirectValidation validation = validate_direct(root, departure_planet, target_planet, encounter_time);
        if (validation.validated) {
            validation_success_count += 1;
            max_abs_new_root_direct_residual =
                std::max(max_abs_new_root_direct_residual, std::abs(validation.direct_slingshot_residual));
        } else {
            validation_failure_count += 1;
        }

        std::cout << "Problem2Scenario2NewRootDirectValidation\n";
        std::cout << "root_index=" << output_root_index << '\n';
        std::cout << "matched=" << (validation.matched ? 1 : 0) << '\n';
        std::cout << "matched_k=" << validation.matched_k << '\n';
        std::cout << "matched_q=" << validation.matched_q << '\n';
        std::cout << "alpha_diff=" << validation.alpha_diff << '\n';
        std::cout << "time_diff_seconds=" << validation.time_diff_seconds << '\n';
        std::cout << "direct_problem1_residual_seconds=" << validation.direct_problem1_residual_seconds << '\n';
        std::cout << "direct_slingshot_residual=" << validation.direct_slingshot_residual << '\n';
        std::cout << "direct_residual_source=" << validation.direct_residual_source << '\n';
        std::cout << "validated=" << (validation.validated ? 1 : 0) << '\n';
        output_root_index += 1;
    }

    std::cout << "Problem2Scenario2AdaptiveScalingComparison\n";
    std::cout << "sample_a=64\n";
    std::cout << "sample_b=128\n";
    std::cout << "root_count_a=" << adaptive64.solutions.size() << '\n';
    std::cout << "root_count_b=" << adaptive128.solutions.size() << '\n';
    std::cout << "matched_root_count=" << adaptive64_vs_128.matched_root_count << '\n';
    std::cout << "new_root_in_b_count=" << adaptive64_vs_128.new_root_in_b_count << '\n';
    std::cout << "missing_root_from_a_count=" << adaptive64_vs_128.missing_root_from_a_count << '\n';
    std::cout << "max_theta_diff_matched=" << adaptive64_vs_128.max_theta_diff_matched << '\n';
    std::cout << "max_alpha_diff_matched=" << adaptive64_vs_128.max_alpha_diff_matched << '\n';

    const bool all_new_roots_validated = validation_success_count == 4 && validation_failure_count == 0;
    const bool scenario_ok = all_new_roots_validated &&
                             adaptive64_vs_128.missing_root_from_a_count == 0;

    std::cout << "Problem2Scenario2AdaptiveNewRootsSummary\n";
    std::cout << "baseline64_solution_count=" << baseline64.solutions.size() << '\n';
    std::cout << "adaptive64_solution_count=" << adaptive64.solutions.size() << '\n';
    std::cout << "adaptive128_solution_count=" << adaptive128.solutions.size() << '\n';
    std::cout << "new_adaptive64_root_count=" << baseline_vs_adaptive64.new_root_in_b_count << '\n';
    std::cout << "new_adaptive64_direct_validation_success_count=" << validation_success_count << '\n';
    std::cout << "new_adaptive64_direct_validation_failure_count=" << validation_failure_count << '\n';
    std::cout << "adaptive64_vs_128_matched_count=" << adaptive64_vs_128.matched_root_count << '\n';
    std::cout << "adaptive64_vs_128_new_count=" << adaptive64_vs_128.new_root_in_b_count << '\n';
    std::cout << "adaptive64_vs_128_missing_count=" << adaptive64_vs_128.missing_root_from_a_count << '\n';
    std::cout << "max_abs_new_root_direct_residual=" << max_abs_new_root_direct_residual << '\n';
    std::cout << "all_new_roots_validated=" << (all_new_roots_validated ? 1 : 0) << '\n';
    std::cout << "scenario2_adaptive_new_roots_ok=" << (scenario_ok ? 1 : 0) << '\n';

    std::cout << "Problem2Scenario2TopologyAdaptiveRegressionSummary\n";
    std::cout << "fast_solution_count=" << baseline64.solutions.size() << '\n';
    std::cout << "adaptive_solution_count=" << adaptive64.solutions.size() << '\n';
    std::cout << "new_topology_change_root_count=" << baseline_vs_adaptive64.new_root_in_b_count << '\n';
    std::cout << "new_roots_direct_validation_success_count=" << validation_success_count << '\n';
    std::cout << "adaptive64_vs_128_missing_count=" << adaptive64_vs_128.missing_root_from_a_count << '\n';
    std::cout << "regression_ok=" << (scenario_ok ? 1 : 0) << '\n';

    return scenario_ok ? 0 : 1;
}
