#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace spaceship_cpp::problem1 {

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::common::orbit_F;
using spaceship_cpp::common::orbit_F_e_derivative;
using spaceship_cpp::common::orbit_F_theta_derivative;

struct Problem1DepartureAnomalyResidualInput {
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;
    double nu_A_depart;
    double nu_B_depart;
    double theta_A;
    double encounter_global_angle;
    int transfer_revolution;
    int target_revolution;
};

struct Problem1DepartureAnomalySolveInput {
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;
    double nu_A_depart;
    double nu_B_depart;
    double theta_A;
    int max_transfer_revolution;
    int max_target_revolution;
    int phi_scan_count;
    double phi_tolerance;
    double residual_tolerance;
    int max_bisection_iterations;
    double max_candidate_relative_residual;
};

struct Problem1DepartureAnomalyResidualResult {
    Problem1ResidualStatus status = Problem1ResidualStatus::InvalidInput;
    double residual = std::numeric_limits<double>::quiet_NaN();
    double encounter_global_angle = std::numeric_limits<double>::quiet_NaN();
    double target_arrival_true_anomaly = std::numeric_limits<double>::quiet_NaN();
    double transfer_e_raw = std::numeric_limits<double>::quiet_NaN();
    double transfer_e = std::numeric_limits<double>::quiet_NaN();
    double transfer_p = std::numeric_limits<double>::quiet_NaN();
    double theta_B = std::numeric_limits<double>::quiet_NaN();
    double deltaF_transfer = std::numeric_limits<double>::quiet_NaN();
    double deltaF_target = std::numeric_limits<double>::quiet_NaN();
    double transfer_time_scale_free = std::numeric_limits<double>::quiet_NaN();
    double target_time_scale_free = std::numeric_limits<double>::quiet_NaN();
    std::string invalid_reason;
};

struct Problem1SolutionBranchCandidate {
    Problem1SolutionBranch branch;
    double relative_residual = std::numeric_limits<double>::quiet_NaN();
    double root_bracket_width = std::numeric_limits<double>::quiet_NaN();
};

struct Problem1TransferOrbitGeometry {
    bool valid = false;
    double e_raw = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double e_sign = 1.0;
    bool negative_e_normalized = false;
    double p = std::numeric_limits<double>::quiet_NaN();
    double theta_A_raw = std::numeric_limits<double>::quiet_NaN();
    double theta_B_raw = std::numeric_limits<double>::quiet_NaN();
    double theta_A_canonical = std::numeric_limits<double>::quiet_NaN();
    double theta_B_canonical = std::numeric_limits<double>::quiet_NaN();
    double theta_B_star_canonical = std::numeric_limits<double>::quiet_NaN();
    double r_A = std::numeric_limits<double>::quiet_NaN();
    double r_B = std::numeric_limits<double>::quiet_NaN();
    double N = std::numeric_limits<double>::quiet_NaN();
    double D = std::numeric_limits<double>::quiet_NaN();
    double lambda_A = std::numeric_limits<double>::quiet_NaN();
    double transfer_perihelion_angle_raw = std::numeric_limits<double>::quiet_NaN();
    double transfer_perihelion_angle_canonical = std::numeric_limits<double>::quiet_NaN();
    std::string invalid_reason;
};

struct Problem1AdaptiveScanPoint {
    double phi = 0.0;
    Problem1DepartureAnomalyResidualResult residual_result;
    bool valid = false;
    double residual = std::numeric_limits<double>::quiet_NaN();
    double abs_residual = std::numeric_limits<double>::quiet_NaN();
    double residual_sq = std::numeric_limits<double>::quiet_NaN();
};

enum class Problem1SuspiciousIntervalReason {
    SignChange,
    NearZero,
    LocalMinimum,
    LocalMaximum,
    ValidBoundary,
    WrapSignChange,
    RapidChange,
};

struct Problem1SuspiciousInterval {
    double left_phi = 0.0;
    double right_phi = 0.0;
    Problem1SuspiciousIntervalReason reason = Problem1SuspiciousIntervalReason::SignChange;
    int left_index = -1;
    int right_index = -1;
    double left_residual = std::numeric_limits<double>::quiet_NaN();
    double right_residual = std::numeric_limits<double>::quiet_NaN();
    double center_residual = std::numeric_limits<double>::quiet_NaN();
};

const char* suspicious_interval_reason_name(Problem1SuspiciousIntervalReason reason) {
    switch (reason) {
    case Problem1SuspiciousIntervalReason::SignChange:
        return "SignChange";
    case Problem1SuspiciousIntervalReason::NearZero:
        return "NearZero";
    case Problem1SuspiciousIntervalReason::LocalMinimum:
        return "LocalMinimum";
    case Problem1SuspiciousIntervalReason::LocalMaximum:
        return "LocalMaximum";
    case Problem1SuspiciousIntervalReason::ValidBoundary:
        return "ValidBoundary";
    case Problem1SuspiciousIntervalReason::WrapSignChange:
        return "WrapSignChange";
    case Problem1SuspiciousIntervalReason::RapidChange:
        return "RapidChange";
    }
    return "Unknown";
}

struct Problem1AdaptiveScanSummary {
    long long coarse_scan_count = 0;
    long long interval_total = 0;
    long long sign_change_interval_count = 0;
    long long near_zero_interval_count = 0;
    long long local_min_interval_count = 0;
    long long local_max_interval_count = 0;
    long long valid_boundary_interval_count = 0;
    long long wrap_interval_count = 0;
    long long rapid_change_interval_count = 0;
    long long residual_eval_count = 0;
    long long refined_interval_count = 0;
    long long bisection_interval_count = 0;
    long long ternary_interval_count = 0;
    long long local_fine_scan_interval_count = 0;
    long long candidate_count_before_dedup = 0;
    long long candidate_count_after_dedup = 0;
    long long ternary_accept_count = 0;
    long long ternary_reject_count = 0;
    long long local_fine_scan_root_count = 0;
    bool fallback_to_fullscan = false;
    double coarse_scan_seconds = 0.0;
    double interval_collection_seconds = 0.0;
    double interval_refine_seconds = 0.0;
    double local_fine_scan_seconds = 0.0;
    double ternary_seconds = 0.0;
    double bisection_seconds = 0.0;
    double candidate_dedup_seconds = 0.0;
    double sorting_seconds = 0.0;
};

class ScopedSecondsAccumulator {
public:
    explicit ScopedSecondsAccumulator(double* target)
        : target_(target), start_(std::chrono::steady_clock::now()) {}

    ScopedSecondsAccumulator(const ScopedSecondsAccumulator&) = delete;
    ScopedSecondsAccumulator& operator=(const ScopedSecondsAccumulator&) = delete;

    ~ScopedSecondsAccumulator() {
        if (target_ != nullptr) {
            *target_ += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_).count();
        }
    }

private:
    double* target_;
    std::chrono::steady_clock::time_point start_;
};

Problem1SolveWithDiagnosticResult solve_problem1_from_departure_anomalies_full_scan_internal(
    const Problem1DepartureAnomalySolveInput& input
);

double scale_free_time_to_seconds(double scale_free_time) {
    const double mu = planet_params::get_solar_system_physical_params().GM_sun;
    if (!is_finite(scale_free_time) || !(mu > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return scale_free_time / std::sqrt(mu);
}

double problem1_residual_seconds_to_scale_free_internal(double residual_seconds) {
    const double mu = planet_params::get_solar_system_physical_params().GM_sun;
    if (!is_finite(residual_seconds) || !(mu > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return residual_seconds * std::sqrt(mu);
}

double compute_transfer_a(double transfer_p, double transfer_e) {
    const double denominator = 1.0 - transfer_e * transfer_e;
    if (!is_finite(transfer_p) || !is_finite(transfer_e) || std::abs(denominator) <= kDefaultEpsilon) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return transfer_p / denominator;
}

Problem1DepartureAnomalyResidualResult make_departure_anomaly_result(Problem1ResidualStatus status) {
    Problem1DepartureAnomalyResidualResult result{};
    result.status = status;
    result.invalid_reason = "invalid departure-anomaly residual";
    return result;
}

void validate_problem1_root_table_config(const Problem1RootTableConfig& config) {
    if (!is_finite(config.nu_A_start) ||
        !is_finite(config.nu_A_step) ||
        config.nu_A_count <= 0 ||
        !is_finite(config.nu_B_depart_start) ||
        !is_finite(config.nu_B_depart_step) ||
        config.nu_B_depart_count <= 0 ||
        !is_finite(config.theta_A_start) ||
        !is_finite(config.theta_A_step) ||
        config.theta_A_count <= 0 ||
        config.max_transfer_revolution < 0 ||
        config.max_target_revolution < 0 ||
        config.schema_version.empty()) {
        throw std::invalid_argument("invalid Problem1RootTableConfig");
    }
}

void validate_departure_anomaly_solve_input(const Problem1DepartureAnomalySolveInput& input) {
    if (!is_finite(input.nu_A_depart) ||
        !is_finite(input.nu_B_depart) ||
        !is_finite(input.theta_A) ||
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
        throw std::invalid_argument("invalid departure-anomaly solve input");
    }
}

bool positive_orbit_denominator(double e, double angle) {
    const double denominator = 1.0 + e * std::cos(angle);
    return is_finite(denominator) && denominator > kDefaultEpsilon;
}

bool is_success_with_finite_residual(const Problem1DepartureAnomalyResidualResult& result) {
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

double compute_problem1_residual_scale(const Problem1DepartureAnomalyResidualResult& result) {
    if (result.status != Problem1ResidualStatus::Success ||
        !is_finite(result.transfer_time_scale_free) ||
        !is_finite(result.target_time_scale_free)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::max(
        {std::abs(result.transfer_time_scale_free), std::abs(result.target_time_scale_free), 1.0});
}

double compute_problem1_relative_residual(const Problem1DepartureAnomalyResidualResult& result) {
    const double scale = compute_problem1_residual_scale(result);
    if (!is_finite(scale) || !(scale > 0.0) || !is_finite(result.residual)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::abs(result.residual) / scale;
}

const char* problem1_root_hessian_method_name(Problem1RootHessianMethod method) {
    switch (method) {
    case Problem1RootHessianMethod::TangentFiniteDifference:
        return "tangent_finite_difference_of_implicit_first_derivatives";
    case Problem1RootHessianMethod::ProjectedTangentFiniteDifference:
        return "projected_tangent_finite_difference_of_implicit_first_derivatives";
    case Problem1RootHessianMethod::NewtonRefinedFiniteDifference:
        return "newton_refined_finite_difference_of_implicit_first_derivatives";
    }
    return "unknown_hessian_method";
}

void update_abs_max_if_finite(double value, double* current_max) {
    if (!is_finite(value)) {
        return;
    }
    const double abs_value = std::abs(value);
    if (!is_finite(*current_max) || abs_value > *current_max) {
        *current_max = abs_value;
    }
}

Problem1SolutionBranch make_solution_branch_from_public_residual(
    const Problem1RootResidualResult& residual,
    int transfer_revolution,
    int target_revolution
) {
    Problem1SolutionBranch branch{};
    branch.valid = true;
    branch.encounter_global_angle = residual.encounter_global_angle;
    branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    branch.transfer_revolution = transfer_revolution;
    branch.target_revolution = target_revolution;
    branch.time_of_flight_seconds = residual.transfer_time_seconds;
    branch.target_time_seconds = residual.target_time_seconds;
    branch.residual_seconds = residual.residual_seconds;
    branch.arrival_time_seconds_since_j2000 = std::numeric_limits<double>::quiet_NaN();
    branch.transfer_e = residual.transfer_e;
    branch.transfer_p = residual.transfer_p;
    branch.transfer_a = residual.transfer_a;
    branch.theta_B = residual.theta_B;
    return branch;
}

struct Problem1ProjectedTangentStencilBranchResult {
    bool valid = false;
    Problem1SolutionBranch branch;
    double residual_before_projection_seconds = std::numeric_limits<double>::quiet_NaN();
    double residual_after_projection_seconds = std::numeric_limits<double>::quiet_NaN();
    double projection_delta_alpha = std::numeric_limits<double>::quiet_NaN();
    double alpha_tangent = std::numeric_limits<double>::quiet_NaN();
    double alpha_projected = std::numeric_limits<double>::quiet_NaN();
    bool projection_derivative_valid = false;
    bool projected_stencil_attach_success = false;
    std::string invalid_reason;
};

Problem1ProjectedTangentStencilBranchResult evaluate_problem1_root_projected_tangent_stencil_branch(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& differentiated_source_branch,
    double stencil_nu_A,
    double stencil_nu_B,
    double stencil_theta_A,
    double residual_tolerance_seconds,
    double projection_derivative_min_abs
) {
    Problem1ProjectedTangentStencilBranchResult result{};
    if (!differentiated_source_branch.valid || !differentiated_source_branch.derivatives_available) {
        result.invalid_reason = "projected_stencil_source_derivatives_unavailable";
        return result;
    }
    if (!is_finite(residual_tolerance_seconds) || !(residual_tolerance_seconds > 0.0)) {
        result.invalid_reason = "projected_stencil_residual_tolerance_seconds_must_be_positive";
        return result;
    }
    if (!is_finite(projection_derivative_min_abs) || !(projection_derivative_min_abs > 0.0)) {
        result.invalid_reason = "projected_stencil_projection_derivative_min_abs_must_be_positive";
        return result;
    }

    const double dx1 = normalize_angle_minus_pi_pi(stencil_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(stencil_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(stencil_theta_A - node_theta_A);
    result.alpha_tangent = normalize_angle_0_2pi(
        differentiated_source_branch.encounter_global_angle +
        differentiated_source_branch.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated_source_branch.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated_source_branch.d_encounter_global_angle_d_theta_A * dx3);
    const Problem1RootResidualResult residual_before = evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        result.alpha_tangent,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    result.residual_before_projection_seconds = residual_before.residual_seconds;
    const Problem1RootResidualDerivatives derivative_before =
        evaluate_problem1_root_residual_derivatives(
            departure_planet,
            target_planet,
            stencil_nu_A,
            stencil_nu_B,
            stencil_theta_A,
            result.alpha_tangent,
            differentiated_source_branch.transfer_revolution,
            differentiated_source_branch.target_revolution);
    result.projection_derivative_valid =
        derivative_before.valid && is_finite(derivative_before.R_alpha);
    if (!residual_before.valid ||
        !is_finite(residual_before.residual_seconds) ||
        !result.projection_derivative_valid) {
        result.invalid_reason = "projected_stencil_residual_or_derivative_invalid";
        return result;
    }
    if (std::abs(derivative_before.R_alpha) <= projection_derivative_min_abs) {
        result.invalid_reason = "projected_stencil_derivative_too_small";
        return result;
    }

    const double residual_scale_free =
        problem1_residual_seconds_to_scale_free(residual_before.residual_seconds);
    result.projection_delta_alpha = -residual_scale_free / derivative_before.R_alpha;
    result.alpha_projected = normalize_angle_0_2pi(result.alpha_tangent + result.projection_delta_alpha);
    const Problem1RootResidualResult residual_after = evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        result.alpha_projected,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    result.residual_after_projection_seconds = residual_after.residual_seconds;
    if (!residual_after.valid || !is_finite(residual_after.residual_seconds)) {
        result.invalid_reason = "projected_stencil_residual_after_invalid";
        return result;
    }
    if (std::abs(residual_after.residual_seconds) > residual_tolerance_seconds) {
        result.invalid_reason = "projected_stencil_residual_after_too_large";
        return result;
    }

    Problem1SolutionBranch projected_branch = make_solution_branch_from_public_residual(
        residual_after,
        differentiated_source_branch.transfer_revolution,
        differentiated_source_branch.target_revolution);
    projected_branch = attach_problem1_root_derivatives(
        departure_planet,
        target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        projected_branch);
    result.projected_stencil_attach_success =
        projected_branch.valid && projected_branch.derivatives_available;
    if (!result.projected_stencil_attach_success) {
        result.invalid_reason = "hessian_attach_failed";
        return result;
    }
    result.valid = true;
    result.branch = projected_branch;
    return result;
}

Problem1RootRefinementResult make_invalid_refinement_result(
    int transfer_revolution,
    int target_revolution,
    double initial_encounter_global_angle,
    const std::string& invalid_reason
) {
    Problem1RootRefinementResult result{};
    result.valid = false;
    result.branch.valid = false;
    result.branch.transfer_revolution = transfer_revolution;
    result.branch.target_revolution = target_revolution;
    result.branch.invalid_reason = invalid_reason;
    result.diagnostic.transfer_revolution = transfer_revolution;
    result.diagnostic.target_revolution = target_revolution;
    result.diagnostic.initial_alpha = initial_encounter_global_angle;
    result.diagnostic.final_alpha = initial_encounter_global_angle;
    result.diagnostic.invalid_reason = invalid_reason;
    return result;
}

Problem1DepartureAnomalyResidualResult evaluate_problem1_residual_from_departure_anomalies(
    const Problem1DepartureAnomalyResidualInput& input
);

Problem1DepartureAnomalyResidualResult evaluate_problem1_residual_from_departure_anomalies_profiled(
    const Problem1DepartureAnomalyResidualInput& input,
    Problem1SolveDiagnostic* diagnostic
) {
    using clock = std::chrono::steady_clock;
    const auto start = clock::now();
    const Problem1DepartureAnomalyResidualResult result =
        evaluate_problem1_residual_from_departure_anomalies(input);
    const auto end = clock::now();
    if (diagnostic != nullptr) {
        diagnostic->residual_evaluations += 1;
        diagnostic->residual_evaluation_seconds += std::chrono::duration<double>(end - start).count();
    }
    return result;
}

Problem1DepartureAnomalyResidualInput make_departure_anomaly_residual_input(
    const Problem1DepartureAnomalySolveInput& input,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    return Problem1DepartureAnomalyResidualInput{
        input.departure_planet,
        input.target_planet,
        input.nu_A_depart,
        input.nu_B_depart,
        input.theta_A,
        encounter_global_angle,
        transfer_revolution,
        target_revolution,
    };
}

Problem1TransferOrbitGeometry compute_problem1_transfer_orbit_geometry(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double theta_A,
    double encounter_global_angle
) {
    Problem1TransferOrbitGeometry geometry{};
    const planet_params::PlanetParams& departure = planet_params::get_planet_params(departure_planet);
    const planet_params::PlanetParams& target = planet_params::get_planet_params(target_planet);

    const double nu_A = normalize_angle_0_2pi(nu_A_depart);
    const double theta_A_wrapped = normalize_angle_0_2pi(theta_A);
    const double encounter_wrapped = normalize_angle_0_2pi(encounter_global_angle);
    const double lambda_A = normalize_angle_0_2pi(departure.orbit.theta_0 + nu_A);
    const double target_arrival_true_anomaly = normalize_angle_0_2pi(encounter_wrapped - target.orbit.theta_0);
    const double r_A = planet_params::planet_radius_at_true_anomaly(departure_planet, nu_A);
    const double r_B = planet_params::planet_radius_at_true_anomaly(target_planet, target_arrival_true_anomaly);
    const double transfer_perihelion_angle_raw = normalize_angle_0_2pi(lambda_A - theta_A_wrapped);
    const double theta_A_raw = theta_A_wrapped;
    const double theta_B_raw = normalize_angle_0_2pi(encounter_wrapped - transfer_perihelion_angle_raw);

    double e_raw = std::numeric_limits<double>::quiet_NaN();
    double p = std::numeric_limits<double>::quiet_NaN();
    try {
        e_raw = compute_transfer_e_from_two_points(r_A, theta_A_raw, r_B, theta_B_raw);
        p = compute_transfer_p_from_departure(r_A, e_raw, theta_A_raw);
    } catch (const std::domain_error&) {
        geometry.invalid_reason = "singular geometry";
        return geometry;
    }

    geometry.e_raw = e_raw;
    geometry.e_sign = e_raw < 0.0 ? -1.0 : 1.0;
    geometry.negative_e_normalized = e_raw < 0.0;
    geometry.e = std::abs(e_raw);
    geometry.p = p;
    geometry.theta_A_raw = theta_A_raw;
    geometry.theta_B_raw = theta_B_raw;
    geometry.theta_A_canonical = theta_A_raw;
    geometry.theta_B_canonical = theta_B_raw;
    geometry.theta_B_star_canonical = theta_B_raw;
    geometry.r_A = r_A;
    geometry.r_B = r_B;
    geometry.lambda_A = lambda_A;
    geometry.transfer_perihelion_angle_raw = transfer_perihelion_angle_raw;
    geometry.transfer_perihelion_angle_canonical = transfer_perihelion_angle_raw;
    if (geometry.negative_e_normalized) {
        geometry.transfer_perihelion_angle_canonical =
            normalize_angle_0_2pi(transfer_perihelion_angle_raw + spaceship_cpp::common::kPi);
        geometry.theta_A_canonical = normalize_angle_0_2pi(lambda_A - geometry.transfer_perihelion_angle_canonical);
        geometry.theta_B_canonical =
            normalize_angle_0_2pi(encounter_wrapped - geometry.transfer_perihelion_angle_canonical);
        geometry.theta_B_star_canonical = geometry.theta_B_canonical;
    }
    geometry.N = r_B - r_A;
    geometry.D =
        r_A * std::cos(geometry.theta_A_canonical) -
        r_B * std::cos(geometry.theta_B_canonical);

    if (!is_finite(geometry.e) || !(geometry.e >= 0.0) ||
        !is_finite(geometry.p) || !(geometry.p > 0.0) ||
        !positive_orbit_denominator(geometry.e, geometry.theta_A_canonical) ||
        !positive_orbit_denominator(geometry.e, geometry.theta_B_canonical)) {
        geometry.invalid_reason = "invalid transfer orbit";
        return geometry;
    }

    geometry.valid = true;
    return geometry;
}

Problem1DepartureAnomalyResidualResult evaluate_problem1_residual_from_departure_anomalies(
    const Problem1DepartureAnomalyResidualInput& input
) {
    if (!is_finite(input.nu_A_depart) ||
        !is_finite(input.nu_B_depart) ||
        !is_finite(input.theta_A) ||
        !is_finite(input.encounter_global_angle) ||
        input.transfer_revolution < 0 ||
        input.target_revolution < 0) {
        Problem1DepartureAnomalyResidualResult result =
            make_departure_anomaly_result(Problem1ResidualStatus::InvalidInput);
        result.invalid_reason = "invalid residual inputs";
        return result;
    }

    const planet_params::PlanetParams& departure = planet_params::get_planet_params(input.departure_planet);
    const planet_params::PlanetParams& target = planet_params::get_planet_params(input.target_planet);

    const double nu_A_depart = normalize_angle_0_2pi(input.nu_A_depart);
    const double nu_B_depart = normalize_angle_0_2pi(input.nu_B_depart);
    const double theta_A = normalize_angle_0_2pi(input.theta_A);

    const double r2_true_anomaly = normalize_angle_0_2pi(input.encounter_global_angle - target.orbit.theta_0);
    const Problem1TransferOrbitGeometry geometry = compute_problem1_transfer_orbit_geometry(
        input.departure_planet,
        input.target_planet,
        nu_A_depart,
        theta_A,
        input.encounter_global_angle);
    if (!geometry.valid) {
        Problem1DepartureAnomalyResidualResult result =
            make_departure_anomaly_result(
                geometry.invalid_reason == "singular geometry"
                    ? Problem1ResidualStatus::SingularGeometry
                    : Problem1ResidualStatus::InvalidTransferOrbit);
        result.invalid_reason = geometry.invalid_reason;
        return result;
    }

    double theta2_end = r2_true_anomaly;
    while (theta2_end <= nu_B_depart) {
        theta2_end += kTwoPi;
    }
    theta2_end += static_cast<double>(input.target_revolution) * kTwoPi;

    double xi1 = std::numeric_limits<double>::quiet_NaN();
    double xi2 = std::numeric_limits<double>::quiet_NaN();
    if (geometry.e < 1.0) {
        xi1 = geometry.theta_A_canonical;
        xi2 = geometry.theta_B_canonical;
        while (xi2 <= xi1) {
            xi2 += kTwoPi;
        }
        xi2 += static_cast<double>(input.transfer_revolution) * kTwoPi;
    } else {
        if (input.transfer_revolution != 0) {
            Problem1DepartureAnomalyResidualResult result =
                make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);
            result.invalid_reason = "hyperbolic transfer does not support k > 0";
            return result;
        }
        xi1 = normalize_angle_minus_pi_pi(geometry.theta_A_canonical);
        xi2 = normalize_angle_minus_pi_pi(geometry.theta_B_canonical);
        if (xi2 <= xi1) {
            Problem1DepartureAnomalyResidualResult result =
                make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);
            result.invalid_reason = "hyperbolic branch is not forward";
            return result;
        }
    }

    double deltaF_transfer = std::numeric_limits<double>::quiet_NaN();
    double deltaF_target = std::numeric_limits<double>::quiet_NaN();
    try {
        deltaF_transfer = orbit_F(geometry.e, xi2) - orbit_F(geometry.e, xi1);
        deltaF_target = orbit_F(target.orbit.e, theta2_end) - orbit_F(target.orbit.e, nu_B_depart);
    } catch (const std::domain_error&) {
        Problem1DepartureAnomalyResidualResult result =
            make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);
        result.invalid_reason = "orbit_F is undefined on requested branch";
        return result;
    }

    if (!is_finite(deltaF_transfer) || !is_finite(deltaF_target) ||
        !(deltaF_transfer > 0.0) || !(deltaF_target > 0.0)) {
        Problem1DepartureAnomalyResidualResult result =
            make_departure_anomaly_result(Problem1ResidualStatus::InvalidTimeOfFlight);
        result.invalid_reason = "non-positive branch time of flight";
        return result;
    }

    const double transfer_time_scale_free = std::pow(geometry.p, 1.5) * deltaF_transfer;
    const double target_time_scale_free = std::pow(target.orbit.p, 1.5) * deltaF_target;
    if (!is_finite(transfer_time_scale_free) || !is_finite(target_time_scale_free)) {
        Problem1DepartureAnomalyResidualResult result =
            make_departure_anomaly_result(Problem1ResidualStatus::InvalidTimeOfFlight);
        result.invalid_reason = "non-finite branch time of flight";
        return result;
    }

    Problem1DepartureAnomalyResidualResult result{};
    result.status = Problem1ResidualStatus::Success;
    result.residual = transfer_time_scale_free - target_time_scale_free;
    result.encounter_global_angle = normalize_angle_0_2pi(input.encounter_global_angle);
    result.target_arrival_true_anomaly = r2_true_anomaly;
    result.transfer_e_raw = geometry.e_raw;
    result.transfer_e = geometry.e;
    result.transfer_p = geometry.p;
    result.theta_B = xi2;
    result.deltaF_transfer = deltaF_transfer;
    result.deltaF_target = deltaF_target;
    result.transfer_time_scale_free = transfer_time_scale_free;
    result.target_time_scale_free = target_time_scale_free;
    result.invalid_reason.clear();
    return result;
}

Problem1RootResidualResult make_public_root_residual_result(const Problem1DepartureAnomalyResidualResult& internal) {
    Problem1RootResidualResult result{};
    result.status = internal.status;
    result.valid = internal.status == Problem1ResidualStatus::Success && is_finite(internal.residual);
    result.residual_scale_free = internal.residual;
    result.residual_seconds = scale_free_time_to_seconds(internal.residual);
    result.transfer_time_scale_free = internal.transfer_time_scale_free;
    result.target_time_scale_free = internal.target_time_scale_free;
    result.transfer_time_seconds = scale_free_time_to_seconds(internal.transfer_time_scale_free);
    result.target_time_seconds = scale_free_time_to_seconds(internal.target_time_scale_free);
    result.encounter_global_angle = internal.encounter_global_angle;
    result.target_arrival_true_anomaly = internal.target_arrival_true_anomaly;
    result.transfer_e_raw = internal.transfer_e_raw;
    result.transfer_e = internal.transfer_e;
    result.transfer_p = internal.transfer_p;
    result.transfer_a = compute_transfer_a(internal.transfer_p, internal.transfer_e);
    result.theta_B = internal.theta_B;
    result.invalid_reason = internal.invalid_reason;
    return result;
}

bool near_unwrap_boundary(double start_angle, double end_angle, double eps = 1e-10) {
    const double wrapped_end = normalize_angle_0_2pi(end_angle);
    return std::abs(wrapped_end - start_angle) <= eps;
}

double wrapped_angle_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

double branch_relative_residual_seconds(const Problem1SolutionBranch& branch) {
    if (!is_finite(branch.residual_seconds)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double scale = std::max(
        {std::abs(branch.time_of_flight_seconds), std::abs(branch.target_time_seconds), 1.0});
    if (!is_finite(scale) || !(scale > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::abs(branch.residual_seconds) / scale;
}

double branch_time_identity_tolerance(
    const Problem1SolutionBranch& lhs,
    const Problem1SolutionBranch& rhs,
    const Problem1BranchCompareOptions& options
) {
    const double scale = std::max(
        {std::abs(lhs.time_of_flight_seconds),
         std::abs(rhs.time_of_flight_seconds),
         std::abs(lhs.target_time_seconds),
         std::abs(rhs.target_time_seconds),
         1.0});
    return std::max(options.time_tolerance_seconds, options.relative_time_tolerance * scale);
}

Problem1BranchCompareDetail make_branch_compare_detail(
    const Problem1SolutionBranch& branch,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A
) {
    Problem1BranchCompareDetail detail{};
    detail.nu_A_depart = nu_A_depart;
    detail.nu_B_depart = nu_B_depart;
    detail.theta_A = theta_A;
    detail.transfer_revolution = branch.transfer_revolution;
    detail.target_revolution = branch.target_revolution;
    detail.encounter_global_angle = branch.encounter_global_angle;
    detail.target_arrival_true_anomaly = branch.target_arrival_true_anomaly;
    detail.time_of_flight_seconds = branch.time_of_flight_seconds;
    detail.target_time_seconds = branch.target_time_seconds;
    detail.residual_seconds = branch.residual_seconds;
    detail.relative_residual = branch_relative_residual_seconds(branch);
    return detail;
}

bool branch_matches_by_identity(
    const Problem1SolutionBranch& candidate,
    const Problem1SolutionBranch& baseline,
    const Problem1BranchCompareOptions& options
) {
    if (candidate.transfer_revolution != baseline.transfer_revolution ||
        candidate.target_revolution != baseline.target_revolution) {
        return false;
    }
    if (wrapped_angle_distance(candidate.encounter_global_angle, baseline.encounter_global_angle) >
        options.angle_tolerance) {
        return false;
    }
    if (wrapped_angle_distance(candidate.target_arrival_true_anomaly, baseline.target_arrival_true_anomaly) >
        options.target_anomaly_tolerance) {
        return false;
    }
    const double time_tol = branch_time_identity_tolerance(candidate, baseline, options);
    return std::abs(candidate.time_of_flight_seconds - baseline.time_of_flight_seconds) <= time_tol ||
        std::abs(candidate.target_time_seconds - baseline.target_time_seconds) <= time_tol;
}

double branch_identity_match_score(
    const Problem1SolutionBranch& candidate,
    const Problem1SolutionBranch& baseline,
    const Problem1BranchCompareOptions& options
) {
    const double time_tol = branch_time_identity_tolerance(candidate, baseline, options);
    const double angle_score =
        wrapped_angle_distance(candidate.encounter_global_angle, baseline.encounter_global_angle) /
        std::max(options.angle_tolerance, kDefaultEpsilon);
    const double target_anomaly_score =
        wrapped_angle_distance(candidate.target_arrival_true_anomaly, baseline.target_arrival_true_anomaly) /
        std::max(options.target_anomaly_tolerance, kDefaultEpsilon);
    const double tof_score =
        std::abs(candidate.time_of_flight_seconds - baseline.time_of_flight_seconds) /
        std::max(time_tol, kDefaultEpsilon);
    const double target_time_score =
        std::abs(candidate.target_time_seconds - baseline.target_time_seconds) /
        std::max(time_tol, kDefaultEpsilon);
    return angle_score + target_anomaly_score + std::min(tof_score, target_time_score);
}

void apply_branch_compare_report_to_profile(
    const Problem1BranchCompareReport& report,
    Problem1SolveModeProfile* profile
) {
    profile->debug_compare_matched_count = report.matched_count;
    profile->debug_compare_missing_from_adaptive = report.missing_count;
    profile->debug_compare_new_in_adaptive = report.extra_count;
    profile->debug_compare_max_alpha_diff = report.max_angle_error;
    profile->debug_compare_max_theta_prime_diff = report.max_target_anomaly_error;
    profile->debug_compare_max_target_time_diff = report.max_time_error;
    profile->debug_compare_max_residual_diff = report.max_residual_error;
}

void apply_adaptive_scan_summary_to_profile(
    const Problem1AdaptiveScanSummary& summary,
    Problem1SolveModeProfile* profile
) {
    profile->adaptive_coarse_scan_count = summary.coarse_scan_count;
    profile->adaptive_interval_total = summary.interval_total;
    profile->adaptive_sign_change_interval_count = summary.sign_change_interval_count;
    profile->adaptive_near_zero_interval_count = summary.near_zero_interval_count;
    profile->adaptive_local_min_interval_count = summary.local_min_interval_count;
    profile->adaptive_local_max_interval_count = summary.local_max_interval_count;
    profile->adaptive_valid_boundary_interval_count = summary.valid_boundary_interval_count;
    profile->adaptive_wrap_interval_count = summary.wrap_interval_count;
    profile->adaptive_rapid_change_interval_count = summary.rapid_change_interval_count;
    profile->adaptive_residual_eval_count = summary.residual_eval_count;
    profile->adaptive_refined_interval_count = summary.refined_interval_count;
    profile->adaptive_bisection_interval_count = summary.bisection_interval_count;
    profile->adaptive_ternary_interval_count = summary.ternary_interval_count;
    profile->adaptive_local_fine_scan_interval_count = summary.local_fine_scan_interval_count;
    profile->adaptive_candidate_count_before_dedup = summary.candidate_count_before_dedup;
    profile->adaptive_candidate_count_after_dedup = summary.candidate_count_after_dedup;
    profile->adaptive_ternary_accept_count = summary.ternary_accept_count;
    profile->adaptive_ternary_reject_count = summary.ternary_reject_count;
    profile->adaptive_local_fine_scan_root_count = summary.local_fine_scan_root_count;
    profile->adaptive_fallback_to_fullscan = summary.fallback_to_fullscan;
    profile->adaptive_coarse_scan_seconds = summary.coarse_scan_seconds;
    profile->adaptive_interval_collection_seconds = summary.interval_collection_seconds;
    profile->adaptive_interval_refine_seconds = summary.interval_refine_seconds;
    profile->adaptive_local_fine_scan_seconds = summary.local_fine_scan_seconds;
    profile->adaptive_ternary_seconds = summary.ternary_seconds;
    profile->adaptive_bisection_seconds = summary.bisection_seconds;
    profile->adaptive_candidate_dedup_seconds = summary.candidate_dedup_seconds;
    profile->adaptive_sorting_seconds = summary.sorting_seconds;
}

std::string branch_identity_failure_reason(
    const Problem1SolutionBranch& candidate,
    const Problem1SolutionBranch& baseline,
    const Problem1BranchCompareOptions& options
) {
    if (candidate.transfer_revolution != baseline.transfer_revolution) {
        return "transfer_revolution_mismatch";
    }
    if (candidate.target_revolution != baseline.target_revolution) {
        return "target_revolution_mismatch";
    }
    if (wrapped_angle_distance(candidate.encounter_global_angle, baseline.encounter_global_angle) >
        options.angle_tolerance) {
        return "encounter_angle_tolerance_exceeded";
    }
    if (wrapped_angle_distance(candidate.target_arrival_true_anomaly, baseline.target_arrival_true_anomaly) >
        options.target_anomaly_tolerance) {
        return "target_anomaly_tolerance_exceeded";
    }
    const double time_tol = branch_time_identity_tolerance(candidate, baseline, options);
    if (std::abs(candidate.time_of_flight_seconds - baseline.time_of_flight_seconds) > time_tol &&
        std::abs(candidate.target_time_seconds - baseline.target_time_seconds) > time_tol) {
        return "time_tolerance_exceeded";
    }
    return "would_match";
}

const Problem1SolutionBranch* find_branch_for_compare_detail(
    const std::vector<Problem1SolutionBranch>& branches,
    const Problem1BranchCompareDetail& detail
) {
    constexpr double kAngleTol = 1e-10;
    constexpr double kTimeTol = 1e-5;
    for (const Problem1SolutionBranch& branch : branches) {
        if (branch.transfer_revolution == detail.transfer_revolution &&
            branch.target_revolution == detail.target_revolution &&
            wrapped_angle_distance(branch.encounter_global_angle, detail.encounter_global_angle) <= kAngleTol &&
            std::abs(branch.time_of_flight_seconds - detail.time_of_flight_seconds) <= kTimeTol) {
            return &branch;
        }
    }
    return nullptr;
}

void fill_local_residual_sweep_diagnostic(
    const Problem1DepartureAnomalySolveInput& input,
    const Problem1SolutionBranch& branch,
    Problem1AdaptiveExtraBranchDiagnostic* diagnostic
) {
    const double full_scan_step = kTwoPi / static_cast<double>(input.phi_scan_count);
    const double window = 2.0 * full_scan_step;
    const int sample_count = 128;
    const double center_phi = branch.encounter_global_angle;
    const double left_phi = center_phi - window;
    const double right_phi = center_phi + window;

    bool has_previous_valid = false;
    double previous_residual = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i <= sample_count; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(sample_count);
        const double phi = left_phi + (right_phi - left_phi) * t;
        const Problem1DepartureAnomalyResidualResult residual =
            evaluate_problem1_residual_from_departure_anomalies(
                make_departure_anomaly_residual_input(
                    input, phi, branch.transfer_revolution, branch.target_revolution));
        if (!is_success_with_finite_residual(residual)) {
            diagnostic->sweep_invalid_sample_count += 1;
            has_previous_valid = false;
            continue;
        }
        diagnostic->sweep_valid_sample_count += 1;
        const double residual_seconds = scale_free_time_to_seconds(residual.residual);
        const double abs_residual = std::abs(residual_seconds);
        if (!is_finite(diagnostic->sweep_min_abs_residual) ||
            abs_residual < diagnostic->sweep_min_abs_residual) {
            diagnostic->sweep_min_abs_residual = abs_residual;
            diagnostic->sweep_phi_at_min_abs_residual = phi;
        }
        if (has_previous_valid && residual_sign_changed(previous_residual, residual.residual)) {
            diagnostic->sweep_sign_change_exists = true;
        }
        if (i == sample_count / 2 - 1) {
            diagnostic->sweep_left_residual = residual_seconds;
        } else if (i == sample_count / 2) {
            diagnostic->sweep_residual_at_adaptive_phi = residual_seconds;
        } else if (i == sample_count / 2 + 1) {
            diagnostic->sweep_right_residual = residual_seconds;
        }
        previous_residual = residual.residual;
        has_previous_valid = true;
    }
}

void fill_nearest_baseline_diagnostic(
    const Problem1SolutionBranch& branch,
    const std::vector<Problem1SolutionBranch>& baseline_branches,
    Problem1AdaptiveExtraBranchDiagnostic* diagnostic
) {
    const Problem1BranchCompareOptions options{};
    const Problem1SolutionBranch* nearest = nullptr;
    double best_score = std::numeric_limits<double>::infinity();
    for (const Problem1SolutionBranch& baseline : baseline_branches) {
        const double score =
            (branch.transfer_revolution == baseline.transfer_revolution ? 0.0 : 1.0e6) +
            (branch.target_revolution == baseline.target_revolution ? 0.0 : 1.0e6) +
            wrapped_angle_distance(branch.encounter_global_angle, baseline.encounter_global_angle) +
            1e-7 * std::abs(branch.time_of_flight_seconds - baseline.time_of_flight_seconds);
        if (score < best_score) {
            best_score = score;
            nearest = &baseline;
        }
    }
    if (nearest == nullptr) {
        diagnostic->nearest_baseline_match_failure_reason = "no_baseline_branches";
        return;
    }
    diagnostic->nearest_baseline_transfer_revolution = nearest->transfer_revolution;
    diagnostic->nearest_baseline_target_revolution = nearest->target_revolution;
    diagnostic->nearest_baseline_angle_diff =
        wrapped_angle_distance(branch.encounter_global_angle, nearest->encounter_global_angle);
    diagnostic->nearest_baseline_target_anomaly_diff =
        wrapped_angle_distance(branch.target_arrival_true_anomaly, nearest->target_arrival_true_anomaly);
    diagnostic->nearest_baseline_time_diff =
        std::min(
            std::abs(branch.time_of_flight_seconds - nearest->time_of_flight_seconds),
            std::abs(branch.target_time_seconds - nearest->target_time_seconds));
    diagnostic->nearest_baseline_residual_diff = std::abs(branch.residual_seconds - nearest->residual_seconds);
    diagnostic->nearest_baseline_match_failure_reason =
        branch_identity_failure_reason(branch, *nearest, options);
}

bool env_flag_enabled(const char* name) {
    const char* raw = std::getenv(name);
    return raw != nullptr && raw[0] != '\0' && raw[0] != '0';
}

int env_positive_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || raw[0] == '\0') {
        return fallback;
    }
    const int parsed = std::atoi(raw);
    return parsed > 0 ? parsed : fallback;
}

int adaptive_coarse_scan_count_from_environment(int full_scan_count) {
    const int fallback = std::max(16, full_scan_count / 2);
    return std::max(3, env_positive_int("PROBLEM1_ADAPTIVE_COARSE_SCAN_COUNT", fallback));
}

void fill_high_res_oracle_diagnostic(
    Problem1DepartureAnomalySolveInput input,
    const Problem1SolutionBranch& branch,
    Problem1AdaptiveExtraBranchDiagnostic* diagnostic
) {
    if (!env_flag_enabled("PROBLEM1_DEBUG_HIGH_RES_ORACLE")) {
        return;
    }
    input.phi_scan_count = env_positive_int("PROBLEM1_DEBUG_HIGH_RES_SCAN_COUNT", input.phi_scan_count * 4);
    diagnostic->high_res_oracle_enabled = true;
    diagnostic->high_res_scan_count = input.phi_scan_count;
    const Problem1SolveWithDiagnosticResult high_res =
        solve_problem1_from_departure_anomalies_full_scan_internal(input);
    diagnostic->high_res_oracle_branch_count = static_cast<int>(high_res.branches.size());
    std::vector<Problem1SolutionBranch> one_branch{branch};
    const Problem1BranchCompareReport report =
        compare_problem1_solution_branches_by_identity(
            one_branch,
            high_res.branches,
            input.nu_A_depart,
            input.nu_B_depart,
            input.theta_A);
    diagnostic->high_res_oracle_matched = report.extra_count == 0;
}

std::vector<Problem1AdaptiveExtraBranchDiagnostic> make_adaptive_extra_branch_diagnostics(
    const Problem1DepartureAnomalySolveInput& input,
    const std::vector<Problem1SolutionBranch>& adaptive_branches,
    const std::vector<Problem1SolutionBranch>& baseline_branches,
    const Problem1BranchCompareReport& compare_report
) {
    std::vector<Problem1AdaptiveExtraBranchDiagnostic> diagnostics;
    diagnostics.reserve(compare_report.extra_branches.size());
    for (const Problem1BranchCompareDetail& extra : compare_report.extra_branches) {
        const Problem1SolutionBranch* branch = find_branch_for_compare_detail(adaptive_branches, extra);
        if (branch == nullptr) {
            continue;
        }
        Problem1AdaptiveExtraBranchDiagnostic diagnostic{};
        diagnostic.branch = extra;
        diagnostic.root_bracket_width = branch->root_bracket_width;
        diagnostic.bisection_iterations = branch->bisection_iterations;
        diagnostic.adaptive_source_reason = branch->adaptive_source_reason;
        diagnostic.ternary_phi_star = branch->adaptive_ternary_phi_star;
        diagnostic.ternary_residual = branch->adaptive_ternary_residual;
        diagnostic.ternary_residual_sq = branch->adaptive_ternary_residual_sq;
        diagnostic.local_interval_left_phi = branch->adaptive_local_interval_left_phi;
        diagnostic.local_interval_right_phi = branch->adaptive_local_interval_right_phi;
        diagnostic.local_subdivision_index = branch->adaptive_local_subdivision_index;
        fill_local_residual_sweep_diagnostic(input, *branch, &diagnostic);
        fill_nearest_baseline_diagnostic(*branch, baseline_branches, &diagnostic);
        fill_high_res_oracle_diagnostic(input, *branch, &diagnostic);
        diagnostics.push_back(std::move(diagnostic));
    }
    return diagnostics;
}

Problem1DepartureAnomalySolveInput make_departure_anomaly_solve_input_from_defaults(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int max_transfer_revolution,
    int max_target_revolution
) {
    const config::Problem1SolveDefaults& defaults = config::global_config().problem1_solve;
    return Problem1DepartureAnomalySolveInput{
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        max_transfer_revolution,
        max_target_revolution,
        defaults.phi_scan_count,
        defaults.phi_tolerance,
        defaults.residual_tolerance,
        defaults.max_bisection_iterations,
        defaults.max_candidate_relative_residual,
    };
}

std::optional<Problem1SolutionBranchCandidate> make_branch_candidate_from_residual_result(
    double encounter_global_angle,
    const Problem1DepartureAnomalyResidualResult& residual_result,
    int transfer_revolution,
    int target_revolution,
    double root_bracket_width,
    double max_candidate_relative_residual
) {
    if (!is_success_with_finite_residual(residual_result)) {
        return std::nullopt;
    }

    const double time_of_flight_seconds = scale_free_time_to_seconds(residual_result.transfer_time_scale_free);
    const double target_time_seconds = scale_free_time_to_seconds(residual_result.target_time_scale_free);
    const double residual_seconds = scale_free_time_to_seconds(residual_result.residual);
    if (!is_finite(time_of_flight_seconds) || !(time_of_flight_seconds > 0.0) ||
        !is_finite(target_time_seconds) || !(target_time_seconds > 0.0) ||
        !is_finite(residual_seconds)) {
        return std::nullopt;
    }

    const double relative_residual = compute_problem1_relative_residual(residual_result);
    if (!is_finite(relative_residual) || relative_residual > max_candidate_relative_residual ||
        !is_finite(root_bracket_width) || root_bracket_width < 0.0) {
        return std::nullopt;
    }

    Problem1SolutionBranch branch{};
    branch.valid = true;
    branch.encounter_global_angle = normalize_angle_0_2pi(encounter_global_angle);
    branch.target_arrival_true_anomaly = residual_result.target_arrival_true_anomaly;
    branch.transfer_revolution = transfer_revolution;
    branch.target_revolution = target_revolution;
    branch.time_of_flight_seconds = time_of_flight_seconds;
    branch.target_time_seconds = target_time_seconds;
    branch.residual_seconds = residual_seconds;
    branch.arrival_time_seconds_since_j2000 = std::numeric_limits<double>::quiet_NaN();
    branch.transfer_e = residual_result.transfer_e;
    branch.transfer_p = residual_result.transfer_p;
    branch.transfer_a = compute_transfer_a(residual_result.transfer_p, residual_result.transfer_e);
    branch.theta_B = residual_result.theta_B;
    branch.root_bracket_width = root_bracket_width;

    return Problem1SolutionBranchCandidate{
        branch,
        relative_residual,
        root_bracket_width,
    };
}

std::optional<Problem1SolutionBranchCandidate> refine_problem1_root_by_bisection_from_departure_anomalies(
    const Problem1DepartureAnomalySolveInput& input,
    double left_phi,
    const Problem1DepartureAnomalyResidualResult& left_result,
    double right_phi,
    const Problem1DepartureAnomalyResidualResult& right_result,
    int transfer_revolution,
    int target_revolution,
    Problem1SolveDiagnostic* diagnostic = nullptr
) {
    if (!is_success_with_finite_residual(left_result) || !is_success_with_finite_residual(right_result)) {
        return std::nullopt;
    }

    if (std::abs(left_result.residual) <= input.residual_tolerance) {
        return make_branch_candidate_from_residual_result(
            left_phi,
            left_result,
            transfer_revolution,
            target_revolution,
            std::abs(right_phi - left_phi),
            input.max_candidate_relative_residual);
    }
    if (std::abs(right_result.residual) <= input.residual_tolerance) {
        return make_branch_candidate_from_residual_result(
            right_phi,
            right_result,
            transfer_revolution,
            target_revolution,
            std::abs(right_phi - left_phi),
            input.max_candidate_relative_residual);
    }
    if (!residual_sign_changed(left_result.residual, right_result.residual)) {
        return std::nullopt;
    }

    double current_left_phi = left_phi;
    double current_right_phi = right_phi;
    Problem1DepartureAnomalyResidualResult current_left_result = left_result;
    Problem1DepartureAnomalyResidualResult current_right_result = right_result;
    double mid_phi = 0.5 * (current_left_phi + current_right_phi);
    Problem1DepartureAnomalyResidualResult mid_result = make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);

    for (int iteration = 0; iteration < input.max_bisection_iterations; ++iteration) {
        mid_phi = 0.5 * (current_left_phi + current_right_phi);
        mid_result = evaluate_problem1_residual_from_departure_anomalies_profiled(
            make_departure_anomaly_residual_input(input, mid_phi, transfer_revolution, target_revolution),
            diagnostic);
        if (!is_success_with_finite_residual(mid_result)) {
            return std::nullopt;
        }
        if (std::abs(mid_result.residual) <= input.residual_tolerance ||
            std::abs(current_right_phi - current_left_phi) <= input.phi_tolerance) {
            return make_branch_candidate_from_residual_result(
                mid_phi,
                mid_result,
                transfer_revolution,
                target_revolution,
                std::abs(current_right_phi - current_left_phi),
                input.max_candidate_relative_residual);
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

    return make_branch_candidate_from_residual_result(
        mid_phi,
        mid_result,
        transfer_revolution,
        target_revolution,
        std::abs(current_right_phi - current_left_phi),
        input.max_candidate_relative_residual);
}

void add_branch_if_not_duplicate(
    std::vector<Problem1SolutionBranchCandidate>* branches,
    Problem1SolutionBranchCandidate candidate
) {
    for (Problem1SolutionBranchCandidate& existing : *branches) {
        if (existing.branch.transfer_revolution != candidate.branch.transfer_revolution ||
            existing.branch.target_revolution != candidate.branch.target_revolution) {
            continue;
        }
        const double angle_diff = normalize_angle_minus_pi_pi(
            existing.branch.encounter_global_angle - candidate.branch.encounter_global_angle);
        if (std::abs(angle_diff) < 1e-8) {
            if (candidate.relative_residual < existing.relative_residual ||
                (candidate.relative_residual == existing.relative_residual &&
                 (candidate.root_bracket_width < existing.root_bracket_width ||
                  (candidate.root_bracket_width == existing.root_bracket_width &&
                   std::abs(candidate.branch.residual_seconds) < std::abs(existing.branch.residual_seconds))))) {
                existing = std::move(candidate);
            }
            return;
        }
    }
    branches->push_back(std::move(candidate));
}

void finalize_problem1_branch_candidates(
    const std::vector<Problem1SolutionBranchCandidate>& branch_candidates,
    std::vector<Problem1SolutionBranch>* branches,
    Problem1SolveDiagnostic* diagnostic
) {
    using clock = std::chrono::steady_clock;
    const auto sorting_start = clock::now();
    branches->reserve(branch_candidates.size());
    for (const Problem1SolutionBranchCandidate& candidate : branch_candidates) {
        branches->push_back(candidate.branch);
    }

    std::sort(
        branches->begin(),
        branches->end(),
        [](const Problem1SolutionBranch& lhs, const Problem1SolutionBranch& rhs) {
            if (lhs.time_of_flight_seconds < rhs.time_of_flight_seconds) {
                return true;
            }
            if (lhs.time_of_flight_seconds > rhs.time_of_flight_seconds) {
                return false;
            }
            if (lhs.transfer_revolution < rhs.transfer_revolution) {
                return true;
            }
            if (lhs.transfer_revolution > rhs.transfer_revolution) {
                return false;
            }
            if (lhs.target_revolution < rhs.target_revolution) {
                return true;
            }
            if (lhs.target_revolution > rhs.target_revolution) {
                return false;
            }
            return lhs.encounter_global_angle < rhs.encounter_global_angle;
        });
    const auto sorting_end = clock::now();
    if (diagnostic != nullptr) {
        diagnostic->sorting_conversion_seconds =
            std::chrono::duration<double>(sorting_end - sorting_start).count();
        diagnostic->final_branch_count = static_cast<int>(branches->size());
    }
}

void increment_adaptive_interval_reason(
    Problem1SuspiciousIntervalReason reason,
    Problem1AdaptiveScanSummary* summary
) {
    if (summary == nullptr) {
        return;
    }
    summary->interval_total += 1;
    switch (reason) {
    case Problem1SuspiciousIntervalReason::SignChange:
        summary->sign_change_interval_count += 1;
        break;
    case Problem1SuspiciousIntervalReason::NearZero:
        summary->near_zero_interval_count += 1;
        break;
    case Problem1SuspiciousIntervalReason::LocalMinimum:
        summary->local_min_interval_count += 1;
        break;
    case Problem1SuspiciousIntervalReason::LocalMaximum:
        summary->local_max_interval_count += 1;
        break;
    case Problem1SuspiciousIntervalReason::ValidBoundary:
        summary->valid_boundary_interval_count += 1;
        break;
    case Problem1SuspiciousIntervalReason::WrapSignChange:
        summary->wrap_interval_count += 1;
        break;
    case Problem1SuspiciousIntervalReason::RapidChange:
        summary->rapid_change_interval_count += 1;
        break;
    }
}

bool duplicate_adaptive_interval(
    const Problem1SuspiciousInterval& lhs,
    const Problem1SuspiciousInterval& rhs
) {
    if (lhs.left_index == rhs.left_index && lhs.right_index == rhs.right_index) {
        return true;
    }
    constexpr double kPhiIntervalDuplicateTolerance = 1e-12;
    return std::abs(lhs.left_phi - rhs.left_phi) <= kPhiIntervalDuplicateTolerance &&
        std::abs(lhs.right_phi - rhs.right_phi) <= kPhiIntervalDuplicateTolerance;
}

bool mark_problem1_adaptive_interval_processed(
    std::vector<std::pair<int, int>>* processed_intervals,
    const Problem1SuspiciousInterval& interval,
    Problem1AdaptiveScanSummary* summary
) {
    const std::pair<int, int> key{interval.left_index, interval.right_index};
    for (const auto& existing : *processed_intervals) {
        if (existing == key) {
            return false;
        }
    }
    processed_intervals->push_back(key);
    increment_adaptive_interval_reason(interval.reason, summary);
    return true;
}

void add_problem1_adaptive_scan_interval_if_new(
    std::vector<Problem1SuspiciousInterval>* intervals,
    Problem1SuspiciousInterval interval,
    Problem1AdaptiveScanSummary* summary
) {
    for (const Problem1SuspiciousInterval& existing : *intervals) {
        if (duplicate_adaptive_interval(existing, interval)) {
            return;
        }
    }
    increment_adaptive_interval_reason(interval.reason, summary);
    intervals->push_back(interval);
}

Problem1AdaptiveScanPoint make_problem1_adaptive_scan_point(
    const Problem1DepartureAnomalySolveInput& input,
    int transfer_revolution,
    int target_revolution,
    int index,
    int coarse_scan_count,
    Problem1AdaptiveScanSummary* summary
) {
    Problem1AdaptiveScanPoint point{};
    point.phi = static_cast<double>(index) * kTwoPi / static_cast<double>(coarse_scan_count);
    point.residual_result = evaluate_problem1_residual_from_departure_anomalies(
        make_departure_anomaly_residual_input(input, point.phi, transfer_revolution, target_revolution));
    if (summary != nullptr) {
        summary->coarse_scan_count += 1;
        summary->residual_eval_count += 1;
    }
    point.valid = is_success_with_finite_residual(point.residual_result);
    if (point.valid) {
        point.residual = point.residual_result.residual;
        point.abs_residual = std::abs(point.residual);
        point.residual_sq = point.residual * point.residual;
    }
    return point;
}

void collect_problem1_adaptive_scan_intervals(
    const Problem1DepartureAnomalySolveInput& input,
    int transfer_revolution,
    int target_revolution,
    int coarse_scan_count,
    double near_zero_tol,
    double rapid_change_tol,
    std::vector<Problem1SuspiciousInterval>* intervals,
    Problem1AdaptiveScanSummary* summary
) {
    if (coarse_scan_count < 3) {
        return;
    }

    std::vector<Problem1AdaptiveScanPoint> points;
    points.reserve(static_cast<std::size_t>(coarse_scan_count));
    for (int i = 0; i < coarse_scan_count; ++i) {
        points.push_back(make_problem1_adaptive_scan_point(
            input, transfer_revolution, target_revolution, i, coarse_scan_count, summary));
    }

    const double phi_step = kTwoPi / static_cast<double>(coarse_scan_count);
    for (int i = 0; i < coarse_scan_count; ++i) {
        const Problem1AdaptiveScanPoint& point = points[static_cast<std::size_t>(i)];
        if (point.valid && is_finite(near_zero_tol) && near_zero_tol >= 0.0 &&
            point.abs_residual <= near_zero_tol) {
            add_problem1_adaptive_scan_interval_if_new(
                intervals,
                Problem1SuspiciousInterval{
                    std::max(0.0, point.phi - phi_step),
                    std::min(kTwoPi, point.phi + phi_step),
                    Problem1SuspiciousIntervalReason::NearZero,
                    i,
                    i,
                    point.residual,
                    point.residual,
                    point.residual,
                },
                summary);
        }

        if (i == 0 || i + 1 >= coarse_scan_count) {
            continue;
        }
        const Problem1AdaptiveScanPoint& left = points[static_cast<std::size_t>(i - 1)];
        const Problem1AdaptiveScanPoint& right = points[static_cast<std::size_t>(i + 1)];
        if (left.valid && point.valid && right.valid &&
            point.residual <= left.residual &&
            point.residual <= right.residual) {
            add_problem1_adaptive_scan_interval_if_new(
                intervals,
                Problem1SuspiciousInterval{
                    left.phi,
                    right.phi,
                    Problem1SuspiciousIntervalReason::LocalMinimum,
                    i - 1,
                    i + 1,
                    left.residual,
                    right.residual,
                    point.residual,
                },
                summary);
        }
        if (left.valid && point.valid && right.valid &&
            point.residual >= left.residual &&
            point.residual >= right.residual) {
            add_problem1_adaptive_scan_interval_if_new(
                intervals,
                Problem1SuspiciousInterval{
                    left.phi,
                    right.phi,
                    Problem1SuspiciousIntervalReason::LocalMaximum,
                    i - 1,
                    i + 1,
                    left.residual,
                    right.residual,
                    point.residual,
                },
                summary);
        }
    }

    for (int i = 0; i + 1 < coarse_scan_count; ++i) {
        const Problem1AdaptiveScanPoint& left = points[static_cast<std::size_t>(i)];
        const Problem1AdaptiveScanPoint& right = points[static_cast<std::size_t>(i + 1)];
        if (left.valid != right.valid) {
            add_problem1_adaptive_scan_interval_if_new(
                intervals,
                Problem1SuspiciousInterval{
                    left.phi,
                    right.phi,
                    Problem1SuspiciousIntervalReason::ValidBoundary,
                    i,
                    i + 1,
                    left.residual,
                    right.residual,
                    std::numeric_limits<double>::quiet_NaN(),
                },
                summary);
            continue;
        }
        if (!left.valid || !right.valid) {
            continue;
        }
        if (residual_sign_changed(left.residual, right.residual)) {
            add_problem1_adaptive_scan_interval_if_new(
                intervals,
                Problem1SuspiciousInterval{
                    left.phi,
                    right.phi,
                    Problem1SuspiciousIntervalReason::SignChange,
                    i,
                    i + 1,
                    left.residual,
                    right.residual,
                    std::numeric_limits<double>::quiet_NaN(),
                },
                summary);
        }
        if (is_finite(rapid_change_tol) && rapid_change_tol > 0.0 &&
            std::abs(right.residual - left.residual) > rapid_change_tol) {
            add_problem1_adaptive_scan_interval_if_new(
                intervals,
                Problem1SuspiciousInterval{
                    left.phi,
                    right.phi,
                    Problem1SuspiciousIntervalReason::RapidChange,
                    i,
                    i + 1,
                    left.residual,
                    right.residual,
                    std::numeric_limits<double>::quiet_NaN(),
                },
                summary);
        }
    }

    const Problem1AdaptiveScanPoint& last = points.back();
    const Problem1AdaptiveScanPoint& first = points.front();
    if (last.valid != first.valid) {
        add_problem1_adaptive_scan_interval_if_new(
            intervals,
            Problem1SuspiciousInterval{
                last.phi,
                first.phi + kTwoPi,
                Problem1SuspiciousIntervalReason::ValidBoundary,
                coarse_scan_count - 1,
                0,
                last.residual,
                first.residual,
                std::numeric_limits<double>::quiet_NaN(),
            },
            summary);
    } else if (last.valid && first.valid && residual_sign_changed(last.residual, first.residual)) {
        add_problem1_adaptive_scan_interval_if_new(
            intervals,
            Problem1SuspiciousInterval{
                last.phi,
                first.phi + kTwoPi,
                Problem1SuspiciousIntervalReason::WrapSignChange,
                coarse_scan_count - 1,
                0,
                last.residual,
                first.residual,
                std::numeric_limits<double>::quiet_NaN(),
            },
            summary);
    }
}

Problem1DepartureAnomalyResidualResult evaluate_problem1_adaptive_residual_profiled(
    const Problem1DepartureAnomalySolveInput& input,
    double phi,
    int transfer_revolution,
    int target_revolution,
    Problem1SolveDiagnostic* diagnostic
) {
    return evaluate_problem1_residual_from_departure_anomalies_profiled(
        make_departure_anomaly_residual_input(input, phi, transfer_revolution, target_revolution),
        diagnostic);
}

double evaluate_problem1_adaptive_residual_objective(
    const Problem1DepartureAnomalySolveInput& input,
    double phi,
    int transfer_revolution,
    int target_revolution,
    Problem1SolveDiagnostic* diagnostic,
    bool find_minimum
) {
    const Problem1DepartureAnomalyResidualResult result =
        evaluate_problem1_adaptive_residual_profiled(input, phi, transfer_revolution, target_revolution, diagnostic);
    if (!is_success_with_finite_residual(result)) {
        return std::numeric_limits<double>::infinity();
    }
    return find_minimum ? result.residual : -result.residual;
}

std::optional<Problem1SolutionBranchCandidate> try_make_adaptive_candidate(
    const Problem1DepartureAnomalySolveInput& input,
    double phi,
    const Problem1DepartureAnomalyResidualResult& residual_result,
    int transfer_revolution,
    int target_revolution,
    double root_bracket_width,
    const Problem1SuspiciousInterval* interval,
    const char* source_reason,
    double ternary_phi_star,
    double ternary_residual,
    int local_subdivision_index,
    Problem1AdaptiveScanSummary* summary
) {
    auto candidate = make_branch_candidate_from_residual_result(
        phi,
        residual_result,
        transfer_revolution,
        target_revolution,
        root_bracket_width,
        input.max_candidate_relative_residual);
    if (candidate.has_value() && summary != nullptr) {
        summary->candidate_count_before_dedup += 1;
    }
    if (candidate.has_value()) {
        candidate->branch.adaptive_source_reason = source_reason != nullptr ? source_reason : "";
        candidate->branch.adaptive_ternary_phi_star = ternary_phi_star;
        candidate->branch.adaptive_ternary_residual = ternary_residual;
        candidate->branch.adaptive_ternary_residual_sq =
            is_finite(ternary_residual) ? ternary_residual * ternary_residual : std::numeric_limits<double>::quiet_NaN();
        if (interval != nullptr) {
            candidate->branch.adaptive_local_interval_left_phi = interval->left_phi;
            candidate->branch.adaptive_local_interval_right_phi = interval->right_phi;
        }
        candidate->branch.adaptive_local_subdivision_index = local_subdivision_index;
    }
    return candidate;
}

Problem1DepartureAnomalyResidualResult ternary_extremize_problem1_residual(
    const Problem1DepartureAnomalySolveInput& input,
    const Problem1SuspiciousInterval& interval,
    int transfer_revolution,
    int target_revolution,
    int iteration_count,
    Problem1SolveDiagnostic* diagnostic,
    Problem1AdaptiveScanSummary* summary,
    bool find_minimum,
    double* phi_ext
) {
    ScopedSecondsAccumulator timer(summary != nullptr ? &summary->ternary_seconds : nullptr);
    if (summary != nullptr) {
        summary->ternary_interval_count += 1;
    }
    double left = interval.left_phi;
    double right = interval.right_phi;
    if (!is_finite(left) || !is_finite(right) || !(right > left)) {
        if (summary != nullptr) {
            summary->ternary_reject_count += 1;
        }
        if (phi_ext != nullptr) {
            *phi_ext = std::numeric_limits<double>::quiet_NaN();
        }
        return make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);
    }

    for (int iteration = 0; iteration < iteration_count; ++iteration) {
        const double third = (right - left) / 3.0;
        const double m1 = left + third;
        const double m2 = right - third;
        const double f1 =
            evaluate_problem1_adaptive_residual_objective(
                input, m1, transfer_revolution, target_revolution, diagnostic, find_minimum);
        const double f2 =
            evaluate_problem1_adaptive_residual_objective(
                input, m2, transfer_revolution, target_revolution, diagnostic, find_minimum);
        if (f1 < f2) {
            right = m2;
        } else {
            left = m1;
        }
    }

    const double phi = 0.5 * (left + right);
    if (phi_ext != nullptr) {
        *phi_ext = phi;
    }
    return evaluate_problem1_adaptive_residual_profiled(input, phi, transfer_revolution, target_revolution, diagnostic);
}

void local_fine_scan_problem1_interval(
    const Problem1DepartureAnomalySolveInput& input,
    const Problem1SuspiciousInterval& interval,
    int transfer_revolution,
    int target_revolution,
    int local_subdivision_count,
    std::vector<Problem1SolutionBranchCandidate>* branch_candidates,
    Problem1SolveDiagnostic* diagnostic,
    Problem1AdaptiveScanSummary* summary,
    const Problem1AdaptiveScanPoint* cached_left = nullptr,
    const Problem1AdaptiveScanPoint* cached_right = nullptr
) {
    ScopedSecondsAccumulator timer(summary != nullptr ? &summary->local_fine_scan_seconds : nullptr);
    if (summary != nullptr) {
        summary->local_fine_scan_interval_count += 1;
    }
    if (local_subdivision_count <= 0 ||
        !is_finite(interval.left_phi) ||
        !is_finite(interval.right_phi) ||
        !(interval.right_phi > interval.left_phi)) {
        return;
    }

    bool has_previous_valid = false;
    double previous_phi = interval.left_phi;
    Problem1DepartureAnomalyResidualResult previous_result =
        make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);

    for (int i = 0; i <= local_subdivision_count; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(local_subdivision_count);
        const double phi = interval.left_phi + (interval.right_phi - interval.left_phi) * t;
        Problem1DepartureAnomalyResidualResult result;
        if (i == 0 && cached_left != nullptr) {
            result = cached_left->residual_result;
        } else if (i == local_subdivision_count && cached_right != nullptr) {
            result = cached_right->residual_result;
        } else {
            result = evaluate_problem1_adaptive_residual_profiled(
                input, phi, transfer_revolution, target_revolution, diagnostic);
        }
        if (!is_success_with_finite_residual(result)) {
            has_previous_valid = false;
            continue;
        }

        if (std::abs(result.residual) <= input.residual_tolerance) {
            auto candidate = try_make_adaptive_candidate(
                input,
                phi,
                result,
                transfer_revolution,
                target_revolution,
                0.0,
                &interval,
                suspicious_interval_reason_name(interval.reason),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                i,
                summary);
            if (candidate.has_value()) {
                add_branch_if_not_duplicate(branch_candidates, std::move(*candidate));
                if (summary != nullptr) {
                    summary->local_fine_scan_root_count += 1;
                }
            }
        }

        if (has_previous_valid && residual_sign_changed(previous_result.residual, result.residual)) {
            std::optional<Problem1SolutionBranchCandidate> refined;
            {
                ScopedSecondsAccumulator bisection_timer(
                    summary != nullptr ? &summary->bisection_seconds : nullptr);
                refined = refine_problem1_root_by_bisection_from_departure_anomalies(
                    input,
                    previous_phi,
                    previous_result,
                    phi,
                    result,
                    transfer_revolution,
                    target_revolution,
                    diagnostic);
            }
            if (refined.has_value()) {
                if (summary != nullptr) {
                    summary->candidate_count_before_dedup += 1;
                    summary->local_fine_scan_root_count += 1;
                }
                refined->branch.adaptive_source_reason = suspicious_interval_reason_name(interval.reason);
                refined->branch.adaptive_local_interval_left_phi = interval.left_phi;
                refined->branch.adaptive_local_interval_right_phi = interval.right_phi;
                refined->branch.adaptive_local_subdivision_index = i;
                add_branch_if_not_duplicate(branch_candidates, std::move(*refined));
            }
        }

        previous_phi = phi;
        previous_result = result;
        has_previous_valid = true;
    }
}

int refine_problem1_adaptive_extremum_subintervals(
    const Problem1DepartureAnomalySolveInput& input,
    const Problem1SuspiciousInterval& interval,
    double phi_ext,
    const Problem1DepartureAnomalyResidualResult& ext_result,
    int transfer_revolution,
    int target_revolution,
    std::vector<Problem1SolutionBranchCandidate>* branch_candidates,
    Problem1SolveDiagnostic* diagnostic,
    Problem1AdaptiveScanSummary* summary,
    const Problem1AdaptiveScanPoint* cached_left = nullptr,
    const Problem1AdaptiveScanPoint* cached_right = nullptr
) {
    if (!is_success_with_finite_residual(ext_result)) {
        return 0;
    }
    const Problem1DepartureAnomalyResidualResult left_result =
        cached_left != nullptr
            ? cached_left->residual_result
            : evaluate_problem1_adaptive_residual_profiled(
                  input, interval.left_phi, transfer_revolution, target_revolution, diagnostic);
    const Problem1DepartureAnomalyResidualResult right_result =
        cached_right != nullptr
            ? cached_right->residual_result
            : evaluate_problem1_adaptive_residual_profiled(
                  input, interval.right_phi, transfer_revolution, target_revolution, diagnostic);

    int accepted_count = 0;
    auto try_refine = [&](double left_phi,
                          const Problem1DepartureAnomalyResidualResult& left,
                          double right_phi,
                          const Problem1DepartureAnomalyResidualResult& right) {
        if (!is_success_with_finite_residual(left) ||
            !is_success_with_finite_residual(right) ||
            !residual_sign_changed(left.residual, right.residual)) {
            return;
        }
        ScopedSecondsAccumulator bisection_timer(
            summary != nullptr ? &summary->bisection_seconds : nullptr);
        auto refined = refine_problem1_root_by_bisection_from_departure_anomalies(
            input,
            left_phi,
            left,
            right_phi,
            right,
            transfer_revolution,
            target_revolution,
            diagnostic);
        if (!refined.has_value()) {
            return;
        }
        if (summary != nullptr) {
            summary->candidate_count_before_dedup += 1;
        }
        refined->branch.adaptive_source_reason = suspicious_interval_reason_name(interval.reason);
        refined->branch.adaptive_ternary_phi_star = phi_ext;
        refined->branch.adaptive_ternary_residual = ext_result.residual;
        refined->branch.adaptive_ternary_residual_sq = ext_result.residual * ext_result.residual;
        refined->branch.adaptive_local_interval_left_phi = interval.left_phi;
        refined->branch.adaptive_local_interval_right_phi = interval.right_phi;
        add_branch_if_not_duplicate(branch_candidates, std::move(*refined));
        accepted_count += 1;
    };

    try_refine(interval.left_phi, left_result, phi_ext, ext_result);
    try_refine(phi_ext, ext_result, interval.right_phi, right_result);
    return accepted_count;
}

void refine_problem1_adaptive_interval(
    const Problem1DepartureAnomalySolveInput& input,
    const Problem1SuspiciousInterval& interval,
    int transfer_revolution,
    int target_revolution,
    int local_subdivision_count,
    std::vector<Problem1SolutionBranchCandidate>* branch_candidates,
    Problem1SolveDiagnostic* diagnostic,
    Problem1AdaptiveScanSummary* summary,
    const Problem1AdaptiveScanPoint* cached_left = nullptr,
    const Problem1AdaptiveScanPoint* cached_mid = nullptr,
    const Problem1AdaptiveScanPoint* cached_right = nullptr
) {
    (void)cached_mid;
    if (summary != nullptr) {
        summary->refined_interval_count += 1;
    }
    ScopedSecondsAccumulator refine_timer(summary != nullptr ? &summary->interval_refine_seconds : nullptr);
    switch (interval.reason) {
    case Problem1SuspiciousIntervalReason::SignChange:
    case Problem1SuspiciousIntervalReason::WrapSignChange: {
        if (summary != nullptr) {
            summary->bisection_interval_count += 1;
        }
        const Problem1DepartureAnomalyResidualResult left_result =
            cached_left != nullptr
                ? cached_left->residual_result
                : evaluate_problem1_adaptive_residual_profiled(
                      input, interval.left_phi, transfer_revolution, target_revolution, diagnostic);
        const Problem1DepartureAnomalyResidualResult right_result =
            cached_right != nullptr
                ? cached_right->residual_result
                : evaluate_problem1_adaptive_residual_profiled(
                      input, interval.right_phi, transfer_revolution, target_revolution, diagnostic);
        std::optional<Problem1SolutionBranchCandidate> refined;
        {
            ScopedSecondsAccumulator bisection_timer(
                summary != nullptr ? &summary->bisection_seconds : nullptr);
            refined = refine_problem1_root_by_bisection_from_departure_anomalies(
                input,
                interval.left_phi,
                left_result,
                interval.right_phi,
                right_result,
                transfer_revolution,
                target_revolution,
                diagnostic);
        }
        if (refined.has_value()) {
            if (summary != nullptr) {
                summary->candidate_count_before_dedup += 1;
            }
            refined->branch.adaptive_source_reason = suspicious_interval_reason_name(interval.reason);
            refined->branch.adaptive_local_interval_left_phi = interval.left_phi;
            refined->branch.adaptive_local_interval_right_phi = interval.right_phi;
            add_branch_if_not_duplicate(branch_candidates, std::move(*refined));
        }
        break;
    }
    case Problem1SuspiciousIntervalReason::LocalMinimum:
    case Problem1SuspiciousIntervalReason::LocalMaximum: {
        double phi_ext = std::numeric_limits<double>::quiet_NaN();
        const bool find_minimum = interval.reason == Problem1SuspiciousIntervalReason::LocalMinimum;
        const Problem1DepartureAnomalyResidualResult ext_result =
            ternary_extremize_problem1_residual(
            input,
            interval,
            transfer_revolution,
            target_revolution,
            40,
            diagnostic,
            summary,
            find_minimum,
            &phi_ext);
        const int accepted_count = refine_problem1_adaptive_extremum_subintervals(
            input,
            interval,
            phi_ext,
            ext_result,
            transfer_revolution,
            target_revolution,
            branch_candidates,
            diagnostic,
            summary,
            cached_left,
            cached_right);
        if (summary != nullptr) {
            if (accepted_count > 0) {
                summary->ternary_accept_count += accepted_count;
            } else {
                summary->ternary_reject_count += 1;
            }
        }
        if (accepted_count == 0) {
            local_fine_scan_problem1_interval(
                input,
                interval,
                transfer_revolution,
                target_revolution,
                local_subdivision_count,
                branch_candidates,
                diagnostic,
                summary,
                cached_left,
                cached_right);
        }
        break;
    }
    case Problem1SuspiciousIntervalReason::ValidBoundary:
    case Problem1SuspiciousIntervalReason::NearZero:
    case Problem1SuspiciousIntervalReason::RapidChange:
        local_fine_scan_problem1_interval(
            input,
            interval,
            transfer_revolution,
            target_revolution,
            local_subdivision_count,
            branch_candidates,
            diagnostic,
            summary,
            cached_left,
            cached_right);
        break;
    }
}

void solve_problem1_adaptive_streaming_for_revolutions(
    const Problem1DepartureAnomalySolveInput& input,
    int transfer_revolution,
    int target_revolution,
    int coarse_scan_count,
    double near_zero_tol,
    double rapid_change_tol,
    int local_subdivision_count,
    std::vector<Problem1SolutionBranchCandidate>* branch_candidates,
    Problem1SolveDiagnostic* diagnostic,
    Problem1AdaptiveScanSummary* summary
) {
    if (coarse_scan_count < 3) {
        return;
    }

    const double phi_step = kTwoPi / static_cast<double>(coarse_scan_count);
    std::vector<std::pair<int, int>> processed_intervals;
    processed_intervals.reserve(static_cast<std::size_t>(std::max(8, coarse_scan_count / 4)));

    std::optional<Problem1AdaptiveScanPoint> first_point;
    std::optional<Problem1AdaptiveScanPoint> previous2;
    std::optional<Problem1AdaptiveScanPoint> previous1;
    const auto streaming_start = std::chrono::steady_clock::now();
    const double interval_refine_seconds_before =
        summary != nullptr ? summary->interval_refine_seconds : 0.0;

    auto process_interval = [&](const Problem1SuspiciousInterval& interval,
                                const Problem1AdaptiveScanPoint* left,
                                const Problem1AdaptiveScanPoint* mid,
                                const Problem1AdaptiveScanPoint* right) {
        if (!mark_problem1_adaptive_interval_processed(&processed_intervals, interval, summary)) {
            return;
        }
        refine_problem1_adaptive_interval(
            input,
            interval,
            transfer_revolution,
            target_revolution,
            local_subdivision_count,
            branch_candidates,
            diagnostic,
            summary,
            left,
            mid,
            right);
    };

    auto process_near_zero = [&](const Problem1AdaptiveScanPoint& point) {
        if (!point.valid || !is_finite(near_zero_tol) || near_zero_tol < 0.0 ||
            point.abs_residual > near_zero_tol) {
            return;
        }
        const Problem1SuspiciousInterval interval{
            std::max(0.0, point.phi - phi_step),
            std::min(kTwoPi, point.phi + phi_step),
            Problem1SuspiciousIntervalReason::NearZero,
            static_cast<int>(std::llround(point.phi / phi_step)),
            static_cast<int>(std::llround(point.phi / phi_step)),
            point.residual,
            point.residual,
            point.residual,
        };
        if (!mark_problem1_adaptive_interval_processed(&processed_intervals, interval, summary)) {
            return;
        }
        if (summary != nullptr) {
            summary->refined_interval_count += 1;
        }
        auto candidate = try_make_adaptive_candidate(
            input,
            point.phi,
            point.residual_result,
            transfer_revolution,
            target_revolution,
            0.0,
            &interval,
            suspicious_interval_reason_name(interval.reason),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            -1,
            summary);
        if (candidate.has_value()) {
            add_branch_if_not_duplicate(branch_candidates, std::move(*candidate));
        }
    };

    for (int i = 0; i < coarse_scan_count; ++i) {
        Problem1AdaptiveScanPoint current = make_problem1_adaptive_scan_point(
            input, transfer_revolution, target_revolution, i, coarse_scan_count, summary);
            if (!first_point.has_value()) {
                first_point = current;
            }

            process_near_zero(current);

            if (previous1.has_value()) {
                const Problem1AdaptiveScanPoint& left = *previous1;
                if (left.valid != current.valid) {
                    process_interval(
                        Problem1SuspiciousInterval{
                            left.phi,
                            current.phi,
                            Problem1SuspiciousIntervalReason::ValidBoundary,
                            i - 1,
                            i,
                            left.residual,
                            current.residual,
                            std::numeric_limits<double>::quiet_NaN(),
                        },
                        &left,
                        nullptr,
                        &current);
                } else if (left.valid && current.valid) {
                    if (residual_sign_changed(left.residual, current.residual)) {
                        process_interval(
                            Problem1SuspiciousInterval{
                                left.phi,
                                current.phi,
                                Problem1SuspiciousIntervalReason::SignChange,
                                i - 1,
                                i,
                                left.residual,
                                current.residual,
                                std::numeric_limits<double>::quiet_NaN(),
                            },
                            &left,
                            nullptr,
                            &current);
                    }
                    if (is_finite(rapid_change_tol) && rapid_change_tol > 0.0 &&
                        std::abs(current.residual - left.residual) > rapid_change_tol) {
                        process_interval(
                            Problem1SuspiciousInterval{
                                left.phi,
                                current.phi,
                                Problem1SuspiciousIntervalReason::RapidChange,
                                i - 1,
                                i,
                                left.residual,
                                current.residual,
                                std::numeric_limits<double>::quiet_NaN(),
                            },
                            &left,
                            nullptr,
                            &current);
                    }
                }
            }

            if (previous2.has_value() && previous1.has_value()) {
                const Problem1AdaptiveScanPoint& left = *previous2;
                const Problem1AdaptiveScanPoint& mid = *previous1;
                const Problem1AdaptiveScanPoint& right = current;
                if (left.valid && mid.valid && right.valid &&
                    mid.residual <= left.residual &&
                    mid.residual <= right.residual) {
                    process_interval(
                        Problem1SuspiciousInterval{
                            left.phi,
                            right.phi,
                            Problem1SuspiciousIntervalReason::LocalMinimum,
                            i - 2,
                            i,
                            left.residual,
                            right.residual,
                            mid.residual,
                        },
                        &left,
                        &mid,
                        &right);
                }
                if (left.valid && mid.valid && right.valid &&
                    mid.residual >= left.residual &&
                    mid.residual >= right.residual) {
                    process_interval(
                        Problem1SuspiciousInterval{
                            left.phi,
                            right.phi,
                            Problem1SuspiciousIntervalReason::LocalMaximum,
                            i - 2,
                            i,
                            left.residual,
                            right.residual,
                            mid.residual,
                        },
                        &left,
                        &mid,
                        &right);
                }
            }

            previous2 = previous1;
            previous1 = current;
        }

    if (first_point.has_value() && previous1.has_value()) {
        const Problem1AdaptiveScanPoint& first = *first_point;
        const Problem1AdaptiveScanPoint& last = *previous1;
        if (last.valid != first.valid) {
            process_interval(
                Problem1SuspiciousInterval{
                    last.phi,
                    first.phi + kTwoPi,
                    Problem1SuspiciousIntervalReason::ValidBoundary,
                    coarse_scan_count - 1,
                    0,
                    last.residual,
                    first.residual,
                    std::numeric_limits<double>::quiet_NaN(),
                },
                &last,
                nullptr,
                &first);
        } else if (last.valid && first.valid && residual_sign_changed(last.residual, first.residual)) {
            process_interval(
                Problem1SuspiciousInterval{
                    last.phi,
                    first.phi + kTwoPi,
                    Problem1SuspiciousIntervalReason::WrapSignChange,
                    coarse_scan_count - 1,
                    0,
                    last.residual,
                    first.residual,
                    std::numeric_limits<double>::quiet_NaN(),
                },
                &last,
                nullptr,
                &first);
        }
    }

    if (summary != nullptr) {
        const double streaming_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - streaming_start).count();
        const double interval_refine_seconds_delta =
            summary->interval_refine_seconds - interval_refine_seconds_before;
        const double scan_detection_seconds = std::max(0.0, streaming_seconds - interval_refine_seconds_delta);
        summary->coarse_scan_seconds += scan_detection_seconds;
        summary->interval_collection_seconds += scan_detection_seconds;
    }
}

Problem1SolveWithDiagnosticResult solve_problem1_from_departure_anomalies_full_scan_internal(
    const Problem1DepartureAnomalySolveInput& input
) {
    using clock = std::chrono::steady_clock;
    validate_departure_anomaly_solve_input(input);

    Problem1SolveWithDiagnosticResult result{};
    result.diagnostic.valid = true;
    result.diagnostic.enumerated_branch_pairs =
        (input.max_transfer_revolution + 1) * (input.max_target_revolution + 1);

    std::vector<Problem1SolutionBranchCandidate> branch_candidates;
    const auto root_scanning_start = clock::now();
    for (int transfer_revolution = 0; transfer_revolution <= input.max_transfer_revolution; ++transfer_revolution) {
        for (int target_revolution = 0; target_revolution <= input.max_target_revolution; ++target_revolution) {
            bool has_previous_valid = false;
            bool has_first_valid = false;
            double previous_phi = 0.0;
            Problem1DepartureAnomalyResidualResult previous_result =
                make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);
            double first_phi = 0.0;
            Problem1DepartureAnomalyResidualResult first_result =
                make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);
            double last_valid_phi = 0.0;
            Problem1DepartureAnomalyResidualResult last_valid_result =
                make_departure_anomaly_result(Problem1ResidualStatus::InvalidBranch);

            for (int i = 0; i < input.phi_scan_count; ++i) {
                const double phi = static_cast<double>(i) * kTwoPi / static_cast<double>(input.phi_scan_count);
                result.diagnostic.alpha_scan_samples += 1;
                const Problem1DepartureAnomalyResidualResult residual_result =
                    evaluate_problem1_residual_from_departure_anomalies_profiled(
                        make_departure_anomaly_residual_input(input, phi, transfer_revolution, target_revolution),
                        &result.diagnostic);

                if (!is_success_with_finite_residual(residual_result)) {
                    has_previous_valid = false;
                    continue;
                }

                if (!has_first_valid) {
                    has_first_valid = true;
                    first_phi = phi;
                    first_result = residual_result;
                }

                if (std::abs(residual_result.residual) <= input.residual_tolerance) {
                    auto branch = make_branch_candidate_from_residual_result(
                        phi,
                        residual_result,
                        transfer_revolution,
                        target_revolution,
                        0.0,
                        input.max_candidate_relative_residual);
                    if (branch.has_value()) {
                        add_branch_if_not_duplicate(&branch_candidates, std::move(*branch));
                    }
                }

                if (has_previous_valid && residual_sign_changed(previous_result.residual, residual_result.residual)) {
                    result.diagnostic.bisection_refinements += 1;
                    auto refined = refine_problem1_root_by_bisection_from_departure_anomalies(
                        input,
                        previous_phi,
                        previous_result,
                        phi,
                        residual_result,
                        transfer_revolution,
                        target_revolution,
                        &result.diagnostic);
                    if (refined.has_value()) {
                        add_branch_if_not_duplicate(&branch_candidates, std::move(*refined));
                    }
                }

                previous_phi = phi;
                previous_result = residual_result;
                last_valid_phi = phi;
                last_valid_result = residual_result;
                has_previous_valid = true;
            }

            if (has_first_valid && has_previous_valid &&
                residual_sign_changed(last_valid_result.residual, first_result.residual)) {
                result.diagnostic.bisection_refinements += 1;
                auto refined = refine_problem1_root_by_bisection_from_departure_anomalies(
                    input,
                    last_valid_phi,
                    last_valid_result,
                    first_phi + kTwoPi,
                    first_result,
                    transfer_revolution,
                    target_revolution,
                    &result.diagnostic);
                if (refined.has_value()) {
                    add_branch_if_not_duplicate(&branch_candidates, std::move(*refined));
                }
            }
        }
    }
    const auto root_scanning_end = clock::now();
    result.diagnostic.root_scanning_seconds =
        std::chrono::duration<double>(root_scanning_end - root_scanning_start).count();

    finalize_problem1_branch_candidates(branch_candidates, &result.branches, &result.diagnostic);
    return result;
}

Problem1SolveWithDiagnosticResult solve_problem1_from_departure_anomalies_adaptive_internal(
    const Problem1DepartureAnomalySolveInput& input,
    const Problem1AdaptiveScanOptions& options,
    Problem1AdaptiveScanSummary* summary = nullptr
) {
    validate_departure_anomaly_solve_input(input);
    const int coarse_scan_count = adaptive_coarse_scan_count_from_environment(input.phi_scan_count);
    const double near_zero_tol = std::max(input.residual_tolerance, 1e-12);
    const double rapid_change_tol = std::numeric_limits<double>::infinity();
    Problem1SolveWithDiagnosticResult result{};
    result.diagnostic.valid = true;
    result.diagnostic.enumerated_branch_pairs =
        (input.max_transfer_revolution + 1) * (input.max_target_revolution + 1);
    std::vector<Problem1SolutionBranchCandidate> branch_candidates;
    for (int transfer_revolution = 0; transfer_revolution <= input.max_transfer_revolution; ++transfer_revolution) {
        for (int target_revolution = 0; target_revolution <= input.max_target_revolution; ++target_revolution) {
            solve_problem1_adaptive_streaming_for_revolutions(
                input,
                transfer_revolution,
                target_revolution,
                coarse_scan_count,
                near_zero_tol,
                rapid_change_tol,
                options.fine_subdivision_count,
                &branch_candidates,
                &result.diagnostic,
                summary);
        }
    }
    // TODO:
    // 1. coarse scan
    // 2. collect suspicious intervals
    // 3. local fine scan
    // 4. ternary search on residual extrema
    // 5. bisection refine sign-change intervals
    // 6. deduplicate candidates
    // 7. fallback to full scan if debug compare detects missing branches
    {
        ScopedSecondsAccumulator dedup_timer(
            summary != nullptr ? &summary->candidate_dedup_seconds : nullptr);
        ScopedSecondsAccumulator sorting_timer(
            summary != nullptr ? &summary->sorting_seconds : nullptr);
        finalize_problem1_branch_candidates(branch_candidates, &result.branches, &result.diagnostic);
    }
    if (summary != nullptr) {
        summary->candidate_count_after_dedup = static_cast<long long>(result.branches.size());
    }
    return result;
}

void compare_problem1_solve_branches_in_order(
    const std::vector<Problem1SolutionBranch>& adaptive_branches,
    const std::vector<Problem1SolutionBranch>& baseline_branches,
    Problem1SolveModeProfile* profile
) {
    const int common_count = static_cast<int>(std::min(adaptive_branches.size(), baseline_branches.size()));
    profile->debug_compare_matched_count = common_count;
    profile->debug_compare_missing_from_adaptive =
        std::max(0, static_cast<int>(baseline_branches.size()) - common_count);
    profile->debug_compare_new_in_adaptive =
        std::max(0, static_cast<int>(adaptive_branches.size()) - common_count);
    for (int i = 0; i < common_count; ++i) {
        const auto& adaptive = adaptive_branches[static_cast<std::size_t>(i)];
        const auto& baseline = baseline_branches[static_cast<std::size_t>(i)];
        if (adaptive.transfer_revolution != baseline.transfer_revolution ||
            adaptive.target_revolution != baseline.target_revolution) {
            profile->debug_compare_missing_from_adaptive += 1;
            profile->debug_compare_new_in_adaptive += 1;
            continue;
        }
        profile->debug_compare_max_alpha_diff = std::max(
            profile->debug_compare_max_alpha_diff,
            std::abs(normalize_angle_minus_pi_pi(
                adaptive.encounter_global_angle - baseline.encounter_global_angle)));
        profile->debug_compare_max_theta_prime_diff = std::max(
            profile->debug_compare_max_theta_prime_diff,
            std::abs(normalize_angle_minus_pi_pi(
                adaptive.target_arrival_true_anomaly - baseline.target_arrival_true_anomaly)));
        profile->debug_compare_max_target_time_diff = std::max(
            profile->debug_compare_max_target_time_diff,
            std::abs(adaptive.target_time_seconds - baseline.target_time_seconds));
        profile->debug_compare_max_residual_diff = std::max(
            profile->debug_compare_max_residual_diff,
            std::abs(adaptive.residual_seconds - baseline.residual_seconds));
    }
}

void add_query_branch_if_not_duplicate(
    std::vector<Problem1SolutionBranch>* branches,
    Problem1SolutionBranch branch
) {
    for (Problem1SolutionBranch& existing : *branches) {
        if (existing.transfer_revolution != branch.transfer_revolution ||
            existing.target_revolution != branch.target_revolution) {
            continue;
        }
        const double angle_diff = wrapped_angle_distance(
            existing.encounter_global_angle,
            branch.encounter_global_angle);
        const double time_diff = std::abs(existing.time_of_flight_seconds - branch.time_of_flight_seconds);
        if (angle_diff <= 1e-8 && time_diff <= 1e-8) {
            if (std::abs(branch.residual_seconds) < std::abs(existing.residual_seconds)) {
                existing = std::move(branch);
            }
            return;
        }
    }
    branches->push_back(std::move(branch));
}

bool same_route_b_approximation_identity(
    const Problem1RootApproximationResult& lhs,
    const Problem1RootApproximationResult& rhs
) {
    if (lhs.transfer_revolution != rhs.transfer_revolution) {
        return false;
    }
    const double time_tol = std::max(
        1e-4,
        1e-10 * std::max(std::abs(lhs.transfer_time_seconds), std::abs(rhs.transfer_time_seconds)));
    const double alpha_tol = 1e-4;
    const double time_diff = std::abs(lhs.transfer_time_seconds - rhs.transfer_time_seconds);
    const double angle_diff = wrapped_angle_distance(
        lhs.predicted_encounter_global_angle,
        rhs.predicted_encounter_global_angle);
    return time_diff <= time_tol && angle_diff <= alpha_tol;
}

void add_route_b_approximation_if_not_duplicate(
    std::vector<Problem1RootApproximationResult>* results,
    Problem1RootApproximationResult candidate
) {
    for (Problem1RootApproximationResult& existing : *results) {
        if (!same_route_b_approximation_identity(existing, candidate)) {
            continue;
        }
        if (std::abs(candidate.residual_seconds) < std::abs(existing.residual_seconds)) {
            existing = std::move(candidate);
        }
        return;
    }
    results->push_back(std::move(candidate));
}

int lower_corner_index_from_nearest(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
}

void populate_route_b_branch_count_diagnostics(
    const Problem1RootCellAdmissibilityResult& cell_admissibility,
    const std::vector<Problem1RootApproximationResult>& approximations,
    Problem1RouteBSafeQueryResult* result
) {
    result->expected_count_by_k = cell_admissibility.reference_root_count_by_k;
    result->candidate_count_by_k.clear();
    result->missing_count_by_k.clear();
    result->extra_count_by_k.clear();
    result->incomplete_by_k.clear();

    for (const Problem1RootApproximationResult& approximation : approximations) {
        result->candidate_count_by_k[approximation.transfer_revolution] += 1;
    }
    result->branch_count_complete = true;
    for (const auto& [k, expected] : result->expected_count_by_k) {
        const int actual = result->candidate_count_by_k.count(k) > 0 ? result->candidate_count_by_k.at(k) : 0;
        result->missing_count_by_k[k] = std::max(0, expected - actual);
        result->extra_count_by_k[k] = std::max(0, actual - expected);
        result->incomplete_by_k[k] = actual != expected;
        if (actual != expected) {
            result->branch_count_complete = false;
        }
    }
}

}  // namespace

Problem1RootTable::Problem1RootTable(Problem1RootTableConfig config)
    : config_(std::move(config)),
      cells_(static_cast<std::size_t>(config_.nu_A_count) *
             static_cast<std::size_t>(config_.nu_B_depart_count) *
             static_cast<std::size_t>(config_.theta_A_count)) {}

const Problem1RootTableConfig& Problem1RootTable::config() const {
    return config_;
}

std::size_t Problem1RootTable::flat_index(int nu_A_index, int nu_B_depart_index, int theta_A_index) const {
    if (nu_A_index < 0 || nu_A_index >= config_.nu_A_count ||
        nu_B_depart_index < 0 || nu_B_depart_index >= config_.nu_B_depart_count ||
        theta_A_index < 0 || theta_A_index >= config_.theta_A_count) {
        throw std::out_of_range("Problem1RootTable index out of range");
    }

    return (static_cast<std::size_t>(nu_A_index) *
            static_cast<std::size_t>(config_.nu_B_depart_count) +
            static_cast<std::size_t>(nu_B_depart_index)) *
               static_cast<std::size_t>(config_.theta_A_count) +
        static_cast<std::size_t>(theta_A_index);
}

const Problem1RootTableCell& Problem1RootTable::at(int nu_A_index, int nu_B_depart_index, int theta_A_index) const {
    return cells_.at(flat_index(nu_A_index, nu_B_depart_index, theta_A_index));
}

const std::vector<Problem1RootTableCell>& Problem1RootTable::cells() const {
    return cells_;
}

std::vector<Problem1RootTableCell>& Problem1RootTable::mutable_cells() {
    return cells_;
}

Problem1SolutionBranch problem1_solution_branch_from_candidate(
    planet_params::PlanetId target_planet,
    const Problem1Candidate& candidate
) {
    Problem1SolutionBranch branch{};
    branch.valid = candidate.residual_result.status == Problem1ResidualStatus::Success &&
        is_finite(candidate.encounter_global_angle) &&
        is_finite(candidate.time_of_flight_seconds);
    branch.encounter_global_angle = normalize_angle_0_2pi(candidate.encounter_global_angle);
    branch.target_arrival_true_anomaly = normalize_angle_0_2pi(
        candidate.encounter_global_angle -
        planet_params::get_planet_params(target_planet).orbit.theta_0);
    branch.transfer_revolution = candidate.transfer_revolution;
    branch.target_revolution = candidate.target_revolution;
    branch.time_of_flight_seconds = candidate.time_of_flight_seconds;
    branch.target_time_seconds = scale_free_time_to_seconds(candidate.residual_result.target_time_scale_free);
    branch.residual_seconds = scale_free_time_to_seconds(candidate.residual_result.residual);
    branch.arrival_time_seconds_since_j2000 = candidate.arrival_time_seconds_since_j2000;
    branch.transfer_e = candidate.residual_result.transfer_e;
    branch.transfer_p = candidate.residual_result.transfer_p;
    branch.transfer_a = compute_transfer_a(candidate.residual_result.transfer_p, candidate.residual_result.transfer_e);
    branch.theta_B = candidate.residual_result.xi2;
    if (!branch.valid) {
        branch.invalid_reason = "candidate is not a finite success root";
    }
    return branch;
}

Problem1BranchCompareReport compare_problem1_solution_branches_by_identity(
    const std::vector<Problem1SolutionBranch>& candidate_branches,
    const std::vector<Problem1SolutionBranch>& baseline_branches,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1BranchCompareOptions& options
) {
    Problem1BranchCompareReport report{};
    const int max_detail_count = std::max(0, options.max_detail_count);
    std::vector<bool> candidate_matched(candidate_branches.size(), false);

    for (const Problem1SolutionBranch& baseline : baseline_branches) {
        int best_index = -1;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < candidate_branches.size(); ++i) {
            if (candidate_matched[i]) {
                continue;
            }
            const Problem1SolutionBranch& candidate = candidate_branches[i];
            if (!branch_matches_by_identity(candidate, baseline, options)) {
                continue;
            }
            const double score = branch_identity_match_score(candidate, baseline, options);
            if (score < best_score) {
                best_score = score;
                best_index = static_cast<int>(i);
            }
        }

        if (best_index < 0) {
            report.missing_count += 1;
            if (static_cast<int>(report.missing_branches.size()) < max_detail_count) {
                report.missing_branches.push_back(
                    make_branch_compare_detail(baseline, nu_A_depart, nu_B_depart, theta_A));
            }
            continue;
        }

        candidate_matched[static_cast<std::size_t>(best_index)] = true;
        const Problem1SolutionBranch& candidate = candidate_branches[static_cast<std::size_t>(best_index)];
        report.matched_count += 1;
        report.max_angle_error = std::max(
            report.max_angle_error,
            wrapped_angle_distance(candidate.encounter_global_angle, baseline.encounter_global_angle));
        report.max_target_anomaly_error = std::max(
            report.max_target_anomaly_error,
            wrapped_angle_distance(candidate.target_arrival_true_anomaly, baseline.target_arrival_true_anomaly));
        report.max_time_error = std::max(
            report.max_time_error,
            std::min(
                std::abs(candidate.time_of_flight_seconds - baseline.time_of_flight_seconds),
                std::abs(candidate.target_time_seconds - baseline.target_time_seconds)));
        report.max_residual_error = std::max(
            report.max_residual_error,
            std::abs(candidate.residual_seconds - baseline.residual_seconds));
    }

    for (std::size_t i = 0; i < candidate_branches.size(); ++i) {
        if (candidate_matched[i]) {
            continue;
        }
        report.extra_count += 1;
        if (static_cast<int>(report.extra_branches.size()) < max_detail_count) {
            report.extra_branches.push_back(
                make_branch_compare_detail(candidate_branches[i], nu_A_depart, nu_B_depart, theta_A));
        }
    }

    return report;
}

Problem1RootResidualResult evaluate_problem1_root_residual(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    const Problem1DepartureAnomalyResidualResult internal =
        evaluate_problem1_residual_from_departure_anomalies({
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            encounter_global_angle,
            transfer_revolution,
            target_revolution,
        });
    return make_public_root_residual_result(internal);
}

double problem1_residual_seconds_to_scale_free(double residual_seconds) {
    return problem1_residual_seconds_to_scale_free_internal(residual_seconds);
}

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
) {
    Problem1RootQSheetSelectionResult result{};
    result.selected_q = source_branch.target_revolution;
    if (!is_finite(alpha_linear) || max_target_revolution < 0) {
        result.selection_failed = true;
        result.selected_continuity_error = std::numeric_limits<double>::infinity();
        result.source_q_continuity_error = std::numeric_limits<double>::infinity();
        return result;
    }

    const double source_time_reference =
        is_finite(source_branch.target_time_seconds)
            ? source_branch.target_time_seconds
            : source_branch.time_of_flight_seconds;

    double best_error = std::numeric_limits<double>::infinity();
    bool any_valid = false;
    for (int q = 0; q <= max_target_revolution; ++q) {
        const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_linear,
            transfer_revolution,
            q);
        if (!residual.valid || !is_finite(residual.target_time_seconds)) {
            continue;
        }
        const double continuity_error = std::abs(residual.target_time_seconds - source_time_reference);
        if (q == source_branch.target_revolution) {
            result.source_q_continuity_error = continuity_error;
        }
        if (!any_valid || continuity_error < best_error) {
            any_valid = true;
            best_error = continuity_error;
            result.selected_q = q;
            result.selected_continuity_error = continuity_error;
        }
    }

    if (!any_valid) {
        result.selection_failed = true;
        result.selected_q = source_branch.target_revolution;
        result.q_changed = false;
        result.selected_continuity_error = std::numeric_limits<double>::infinity();
        if (!is_finite(result.source_q_continuity_error)) {
            result.source_q_continuity_error = std::numeric_limits<double>::infinity();
        }
        return result;
    }

    result.q_changed = result.selected_q != source_branch.target_revolution;
    return result;
}

Problem1RootResidualDerivatives evaluate_problem1_root_residual_derivatives_analytic(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    Problem1RootResidualDerivatives derivatives{};

    const Problem1DepartureAnomalyResidualResult base =
        evaluate_problem1_residual_from_departure_anomalies({
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            encounter_global_angle,
            transfer_revolution,
            target_revolution,
        });
    if (!is_success_with_finite_residual(base)) {
        derivatives.invalid_reason = base.invalid_reason.empty() ? "invalid base residual" : base.invalid_reason;
        return derivatives;
    }
    if (std::abs(base.transfer_e - 1.0) <= 1e-10) {
        derivatives.invalid_reason = "analytic derivatives are unsupported near e = 1";
        return derivatives;
    }
    if (near_unwrap_boundary(normalize_angle_0_2pi(nu_B_depart), base.target_arrival_true_anomaly)) {
        derivatives.invalid_reason = "analytic derivatives are unsupported on target unwrap boundary";
        return derivatives;
    }

    const planet_params::PlanetParams& departure = planet_params::get_planet_params(departure_planet);
    const planet_params::PlanetParams& target = planet_params::get_planet_params(target_planet);

    const double nu_A = normalize_angle_0_2pi(nu_A_depart);
    const double nu_B = normalize_angle_0_2pi(nu_B_depart);
    const double theta_A_wrapped = normalize_angle_0_2pi(theta_A);
    const double target_arrival_true_anomaly = base.target_arrival_true_anomaly;
    const double p_A = departure.orbit.p;
    const double e_A = departure.orbit.e;
    const double p_B = target.orbit.p;
    const double e_B = target.orbit.e;
    const Problem1TransferOrbitGeometry geometry = compute_problem1_transfer_orbit_geometry(
        departure_planet,
        target_planet,
        nu_A_depart,
        theta_A,
        encounter_global_angle);
    if (!geometry.valid) {
        derivatives.invalid_reason = geometry.invalid_reason.empty()
            ? "analytic derivatives geometry reconstruction failed"
            : geometry.invalid_reason;
        return derivatives;
    }
    if (std::abs(geometry.e_raw) <= 1e-12) {
        derivatives.invalid_reason = "analytic derivatives are unsupported near |e_raw| = 0";
        return derivatives;
    }

    const double r_A = geometry.r_A;
    const double r_B = geometry.r_B;
    const double theta_B = base.theta_B;
    const double e_t = geometry.e;
    const double p_t = geometry.p;
    const double theta_A_eval =
        e_t < 1.0 ? geometry.theta_A_canonical : normalize_angle_minus_pi_pi(geometry.theta_A_canonical);
    if (std::abs(std::abs(theta_A_eval) - spaceship_cpp::common::kPi) <= 1e-10) {
        derivatives.invalid_reason = "analytic derivatives are unsupported on transfer unwrap boundary";
        return derivatives;
    }

    const double denominator_A = 1.0 + e_A * std::cos(nu_A);
    const double denominator_B = 1.0 + e_B * std::cos(target_arrival_true_anomaly);
    const double r_A_nu_A = p_A * e_A * std::sin(nu_A) / (denominator_A * denominator_A);
    const double r_B_alpha =
        p_B * e_B * std::sin(target_arrival_true_anomaly) / (denominator_B * denominator_B);

    const double theta_B_alpha = 1.0;
    const double theta_B_nu_A = -1.0;
    const double theta_B_nu_B = 0.0;
    const double theta_B_theta_A = 1.0;
    const double theta_A_alpha = 0.0;
    const double theta_A_nu_A = 0.0;
    const double theta_A_nu_B = 0.0;
    const double theta_A_theta_A = 1.0;

    const double N = geometry.N;
    const double D = geometry.D;
    if (!is_finite(D) || std::abs(D) <= 1e-12) {
        derivatives.invalid_reason = "analytic derivatives are unsupported for singular e denominator";
        return derivatives;
    }

    const auto compute_e_x = [&](double r_A_x, double r_B_x, double theta_B_x, double theta_A_x) {
        const double N_x = r_B_x - r_A_x;
        const double D_x =
            r_A_x * std::cos(theta_A_eval) -
            r_A * std::sin(theta_A_eval) * theta_A_x -
            r_B_x * std::cos(theta_B) +
            r_B * std::sin(theta_B) * theta_B_x;
        return (N_x * D - N * D_x) / (D * D);
    };

    const double e_alpha = compute_e_x(0.0, r_B_alpha, theta_B_alpha, theta_A_alpha);
    const double e_nu_A = compute_e_x(r_A_nu_A, 0.0, theta_B_nu_A, theta_A_nu_A);
    const double e_nu_B = 0.0;
    const double e_theta_A = compute_e_x(0.0, 0.0, theta_B_theta_A, theta_A_theta_A);

    const auto compute_p_x = [&](double r_A_x, double e_x, double theta_A_x) {
        return
            r_A_x * (1.0 + e_t * std::cos(theta_A_eval)) +
            r_A * (e_x * std::cos(theta_A_eval) - e_t * std::sin(theta_A_eval) * theta_A_x);
    };
    const double p_alpha = compute_p_x(0.0, e_alpha, theta_A_alpha);
    const double p_nu_A = compute_p_x(r_A_nu_A, e_nu_A, theta_A_nu_A);
    const double p_nu_B = 0.0;
    const double p_theta_A = compute_p_x(0.0, e_theta_A, theta_A_theta_A);

    const double F_e_theta_B = orbit_F_e_derivative(e_t, theta_B);
    const double F_e_theta_A = orbit_F_e_derivative(e_t, theta_A_eval);
    const double F_theta_theta_B = orbit_F_theta_derivative(e_t, theta_B);
    const double F_theta_theta_A = orbit_F_theta_derivative(e_t, theta_A_eval);
    const double deltaF_transfer = base.deltaF_transfer;

    const auto compute_transfer_time_x = [&](double p_x, double e_x, double theta_B_x, double theta_A_x) {
        const double deltaF_x =
            (F_e_theta_B - F_e_theta_A) * e_x +
            F_theta_theta_B * theta_B_x -
            F_theta_theta_A * theta_A_x;
        return
            1.5 * std::sqrt(p_t) * p_x * deltaF_transfer +
            std::pow(p_t, 1.5) * deltaF_x;
    };

    const double T_transfer_alpha =
        compute_transfer_time_x(p_alpha, e_alpha, theta_B_alpha, theta_A_alpha);
    const double T_transfer_nu_A =
        compute_transfer_time_x(p_nu_A, e_nu_A, theta_B_nu_A, theta_A_nu_A);
    const double T_transfer_nu_B =
        compute_transfer_time_x(p_nu_B, e_nu_B, theta_B_nu_B, theta_A_nu_B);
    const double T_transfer_theta_A =
        compute_transfer_time_x(p_theta_A, e_theta_A, theta_B_theta_A, theta_A_theta_A);

    const double target_scale = std::pow(target.orbit.p, 1.5);
    // 中文注释：为了避免重新推导 q 对应的 unwrap 偏移，这里直接从 residual evaluator 的
    // scale-free target 时间恢复分支末端 F 值，再配合 F_theta 的局部导数使用。
    const double target_end_unwrapped_F = orbit_F(e_B, nu_B) + base.target_time_scale_free / target_scale;
    double target_end_unwrapped = target_arrival_true_anomaly;
    while (orbit_F(e_B, target_end_unwrapped) < target_end_unwrapped_F - 1e-12) {
        target_end_unwrapped += kTwoPi;
    }
    const double T_target_alpha_checked =
        target_scale * orbit_F_theta_derivative(e_B, target_end_unwrapped);
    const double T_target_nu_B =
        -target_scale * orbit_F_theta_derivative(e_B, nu_B);

    derivatives.R_alpha = T_transfer_alpha - T_target_alpha_checked;
    derivatives.R_nu_A = T_transfer_nu_A;
    derivatives.R_nu_B = T_transfer_nu_B - T_target_nu_B;
    derivatives.R_theta_A = T_transfer_theta_A;
    derivatives.F_alpha = derivatives.R_alpha;
    if (!is_finite(derivatives.R_alpha) || std::abs(derivatives.R_alpha) <= 1e-12) {
        derivatives.invalid_reason = "analytic derivatives are unsupported for near-zero R_alpha";
        return derivatives;
    }

    derivatives.d_alpha_d_nu_A = -derivatives.R_nu_A / derivatives.R_alpha;
    derivatives.d_alpha_d_nu_B = -derivatives.R_nu_B / derivatives.R_alpha;
    derivatives.d_alpha_d_theta_A = -derivatives.R_theta_A / derivatives.R_alpha;
    if (!is_finite(derivatives.d_alpha_d_nu_A) ||
        !is_finite(derivatives.d_alpha_d_nu_B) ||
        !is_finite(derivatives.d_alpha_d_theta_A)) {
        derivatives.invalid_reason = "implicit derivatives produced non-finite values";
        return derivatives;
    }
    derivatives.valid = true;
    derivatives.invalid_reason.clear();
    return derivatives;
}

Problem1RootResidualDerivatives evaluate_problem1_root_residual_derivatives(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    double encounter_global_angle,
    int transfer_revolution,
    int target_revolution
) {
    return evaluate_problem1_root_residual_derivatives_analytic(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        encounter_global_angle,
        transfer_revolution,
        target_revolution);
}

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
) {
    Problem1RootResidualDerivatives derivatives{};
    if (!is_finite(step) || !(step > 0.0)) {
        derivatives.invalid_reason = "finite-difference derivative step must be positive";
        return derivatives;
    }

    const auto central_difference = [&](const char* axis_name,
                                        double minus_nu_A,
                                        double plus_nu_A,
                                        double minus_nu_B,
                                        double plus_nu_B,
                                        double minus_theta_A,
                                        double plus_theta_A,
                                        double minus_alpha,
                                        double plus_alpha,
                                        double* output) -> bool {
        const Problem1RootResidualResult minus_residual = evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            minus_nu_A,
            minus_nu_B,
            minus_theta_A,
            minus_alpha,
            transfer_revolution,
            target_revolution);
        if (!minus_residual.valid || !is_finite(minus_residual.residual_scale_free)) {
            derivatives.invalid_reason = std::string("finite-difference residual failed on ") + axis_name + " minus side: " +
                (minus_residual.invalid_reason.empty() ? "invalid residual" : minus_residual.invalid_reason);
            return false;
        }
        const Problem1RootResidualResult plus_residual = evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            plus_nu_A,
            plus_nu_B,
            plus_theta_A,
            plus_alpha,
            transfer_revolution,
            target_revolution);
        if (!plus_residual.valid || !is_finite(plus_residual.residual_scale_free)) {
            derivatives.invalid_reason = std::string("finite-difference residual failed on ") + axis_name + " plus side: " +
                (plus_residual.invalid_reason.empty() ? "invalid residual" : plus_residual.invalid_reason);
            return false;
        }
        *output = (plus_residual.residual_scale_free - minus_residual.residual_scale_free) / (2.0 * step);
        return is_finite(*output);
    };

    if (!central_difference(
            "alpha",
            normalize_angle_0_2pi(nu_A_depart),
            normalize_angle_0_2pi(nu_A_depart),
            normalize_angle_0_2pi(nu_B_depart),
            normalize_angle_0_2pi(nu_B_depart),
            normalize_angle_0_2pi(theta_A),
            normalize_angle_0_2pi(theta_A),
            normalize_angle_0_2pi(encounter_global_angle - step),
            normalize_angle_0_2pi(encounter_global_angle + step),
            &derivatives.R_alpha)) {
        if (derivatives.invalid_reason.empty()) {
            derivatives.invalid_reason = "finite-difference R_alpha is non-finite";
        }
        return derivatives;
    }
    if (!central_difference(
            "nu_A",
            normalize_angle_0_2pi(nu_A_depart - step),
            normalize_angle_0_2pi(nu_A_depart + step),
            normalize_angle_0_2pi(nu_B_depart),
            normalize_angle_0_2pi(nu_B_depart),
            normalize_angle_0_2pi(theta_A),
            normalize_angle_0_2pi(theta_A),
            normalize_angle_0_2pi(encounter_global_angle),
            normalize_angle_0_2pi(encounter_global_angle),
            &derivatives.R_nu_A)) {
        if (derivatives.invalid_reason.empty()) {
            derivatives.invalid_reason = "finite-difference R_nu_A is non-finite";
        }
        return derivatives;
    }
    if (!central_difference(
            "nu_B",
            normalize_angle_0_2pi(nu_A_depart),
            normalize_angle_0_2pi(nu_A_depart),
            normalize_angle_0_2pi(nu_B_depart - step),
            normalize_angle_0_2pi(nu_B_depart + step),
            normalize_angle_0_2pi(theta_A),
            normalize_angle_0_2pi(theta_A),
            normalize_angle_0_2pi(encounter_global_angle),
            normalize_angle_0_2pi(encounter_global_angle),
            &derivatives.R_nu_B)) {
        if (derivatives.invalid_reason.empty()) {
            derivatives.invalid_reason = "finite-difference R_nu_B is non-finite";
        }
        return derivatives;
    }
    if (!central_difference(
            "theta_A",
            normalize_angle_0_2pi(nu_A_depart),
            normalize_angle_0_2pi(nu_A_depart),
            normalize_angle_0_2pi(nu_B_depart),
            normalize_angle_0_2pi(nu_B_depart),
            normalize_angle_0_2pi(theta_A - step),
            normalize_angle_0_2pi(theta_A + step),
            normalize_angle_0_2pi(encounter_global_angle),
            normalize_angle_0_2pi(encounter_global_angle),
            &derivatives.R_theta_A)) {
        if (derivatives.invalid_reason.empty()) {
            derivatives.invalid_reason = "finite-difference R_theta_A is non-finite";
        }
        return derivatives;
    }

    derivatives.F_alpha = derivatives.R_alpha;
    if (!is_finite(derivatives.R_alpha) || std::abs(derivatives.R_alpha) <= 1e-12) {
        derivatives.invalid_reason = "finite-difference derivatives produced near-zero R_alpha";
        return derivatives;
    }
    derivatives.d_alpha_d_nu_A = -derivatives.R_nu_A / derivatives.R_alpha;
    derivatives.d_alpha_d_nu_B = -derivatives.R_nu_B / derivatives.R_alpha;
    derivatives.d_alpha_d_theta_A = -derivatives.R_theta_A / derivatives.R_alpha;
    if (!is_finite(derivatives.d_alpha_d_nu_A) ||
        !is_finite(derivatives.d_alpha_d_nu_B) ||
        !is_finite(derivatives.d_alpha_d_theta_A)) {
        derivatives.invalid_reason = "finite-difference implicit derivatives produced non-finite values";
        return derivatives;
    }
    derivatives.valid = true;
    derivatives.invalid_reason.clear();
    return derivatives;
}

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
) {
    switch (mode) {
    case Problem1RootDerivativeMode::AnalyticOnly:
        return evaluate_problem1_root_residual_derivatives_analytic(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            encounter_global_angle,
            transfer_revolution,
            target_revolution);
    case Problem1RootDerivativeMode::FiniteDifferenceOnly:
        return estimate_problem1_root_residual_derivatives_finite_difference(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            encounter_global_angle,
            transfer_revolution,
            target_revolution,
            finite_difference_step);
    case Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback: {
        Problem1RootResidualDerivatives analytic = evaluate_problem1_root_residual_derivatives_analytic(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            encounter_global_angle,
            transfer_revolution,
            target_revolution);
        if (analytic.valid) {
            return analytic;
        }
        Problem1RootResidualDerivatives fd = estimate_problem1_root_residual_derivatives_finite_difference(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            encounter_global_angle,
            transfer_revolution,
            target_revolution,
            finite_difference_step);
        if (!fd.valid && !analytic.invalid_reason.empty() && !fd.invalid_reason.empty()) {
            fd.invalid_reason = "analytic failed: " + analytic.invalid_reason + "; finite-difference failed: " +
                fd.invalid_reason;
        }
        return fd;
    }
    }
    Problem1RootResidualDerivatives invalid{};
    invalid.invalid_reason = "unknown derivative mode";
    return invalid;
}

Problem1SolutionBranch attach_problem1_root_derivatives(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1SolutionBranch& branch
) {
    return attach_problem1_root_derivatives_with_mode(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        branch,
        Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
}

Problem1SolutionBranch attach_problem1_root_derivatives_with_mode(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1SolutionBranch& branch,
    Problem1RootDerivativeMode mode,
    double finite_difference_step
) {
    Problem1SolutionBranch enriched = branch;
    const Problem1RootResidualDerivatives derivatives = evaluate_problem1_root_residual_derivatives_with_mode(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        branch.encounter_global_angle,
        branch.transfer_revolution,
        branch.target_revolution,
        mode,
        finite_difference_step);
    if (!derivatives.valid) {
        enriched.derivatives_available = false;
        enriched.invalid_reason = derivatives.invalid_reason;
        return enriched;
    }
    enriched.derivatives_available = true;
    enriched.d_encounter_global_angle_d_nu_A = derivatives.d_alpha_d_nu_A;
    enriched.d_encounter_global_angle_d_nu_B = derivatives.d_alpha_d_nu_B;
    enriched.d_encounter_global_angle_d_theta_A = derivatives.d_alpha_d_theta_A;
    return enriched;
}

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
    double residual_tolerance,
    double alpha_step_tolerance
) {
    return refine_problem1_root_branch_newton_diagnostic(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        transfer_revolution,
        target_revolution,
        initial_encounter_global_angle,
        max_iterations,
        residual_tolerance,
        alpha_step_tolerance).branch;
}

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
) {
    return refine_problem1_root_branch_newton(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        transfer_revolution,
        target_revolution,
        initial_encounter_global_angle,
        max_iterations,
        problem1_residual_seconds_to_scale_free(residual_tolerance_seconds),
        alpha_step_tolerance);
}

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
    double residual_tolerance,
    double alpha_step_tolerance,
    double max_alpha_step,
    bool enable_backtracking,
    int max_backtracking_steps,
    Problem1RootDerivativeMode derivative_mode,
    double finite_difference_step
) {
    if (max_iterations <= 0 || !is_finite(residual_tolerance) || residual_tolerance < 0.0 ||
        !is_finite(alpha_step_tolerance) || !(alpha_step_tolerance > 0.0) ||
        !is_finite(initial_encounter_global_angle) ||
        (!std::isinf(max_alpha_step) && (!is_finite(max_alpha_step) || !(max_alpha_step > 0.0))) ||
        max_backtracking_steps < 0 ||
        ((derivative_mode == Problem1RootDerivativeMode::FiniteDifferenceOnly ||
          derivative_mode == Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback) &&
         (!is_finite(finite_difference_step) || !(finite_difference_step > 0.0)))) {
        return make_invalid_refinement_result(
            transfer_revolution,
            target_revolution,
            initial_encounter_global_angle,
            "invalid Newton refinement inputs");
    }

    Problem1RootRefinementResult result{};
    result.diagnostic.valid = true;
    result.diagnostic.transfer_revolution = transfer_revolution;
    result.diagnostic.target_revolution = target_revolution;
    result.diagnostic.initial_alpha = normalize_angle_0_2pi(initial_encounter_global_angle);

    double alpha_unwrapped = initial_encounter_global_angle;
    double previous_abs_residual = std::numeric_limits<double>::quiet_NaN();

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Problem1RootNewtonTraceStep trace{};
        trace.iteration = iteration;
        trace.alpha = normalize_angle_0_2pi(alpha_unwrapped);
        trace.alpha_normalized = std::abs(alpha_unwrapped - trace.alpha) > 1e-12;

        const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trace.alpha,
            transfer_revolution,
            target_revolution);
        trace.residual_scale_free = residual.residual_scale_free;
        trace.residual_seconds = residual.residual_seconds;
        if (iteration == 0) {
            result.diagnostic.initial_residual_seconds = residual.residual_seconds;
        }
        if (!residual.valid || !is_finite(residual.residual_scale_free)) {
            trace.reason = residual.invalid_reason.empty()
                ? "Newton residual evaluation failed"
                : residual.invalid_reason;
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }

        const double abs_residual = std::abs(residual.residual_scale_free);
        if (is_finite(previous_abs_residual) && abs_residual > previous_abs_residual) {
            trace.residual_increased = true;
            result.diagnostic.residual_increased = true;
        }
        previous_abs_residual = abs_residual;

        if (abs_residual <= residual_tolerance) {
            trace.reason = "residual tolerance reached";
            result.diagnostic.trace.push_back(trace);
            result.branch = attach_problem1_root_derivatives(
                departure_planet,
                target_planet,
                nu_A_depart,
                nu_B_depart,
                theta_A,
                make_solution_branch_from_public_residual(residual, transfer_revolution, target_revolution));
            result.valid = result.branch.valid;
            result.diagnostic.converged = result.branch.valid;
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = result.branch.encounter_global_angle;
            result.diagnostic.final_residual_seconds = result.branch.residual_seconds;
            result.diagnostic.likely_wrong_root =
                std::abs(normalize_angle_minus_pi_pi(
                    result.branch.encounter_global_angle - result.diagnostic.initial_alpha)) > 0.2;
            return result;
        }

        const Problem1RootResidualDerivatives analytic_probe = evaluate_problem1_root_residual_derivatives_analytic(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trace.alpha,
            transfer_revolution,
            target_revolution);
        const Problem1RootResidualDerivatives derivatives = evaluate_problem1_root_residual_derivatives_with_mode(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trace.alpha,
            transfer_revolution,
            target_revolution,
            derivative_mode,
            finite_difference_step);
        trace.derivative_valid =
            derivatives.valid && is_finite(derivatives.R_alpha) && std::abs(derivatives.R_alpha) > 1e-12;
        if (derivative_mode == Problem1RootDerivativeMode::FiniteDifferenceOnly) {
            trace.derivative_source = derivatives.valid ? "finite_difference" : "none";
        } else if (analytic_probe.valid) {
            trace.derivative_source = "analytic";
        } else if (derivatives.valid) {
            trace.derivative_source = "finite_difference_fallback";
        } else {
            trace.derivative_source = "none";
        }
        if (derivative_mode == Problem1RootDerivativeMode::FiniteDifferenceOnly && derivatives.valid) {
            result.diagnostic.finite_difference_success_count += 1;
        } else if (!analytic_probe.valid &&
                   derivative_mode == Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback &&
                   derivatives.valid) {
            result.diagnostic.fallback_used_count += 1;
            result.diagnostic.finite_difference_success_count += 1;
        }
        trace.R_alpha = derivatives.R_alpha;
        if (!trace.derivative_valid) {
            trace.reason = derivatives.invalid_reason.empty()
                ? "Newton derivative evaluation failed"
                : derivatives.invalid_reason;
            result.diagnostic.derivative_failed = true;
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }

        double delta_alpha = -residual.residual_scale_free / derivatives.R_alpha;
        if (!is_finite(delta_alpha)) {
            trace.reason = "Newton step is non-finite";
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }

        if (is_finite(max_alpha_step) && std::abs(delta_alpha) > max_alpha_step) {
            delta_alpha = std::copysign(max_alpha_step, delta_alpha);
            trace.step_clamped = true;
            result.diagnostic.step_clamped = true;
        }

        double accepted_delta = delta_alpha;
        Problem1RootResidualResult accepted_trial_residual{};
        bool accepted_trial_valid = false;
        bool residual_increased = false;
        if (enable_backtracking) {
            const double current_abs_residual = std::abs(residual.residual_scale_free);
            bool found_reducing_step = false;
            for (int backtracking_iteration = 0;
                 backtracking_iteration <= std::max(0, max_backtracking_steps);
                 ++backtracking_iteration) {
                const double scale = std::ldexp(1.0, -backtracking_iteration);
                const double trial_delta = delta_alpha * scale;
                const double trial_alpha = normalize_angle_0_2pi(alpha_unwrapped + trial_delta);
                const Problem1RootResidualResult trial_residual = evaluate_problem1_root_residual(
                    departure_planet,
                    target_planet,
                    nu_A_depart,
                    nu_B_depart,
                    theta_A,
                    trial_alpha,
                    transfer_revolution,
                    target_revolution);
                if (!trial_residual.valid || !is_finite(trial_residual.residual_scale_free)) {
                    continue;
                }
                accepted_delta = trial_delta;
                accepted_trial_residual = trial_residual;
                accepted_trial_valid = true;
                if (std::abs(trial_residual.residual_scale_free) <= current_abs_residual) {
                    found_reducing_step = true;
                    break;
                }
            }
            if (accepted_trial_valid) {
                residual_increased = std::abs(accepted_trial_residual.residual_scale_free) >
                    std::abs(residual.residual_scale_free);
            }
            if (!found_reducing_step && accepted_trial_valid) {
                result.diagnostic.residual_increased = true;
            }
        }

        trace.delta_alpha = accepted_delta;
        trace.residual_increased = residual_increased;
        if (trace.residual_increased) {
            result.diagnostic.residual_increased = true;
        }

        alpha_unwrapped += accepted_delta;
        if (std::abs(accepted_delta) <= alpha_step_tolerance) {
            const Problem1RootResidualResult final_residual =
                accepted_trial_valid
                ? accepted_trial_residual
                : evaluate_problem1_root_residual(
                      departure_planet,
                      target_planet,
                      nu_A_depart,
                      nu_B_depart,
                      theta_A,
                      normalize_angle_0_2pi(alpha_unwrapped),
                      transfer_revolution,
                      target_revolution);
            if (!final_residual.valid || std::abs(final_residual.residual_scale_free) > residual_tolerance) {
                trace.reason = "Newton step tolerance reached without acceptable residual";
                result.diagnostic.trace.push_back(trace);
                result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
                result.diagnostic.final_alpha = normalize_angle_0_2pi(alpha_unwrapped);
                result.diagnostic.final_residual_seconds = final_residual.residual_seconds;
                result.diagnostic.invalid_reason = trace.reason;
                Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                    transfer_revolution,
                    target_revolution,
                    result.diagnostic.initial_alpha,
                    trace.reason);
                invalid.diagnostic = result.diagnostic;
                return invalid;
            }
            trace.reason = "step tolerance reached with acceptable residual";
            result.diagnostic.trace.push_back(trace);
            result.branch = attach_problem1_root_derivatives(
                departure_planet,
                target_planet,
                nu_A_depart,
                nu_B_depart,
                theta_A,
                make_solution_branch_from_public_residual(final_residual, transfer_revolution, target_revolution));
            result.valid = result.branch.valid;
            result.diagnostic.converged = result.branch.valid;
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = result.branch.encounter_global_angle;
            result.diagnostic.final_residual_seconds = result.branch.residual_seconds;
            result.diagnostic.likely_wrong_root =
                std::abs(normalize_angle_minus_pi_pi(
                    result.branch.encounter_global_angle - result.diagnostic.initial_alpha)) > 0.2;
            return result;
        }

        trace.reason = "continue";
        result.diagnostic.trace.push_back(trace);
    }

    result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
    result.diagnostic.final_alpha = normalize_angle_0_2pi(alpha_unwrapped);
    result.diagnostic.final_residual_seconds = result.diagnostic.trace.empty()
        ? std::numeric_limits<double>::quiet_NaN()
        : result.diagnostic.trace.back().residual_seconds;
    result.diagnostic.invalid_reason = "Newton refinement exceeded maximum iterations";
    Problem1RootRefinementResult invalid = make_invalid_refinement_result(
        transfer_revolution,
        target_revolution,
        result.diagnostic.initial_alpha,
        result.diagnostic.invalid_reason);
    invalid.diagnostic = result.diagnostic;
    return invalid;
}

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
    double max_alpha_step,
    bool enable_backtracking,
    int max_backtracking_steps,
    Problem1RootDerivativeMode derivative_mode,
    double finite_difference_step
) {
    return refine_problem1_root_branch_newton_diagnostic(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        transfer_revolution,
        target_revolution,
        initial_encounter_global_angle,
        max_iterations,
        problem1_residual_seconds_to_scale_free(residual_tolerance_seconds),
        alpha_step_tolerance,
        max_alpha_step,
        enable_backtracking,
        max_backtracking_steps,
        derivative_mode,
        finite_difference_step);
}

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
    Problem1RootDerivativeMode derivative_mode,
    double finite_difference_step
) {
    if (max_iterations <= 0 ||
        !is_finite(residual_tolerance_seconds) || residual_tolerance_seconds < 0.0 ||
        !is_finite(alpha_step_tolerance) || !(alpha_step_tolerance > 0.0) ||
        !is_finite(initial_encounter_global_angle) ||
        ((derivative_mode == Problem1RootDerivativeMode::FiniteDifferenceOnly ||
          derivative_mode == Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback) &&
         (!is_finite(finite_difference_step) || !(finite_difference_step > 0.0)))) {
        return make_invalid_refinement_result(
            transfer_revolution,
            target_revolution,
            initial_encounter_global_angle,
            "invalid residual-first Newton inputs");
    }

    Problem1RootRefinementResult result{};
    result.diagnostic.valid = true;
    result.diagnostic.transfer_revolution = transfer_revolution;
    result.diagnostic.target_revolution = target_revolution;
    result.diagnostic.initial_alpha = normalize_angle_0_2pi(initial_encounter_global_angle);

    double alpha_unwrapped = initial_encounter_global_angle;
    double previous_abs_residual_scale_free = std::numeric_limits<double>::infinity();

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Problem1RootNewtonTraceStep trace{};
        trace.iteration = iteration;
        trace.alpha = normalize_angle_0_2pi(alpha_unwrapped);
        trace.alpha_normalized = std::abs(alpha_unwrapped - trace.alpha) > 1e-12;

        const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trace.alpha,
            transfer_revolution,
            target_revolution);
        trace.residual_scale_free = residual.residual_scale_free;
        trace.residual_seconds = residual.residual_seconds;
        if (iteration == 0) {
            result.diagnostic.initial_residual_seconds = residual.residual_seconds;
        }
        if (!residual.valid || !is_finite(residual.residual_seconds) || !is_finite(residual.residual_scale_free)) {
            trace.reason = residual.invalid_reason.empty()
                ? "Newton residual evaluation failed"
                : residual.invalid_reason;
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }

        if (std::abs(residual.residual_seconds) <= residual_tolerance_seconds) {
            trace.reason = "residual tolerance reached";
            result.diagnostic.trace.push_back(trace);
            const Problem1SolutionBranch base_branch =
                make_solution_branch_from_public_residual(residual, transfer_revolution, target_revolution);
            Problem1SolutionBranch attached_branch = attach_problem1_root_derivatives(
                departure_planet,
                target_planet,
                nu_A_depart,
                nu_B_depart,
                theta_A,
                base_branch);
            if (attached_branch.valid && attached_branch.derivatives_available) {
                result.branch = attached_branch;
            } else {
                result.branch = base_branch;
                result.branch.derivatives_available = false;
                result.branch.invalid_reason = attached_branch.invalid_reason;
                result.diagnostic.derivative_attach_failed_after_convergence = true;
            }
            result.valid = result.branch.valid;
            result.diagnostic.converged = result.branch.valid;
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = result.branch.encounter_global_angle;
            result.diagnostic.final_residual_seconds = result.branch.residual_seconds;
            result.diagnostic.likely_wrong_root =
                std::abs(normalize_angle_minus_pi_pi(
                    result.branch.encounter_global_angle - result.diagnostic.initial_alpha)) > 0.2;
            return result;
        }

        const double abs_residual_scale_free = std::abs(residual.residual_scale_free);
        if (abs_residual_scale_free > previous_abs_residual_scale_free) {
            trace.residual_increased = true;
            trace.reason = "residual_increase_stop";
            result.diagnostic.residual_increased = true;
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }
        previous_abs_residual_scale_free = abs_residual_scale_free;

        const Problem1RootResidualDerivatives analytic_probe = evaluate_problem1_root_residual_derivatives_analytic(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trace.alpha,
            transfer_revolution,
            target_revolution);
        const Problem1RootResidualDerivatives derivatives = evaluate_problem1_root_residual_derivatives_with_mode(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            trace.alpha,
            transfer_revolution,
            target_revolution,
            derivative_mode,
            finite_difference_step);
        trace.derivative_valid =
            derivatives.valid && is_finite(derivatives.R_alpha) && std::abs(derivatives.R_alpha) > 1e-12;
        if (derivative_mode == Problem1RootDerivativeMode::FiniteDifferenceOnly) {
            trace.derivative_source = derivatives.valid ? "finite_difference" : "none";
        } else if (analytic_probe.valid) {
            trace.derivative_source = "analytic";
        } else if (derivatives.valid) {
            trace.derivative_source = "finite_difference_fallback";
        } else {
            trace.derivative_source = "none";
        }
        if (derivative_mode == Problem1RootDerivativeMode::FiniteDifferenceOnly && derivatives.valid) {
            result.diagnostic.finite_difference_success_count += 1;
        } else if (!analytic_probe.valid &&
                   derivative_mode == Problem1RootDerivativeMode::AnalyticWithFiniteDifferenceFallback &&
                   derivatives.valid) {
            result.diagnostic.fallback_used_count += 1;
            result.diagnostic.finite_difference_success_count += 1;
        }
        trace.R_alpha = derivatives.R_alpha;
        if (!trace.derivative_valid) {
            trace.reason = derivatives.invalid_reason.empty()
                ? "Newton derivative evaluation failed"
                : derivatives.invalid_reason;
            result.diagnostic.derivative_failed = true;
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }

        const double delta_alpha = -residual.residual_scale_free / derivatives.R_alpha;
        trace.delta_alpha = delta_alpha;
        if (!is_finite(delta_alpha)) {
            trace.reason = "Newton step is non-finite";
            result.diagnostic.trace.push_back(trace);
            result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
            result.diagnostic.final_alpha = trace.alpha;
            result.diagnostic.final_residual_seconds = residual.residual_seconds;
            result.diagnostic.invalid_reason = trace.reason;
            Problem1RootRefinementResult invalid = make_invalid_refinement_result(
                transfer_revolution,
                target_revolution,
                result.diagnostic.initial_alpha,
                trace.reason);
            invalid.diagnostic = result.diagnostic;
            return invalid;
        }

        alpha_unwrapped += delta_alpha;
        trace.reason = "continue";
        result.diagnostic.trace.push_back(trace);
    }

    result.diagnostic.iterations = static_cast<int>(result.diagnostic.trace.size());
    result.diagnostic.final_alpha = normalize_angle_0_2pi(alpha_unwrapped);
    result.diagnostic.final_residual_seconds = result.diagnostic.trace.empty()
        ? std::numeric_limits<double>::quiet_NaN()
        : result.diagnostic.trace.back().residual_seconds;
    result.diagnostic.invalid_reason = "Newton refinement exceeded maximum iterations";
    Problem1RootRefinementResult invalid = make_invalid_refinement_result(
        transfer_revolution,
        target_revolution,
        result.diagnostic.initial_alpha,
        result.diagnostic.invalid_reason);
    invalid.diagnostic = result.diagnostic;
    return invalid;
}

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
    Problem1RootDerivativeMode derivative_mode,
    double finite_difference_step
) {
    return refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        transfer_revolution,
        target_revolution,
        initial_encounter_global_angle,
        max_iterations,
        residual_tolerance_seconds,
        alpha_step_tolerance,
        derivative_mode,
        finite_difference_step).branch;
}

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
) {
    Problem1RootLinearPrediction prediction{};
    prediction.method = "linear_nearest";
    prediction.transfer_revolution = node_branch.transfer_revolution;
    prediction.target_revolution = node_branch.target_revolution;

    if (!node_branch.valid) {
        prediction.invalid_reason = "node branch is invalid";
        return prediction;
    }

    Problem1SolutionBranch differentiated = node_branch;
    if (!differentiated.derivatives_available) {
        differentiated = attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated.derivatives_available) {
        prediction.invalid_reason = "node branch derivatives are unavailable";
        return prediction;
    }

    const double d_nu_A = normalize_angle_minus_pi_pi(query_nu_A - node_nu_A);
    const double d_nu_B = normalize_angle_minus_pi_pi(query_nu_B - node_nu_B);
    const double d_theta_A = normalize_angle_minus_pi_pi(query_theta_A - node_theta_A);
    const double alpha_pred =
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * d_nu_A +
        differentiated.d_encounter_global_angle_d_nu_B * d_nu_B +
        differentiated.d_encounter_global_angle_d_theta_A * d_theta_A;
    if (!is_finite(alpha_pred)) {
        prediction.invalid_reason = "linear nearest prediction produced non-finite alpha";
        return prediction;
    }

    prediction.valid = true;
    prediction.invalid_reason.clear();
    prediction.predicted_encounter_global_angle = normalize_angle_0_2pi(alpha_pred);
    return prediction;
}

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
) {
    const Problem1RootLinearPrediction prediction = predict_problem1_root_branch_linear_from_node(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        node_branch,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!prediction.valid) {
        Problem1SolutionBranch invalid{};
        invalid.valid = false;
        invalid.transfer_revolution = node_branch.transfer_revolution;
        invalid.target_revolution = node_branch.target_revolution;
        invalid.invalid_reason = prediction.invalid_reason.empty()
            ? "linear nearest prediction failed"
            : prediction.invalid_reason;
        return invalid;
    }

    return refine_problem1_root_branch_newton(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        prediction.transfer_revolution,
        prediction.target_revolution,
        prediction.predicted_encounter_global_angle,
        max_iterations,
        residual_tolerance,
        alpha_step_tolerance);
}

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
) {
    return refine_problem1_root_branch_from_linear_prediction(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        node_branch,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        max_iterations,
        problem1_residual_seconds_to_scale_free(residual_tolerance_seconds),
        alpha_step_tolerance);
}

Problem1RootNearestNode find_nearest_problem1_root_table_node(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    Problem1RootNearestNode nearest{};
    const Problem1RootTableConfig& config = table.config();
    if (table.cells().empty() ||
        config.nu_A_count <= 0 ||
        config.nu_B_depart_count <= 0 ||
        config.theta_A_count <= 0) {
        nearest.invalid_reason = "root table is empty or has invalid dimensions";
        return nearest;
    }
    if (!is_finite(query_nu_A) || !is_finite(query_nu_B) || !is_finite(query_theta_A)) {
        nearest.invalid_reason = "query angles are not finite";
        return nearest;
    }

    const double wrapped_query_nu_A = normalize_angle_0_2pi(query_nu_A);
    const double wrapped_query_nu_B = normalize_angle_0_2pi(query_nu_B);
    const double wrapped_query_theta_A = normalize_angle_0_2pi(query_theta_A);
    double best_distance_squared = std::numeric_limits<double>::infinity();

    for (int i = 0; i < config.nu_A_count; ++i) {
        for (int j = 0; j < config.nu_B_depart_count; ++j) {
            for (int k = 0; k < config.theta_A_count; ++k) {
                const Problem1RootTableCell& cell = table.at(i, j, k);
                const double d_nu_A = normalize_angle_minus_pi_pi(wrapped_query_nu_A - cell.nu_A_depart);
                const double d_nu_B = normalize_angle_minus_pi_pi(wrapped_query_nu_B - cell.nu_B_depart);
                const double d_theta_A = normalize_angle_minus_pi_pi(wrapped_query_theta_A - cell.theta_A);
                const double distance_squared =
                    d_nu_A * d_nu_A + d_nu_B * d_nu_B + d_theta_A * d_theta_A;
                if (distance_squared < best_distance_squared) {
                    best_distance_squared = distance_squared;
                    nearest.valid = true;
                    nearest.i = i;
                    nearest.j = j;
                    nearest.k = k;
                    nearest.node_nu_A = cell.nu_A_depart;
                    nearest.node_nu_B = cell.nu_B_depart;
                    nearest.node_theta_A = cell.theta_A;
                    nearest.cell = &cell;
                    nearest.invalid_reason.clear();
                }
            }
        }
    }

    if (!nearest.valid) {
        nearest.invalid_reason = "failed to locate nearest root-table node";
    }
    return nearest;
}

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
) {
    Problem1RootApproximationResult result{};
    result.method = "linear_nearest_raw";
    result.transfer_revolution = node_branch.transfer_revolution;
    result.target_revolution = node_branch.target_revolution;

    const Problem1RootLinearPrediction prediction = predict_problem1_root_branch_linear_from_node(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        node_branch,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!prediction.valid) {
        result.invalid_reason = prediction.invalid_reason.empty()
            ? "linear nearest prediction failed"
            : prediction.invalid_reason;
        return result;
    }

    const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        prediction.predicted_encounter_global_angle,
        prediction.transfer_revolution,
        prediction.target_revolution);
    if (!residual.valid) {
        result.invalid_reason = residual.invalid_reason.empty()
            ? "linear nearest residual evaluation failed"
            : residual.invalid_reason;
        return result;
    }

    result.valid = true;
    result.predicted_encounter_global_angle = residual.encounter_global_angle;
    result.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    result.residual_scale_free = residual.residual_scale_free;
    result.residual_seconds = residual.residual_seconds;
    result.transfer_time_seconds = residual.transfer_time_seconds;
    result.target_time_seconds = residual.target_time_seconds;
    result.transfer_e = residual.transfer_e;
    result.transfer_p = residual.transfer_p;
    result.transfer_a = residual.transfer_a;
    result.theta_B = residual.theta_B;
    result.diagnostics.raw_residual_scale_free = residual.residual_scale_free;
    result.diagnostics.raw_residual_seconds = residual.residual_seconds;
    result.diagnostics.admissibility_reason = "Route A linear raw diagnostics not used for fast approximation";
    return result;
}

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
) {
    Problem1RootApproximationResult result{};
    result.method = "route_b_linear_no_newton";
    result.transfer_revolution = node_branch.transfer_revolution;
    result.target_revolution = node_branch.target_revolution;
    result.diagnostics.source_target_revolution = node_branch.target_revolution;
    result.diagnostics.selected_target_revolution = node_branch.target_revolution;
    result.diagnostics.hessian_method = "none";

    const Problem1RootLinearPrediction prediction = predict_problem1_root_branch_linear_from_node(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        node_branch,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!prediction.valid) {
        result.invalid_reason = prediction.invalid_reason.empty()
            ? "route_b linear prediction failed"
            : prediction.invalid_reason;
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }
    result.diagnostics.alpha_linear = prediction.predicted_encounter_global_angle;
    result.diagnostics.alpha_quadratic = prediction.predicted_encounter_global_angle;

    const Problem1RootQSheetSelectionResult q_selection = select_q_by_target_time_sheet_continuity(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        prediction.transfer_revolution,
        prediction.predicted_encounter_global_angle,
        node_branch,
        max_target_revolution);
    const int selected_q = q_selection.selection_failed
        ? node_branch.target_revolution
        : q_selection.selected_q;
    result.diagnostics.selected_target_revolution = selected_q;
    result.diagnostics.q_sheet_selection_changed = q_selection.q_changed;
    result.diagnostics.q_sheet_selection_failed = q_selection.selection_failed;
    result.diagnostics.selected_q_continuity_error = q_selection.selected_continuity_error;
    result.diagnostics.source_q_continuity_error = q_selection.source_q_continuity_error;

    const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        prediction.predicted_encounter_global_angle,
        prediction.transfer_revolution,
        selected_q);
    if (!residual.valid) {
        result.invalid_reason = residual.invalid_reason.empty()
            ? "route_b linear residual evaluation failed"
            : residual.invalid_reason;
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    result.valid = true;
    result.target_revolution = selected_q;
    result.predicted_encounter_global_angle = residual.encounter_global_angle;
    result.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    result.residual_scale_free = residual.residual_scale_free;
    result.residual_seconds = residual.residual_seconds;
    result.transfer_time_seconds = residual.transfer_time_seconds;
    result.target_time_seconds = residual.target_time_seconds;
    result.transfer_e = residual.transfer_e;
    result.transfer_p = residual.transfer_p;
    result.transfer_a = residual.transfer_a;
    result.theta_B = residual.theta_B;
    result.diagnostics.raw_residual_scale_free = residual.residual_scale_free;
    result.diagnostics.raw_residual_seconds = residual.residual_seconds;
    result.diagnostics.admissible_for_fast_approximation =
        is_finite(residual.residual_seconds) && is_finite(residual.residual_scale_free);
    result.diagnostics.admissibility_reason = result.diagnostics.admissible_for_fast_approximation
        ? "admissible"
        : "route_b linear residual is non-finite";
    return result;
}

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
    Problem1RootHessianMethod hessian_method,
    double tangent_residual_tolerance
) {
    Problem1RootApproximationResult result{};
    result.method = "route_b_quadratic_no_newton";
    result.transfer_revolution = node_branch.transfer_revolution;
    result.target_revolution = node_branch.target_revolution;
    result.diagnostics.source_target_revolution = node_branch.target_revolution;
    result.diagnostics.selected_target_revolution = node_branch.target_revolution;
    result.diagnostics.hessian_method = problem1_root_hessian_method_name(hessian_method);
    result.diagnostics.hessian_step = hessian_step;

    if (hessian_method == Problem1RootHessianMethod::NewtonRefinedFiniteDifference) {
        result.invalid_reason = "route_b_quadratic_newton_hessian_method_disabled";
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }
    if (!node_branch.valid) {
        result.invalid_reason = "node branch is invalid";
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    Problem1SolutionBranch differentiated = node_branch;
    if (!differentiated.derivatives_available) {
        differentiated = attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated.derivatives_available) {
        result.invalid_reason = "node branch derivatives are unavailable";
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - node_theta_A);
    const double alpha_linear =
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated.d_encounter_global_angle_d_theta_A * dx3;
    result.diagnostics.alpha_linear = alpha_linear;
    if (!is_finite(alpha_linear)) {
        result.invalid_reason = "route_b quadratic linear part produced non-finite alpha";
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    Problem1RootHessian hessian{};
    if (hessian_method == Problem1RootHessianMethod::TangentFiniteDifference) {
        hessian = estimate_problem1_root_hessian_tangent_finite_difference(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            differentiated,
            hessian_step,
            tangent_residual_tolerance);
    } else if (hessian_method == Problem1RootHessianMethod::ProjectedTangentFiniteDifference) {
        hessian = estimate_problem1_root_hessian_projected_tangent_finite_difference(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            differentiated,
            hessian_step,
            1e-2,
            1e-12);
    } else {
        result.invalid_reason = "route_b_quadratic_unsupported_hessian_method";
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }
    result.diagnostics.hessian_valid = hessian.valid;
    result.diagnostics.tangent_residual_max_scale_free = hessian.tangent_residual_max_scale_free;
    result.diagnostics.tangent_residual_max_seconds = hessian.tangent_residual_max_seconds;
    if (!hessian.valid) {
        result.invalid_reason = hessian.invalid_reason.empty()
            ? "route_b quadratic Hessian estimate failed"
            : hessian.invalid_reason;
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    const double alpha_quadratic_unwrapped =
        alpha_linear +
        0.5 * (
            hessian.H_nu_A_nu_A * dx1 * dx1 +
            hessian.H_nu_B_nu_B * dx2 * dx2 +
            hessian.H_theta_A_theta_A * dx3 * dx3 +
            2.0 * hessian.H_nu_A_nu_B * dx1 * dx2 +
            2.0 * hessian.H_nu_A_theta_A * dx1 * dx3 +
            2.0 * hessian.H_nu_B_theta_A * dx2 * dx3);
    if (!is_finite(alpha_quadratic_unwrapped)) {
        result.invalid_reason = "route_b quadratic prediction produced non-finite alpha";
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }
    const double alpha_quadratic = normalize_angle_0_2pi(alpha_quadratic_unwrapped);
    result.diagnostics.alpha_quadratic = alpha_quadratic;

    const Problem1RootQSheetSelectionResult q_selection = select_q_by_target_time_sheet_continuity(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        node_branch.transfer_revolution,
        alpha_quadratic,
        node_branch,
        max_target_revolution);
    const int selected_q = q_selection.selection_failed
        ? node_branch.target_revolution
        : q_selection.selected_q;
    result.diagnostics.selected_target_revolution = selected_q;
    result.diagnostics.q_sheet_selection_changed = q_selection.q_changed;
    result.diagnostics.q_sheet_selection_failed = q_selection.selection_failed;
    result.diagnostics.selected_q_continuity_error = q_selection.selected_continuity_error;
    result.diagnostics.source_q_continuity_error = q_selection.source_q_continuity_error;

    const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        alpha_quadratic,
        node_branch.transfer_revolution,
        selected_q);
    if (!residual.valid) {
        result.invalid_reason = residual.invalid_reason.empty()
            ? "route_b quadratic residual evaluation failed"
            : residual.invalid_reason;
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    result.valid = true;
    result.target_revolution = selected_q;
    result.predicted_encounter_global_angle = residual.encounter_global_angle;
    result.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    result.residual_scale_free = residual.residual_scale_free;
    result.residual_seconds = residual.residual_seconds;
    result.transfer_time_seconds = residual.transfer_time_seconds;
    result.target_time_seconds = residual.target_time_seconds;
    result.transfer_e = residual.transfer_e;
    result.transfer_p = residual.transfer_p;
    result.transfer_a = residual.transfer_a;
    result.theta_B = residual.theta_B;
    result.diagnostics.raw_residual_scale_free = residual.residual_scale_free;
    result.diagnostics.raw_residual_seconds = residual.residual_seconds;
    result.diagnostics.admissible_for_fast_approximation =
        is_finite(residual.residual_seconds) && is_finite(residual.residual_scale_free);
    result.diagnostics.admissibility_reason = result.diagnostics.admissible_for_fast_approximation
        ? "admissible"
        : "route_b quadratic residual is non-finite";
    return result;
}

Problem1RootCellAdmissibilityResult evaluate_problem1_root_cell_admissibility(
    const Problem1RootTable& table,
    int nu_A_index,
    int nu_B_depart_index,
    int theta_A_index,
    int max_transfer_revolution,
    int max_target_revolution
) {
    Problem1RootCellAdmissibilityResult result{};
    result.corner_count = 8;
    const Problem1RootTableConfig& config = table.config();
    if (max_transfer_revolution < 0 || max_target_revolution < 0) {
        result.reason = "invalid_max_revolution_inputs";
        return result;
    }
    if (nu_A_index < 0 || nu_B_depart_index < 0 || theta_A_index < 0 ||
        nu_A_index + 1 >= config.nu_A_count ||
        nu_B_depart_index + 1 >= config.nu_B_depart_count ||
        theta_A_index + 1 >= config.theta_A_count) {
        result.reason = "cell_index_out_of_range";
        return result;
    }

    for (int k = 0; k <= max_transfer_revolution; ++k) {
        result.reference_root_count_by_k[k] = 0;
        result.min_root_count_by_k[k] = std::numeric_limits<int>::max();
        result.max_root_count_by_k[k] = std::numeric_limits<int>::lowest();
    }

    for (int dnu_A = 0; dnu_A <= 1; ++dnu_A) {
        for (int dnu_B = 0; dnu_B <= 1; ++dnu_B) {
            for (int dtheta_A = 0; dtheta_A <= 1; ++dtheta_A) {
                const Problem1RootTableCell& corner = table.at(
                    nu_A_index + dnu_A,
                    nu_B_depart_index + dnu_B,
                    theta_A_index + dtheta_A);
                std::map<int, int> count_by_k;
                for (const Problem1SolutionBranch& branch : corner.solutions_sorted_by_time_of_flight) {
                    if (!branch.valid) {
                        continue;
                    }
                    if (branch.transfer_revolution < 0 || branch.transfer_revolution > max_transfer_revolution) {
                        continue;
                    }
                    if (branch.target_revolution < 0 || branch.target_revolution > max_target_revolution) {
                        continue;
                    }
                    count_by_k[branch.transfer_revolution] += 1;
                }
                for (int k = 0; k <= max_transfer_revolution; ++k) {
                    const int count = count_by_k.count(k) > 0 ? count_by_k.at(k) : 0;
                    if (dnu_A == 0 && dnu_B == 0 && dtheta_A == 0) {
                        result.reference_root_count_by_k[k] = count;
                    }
                    result.min_root_count_by_k[k] = std::min(result.min_root_count_by_k[k], count);
                    result.max_root_count_by_k[k] = std::max(result.max_root_count_by_k[k], count);
                }
            }
        }
    }

    for (int k = 0; k <= max_transfer_revolution; ++k) {
        if (result.min_root_count_by_k[k] != result.max_root_count_by_k[k]) {
            result.admissible = false;
            result.reason = "branch_count_inconsistent_across_corners";
            return result;
        }
    }

    result.admissible = true;
    result.reason = "corner_branch_count_consistent";
    return result;
}

Problem1RootHessian estimate_problem1_root_hessian_finite_difference(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double step
) {
    Problem1RootHessian hessian{};
    hessian.method = "newton_refined_finite_difference_of_implicit_first_derivatives";
    if (!is_finite(step) || !(step > 0.0)) {
        hessian.invalid_reason = "finite-difference Hessian step must be positive";
        return hessian;
    }
    if (!node_branch.valid) {
        hessian.invalid_reason = "node branch is invalid";
        return hessian;
    }

    Problem1SolutionBranch differentiated_node = node_branch;
    if (!differentiated_node.derivatives_available) {
        differentiated_node = attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated_node.derivatives_available) {
        hessian.invalid_reason = "node branch derivatives are unavailable";
        return hessian;
    }

    const auto refine_and_attach = [&](double query_nu_A, double query_nu_B, double query_theta_A) {
        Problem1SolutionBranch refined = refine_problem1_root_branch_from_linear_prediction(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            differentiated_node,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            30,
            1e-10,
            1e-12);
        if (!refined.valid) {
            return refined;
        }
        return attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            refined);
    };

    const auto plus_minus = [&](int axis) {
        double plus_nu_A = node_nu_A;
        double plus_nu_B = node_nu_B;
        double plus_theta_A = node_theta_A;
        double minus_nu_A = node_nu_A;
        double minus_nu_B = node_nu_B;
        double minus_theta_A = node_theta_A;
        if (axis == 0) {
            plus_nu_A = normalize_angle_0_2pi(node_nu_A + step);
            minus_nu_A = normalize_angle_0_2pi(node_nu_A - step);
        } else if (axis == 1) {
            plus_nu_B = normalize_angle_0_2pi(node_nu_B + step);
            minus_nu_B = normalize_angle_0_2pi(node_nu_B - step);
        } else {
            plus_theta_A = normalize_angle_0_2pi(node_theta_A + step);
            minus_theta_A = normalize_angle_0_2pi(node_theta_A - step);
        }
        return std::pair{
            refine_and_attach(plus_nu_A, plus_nu_B, plus_theta_A),
            refine_and_attach(minus_nu_A, minus_nu_B, minus_theta_A)};
    };

    const auto [nu_A_plus, nu_A_minus] = plus_minus(0);
    const auto [nu_B_plus, nu_B_minus] = plus_minus(1);
    const auto [theta_A_plus, theta_A_minus] = plus_minus(2);
    const auto all_valid = [&](const Problem1SolutionBranch& a, const Problem1SolutionBranch& b) {
        return a.valid && b.valid && a.derivatives_available && b.derivatives_available;
    };
    if (!all_valid(nu_A_plus, nu_A_minus) ||
        !all_valid(nu_B_plus, nu_B_minus) ||
        !all_valid(theta_A_plus, theta_A_minus)) {
        hessian.invalid_reason = "failed to refine/attach derivatives at Hessian stencil points";
        return hessian;
    }

    const double inv_2h = 1.0 / (2.0 * step);
    hessian.H_nu_A_nu_A =
        (nu_A_plus.d_encounter_global_angle_d_nu_A - nu_A_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    hessian.H_nu_B_nu_B =
        (nu_B_plus.d_encounter_global_angle_d_nu_B - nu_B_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    hessian.H_theta_A_theta_A =
        (theta_A_plus.d_encounter_global_angle_d_theta_A -
         theta_A_minus.d_encounter_global_angle_d_theta_A) * inv_2h;

    const double H12_from_g1 =
        (nu_B_plus.d_encounter_global_angle_d_nu_A - nu_B_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 =
        (nu_A_plus.d_encounter_global_angle_d_nu_B - nu_A_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 =
        (theta_A_plus.d_encounter_global_angle_d_nu_A - theta_A_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 =
        (nu_A_plus.d_encounter_global_angle_d_theta_A - nu_A_minus.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H23_from_g2 =
        (theta_A_plus.d_encounter_global_angle_d_nu_B - theta_A_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H32_from_g3 =
        (nu_B_plus.d_encounter_global_angle_d_theta_A - nu_B_minus.d_encounter_global_angle_d_theta_A) * inv_2h;

    hessian.H_nu_A_nu_B = 0.5 * (H12_from_g1 + H21_from_g2);
    hessian.H_nu_A_theta_A = 0.5 * (H13_from_g1 + H31_from_g3);
    hessian.H_nu_B_theta_A = 0.5 * (H23_from_g2 + H32_from_g3);

    if (!is_finite(hessian.H_nu_A_nu_A) ||
        !is_finite(hessian.H_nu_B_nu_B) ||
        !is_finite(hessian.H_theta_A_theta_A) ||
        !is_finite(hessian.H_nu_A_nu_B) ||
        !is_finite(hessian.H_nu_A_theta_A) ||
        !is_finite(hessian.H_nu_B_theta_A)) {
        hessian.invalid_reason = "finite-difference Hessian contains non-finite entries";
        return hessian;
    }

    hessian.valid = true;
    return hessian;
}

Problem1RootHessian estimate_problem1_root_hessian_tangent_finite_difference(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double node_nu_A,
    double node_nu_B,
    double node_theta_A,
    const Problem1SolutionBranch& node_branch,
    double step,
    double tangent_residual_tolerance
) {
    Problem1RootHessian hessian{};
    hessian.method = "tangent_finite_difference_of_implicit_first_derivatives";
    if (!is_finite(step) || !(step > 0.0)) {
        hessian.invalid_reason = "tangent finite-difference Hessian step must be positive";
        return hessian;
    }
    if (!is_finite(tangent_residual_tolerance) || !(tangent_residual_tolerance > 0.0)) {
        hessian.invalid_reason = "tangent residual tolerance must be positive";
        return hessian;
    }
    if (!node_branch.valid) {
        hessian.invalid_reason = "node branch is invalid";
        return hessian;
    }

    Problem1SolutionBranch differentiated_node = node_branch;
    if (!differentiated_node.derivatives_available) {
        differentiated_node = attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated_node.derivatives_available) {
        hessian.invalid_reason = "node branch derivatives are unavailable";
        return hessian;
    }

    const double g1 = differentiated_node.d_encounter_global_angle_d_nu_A;
    const double g2 = differentiated_node.d_encounter_global_angle_d_nu_B;
    const double g3 = differentiated_node.d_encounter_global_angle_d_theta_A;
    double tangent_residual_max_scale_free = std::numeric_limits<double>::quiet_NaN();
    double tangent_residual_max_seconds = std::numeric_limits<double>::quiet_NaN();

    const auto evaluate_tangent_point = [&](double query_nu_A, double query_nu_B, double query_theta_A) {
        Problem1SolutionBranch invalid{};
        invalid.valid = false;
        invalid.transfer_revolution = differentiated_node.transfer_revolution;
        invalid.target_revolution = differentiated_node.target_revolution;

        const double d_nu_A = normalize_angle_minus_pi_pi(query_nu_A - node_nu_A);
        const double d_nu_B = normalize_angle_minus_pi_pi(query_nu_B - node_nu_B);
        const double d_theta_A = normalize_angle_minus_pi_pi(query_theta_A - node_theta_A);
        const double alpha_tangent = normalize_angle_0_2pi(
            differentiated_node.encounter_global_angle +
            g1 * d_nu_A + g2 * d_nu_B + g3 * d_theta_A);
        const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_tangent,
            differentiated_node.transfer_revolution,
            differentiated_node.target_revolution);
        if (residual.valid) {
            update_abs_max_if_finite(residual.residual_scale_free, &tangent_residual_max_scale_free);
            update_abs_max_if_finite(residual.residual_seconds, &tangent_residual_max_seconds);
        }
        if (!residual.valid || !is_finite(residual.residual_scale_free) ||
            std::abs(residual.residual_scale_free) > tangent_residual_tolerance) {
            invalid.invalid_reason = residual.invalid_reason.empty()
                ? "tangent residual is too large"
                : residual.invalid_reason;
            return invalid;
        }

        Problem1SolutionBranch tangent_branch{};
        tangent_branch.valid = true;
        tangent_branch.encounter_global_angle = residual.encounter_global_angle;
        tangent_branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
        tangent_branch.transfer_revolution = differentiated_node.transfer_revolution;
        tangent_branch.target_revolution = differentiated_node.target_revolution;
        tangent_branch.time_of_flight_seconds = residual.transfer_time_seconds;
        tangent_branch.target_time_seconds = residual.target_time_seconds;
        tangent_branch.residual_seconds = residual.residual_seconds;
        tangent_branch.arrival_time_seconds_since_j2000 = std::numeric_limits<double>::quiet_NaN();
        tangent_branch.transfer_e = residual.transfer_e;
        tangent_branch.transfer_p = residual.transfer_p;
        tangent_branch.transfer_a = residual.transfer_a;
        tangent_branch.theta_B = residual.theta_B;
        return attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            tangent_branch);
    };

    const auto plus_minus = [&](int axis) {
        double plus_nu_A = node_nu_A;
        double plus_nu_B = node_nu_B;
        double plus_theta_A = node_theta_A;
        double minus_nu_A = node_nu_A;
        double minus_nu_B = node_nu_B;
        double minus_theta_A = node_theta_A;
        if (axis == 0) {
            plus_nu_A = normalize_angle_0_2pi(node_nu_A + step);
            minus_nu_A = normalize_angle_0_2pi(node_nu_A - step);
        } else if (axis == 1) {
            plus_nu_B = normalize_angle_0_2pi(node_nu_B + step);
            minus_nu_B = normalize_angle_0_2pi(node_nu_B - step);
        } else {
            plus_theta_A = normalize_angle_0_2pi(node_theta_A + step);
            minus_theta_A = normalize_angle_0_2pi(node_theta_A - step);
        }
        return std::pair{
            evaluate_tangent_point(plus_nu_A, plus_nu_B, plus_theta_A),
            evaluate_tangent_point(minus_nu_A, minus_nu_B, minus_theta_A)};
    };

    const auto [nu_A_plus, nu_A_minus] = plus_minus(0);
    const auto [nu_B_plus, nu_B_minus] = plus_minus(1);
    const auto [theta_A_plus, theta_A_minus] = plus_minus(2);
    hessian.tangent_residual_max_scale_free = tangent_residual_max_scale_free;
    hessian.tangent_residual_max_seconds = tangent_residual_max_seconds;

    const auto all_valid = [&](const Problem1SolutionBranch& a, const Problem1SolutionBranch& b) {
        return a.valid && b.valid && a.derivatives_available && b.derivatives_available;
    };
    if (!all_valid(nu_A_plus, nu_A_minus) ||
        !all_valid(nu_B_plus, nu_B_minus) ||
        !all_valid(theta_A_plus, theta_A_minus)) {
        hessian.invalid_reason = "failed to attach derivatives at tangent Hessian stencil points";
        return hessian;
    }

    const double inv_2h = 1.0 / (2.0 * step);
    hessian.H_nu_A_nu_A =
        (nu_A_plus.d_encounter_global_angle_d_nu_A - nu_A_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    hessian.H_nu_B_nu_B =
        (nu_B_plus.d_encounter_global_angle_d_nu_B - nu_B_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    hessian.H_theta_A_theta_A =
        (theta_A_plus.d_encounter_global_angle_d_theta_A -
         theta_A_minus.d_encounter_global_angle_d_theta_A) * inv_2h;

    const double H12_from_g1 =
        (nu_B_plus.d_encounter_global_angle_d_nu_A - nu_B_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 =
        (nu_A_plus.d_encounter_global_angle_d_nu_B - nu_A_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 =
        (theta_A_plus.d_encounter_global_angle_d_nu_A - theta_A_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 =
        (nu_A_plus.d_encounter_global_angle_d_theta_A - nu_A_minus.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H23_from_g2 =
        (theta_A_plus.d_encounter_global_angle_d_nu_B - theta_A_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H32_from_g3 =
        (nu_B_plus.d_encounter_global_angle_d_theta_A - nu_B_minus.d_encounter_global_angle_d_theta_A) * inv_2h;

    hessian.H_nu_A_nu_B = 0.5 * (H12_from_g1 + H21_from_g2);
    hessian.H_nu_A_theta_A = 0.5 * (H13_from_g1 + H31_from_g3);
    hessian.H_nu_B_theta_A = 0.5 * (H23_from_g2 + H32_from_g3);

    if (!is_finite(hessian.H_nu_A_nu_A) ||
        !is_finite(hessian.H_nu_B_nu_B) ||
        !is_finite(hessian.H_theta_A_theta_A) ||
        !is_finite(hessian.H_nu_A_nu_B) ||
        !is_finite(hessian.H_nu_A_theta_A) ||
        !is_finite(hessian.H_nu_B_theta_A)) {
        hessian.invalid_reason = "tangent finite-difference Hessian contains non-finite entries";
        return hessian;
    }

    hessian.valid = true;
    return hessian;
}

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
) {
    Problem1RootHessian hessian{};
    hessian.method = "projected_tangent_finite_difference_of_implicit_first_derivatives";
    if (!is_finite(step) || !(step > 0.0)) {
        hessian.invalid_reason = "projected tangent finite-difference Hessian step must be positive";
        return hessian;
    }
    if (!is_finite(residual_tolerance_seconds) || !(residual_tolerance_seconds > 0.0)) {
        hessian.invalid_reason = "projected stencil residual tolerance seconds must be positive";
        return hessian;
    }
    if (!is_finite(projection_derivative_min_abs) || !(projection_derivative_min_abs > 0.0)) {
        hessian.invalid_reason = "projected stencil derivative minimum absolute value must be positive";
        return hessian;
    }
    if (!node_branch.valid) {
        hessian.invalid_reason = "node branch is invalid";
        return hessian;
    }

    Problem1SolutionBranch differentiated_node = node_branch;
    if (!differentiated_node.derivatives_available) {
        differentiated_node = attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated_node.derivatives_available) {
        hessian.invalid_reason = "node branch derivatives are unavailable";
        return hessian;
    }

    double projected_residual_max_scale_free = std::numeric_limits<double>::quiet_NaN();
    double projected_residual_max_seconds = std::numeric_limits<double>::quiet_NaN();

    const auto evaluate_projected_point = [&](double query_nu_A, double query_nu_B, double query_theta_A) {
        const Problem1ProjectedTangentStencilBranchResult projected =
            evaluate_problem1_root_projected_tangent_stencil_branch(
                departure_planet,
                target_planet,
                node_nu_A,
                node_nu_B,
                node_theta_A,
                differentiated_node,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                residual_tolerance_seconds,
                projection_derivative_min_abs);
        if (is_finite(projected.residual_after_projection_seconds)) {
            update_abs_max_if_finite(
                problem1_residual_seconds_to_scale_free(projected.residual_after_projection_seconds),
                &projected_residual_max_scale_free);
            update_abs_max_if_finite(
                projected.residual_after_projection_seconds,
                &projected_residual_max_seconds);
        }
        return projected.branch;
    };

    const auto plus_minus = [&](int axis) {
        double plus_nu_A = node_nu_A;
        double plus_nu_B = node_nu_B;
        double plus_theta_A = node_theta_A;
        double minus_nu_A = node_nu_A;
        double minus_nu_B = node_nu_B;
        double minus_theta_A = node_theta_A;
        if (axis == 0) {
            plus_nu_A = normalize_angle_0_2pi(node_nu_A + step);
            minus_nu_A = normalize_angle_0_2pi(node_nu_A - step);
        } else if (axis == 1) {
            plus_nu_B = normalize_angle_0_2pi(node_nu_B + step);
            minus_nu_B = normalize_angle_0_2pi(node_nu_B - step);
        } else {
            plus_theta_A = normalize_angle_0_2pi(node_theta_A + step);
            minus_theta_A = normalize_angle_0_2pi(node_theta_A - step);
        }
        return std::pair{
            evaluate_projected_point(plus_nu_A, plus_nu_B, plus_theta_A),
            evaluate_projected_point(minus_nu_A, minus_nu_B, minus_theta_A)};
    };

    const auto [nu_A_plus, nu_A_minus] = plus_minus(0);
    const auto [nu_B_plus, nu_B_minus] = plus_minus(1);
    const auto [theta_A_plus, theta_A_minus] = plus_minus(2);
    hessian.tangent_residual_max_scale_free = projected_residual_max_scale_free;
    hessian.tangent_residual_max_seconds = projected_residual_max_seconds;

    const auto all_valid = [&](const Problem1SolutionBranch& a, const Problem1SolutionBranch& b) {
        return a.valid && b.valid && a.derivatives_available && b.derivatives_available;
    };
    if (!all_valid(nu_A_plus, nu_A_minus) ||
        !all_valid(nu_B_plus, nu_B_minus) ||
        !all_valid(theta_A_plus, theta_A_minus)) {
        hessian.invalid_reason = "failed to attach derivatives at projected tangent Hessian stencil points";
        return hessian;
    }

    const double inv_2h = 1.0 / (2.0 * step);
    hessian.H_nu_A_nu_A =
        (nu_A_plus.d_encounter_global_angle_d_nu_A - nu_A_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    hessian.H_nu_B_nu_B =
        (nu_B_plus.d_encounter_global_angle_d_nu_B - nu_B_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    hessian.H_theta_A_theta_A =
        (theta_A_plus.d_encounter_global_angle_d_theta_A -
         theta_A_minus.d_encounter_global_angle_d_theta_A) * inv_2h;

    const double H12_from_g1 =
        (nu_B_plus.d_encounter_global_angle_d_nu_A - nu_B_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 =
        (nu_A_plus.d_encounter_global_angle_d_nu_B - nu_A_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 =
        (theta_A_plus.d_encounter_global_angle_d_nu_A -
         theta_A_minus.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 =
        (nu_A_plus.d_encounter_global_angle_d_theta_A - nu_A_minus.d_encounter_global_angle_d_theta_A) * inv_2h;
    const double H23_from_g2 =
        (theta_A_plus.d_encounter_global_angle_d_nu_B -
         theta_A_minus.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H32_from_g3 =
        (nu_B_plus.d_encounter_global_angle_d_theta_A - nu_B_minus.d_encounter_global_angle_d_theta_A) * inv_2h;

    hessian.H_nu_A_nu_B = 0.5 * (H12_from_g1 + H21_from_g2);
    hessian.H_nu_A_theta_A = 0.5 * (H13_from_g1 + H31_from_g3);
    hessian.H_nu_B_theta_A = 0.5 * (H23_from_g2 + H32_from_g3);

    if (!is_finite(hessian.H_nu_A_nu_A) ||
        !is_finite(hessian.H_nu_B_nu_B) ||
        !is_finite(hessian.H_theta_A_theta_A) ||
        !is_finite(hessian.H_nu_A_nu_B) ||
        !is_finite(hessian.H_nu_A_theta_A) ||
        !is_finite(hessian.H_nu_B_theta_A)) {
        hessian.invalid_reason = "finite-difference Hessian contains non-finite entries";
        return hessian;
    }

    hessian.valid = true;
    return hessian;
}

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
    Problem1RootHessianMethod hessian_method,
    double tangent_residual_tolerance
) {
    Problem1RootQuadraticPrediction prediction{};
    prediction.method = "quadratic_nearest_raw";
    prediction.transfer_revolution = node_branch.transfer_revolution;
    prediction.target_revolution = node_branch.target_revolution;
    prediction.hessian_method = problem1_root_hessian_method_name(hessian_method);
    prediction.hessian_step = hessian_step;
    if (!is_finite(tangent_residual_tolerance) || !(tangent_residual_tolerance > 0.0)) {
        prediction.invalid_reason = "quadratic nearest tangent residual tolerance must be positive";
        return prediction;
    }
    if (!node_branch.valid) {
        prediction.invalid_reason = "node branch is invalid";
        return prediction;
    }

    Problem1SolutionBranch differentiated = node_branch;
    if (!differentiated.derivatives_available) {
        differentiated = attach_problem1_root_derivatives(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            node_branch);
    }
    if (!differentiated.derivatives_available) {
        prediction.invalid_reason = "node branch derivatives are unavailable";
        return prediction;
    }

    const double dx1 = normalize_angle_minus_pi_pi(query_nu_A - node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(query_nu_B - node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(query_theta_A - node_theta_A);

    Problem1RootHessian hessian{};
    if (hessian_method == Problem1RootHessianMethod::TangentFiniteDifference) {
        hessian = estimate_problem1_root_hessian_tangent_finite_difference(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            differentiated,
            hessian_step,
            tangent_residual_tolerance);
    } else if (hessian_method == Problem1RootHessianMethod::ProjectedTangentFiniteDifference) {
        hessian = estimate_problem1_root_hessian_projected_tangent_finite_difference(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            differentiated,
            hessian_step,
            1e-2,
            1e-12);
    } else {
        hessian = estimate_problem1_root_hessian_finite_difference(
            departure_planet,
            target_planet,
            node_nu_A,
            node_nu_B,
            node_theta_A,
            differentiated,
            hessian_step);
    }
    prediction.hessian_valid = hessian.valid;
    prediction.tangent_residual_max_scale_free = hessian.tangent_residual_max_scale_free;
    prediction.tangent_residual_max_seconds = hessian.tangent_residual_max_seconds;
    if (!hessian.valid) {
        prediction.invalid_reason = hessian.invalid_reason.empty()
            ? "quadratic nearest Hessian estimate failed"
            : hessian.invalid_reason;
        return prediction;
    }

    const double alpha_quad =
        differentiated.encounter_global_angle +
        differentiated.d_encounter_global_angle_d_nu_A * dx1 +
        differentiated.d_encounter_global_angle_d_nu_B * dx2 +
        differentiated.d_encounter_global_angle_d_theta_A * dx3 +
        0.5 * (
            hessian.H_nu_A_nu_A * dx1 * dx1 +
            hessian.H_nu_B_nu_B * dx2 * dx2 +
            hessian.H_theta_A_theta_A * dx3 * dx3 +
            2.0 * hessian.H_nu_A_nu_B * dx1 * dx2 +
            2.0 * hessian.H_nu_A_theta_A * dx1 * dx3 +
            2.0 * hessian.H_nu_B_theta_A * dx2 * dx3);
    if (!is_finite(alpha_quad)) {
        prediction.invalid_reason = "quadratic nearest prediction produced non-finite alpha";
        return prediction;
    }

    prediction.valid = true;
    prediction.predicted_encounter_global_angle = normalize_angle_0_2pi(alpha_quad);
    return prediction;
}

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
    Problem1RootHessianMethod hessian_method,
    double tangent_residual_tolerance,
    int max_target_revolution
) {
    Problem1RootApproximationResult result{};
    result.method = "quadratic_nearest_raw";
    result.transfer_revolution = node_branch.transfer_revolution;
    result.target_revolution = node_branch.target_revolution;
    result.diagnostics.source_target_revolution = node_branch.target_revolution;
    result.diagnostics.selected_target_revolution = node_branch.target_revolution;
    result.diagnostics.hessian_method = problem1_root_hessian_method_name(hessian_method);
    result.diagnostics.hessian_step = hessian_step;

    const Problem1RootQuadraticPrediction prediction = predict_problem1_root_branch_quadratic_from_node(
        departure_planet,
        target_planet,
        node_nu_A,
        node_nu_B,
        node_theta_A,
        node_branch,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        hessian_step,
        hessian_method,
        tangent_residual_tolerance);
    result.diagnostics.hessian_method = prediction.hessian_method;
    result.diagnostics.hessian_valid = prediction.hessian_valid;
    result.diagnostics.tangent_residual_max_scale_free = prediction.tangent_residual_max_scale_free;
    result.diagnostics.tangent_residual_max_seconds = prediction.tangent_residual_max_seconds;
    result.diagnostics.hessian_step = prediction.hessian_step;
    if (!prediction.valid) {
        result.invalid_reason = prediction.invalid_reason.empty()
            ? "quadratic nearest prediction failed"
            : prediction.invalid_reason;
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    const Problem1RootQSheetSelectionResult q_selection = select_q_by_target_time_sheet_continuity(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        prediction.transfer_revolution,
        prediction.predicted_encounter_global_angle,
        node_branch,
        std::max(node_branch.target_revolution, max_target_revolution));
    const int selected_q = q_selection.selection_failed
        ? node_branch.target_revolution
        : q_selection.selected_q;
    result.diagnostics.selected_target_revolution = selected_q;
    result.diagnostics.q_sheet_selection_changed = q_selection.q_changed;
    result.diagnostics.q_sheet_selection_failed = q_selection.selection_failed;
    result.diagnostics.selected_q_continuity_error = q_selection.selected_continuity_error;
    result.diagnostics.source_q_continuity_error = q_selection.source_q_continuity_error;

    const Problem1RootResidualResult residual = evaluate_problem1_root_residual(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        prediction.predicted_encounter_global_angle,
        prediction.transfer_revolution,
        selected_q);
    if (!residual.valid) {
        result.invalid_reason = residual.invalid_reason.empty()
            ? "quadratic nearest residual evaluation failed"
            : residual.invalid_reason;
        result.diagnostics.admissibility_reason = result.invalid_reason;
        return result;
    }

    result.valid = true;
    result.target_revolution = selected_q;
    result.predicted_encounter_global_angle = residual.encounter_global_angle;
    result.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    result.residual_scale_free = residual.residual_scale_free;
    result.residual_seconds = residual.residual_seconds;
    result.transfer_time_seconds = residual.transfer_time_seconds;
    result.target_time_seconds = residual.target_time_seconds;
    result.transfer_e = residual.transfer_e;
    result.transfer_p = residual.transfer_p;
    result.transfer_a = residual.transfer_a;
    result.theta_B = residual.theta_B;
    result.diagnostics.raw_residual_scale_free = residual.residual_scale_free;
    result.diagnostics.raw_residual_seconds = residual.residual_seconds;

    const bool raw_residual_is_finite =
        is_finite(residual.residual_scale_free) && is_finite(residual.residual_seconds);
    const bool tangent_hessian_gate_passed =
        hessian_method != Problem1RootHessianMethod::TangentFiniteDifference ||
        (result.diagnostics.hessian_valid &&
         is_finite(result.diagnostics.tangent_residual_max_scale_free) &&
         result.diagnostics.tangent_residual_max_scale_free <= tangent_residual_tolerance);
    if (!result.diagnostics.hessian_valid) {
        result.diagnostics.admissibility_reason = "quadratic Hessian estimate is invalid";
    } else if (!raw_residual_is_finite) {
        result.diagnostics.admissibility_reason = "quadratic raw residual is non-finite";
    } else if (!tangent_hessian_gate_passed) {
        result.diagnostics.admissibility_reason =
            "tangent Hessian stencil residual exceeds experimental tolerance";
    } else if (std::abs(residual.residual_scale_free) > tangent_residual_tolerance) {
        result.diagnostics.admissibility_reason =
            "quadratic raw residual exceeds experimental fast-approx tolerance";
    } else {
        result.diagnostics.admissible_for_fast_approximation = true;
        result.diagnostics.admissibility_reason = "admissible";
    }
    return result;
}

Problem1RootTableQueryResult query_problem1_root_table_linear_newton(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_iterations,
    double residual_tolerance,
    double alpha_step_tolerance
) {
    Problem1RootTableQueryResult result{};
    result.method = "linear_nearest_newton";
    const Problem1RootNearestNode nearest = find_nearest_problem1_root_table_node(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!nearest.valid || nearest.cell == nullptr) {
        result.invalid_reason = nearest.invalid_reason.empty()
            ? "failed to locate nearest root-table node"
            : nearest.invalid_reason;
        return result;
    }

    std::vector<Problem1SolutionBranch> branches;
    for (const Problem1SolutionBranch& node_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
        if (!node_branch.valid) {
            continue;
        }
        // 中文注释：第一版 query 层直接从最近节点的每条 branch 出发，用一阶预测给 Newton 提供 seed。
        Problem1SolutionBranch refined = refine_problem1_root_branch_from_linear_prediction(
            table.config().departure_planet,
            table.config().target_planet,
            nearest.node_nu_A,
            nearest.node_nu_B,
            nearest.node_theta_A,
            node_branch,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            max_iterations,
            residual_tolerance,
            alpha_step_tolerance);
        if (!refined.valid) {
            continue;
        }
        add_query_branch_if_not_duplicate(&branches, std::move(refined));
    }

    std::sort(
        branches.begin(),
        branches.end(),
        [](const Problem1SolutionBranch& lhs, const Problem1SolutionBranch& rhs) {
            if (lhs.time_of_flight_seconds < rhs.time_of_flight_seconds) {
                return true;
            }
            if (lhs.time_of_flight_seconds > rhs.time_of_flight_seconds) {
                return false;
            }
            if (lhs.transfer_revolution < rhs.transfer_revolution) {
                return true;
            }
            if (lhs.transfer_revolution > rhs.transfer_revolution) {
                return false;
            }
            if (lhs.target_revolution < rhs.target_revolution) {
                return true;
            }
            if (lhs.target_revolution > rhs.target_revolution) {
                return false;
            }
            return lhs.encounter_global_angle < rhs.encounter_global_angle;
        });

    result.valid = true;
    result.invalid_reason.clear();
    result.branches = std::move(branches);
    return result;
}

Problem1RootTableQueryResult query_problem1_root_table_linear_newton_seconds(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance
) {
    return query_problem1_root_table_linear_newton(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        max_iterations,
        problem1_residual_seconds_to_scale_free(residual_tolerance_seconds),
        alpha_step_tolerance);
}

std::vector<Problem1RootApproximationResult> query_problem1_root_table_quadratic_raw(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method,
    double tangent_residual_tolerance
) {
    std::vector<Problem1RootApproximationResult> results;
    const Problem1RootNearestNode nearest = find_nearest_problem1_root_table_node(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!nearest.valid || nearest.cell == nullptr) {
        return results;
    }

    for (const Problem1SolutionBranch& node_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
        if (!node_branch.valid) {
            continue;
        }
        Problem1RootApproximationResult approximation =
            evaluate_problem1_root_quadratic_approximation_from_node(
                table.config().departure_planet,
                table.config().target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                hessian_step,
                hessian_method,
                tangent_residual_tolerance,
                table.config().max_target_revolution);
        if (!approximation.valid) {
            continue;
        }
        add_route_b_approximation_if_not_duplicate(&results, std::move(approximation));
    }

    std::sort(
        results.begin(),
        results.end(),
        [](const Problem1RootApproximationResult& lhs, const Problem1RootApproximationResult& rhs) {
            if (lhs.transfer_time_seconds < rhs.transfer_time_seconds) {
                return true;
            }
            if (lhs.transfer_time_seconds > rhs.transfer_time_seconds) {
                return false;
            }
            if (lhs.transfer_revolution < rhs.transfer_revolution) {
                return true;
            }
            if (lhs.transfer_revolution > rhs.transfer_revolution) {
                return false;
            }
            return lhs.predicted_encounter_global_angle < rhs.predicted_encounter_global_angle;
        });
    return results;
}

Problem1RouteBSafeQueryResult query_problem1_root_table_route_b_safe(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method,
    double tangent_residual_tolerance
) {
    Problem1RouteBSafeQueryResult result{};
    result.method = "route_b_safe_quadratic_raw";

    if (hessian_method == Problem1RootHessianMethod::NewtonRefinedFiniteDifference) {
        result.fallback_required = true;
        result.reason = "route_b_newton_hessian_method_disabled";
        return result;
    }

    const Problem1RootNearestNode nearest = find_nearest_problem1_root_table_node(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!nearest.valid || nearest.cell == nullptr) {
        result.fallback_required = true;
        result.reason = nearest.invalid_reason.empty()
            ? "route_b_failed_to_locate_nearest_node"
            : nearest.invalid_reason;
        return result;
    }

    const int cell_i = lower_corner_index_from_nearest(
        nearest.i,
        normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A));
    const int cell_j = lower_corner_index_from_nearest(
        nearest.j,
        normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B));
    const int cell_k = lower_corner_index_from_nearest(
        nearest.k,
        normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A));
    result.cell_admissibility = evaluate_problem1_root_cell_admissibility(
        table,
        cell_i,
        cell_j,
        cell_k,
        table.config().max_transfer_revolution,
        table.config().max_target_revolution);
    if (!result.cell_admissibility.admissible) {
        result.fallback_required = true;
        result.reason = "cell_non_admissible_branch_count_inconsistent";
        return result;
    }

    result.approximations = query_problem1_root_table_quadratic_raw(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        hessian_step,
        hessian_method,
        tangent_residual_tolerance);
    result.valid = true;
    result.fallback_required = false;
    result.reason = "route_b_safe_query_ok";
    return result;
}

Problem1RouteBSafeQueryResult query_problem1_root_table_route_b_linear_safe(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    Problem1RouteBSafeQueryResult result{};
    result.method = "route_b_linear_safe";

    const Problem1RootNearestNode nearest = find_nearest_problem1_root_table_node(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!nearest.valid || nearest.cell == nullptr) {
        result.fallback_required = true;
        result.reason = nearest.invalid_reason.empty()
            ? "route_b_failed_to_locate_nearest_node"
            : nearest.invalid_reason;
        return result;
    }

    const int cell_i = lower_corner_index_from_nearest(
        nearest.i,
        normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A));
    const int cell_j = lower_corner_index_from_nearest(
        nearest.j,
        normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B));
    const int cell_k = lower_corner_index_from_nearest(
        nearest.k,
        normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A));
    result.cell_admissibility = evaluate_problem1_root_cell_admissibility(
        table,
        cell_i,
        cell_j,
        cell_k,
        table.config().max_transfer_revolution,
        table.config().max_target_revolution);
    if (!result.cell_admissibility.admissible) {
        result.fallback_required = true;
        result.reason = "cell_non_admissible_branch_count_inconsistent";
        return result;
    }

    std::vector<Problem1RootApproximationResult> approximations;
    for (const Problem1SolutionBranch& node_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
        if (!node_branch.valid) {
            continue;
        }
        Problem1RootApproximationResult approximation =
            evaluate_problem1_root_linear_route_b_approximation_from_node(
                table.config().departure_planet,
                table.config().target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                table.config().max_target_revolution);
        if (!approximation.valid) {
            continue;
        }
        add_route_b_approximation_if_not_duplicate(&approximations, std::move(approximation));
    }

    populate_route_b_branch_count_diagnostics(result.cell_admissibility, approximations, &result);

    std::sort(
        approximations.begin(),
        approximations.end(),
        [](const Problem1RootApproximationResult& lhs, const Problem1RootApproximationResult& rhs) {
            if (lhs.transfer_revolution < rhs.transfer_revolution) {
                return true;
            }
            if (lhs.transfer_revolution > rhs.transfer_revolution) {
                return false;
            }
            if (lhs.transfer_time_seconds < rhs.transfer_time_seconds) {
                return true;
            }
            if (lhs.transfer_time_seconds > rhs.transfer_time_seconds) {
                return false;
            }
            return lhs.predicted_encounter_global_angle < rhs.predicted_encounter_global_angle;
        });

    if (approximations.empty()) {
        result.valid = false;
        result.fallback_required = true;
        result.reason = "route_b_linear_no_valid_approximations";
        return result;
    }

    result.approximations = approximations;
    result.valid = true;
    if (!result.branch_count_complete) {
        result.fallback_required = true;
        result.reason = "route_b_linear_incomplete_branch_count";
        return result;
    }

    result.fallback_required = false;
    result.reason = "route_b_linear_safe_query_ok";
    return result;
}

Problem1RouteBSafeQueryResult query_problem1_root_table_route_b_quadratic_safe(
    const Problem1RootTable& table,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    double hessian_step,
    Problem1RootHessianMethod hessian_method,
    double tangent_residual_tolerance
) {
    Problem1RouteBSafeQueryResult result{};
    result.method = "route_b_quadratic_safe";

    if (hessian_method == Problem1RootHessianMethod::NewtonRefinedFiniteDifference) {
        result.fallback_required = true;
        result.reason = "route_b_quadratic_newton_hessian_method_disabled";
        return result;
    }

    const Problem1RootNearestNode nearest = find_nearest_problem1_root_table_node(
        table,
        query_nu_A,
        query_nu_B,
        query_theta_A);
    if (!nearest.valid || nearest.cell == nullptr) {
        result.fallback_required = true;
        result.reason = nearest.invalid_reason.empty()
            ? "route_b_failed_to_locate_nearest_node"
            : nearest.invalid_reason;
        return result;
    }

    const int cell_i = lower_corner_index_from_nearest(
        nearest.i,
        normalize_angle_minus_pi_pi(query_nu_A - nearest.node_nu_A));
    const int cell_j = lower_corner_index_from_nearest(
        nearest.j,
        normalize_angle_minus_pi_pi(query_nu_B - nearest.node_nu_B));
    const int cell_k = lower_corner_index_from_nearest(
        nearest.k,
        normalize_angle_minus_pi_pi(query_theta_A - nearest.node_theta_A));
    result.cell_admissibility = evaluate_problem1_root_cell_admissibility(
        table,
        cell_i,
        cell_j,
        cell_k,
        table.config().max_transfer_revolution,
        table.config().max_target_revolution);
    if (!result.cell_admissibility.admissible) {
        result.fallback_required = true;
        result.reason = "cell_non_admissible_branch_count_inconsistent";
        return result;
    }

    std::vector<Problem1RootApproximationResult> approximations;
    for (const Problem1SolutionBranch& node_branch : nearest.cell->solutions_sorted_by_time_of_flight) {
        if (!node_branch.valid) {
            continue;
        }
        Problem1RootApproximationResult approximation =
            evaluate_problem1_root_quadratic_route_b_approximation_from_node(
                table.config().departure_planet,
                table.config().target_planet,
                nearest.node_nu_A,
                nearest.node_nu_B,
                nearest.node_theta_A,
                node_branch,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                table.config().max_target_revolution,
                hessian_step,
                hessian_method,
                tangent_residual_tolerance);
        if (!approximation.valid) {
            continue;
        }
        add_route_b_approximation_if_not_duplicate(&approximations, std::move(approximation));
    }

    std::sort(
        approximations.begin(),
        approximations.end(),
        [](const Problem1RootApproximationResult& lhs, const Problem1RootApproximationResult& rhs) {
            if (lhs.transfer_revolution < rhs.transfer_revolution) {
                return true;
            }
            if (lhs.transfer_revolution > rhs.transfer_revolution) {
                return false;
            }
            if (lhs.transfer_time_seconds < rhs.transfer_time_seconds) {
                return true;
            }
            if (lhs.transfer_time_seconds > rhs.transfer_time_seconds) {
                return false;
            }
            return lhs.predicted_encounter_global_angle < rhs.predicted_encounter_global_angle;
        });

    populate_route_b_branch_count_diagnostics(result.cell_admissibility, approximations, &result);
    if (approximations.empty()) {
        result.valid = false;
        result.fallback_required = true;
        result.reason = "route_b_quadratic_no_valid_approximations";
        return result;
    }

    result.approximations = approximations;
    result.valid = true;
    if (!result.branch_count_complete) {
        result.fallback_required = true;
        result.reason = "route_b_quadratic_incomplete_branch_count";
        return result;
    }

    result.fallback_required = false;
    result.reason = "route_b_quadratic_safe_query_ok";
    return result;
}

std::vector<Problem1SolutionBranch> convert_problem1_candidates_to_solution_branches(
    planet_params::PlanetId target_planet,
    const std::vector<Problem1Candidate>& candidates
) {
    std::vector<Problem1SolutionBranch> branches;
    branches.reserve(candidates.size());
    for (const Problem1Candidate& candidate : candidates) {
        branches.push_back(problem1_solution_branch_from_candidate(target_planet, candidate));
    }

    std::sort(
        branches.begin(),
        branches.end(),
        [](const Problem1SolutionBranch& lhs, const Problem1SolutionBranch& rhs) {
            // 中文注释：对 root table 来说，排序应按相对飞行时间；只有 launch_time 固定时，
            // 它才与绝对 arrival_time 排序等价。
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    return branches;
}

std::vector<Problem1SolutionBranch> solve_problem1_from_departure_anomalies(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int max_transfer_revolution,
    int max_target_revolution
) {
    return solve_problem1_from_departure_anomalies_diagnostic(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        max_transfer_revolution,
        max_target_revolution).branches;
}

Problem1SolveWithModeResult solve_problem1_from_departure_anomalies_with_mode(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    const Problem1SolveOptions& options
) {
    using clock = std::chrono::steady_clock;
    const auto start = clock::now();

    Problem1SolveWithModeResult result{};
    result.mode_profile.mode_requested = options.mode;
    result.mode_profile.mode_used = Problem1SolveMode::FullScan2880;
    result.mode_profile.adaptive_attempted = options.mode != Problem1SolveMode::FullScan2880;
    result.mode_profile.fallback_used = false;
    result.mode_profile.debug_compare_enabled = options.mode == Problem1SolveMode::AdaptiveScanDebugCompare;

    const Problem1DepartureAnomalySolveInput input =
        make_departure_anomaly_solve_input_from_defaults(
            departure_planet,
            target_planet,
            nu_A_depart,
            nu_B_depart,
            theta_A,
            options.max_transfer_revolution,
            options.max_target_revolution);

    if (options.mode == Problem1SolveMode::AdaptiveScanDebugCompare) {
        Problem1AdaptiveScanSummary adaptive_summary{};
        const auto adaptive_placeholder =
            solve_problem1_from_departure_anomalies_adaptive_internal(input, options.adaptive, &adaptive_summary);
        const auto full_baseline = solve_problem1_from_departure_anomalies_full_scan_internal(input);
        const Problem1BranchCompareReport compare_report =
            compare_problem1_solution_branches_by_identity(
                adaptive_placeholder.branches,
                full_baseline.branches,
                nu_A_depart,
                nu_B_depart,
                theta_A);
        apply_branch_compare_report_to_profile(compare_report, &result.mode_profile);
        result.mode_profile.debug_compare_missing_branches = compare_report.missing_branches;
        result.mode_profile.debug_compare_extra_branches = compare_report.extra_branches;
        result.mode_profile.debug_compare_extra_diagnostics =
            make_adaptive_extra_branch_diagnostics(input, adaptive_placeholder.branches, full_baseline.branches, compare_report);
        apply_adaptive_scan_summary_to_profile(adaptive_summary, &result.mode_profile);
        result.branches = full_baseline.branches;
        result.diagnostic = full_baseline.diagnostic;
        result.mode_profile.residual_evaluation_count =
            static_cast<int>(std::min<long long>(
                full_baseline.diagnostic.residual_evaluations,
                static_cast<long long>(std::numeric_limits<int>::max())));
        result.mode_profile.branch_count = static_cast<int>(result.branches.size());
        result.mode_profile.total_ms = std::chrono::duration<double, std::milli>(clock::now() - start).count();
        return result;
    }

    Problem1AdaptiveScanSummary adaptive_summary{};
    Problem1SolveWithDiagnosticResult solved{};
    if (options.mode == Problem1SolveMode::FullScan2880) {
        solved = solve_problem1_from_departure_anomalies_full_scan_internal(input);
    } else {
        if (env_flag_enabled("PROBLEM1_ADAPTIVE_FORCE_FULLSCAN")) {
            adaptive_summary.fallback_to_fullscan = true;
            result.mode_profile.fallback_used = true;
            solved = solve_problem1_from_departure_anomalies_full_scan_internal(input);
        } else {
            solved = solve_problem1_from_departure_anomalies_adaptive_internal(
                input, options.adaptive, &adaptive_summary);
            if (!solved.diagnostic.valid) {
                adaptive_summary.fallback_to_fullscan = true;
                result.mode_profile.fallback_used = true;
                solved = solve_problem1_from_departure_anomalies_full_scan_internal(input);
            }
        }
    }
    if (options.mode != Problem1SolveMode::FullScan2880) {
        apply_adaptive_scan_summary_to_profile(adaptive_summary, &result.mode_profile);
    }
    result.branches = solved.branches;
    result.diagnostic = solved.diagnostic;
    result.mode_profile.residual_evaluation_count =
        static_cast<int>(std::min<long long>(
            solved.diagnostic.residual_evaluations,
            static_cast<long long>(std::numeric_limits<int>::max())));
    result.mode_profile.branch_count = static_cast<int>(result.branches.size());
    result.mode_profile.total_ms = std::chrono::duration<double, std::milli>(clock::now() - start).count();
    return result;
}

Problem1SolveWithDiagnosticResult solve_problem1_from_departure_anomalies_diagnostic(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int max_transfer_revolution,
    int max_target_revolution
) {
    const Problem1DepartureAnomalySolveInput input =
        make_departure_anomaly_solve_input_from_defaults(
        departure_planet,
        target_planet,
        nu_A_depart,
        nu_B_depart,
        theta_A,
        max_transfer_revolution,
        max_target_revolution);
    return solve_problem1_from_departure_anomalies_full_scan_internal(input);
}

Problem1RootTable build_problem1_root_table(const Problem1RootTableConfig& config) {
    validate_problem1_root_table_config(config);

    Problem1RootTable table(config);
    for (int nu_A_index = 0; nu_A_index < config.nu_A_count; ++nu_A_index) {
        for (int nu_B_depart_index = 0; nu_B_depart_index < config.nu_B_depart_count; ++nu_B_depart_index) {
            for (int theta_A_index = 0; theta_A_index < config.theta_A_count; ++theta_A_index) {
                Problem1RootTableCell cell{};
                cell.nu_A_depart = normalize_angle_0_2pi(
                    config.nu_A_start + static_cast<double>(nu_A_index) * config.nu_A_step);
                cell.nu_B_depart = normalize_angle_0_2pi(
                    config.nu_B_depart_start + static_cast<double>(nu_B_depart_index) * config.nu_B_depart_step);
                cell.theta_A = normalize_angle_0_2pi(
                    config.theta_A_start + static_cast<double>(theta_A_index) * config.theta_A_step);
                cell.solutions_sorted_by_time_of_flight = solve_problem1_from_departure_anomalies(
                    config.departure_planet,
                    config.target_planet,
                    cell.nu_A_depart,
                    cell.nu_B_depart,
                    cell.theta_A,
                    config.max_transfer_revolution,
                    config.max_target_revolution);
                cell.solved = true;
                cell.invalid_reason.clear();

                table.mutable_cells().at(table.flat_index(nu_A_index, nu_B_depart_index, theta_A_index)) =
                    std::move(cell);
            }
        }
    }

    return table;
}

}  // namespace spaceship_cpp::problem1
