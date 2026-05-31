#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <array>
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
constexpr double kNearNodeOffsetDegrees = 0.1;
constexpr int kFixedSampleIndex = 0;
constexpr int kFixedSourceBranchIndex = 0;
constexpr double kSecondsTolerance1e2 = 1e-2;
constexpr std::array<double, 4> kHessianSteps{{1e-5, 5e-6, 2e-6, 1e-6}};

struct QuerySample {
    std::string group_name;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

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

struct FixedCase {
    QuerySample sample;
    spaceship_cpp::problem1::Problem1RootTable table;
    spaceship_cpp::problem1::Problem1RootNearestNode nearest;
    spaceship_cpp::problem1::Problem1SolutionBranch source_branch;
    spaceship_cpp::problem1::Problem1SolutionBranch differentiated_source_branch;
    spaceship_cpp::problem1::Problem1SolutionBranch exact_branch;
    int source_branch_index = -1;

    FixedCase(
        QuerySample sample_in,
        spaceship_cpp::problem1::Problem1RootTable table_in,
        spaceship_cpp::problem1::Problem1RootNearestNode nearest_in,
        spaceship_cpp::problem1::Problem1SolutionBranch source_branch_in,
        spaceship_cpp::problem1::Problem1SolutionBranch differentiated_source_branch_in,
        spaceship_cpp::problem1::Problem1SolutionBranch exact_branch_in,
        int source_branch_index_in
    )
        : sample(std::move(sample_in)),
          table(std::move(table_in)),
          nearest(nearest_in),
          source_branch(std::move(source_branch_in)),
          differentiated_source_branch(std::move(differentiated_source_branch_in)),
          exact_branch(std::move(exact_branch_in)),
          source_branch_index(source_branch_index_in) {}
};

struct StencilAttachedBranch {
    bool valid = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct HessianBuildResult {
    bool hessian_valid = false;
    std::string hessian_invalid_reason;
    double hessian_step = std::numeric_limits<double>::quiet_NaN();
    double H_aa = std::numeric_limits<double>::quiet_NaN();
    double H_bb = std::numeric_limits<double>::quiet_NaN();
    double H_cc = std::numeric_limits<double>::quiet_NaN();
    double H_ab = std::numeric_limits<double>::quiet_NaN();
    double H_ac = std::numeric_limits<double>::quiet_NaN();
    double H_bc = std::numeric_limits<double>::quiet_NaN();
};

struct CorrectionStepResult {
    bool residual_valid = false;
    double alpha = std::numeric_limits<double>::quiet_NaN();
    double residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double transfer_time_seconds = std::numeric_limits<double>::quiet_NaN();
    double target_time_seconds = std::numeric_limits<double>::quiet_NaN();
    double alpha_error = std::numeric_limits<double>::quiet_NaN();
    double time_error = std::numeric_limits<double>::quiet_NaN();
    bool derivative_valid = false;
    double dF_dalpha = std::numeric_limits<double>::quiet_NaN();
    double delta_alpha_correction = std::numeric_limits<double>::quiet_NaN();
    double alpha_corrected = std::numeric_limits<double>::quiet_NaN();
};

struct CorrectionExperiment {
    double hessian_step = std::numeric_limits<double>::quiet_NaN();
    double alpha_linear = std::numeric_limits<double>::quiet_NaN();
    double alpha_quadratic = std::numeric_limits<double>::quiet_NaN();
    HessianBuildResult hessian;
    CorrectionStepResult zero;
    CorrectionStepResult one;
    CorrectionStepResult two;
    bool correction_success = false;
};

double wrapped_alpha_distance(double lhs, double rhs) {
    return std::abs(normalize_angle_minus_pi_pi(lhs - rhs));
}

Problem1RootVirtualGridNode find_nearest_virtual_root_table_node(
    double query_nu_A,
    double query_nu_B,
    double query_theta_A
) {
    const int axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));
    const auto nearest_index = [&](double angle) {
        long long index = std::llround(normalize_angle_0_2pi(angle) / kVirtualGridStepRadians);
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
    node.nu_A_node = normalize_angle_0_2pi(static_cast<double>(node.nu_A_index) * kVirtualGridStepRadians);
    node.nu_B_node = normalize_angle_0_2pi(static_cast<double>(node.nu_B_index) * kVirtualGridStepRadians);
    node.theta_A_node = normalize_angle_0_2pi(static_cast<double>(node.theta_A_index) * kVirtualGridStepRadians);
    node.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

int lower_corner_index(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
}

std::vector<QuerySample> build_samples(int samples_per_group) {
    std::vector<QuerySample> samples;
    const double near_offset = kNearNodeOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < samples_per_group; ++i) {
        const double base_nu_A = normalize_angle_0_2pi(0.37 + static_cast<double>(i) * 2.3999632297);
        const double base_nu_B = normalize_angle_0_2pi(1.11 + static_cast<double>(i) * 1.7548776662);
        const double base_theta_A = normalize_angle_0_2pi(0.23 + static_cast<double>(i) * 0.9182736455);
        const auto near_node = find_nearest_virtual_root_table_node(base_nu_A, base_nu_B, base_theta_A);
        samples.push_back({"near_node",
                           normalize_angle_0_2pi(near_node.nu_A_node + near_offset),
                           normalize_angle_0_2pi(near_node.nu_B_node + near_offset),
                           normalize_angle_0_2pi(near_node.theta_A_node + near_offset)});
    }
    return samples;
}

spaceship_cpp::problem1::Problem1RootTable build_local_cell_table_for_query(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    int max_transfer_revolution,
    int max_target_revolution
) {
    const auto nearest = find_nearest_virtual_root_table_node(
        sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
    const int lower_i = lower_corner_index(nearest.nu_A_index, nearest.dnu_A);
    const int lower_j = lower_corner_index(nearest.nu_B_index, nearest.dnu_B);
    const int lower_k = lower_corner_index(nearest.theta_A_index, nearest.dtheta_A);

    spaceship_cpp::problem1::Problem1RootTableConfig config{};
    config.departure_planet = departure_planet;
    config.target_planet = target_planet;
    config.nu_A_start = static_cast<double>(lower_i) * kVirtualGridStepRadians;
    config.nu_A_step = kVirtualGridStepRadians;
    config.nu_A_count = 2;
    config.nu_B_depart_start = static_cast<double>(lower_j) * kVirtualGridStepRadians;
    config.nu_B_depart_step = kVirtualGridStepRadians;
    config.nu_B_depart_count = 2;
    config.theta_A_start = static_cast<double>(lower_k) * kVirtualGridStepRadians;
    config.theta_A_step = kVirtualGridStepRadians;
    config.theta_A_count = 2;
    config.max_transfer_revolution = max_transfer_revolution;
    config.max_target_revolution = max_target_revolution;
    return spaceship_cpp::problem1::build_problem1_root_table(config);
}

FixedCase build_fixed_case() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const auto samples = build_samples(12);
    const QuerySample sample = samples[static_cast<std::size_t>(kFixedSampleIndex)];
    auto table = build_local_cell_table_for_query(departure_planet, target_planet, sample, 1, 1);
    const auto nearest = problem1::find_nearest_problem1_root_table_node(
        table, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
    const auto source_branch =
        nearest.cell->solutions_sorted_by_time_of_flight[static_cast<std::size_t>(kFixedSourceBranchIndex)];
    const auto differentiated_source_branch = problem1::attach_problem1_root_derivatives(
        departure_planet,
        target_planet,
        nearest.node_nu_A,
        nearest.node_nu_B,
        nearest.node_theta_A,
        source_branch);
    auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        sample.query_nu_A,
        sample.query_nu_B,
        sample.query_theta_A,
        1,
        1);
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> exact_same_k;
    for (const auto& branch : exact_branches) {
        if (branch.valid && branch.transfer_revolution == source_branch.transfer_revolution) {
            exact_same_k.push_back(branch);
        }
    }
    std::sort(exact_same_k.begin(), exact_same_k.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
    });
    const auto exact_branch = exact_same_k[static_cast<std::size_t>(kFixedSourceBranchIndex)];
    return FixedCase(
        sample,
        std::move(table),
        nearest,
        source_branch,
        differentiated_source_branch,
        exact_branch,
        kFixedSourceBranchIndex);
}

StencilAttachedBranch evaluate_stencil_branch_under_seconds_gate(
    const FixedCase& fixed_case,
    double stencil_nu_A,
    double stencil_nu_B,
    double stencil_theta_A
) {
    namespace problem1 = spaceship_cpp::problem1;
    const auto& source = fixed_case.differentiated_source_branch;
    const double dx1 = normalize_angle_minus_pi_pi(stencil_nu_A - fixed_case.nearest.node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(stencil_nu_B - fixed_case.nearest.node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(stencil_theta_A - fixed_case.nearest.node_theta_A);
    const double alpha_tangent = normalize_angle_0_2pi(
        source.encounter_global_angle +
        source.d_encounter_global_angle_d_nu_A * dx1 +
        source.d_encounter_global_angle_d_nu_B * dx2 +
        source.d_encounter_global_angle_d_theta_A * dx3);
    const auto residual = problem1::evaluate_problem1_root_residual(
        fixed_case.table.config().departure_planet,
        fixed_case.table.config().target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        alpha_tangent,
        source.transfer_revolution,
        source.target_revolution);
    StencilAttachedBranch result{};
    if (!(residual.valid &&
          std::isfinite(residual.residual_seconds) &&
          std::abs(residual.residual_seconds) <= kSecondsTolerance1e2)) {
        return result;
    }
    spaceship_cpp::problem1::Problem1SolutionBranch branch{};
    branch.valid = true;
    branch.encounter_global_angle = residual.encounter_global_angle;
    branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
    branch.transfer_revolution = source.transfer_revolution;
    branch.target_revolution = source.target_revolution;
    branch.time_of_flight_seconds = residual.transfer_time_seconds;
    branch.target_time_seconds = residual.target_time_seconds;
    branch.residual_seconds = residual.residual_seconds;
    branch.transfer_e = residual.transfer_e;
    branch.transfer_p = residual.transfer_p;
    branch.transfer_a = residual.transfer_a;
    branch.theta_B = residual.theta_B;
    const auto attached = problem1::attach_problem1_root_derivatives(
        fixed_case.table.config().departure_planet,
        fixed_case.table.config().target_planet,
        stencil_nu_A,
        stencil_nu_B,
        stencil_theta_A,
        branch);
    result.valid = attached.valid && attached.derivatives_available;
    result.branch = attached;
    return result;
}

HessianBuildResult build_hessian_under_seconds_gate(
    const FixedCase& fixed_case,
    double hessian_step
) {
    HessianBuildResult result{};
    result.hessian_step = hessian_step;
    const auto& node = fixed_case.nearest;
    const auto nu_A_plus = evaluate_stencil_branch_under_seconds_gate(
        fixed_case, normalize_angle_0_2pi(node.node_nu_A + hessian_step), node.node_nu_B, node.node_theta_A);
    const auto nu_A_minus = evaluate_stencil_branch_under_seconds_gate(
        fixed_case, normalize_angle_0_2pi(node.node_nu_A - hessian_step), node.node_nu_B, node.node_theta_A);
    const auto nu_B_plus = evaluate_stencil_branch_under_seconds_gate(
        fixed_case, node.node_nu_A, normalize_angle_0_2pi(node.node_nu_B + hessian_step), node.node_theta_A);
    const auto nu_B_minus = evaluate_stencil_branch_under_seconds_gate(
        fixed_case, node.node_nu_A, normalize_angle_0_2pi(node.node_nu_B - hessian_step), node.node_theta_A);
    const auto theta_A_plus = evaluate_stencil_branch_under_seconds_gate(
        fixed_case, node.node_nu_A, node.node_nu_B, normalize_angle_0_2pi(node.node_theta_A + hessian_step));
    const auto theta_A_minus = evaluate_stencil_branch_under_seconds_gate(
        fixed_case, node.node_nu_A, node.node_nu_B, normalize_angle_0_2pi(node.node_theta_A - hessian_step));

    if (!(nu_A_plus.valid && nu_A_minus.valid && nu_B_plus.valid && nu_B_minus.valid && theta_A_plus.valid &&
          theta_A_minus.valid)) {
        result.hessian_invalid_reason = "seconds-gate stencil attach failed";
        return result;
    }

    const double inv_2h = 1.0 / (2.0 * hessian_step);
    result.H_aa =
        (nu_A_plus.branch.d_encounter_global_angle_d_nu_A - nu_A_minus.branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    result.H_bb =
        (nu_B_plus.branch.d_encounter_global_angle_d_nu_B - nu_B_minus.branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    result.H_cc =
        (theta_A_plus.branch.d_encounter_global_angle_d_theta_A -
         theta_A_minus.branch.d_encounter_global_angle_d_theta_A) * inv_2h;

    const double H12_from_g1 =
        (nu_B_plus.branch.d_encounter_global_angle_d_nu_A - nu_B_minus.branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H21_from_g2 =
        (nu_A_plus.branch.d_encounter_global_angle_d_nu_B - nu_A_minus.branch.d_encounter_global_angle_d_nu_B) * inv_2h;
    const double H13_from_g1 =
        (theta_A_plus.branch.d_encounter_global_angle_d_nu_A - theta_A_minus.branch.d_encounter_global_angle_d_nu_A) * inv_2h;
    const double H31_from_g3 =
        (nu_A_plus.branch.d_encounter_global_angle_d_theta_A - nu_A_minus.branch.d_encounter_global_angle_d_theta_A) *
        inv_2h;
    const double H23_from_g2 =
        (theta_A_plus.branch.d_encounter_global_angle_d_nu_B - theta_A_minus.branch.d_encounter_global_angle_d_nu_B) *
        inv_2h;
    const double H32_from_g3 =
        (nu_B_plus.branch.d_encounter_global_angle_d_theta_A - nu_B_minus.branch.d_encounter_global_angle_d_theta_A) *
        inv_2h;

    result.H_ab = 0.5 * (H12_from_g1 + H21_from_g2);
    result.H_ac = 0.5 * (H13_from_g1 + H31_from_g3);
    result.H_bc = 0.5 * (H23_from_g2 + H32_from_g3);
    result.hessian_valid =
        std::isfinite(result.H_aa) &&
        std::isfinite(result.H_bb) &&
        std::isfinite(result.H_cc) &&
        std::isfinite(result.H_ab) &&
        std::isfinite(result.H_ac) &&
        std::isfinite(result.H_bc);
    if (!result.hessian_valid) {
        result.hessian_invalid_reason = "non-finite Hessian entry";
    }
    return result;
}

CorrectionStepResult evaluate_alpha_state(
    const FixedCase& fixed_case,
    double alpha
) {
    namespace problem1 = spaceship_cpp::problem1;
    CorrectionStepResult result{};
    result.alpha = normalize_angle_0_2pi(alpha);
    const auto residual = problem1::evaluate_problem1_root_residual(
        fixed_case.table.config().departure_planet,
        fixed_case.table.config().target_planet,
        fixed_case.sample.query_nu_A,
        fixed_case.sample.query_nu_B,
        fixed_case.sample.query_theta_A,
        result.alpha,
        fixed_case.source_branch.transfer_revolution,
        fixed_case.source_branch.target_revolution);
    result.residual_valid = residual.valid;
    result.residual_seconds = residual.residual_seconds;
    result.transfer_time_seconds = residual.transfer_time_seconds;
    result.target_time_seconds = residual.target_time_seconds;
    result.alpha_error = wrapped_alpha_distance(result.alpha, fixed_case.exact_branch.encounter_global_angle);
    result.time_error = std::abs(result.transfer_time_seconds - fixed_case.exact_branch.time_of_flight_seconds);

    const auto derivatives = problem1::evaluate_problem1_root_residual_derivatives(
        fixed_case.table.config().departure_planet,
        fixed_case.table.config().target_planet,
        fixed_case.sample.query_nu_A,
        fixed_case.sample.query_nu_B,
        fixed_case.sample.query_theta_A,
        result.alpha,
        fixed_case.source_branch.transfer_revolution,
        fixed_case.source_branch.target_revolution);
    result.derivative_valid = derivatives.valid && std::isfinite(derivatives.R_alpha);
    result.dF_dalpha = derivatives.R_alpha;
    if (result.residual_valid && result.derivative_valid && std::abs(derivatives.R_alpha) > 1e-12) {
        const double residual_scale_free =
            spaceship_cpp::problem1::problem1_residual_seconds_to_scale_free(result.residual_seconds);
        result.delta_alpha_correction = -residual_scale_free / derivatives.R_alpha;
        result.alpha_corrected = normalize_angle_0_2pi(result.alpha + result.delta_alpha_correction);
    }
    return result;
}

CorrectionExperiment build_correction_experiment(
    const FixedCase& fixed_case,
    double hessian_step
) {
    CorrectionExperiment result{};
    result.hessian_step = hessian_step;
    result.hessian = build_hessian_under_seconds_gate(fixed_case, hessian_step);

    const auto& source = fixed_case.differentiated_source_branch;
    const double dx1 = normalize_angle_minus_pi_pi(fixed_case.sample.query_nu_A - fixed_case.nearest.node_nu_A);
    const double dx2 = normalize_angle_minus_pi_pi(fixed_case.sample.query_nu_B - fixed_case.nearest.node_nu_B);
    const double dx3 = normalize_angle_minus_pi_pi(fixed_case.sample.query_theta_A - fixed_case.nearest.node_theta_A);
    result.alpha_linear = normalize_angle_0_2pi(
        source.encounter_global_angle +
        source.d_encounter_global_angle_d_nu_A * dx1 +
        source.d_encounter_global_angle_d_nu_B * dx2 +
        source.d_encounter_global_angle_d_theta_A * dx3);
    if (!result.hessian.hessian_valid) {
        return result;
    }
    const double quadratic_correction =
        0.5 * (
            result.hessian.H_aa * dx1 * dx1 +
            result.hessian.H_bb * dx2 * dx2 +
            result.hessian.H_cc * dx3 * dx3 +
            2.0 * result.hessian.H_ab * dx1 * dx2 +
            2.0 * result.hessian.H_ac * dx1 * dx3 +
            2.0 * result.hessian.H_bc * dx2 * dx3);
    result.alpha_quadratic = normalize_angle_0_2pi(result.alpha_linear + quadratic_correction);

    result.zero = evaluate_alpha_state(fixed_case, result.alpha_quadratic);
    if (std::isfinite(result.zero.alpha_corrected)) {
        result.one = evaluate_alpha_state(fixed_case, result.zero.alpha_corrected);
    }
    if (std::isfinite(result.one.alpha_corrected)) {
        result.two = evaluate_alpha_state(fixed_case, result.one.alpha_corrected);
    }
    result.correction_success =
        result.one.residual_valid && std::isfinite(result.one.residual_seconds) &&
        std::abs(result.one.residual_seconds) <= 1e-2;
    return result;
}

void print_fixed_case(const FixedCase& fixed_case) {
    std::cout << "FixedInvalidCase\n";
    std::cout << "group=" << fixed_case.sample.group_name << '\n';
    std::cout << "sample_index=" << kFixedSampleIndex << '\n';
    std::cout << "source_branch_index=" << fixed_case.source_branch_index << '\n';
    std::cout << "k=" << fixed_case.source_branch.transfer_revolution << '\n';
    std::cout << "source_q=" << fixed_case.source_branch.target_revolution << '\n';
    std::cout << "node_nu_A=" << fixed_case.nearest.node_nu_A << '\n';
    std::cout << "node_nu_B=" << fixed_case.nearest.node_nu_B << '\n';
    std::cout << "node_theta_A=" << fixed_case.nearest.node_theta_A << '\n';
    std::cout << "query_nu_A=" << fixed_case.sample.query_nu_A << '\n';
    std::cout << "query_nu_B=" << fixed_case.sample.query_nu_B << '\n';
    std::cout << "query_theta_A=" << fixed_case.sample.query_theta_A << '\n';
    std::cout << "exact_alpha=" << fixed_case.exact_branch.encounter_global_angle << '\n';
    std::cout << "exact_time_of_flight=" << fixed_case.exact_branch.time_of_flight_seconds << '\n';
}

void print_correction_step(const std::string& label, const CorrectionStepResult& step) {
    std::cout << label << "_alpha=" << step.alpha << '\n';
    std::cout << label << "_residual_seconds=" << step.residual_seconds << '\n';
    std::cout << label << "_dF_dalpha=" << step.dF_dalpha << '\n';
    std::cout << label << "_delta_alpha_correction=" << step.delta_alpha_correction << '\n';
    std::cout << label << "_alpha_corrected=" << step.alpha_corrected << '\n';
    std::cout << label << "_alpha_error=" << step.alpha_error << '\n';
    std::cout << label << "_time_error=" << step.time_error << '\n';
}

void print_correction_experiment(const CorrectionExperiment& result) {
    std::cout << "CorrectionExperiment\n";
    std::cout << "hessian_step=" << result.hessian_step << '\n';
    std::cout << "alpha_quadratic=" << result.alpha_quadratic << '\n';
    print_correction_step("zero_correction", result.zero);
    print_correction_step("one_correction", result.one);
    print_correction_step("two_correction", result.two);
    std::cout << "alpha_error_before=" << result.zero.alpha_error << '\n';
    std::cout << "alpha_error_after=" << result.one.alpha_error << '\n';
    std::cout << "time_error_before=" << result.zero.time_error << '\n';
    std::cout << "time_error_after=" << result.one.time_error << '\n';
    std::cout << "residual_seconds_before=" << result.zero.residual_seconds << '\n';
    std::cout << "residual_seconds_after=" << result.one.residual_seconds << '\n';
    std::cout << "correction_success=" << (result.correction_success ? "true" : "false") << '\n';
}

}  // namespace

int main() {
    std::cout << std::setprecision(15) << std::scientific;

    const FixedCase fixed_case = build_fixed_case();
    print_fixed_case(fixed_case);
    for (double hessian_step : kHessianSteps) {
        print_correction_experiment(build_correction_experiment(fixed_case, hessian_step));
    }
    return 0;
}
