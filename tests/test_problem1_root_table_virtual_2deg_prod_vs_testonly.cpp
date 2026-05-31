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
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kMidCellOffsetDegrees = 1.0;
constexpr double kEngineeringResidualThresholdSeconds = 1e-2;
constexpr double kPhysicalTimeErrorThresholdSeconds = 1.0;
constexpr double kPhysicalAlphaErrorThresholdRadians = 1e-3;

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

struct IndexedBranch {
    int original_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

enum class StopReason {
    ResidualToleranceSuccess,
    ResidualIncreaseStop,
    MaxIterationStop,
    DerivativeInvalidStop,
    ResidualInvalidStop,
    NonFiniteStepStop,
};

struct TestOnlyRefinementResult {
    bool valid = false;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
    StopReason stop_reason = StopReason::ResidualInvalidStop;
    double final_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
};

struct PipelineRecord {
    int source_node_branch_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch source_node_branch;
    bool attach_derivative_ok = false;
    bool raw_prediction_finite = false;
    double raw_predicted_alpha = std::numeric_limits<double>::quiet_NaN();
    bool raw_residual_valid = false;
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();

    bool production_valid = false;
    spaceship_cpp::problem1::Problem1RootRefinementResult production_result;

    bool test_only_valid = false;
    TestOnlyRefinementResult test_only_result;
};

double wrapped_alpha_distance(double alpha1, double alpha2) {
    return std::abs(normalize_angle_minus_pi_pi(alpha1 - alpha2));
}

const char* stop_reason_name(StopReason reason) {
    switch (reason) {
        case StopReason::ResidualToleranceSuccess: return "residual_tolerance_success";
        case StopReason::ResidualIncreaseStop: return "residual_increase_stop";
        case StopReason::MaxIterationStop: return "max_iteration_stop";
        case StopReason::DerivativeInvalidStop: return "derivative_invalid_stop";
        case StopReason::ResidualInvalidStop: return "residual_invalid_stop";
        case StopReason::NonFiniteStepStop: return "nonfinite_step_stop";
    }
    return "unknown";
}

double branch_time_duplicate_tolerance_seconds(double t1, double t2) {
    return std::max(1e-4, 1e-10 * std::max(std::abs(t1), std::abs(t2)));
}

bool same_physical_branch_by_k_time(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return lhs.transfer_revolution == rhs.transfer_revolution &&
        std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds) <=
            branch_time_duplicate_tolerance_seconds(lhs.time_of_flight_seconds, rhs.time_of_flight_seconds);
}

bool safe_to_merge(
    const spaceship_cpp::problem1::Problem1SolutionBranch& lhs,
    const spaceship_cpp::problem1::Problem1SolutionBranch& rhs
) {
    return same_physical_branch_by_k_time(lhs, rhs) &&
        wrapped_alpha_distance(lhs.encounter_global_angle, rhs.encounter_global_angle) <= 1e-4;
}

void add_if_not_duplicate(
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>* branches,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    for (auto& existing : *branches) {
        if (!same_physical_branch_by_k_time(existing, branch)) {
            continue;
        }
        if (!safe_to_merge(existing, branch)) {
            continue;
        }
        if (std::abs(branch.residual_seconds) < std::abs(existing.residual_seconds)) {
            existing = branch;
        }
        return;
    }
    branches->push_back(branch);
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

std::vector<QuerySample> build_samples(int samples_per_group) {
    std::vector<QuerySample> samples;
    const double mid_offset = kMidCellOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < samples_per_group; ++i) {
        const double base2_nu_A = normalize_angle_0_2pi(0.91 + static_cast<double>(i) * 1.8123456789);
        const double base2_nu_B = normalize_angle_0_2pi(1.41 + static_cast<double>(i) * 2.2718281828);
        const double base2_theta_A = normalize_angle_0_2pi(0.63 + static_cast<double>(i) * 1.1447298860);
        const auto mid_node = find_nearest_virtual_root_table_node(base2_nu_A, base2_nu_B, base2_theta_A);
        samples.push_back({"mid_cell",
                           normalize_angle_0_2pi(mid_node.nu_A_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.nu_B_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.theta_A_node + mid_offset)});
    }
    return samples;
}

std::vector<IndexedBranch> collect_sorted_by_time(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::vector<IndexedBranch> result;
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (branches[i].valid) {
            result.push_back({static_cast<int>(i), branches[i]});
        }
    }
    std::sort(result.begin(), result.end(), [](const IndexedBranch& a, const IndexedBranch& b) {
        if (a.branch.time_of_flight_seconds < b.branch.time_of_flight_seconds) {
            return true;
        }
        if (a.branch.time_of_flight_seconds > b.branch.time_of_flight_seconds) {
            return false;
        }
        return a.branch.encounter_global_angle < b.branch.encounter_global_angle;
    });
    return result;
}

std::map<int, std::vector<IndexedBranch>> group_by_k_sorted(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::map<int, std::vector<IndexedBranch>> grouped;
    const auto sorted = collect_sorted_by_time(branches);
    for (const auto& entry : sorted) {
        grouped[entry.branch.transfer_revolution].push_back(entry);
    }
    return grouped;
}

bool is_engineering_match(
    const spaceship_cpp::problem1::Problem1SolutionBranch& candidate,
    const spaceship_cpp::problem1::Problem1SolutionBranch& exact
) {
    return wrapped_alpha_distance(candidate.encounter_global_angle, exact.encounter_global_angle) <=
            kPhysicalAlphaErrorThresholdRadians &&
        std::abs(candidate.time_of_flight_seconds - exact.time_of_flight_seconds) <=
            kPhysicalTimeErrorThresholdSeconds &&
        std::abs(candidate.residual_seconds) <= kEngineeringResidualThresholdSeconds;
}

TestOnlyRefinementResult refine_problem1_root_branch_newton_residual_first_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double nu_A_depart,
    double nu_B_depart,
    double theta_A,
    int transfer_revolution,
    int target_revolution,
    double initial_alpha,
    int max_iterations,
    double residual_tolerance_seconds,
    double alpha_step_tolerance
) {
    namespace problem1 = spaceship_cpp::problem1;
    TestOnlyRefinementResult result{};
    double alpha_unwrapped = initial_alpha;
    double previous_abs_residual_scale_free = std::numeric_limits<double>::infinity();
    (void)alpha_step_tolerance;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        const double alpha_wrapped = normalize_angle_0_2pi(alpha_unwrapped);
        const auto residual = problem1::evaluate_problem1_root_residual(
            departure_planet, target_planet, nu_A_depart, nu_B_depart, theta_A, alpha_wrapped,
            transfer_revolution, target_revolution);
        if (!residual.valid || !std::isfinite(residual.residual_seconds) ||
            !std::isfinite(residual.residual_scale_free)) {
            result.stop_reason = StopReason::ResidualInvalidStop;
            result.final_residual_seconds = residual.residual_seconds;
            result.iterations = iteration + 1;
            return result;
        }
        result.final_residual_seconds = residual.residual_seconds;
        if (std::abs(residual.residual_seconds) <= residual_tolerance_seconds) {
            result.valid = true;
            result.stop_reason = StopReason::ResidualToleranceSuccess;
            result.iterations = iteration + 1;
            result.branch.valid = true;
            result.branch.encounter_global_angle = residual.encounter_global_angle;
            result.branch.target_arrival_true_anomaly = residual.target_arrival_true_anomaly;
            result.branch.transfer_revolution = transfer_revolution;
            result.branch.target_revolution = target_revolution;
            result.branch.time_of_flight_seconds = residual.transfer_time_seconds;
            result.branch.target_time_seconds = residual.target_time_seconds;
            result.branch.residual_seconds = residual.residual_seconds;
            result.branch.transfer_e = residual.transfer_e;
            result.branch.transfer_p = residual.transfer_p;
            result.branch.transfer_a = residual.transfer_a;
            result.branch.theta_B = residual.theta_B;
            return result;
        }
        const auto derivatives = problem1::evaluate_problem1_root_residual_derivatives_with_mode(
            departure_planet, target_planet, nu_A_depart, nu_B_depart, theta_A, alpha_wrapped,
            transfer_revolution, target_revolution,
            problem1::Problem1RootDerivativeMode::AnalyticOnly, 1e-6);
        if (!derivatives.valid || !std::isfinite(derivatives.R_alpha) || std::abs(derivatives.R_alpha) <= 1e-12) {
            result.stop_reason = StopReason::DerivativeInvalidStop;
            result.iterations = iteration + 1;
            return result;
        }
        const double delta_alpha = -residual.residual_scale_free / derivatives.R_alpha;
        if (!std::isfinite(delta_alpha)) {
            result.stop_reason = StopReason::NonFiniteStepStop;
            result.iterations = iteration + 1;
            return result;
        }
        const double trial_alpha_wrapped = normalize_angle_0_2pi(alpha_unwrapped + delta_alpha);
        const auto trial_residual = problem1::evaluate_problem1_root_residual(
            departure_planet, target_planet, nu_A_depart, nu_B_depart, theta_A, trial_alpha_wrapped,
            transfer_revolution, target_revolution);
        if (!trial_residual.valid || !std::isfinite(trial_residual.residual_scale_free)) {
            result.stop_reason = StopReason::ResidualInvalidStop;
            result.final_residual_seconds = trial_residual.residual_seconds;
            result.iterations = iteration + 1;
            return result;
        }
        const double current_abs = std::abs(residual.residual_scale_free);
        const double trial_abs = std::abs(trial_residual.residual_scale_free);
        if (trial_abs > current_abs || trial_abs > previous_abs_residual_scale_free) {
            result.stop_reason = StopReason::ResidualIncreaseStop;
            result.iterations = iteration + 1;
            result.final_residual_seconds = residual.residual_seconds;
            return result;
        }
        previous_abs_residual_scale_free = trial_abs;
        alpha_unwrapped += delta_alpha;
    }
    result.stop_reason = StopReason::MaxIterationStop;
    result.iterations = max_iterations;
    return result;
}

std::vector<PipelineRecord> build_midcell_records(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& node_branches
) {
    namespace problem1 = spaceship_cpp::problem1;
    std::vector<PipelineRecord> records;
    for (std::size_t i = 0; i < node_branches.size(); ++i) {
        const auto& source = node_branches[i];
        if (!source.valid) {
            continue;
        }
        PipelineRecord record{};
        record.source_node_branch_index = static_cast<int>(i);
        record.source_node_branch = source;
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet, target_planet, node.nu_A_node, node.nu_B_node, node.theta_A_node, source,
            problem1::Problem1RootDerivativeMode::AnalyticOnly, 1e-6);
        record.attach_derivative_ok = attached.valid && attached.derivatives_available;
        if (!record.attach_derivative_ok) {
            records.push_back(record);
            continue;
        }
        record.raw_predicted_alpha = normalize_angle_0_2pi(
            attached.encounter_global_angle +
            attached.d_encounter_global_angle_d_nu_A * node.dnu_A +
            attached.d_encounter_global_angle_d_nu_B * node.dnu_B +
            attached.d_encounter_global_angle_d_theta_A * node.dtheta_A);
        record.raw_prediction_finite = std::isfinite(record.raw_predicted_alpha);
        if (!record.raw_prediction_finite) {
            records.push_back(record);
            continue;
        }
        const auto raw_residual = problem1::evaluate_problem1_root_residual(
            departure_planet, target_planet, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            record.raw_predicted_alpha, source.transfer_revolution, source.target_revolution);
        record.raw_residual_valid = raw_residual.valid;
        record.raw_residual_seconds = raw_residual.residual_seconds;

        record.production_result = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            departure_planet, target_planet, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            source.transfer_revolution, source.target_revolution, record.raw_predicted_alpha,
            80, 1e-2, 1e-12, problem1::Problem1RootDerivativeMode::AnalyticOnly, 1e-6);
        record.production_valid = record.production_result.valid;

        record.test_only_result = refine_problem1_root_branch_newton_residual_first_test_only(
            departure_planet, target_planet, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            source.transfer_revolution, source.target_revolution, record.raw_predicted_alpha,
            80, 1e-2, 1e-14);
        record.test_only_valid = record.test_only_result.valid;

        records.push_back(record);
    }
    return records;
}

std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> collect_before_dedup(
    const std::vector<PipelineRecord>& records,
    bool production
) {
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> out;
    for (const auto& r : records) {
        if (production) {
            if (r.production_valid) {
                out.push_back(r.production_result.branch);
            }
        } else {
            if (r.test_only_valid) {
                out.push_back(r.test_only_result.branch);
            }
        }
    }
    return out;
}

std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> dedup_candidates(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> out;
    for (const auto& b : branches) {
        add_if_not_duplicate(&out, b);
    }
    return out;
}

int engineering_coverage_for_k(
    const std::vector<IndexedBranch>& exact_group,
    const std::vector<IndexedBranch>& cand_group
) {
    const int n = std::min(exact_group.size(), cand_group.size());
    int covered = 0;
    for (int i = 0; i < n; ++i) {
        if (is_engineering_match(cand_group[static_cast<std::size_t>(i)].branch,
                                 exact_group[static_cast<std::size_t>(i)].branch)) {
            covered += 1;
        }
    }
    return covered;
}

void print_candidate_table(
    const char* title,
    const std::vector<IndexedBranch>& group,
    const std::vector<PipelineRecord>& records,
    bool production
) {
    std::cout << title << "\n";
    for (const auto& entry : group) {
        int src_index = -1;
        std::string stop_reason = "n/a";
        for (const auto& r : records) {
            const auto* branch = production
                ? (r.production_valid ? &r.production_result.branch : nullptr)
                : (r.test_only_valid ? &r.test_only_result.branch : nullptr);
            if (branch == nullptr) {
                continue;
            }
            if (branch->transfer_revolution == entry.branch.transfer_revolution &&
                std::abs(branch->time_of_flight_seconds - entry.branch.time_of_flight_seconds) <= 1e-8 &&
                wrapped_alpha_distance(branch->encounter_global_angle, entry.branch.encounter_global_angle) <= 1e-10) {
                src_index = r.source_node_branch_index;
                stop_reason = production
                    ? (r.production_result.valid ? "valid"
                       : (r.production_result.diagnostic.invalid_reason.empty() ? "invalid" : r.production_result.diagnostic.invalid_reason))
                    : stop_reason_name(r.test_only_result.stop_reason);
                break;
            }
        }
        std::cout << "  source_node_branch_index=" << src_index
                  << " k=" << entry.branch.transfer_revolution
                  << " q=" << entry.branch.target_revolution
                  << " time=" << entry.branch.time_of_flight_seconds
                  << " alpha=" << entry.branch.encounter_global_angle
                  << " residual=" << entry.branch.residual_seconds
                  << " stop_reason=" << stop_reason << "\n";
    }
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const auto samples = build_samples(12);
    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "midcell_prod_vs_testonly_debug\n";
    std::cout << "sample_count=" << samples.size() << "\n";
    std::cout << "shared_sample_generation=1\n";
    std::cout << "shared_nearest_node_logic=1\n";
    std::cout << "shared_max_transfer_revolution=" << max_transfer_revolution << "\n";
    std::cout << "shared_max_target_revolution=" << max_target_revolution << "\n";
    std::cout << "shared_derivative_mode=AnalyticOnly\n";

    bool printed_first_diff = false;
    int total_prod_coverage = 0;
    int total_test_coverage = 0;
    int total_exact = 0;

    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const auto& sample = samples[sample_index];
        const auto node = find_nearest_virtual_root_table_node(
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);

        auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet, target_planet, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            max_transfer_revolution, max_target_revolution);
        for (auto& branch : exact_branches) {
            if (!branch.valid) {
                continue;
            }
            const auto polished = problem1::refine_problem1_root_branch_newton_seconds(
                departure_planet, target_planet, sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                branch.transfer_revolution, branch.target_revolution, branch.encounter_global_angle,
                30, 1e-6, 1e-14);
            if (polished.valid) {
                branch = polished;
            }
        }

        const auto node_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet, target_planet, node.nu_A_node, node.nu_B_node, node.theta_A_node,
            max_transfer_revolution, max_target_revolution);
        const auto records = build_midcell_records(
            departure_planet, target_planet, sample, node, node_branches);

        const auto prod_before = collect_before_dedup(records, true);
        const auto test_before = collect_before_dedup(records, false);
        const auto prod_after = dedup_candidates(prod_before);
        const auto test_after = dedup_candidates(test_before);

        const auto exact_by_k = group_by_k_sorted(exact_branches);
        const auto prod_before_by_k = group_by_k_sorted(prod_before);
        const auto test_before_by_k = group_by_k_sorted(test_before);
        const auto prod_after_by_k = group_by_k_sorted(prod_after);
        const auto test_after_by_k = group_by_k_sorted(test_after);

        int sample_prod_cov = 0;
        int sample_test_cov = 0;
        int sample_exact = 0;
        for (const auto& [k, exact_group] : exact_by_k) {
            sample_exact += static_cast<int>(exact_group.size());
            const auto prod_it = prod_after_by_k.find(k);
            const auto test_it = test_after_by_k.find(k);
            const std::vector<IndexedBranch> empty;
            const auto& prod_group = prod_it == prod_after_by_k.end() ? empty : prod_it->second;
            const auto& test_group = test_it == test_after_by_k.end() ? empty : test_it->second;
            sample_prod_cov += engineering_coverage_for_k(exact_group, prod_group);
            sample_test_cov += engineering_coverage_for_k(exact_group, test_group);
        }
        total_prod_coverage += sample_prod_cov;
        total_test_coverage += sample_test_cov;
        total_exact += sample_exact;

        bool sample_diff = false;
        for (const auto& [k, exact_group] : exact_by_k) {
            const int prod_before_count = prod_before_by_k.count(k) ? static_cast<int>(prod_before_by_k.at(k).size()) : 0;
            const int test_before_count = test_before_by_k.count(k) ? static_cast<int>(test_before_by_k.at(k).size()) : 0;
            const int prod_after_count = prod_after_by_k.count(k) ? static_cast<int>(prod_after_by_k.at(k).size()) : 0;
            const int test_after_count = test_after_by_k.count(k) ? static_cast<int>(test_after_by_k.at(k).size()) : 0;
            if (prod_before_count != test_before_count || prod_after_count != test_after_count) {
                sample_diff = true;
                break;
            }
        }
        if (!sample_diff && sample_prod_cov == sample_test_cov) {
            continue;
        }

        if (!printed_first_diff) {
            printed_first_diff = true;
            std::cout << "first_midcell_diff_case\n";
            std::cout << "sample_index=" << sample_index << "\n";
            std::cout << "query_nu_A=" << sample.query_nu_A << "\n";
            std::cout << "query_nu_B=" << sample.query_nu_B << "\n";
            std::cout << "query_theta_A=" << sample.query_theta_A << "\n";
            std::cout << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
            std::cout << "node_offset_degrees=("
                      << (node.dnu_A * 180.0 / kPi) << ","
                      << (node.dnu_B * 180.0 / kPi) << ","
                      << (node.dtheta_A * 180.0 / kPi) << ")\n";
            std::cout << "exact_total_count=" << sample_exact << "\n";
            std::cout << "production_engineering_coverage=" << sample_prod_cov << "\n";
            std::cout << "test_only_engineering_coverage=" << sample_test_cov << "\n";

            std::cout << "source_node_branches\n";
            for (const auto& r : records) {
                std::cout << "  source_node_branch_index=" << r.source_node_branch_index
                          << " source_k=" << r.source_node_branch.transfer_revolution
                          << " source_q=" << r.source_node_branch.target_revolution
                          << " source_time=" << r.source_node_branch.time_of_flight_seconds
                          << " source_alpha=" << r.source_node_branch.encounter_global_angle
                          << " raw_predicted_alpha=" << r.raw_predicted_alpha
                          << " raw_residual_seconds=" << r.raw_residual_seconds
                          << " production_valid=" << r.production_valid
                          << " production_final_k=" << (r.production_valid ? r.production_result.branch.transfer_revolution : -1)
                          << " production_final_q=" << (r.production_valid ? r.production_result.branch.target_revolution : -1)
                          << " production_final_time=" << (r.production_valid ? r.production_result.branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN())
                          << " production_final_alpha=" << (r.production_valid ? r.production_result.branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN())
                          << " production_final_residual=" << (r.production_valid ? r.production_result.branch.residual_seconds : std::numeric_limits<double>::quiet_NaN())
                          << " production_stop_reason=" << (r.production_valid ? "residual_tolerance_success" : r.production_result.diagnostic.invalid_reason)
                          << " production_derivative_attach_failed_after_convergence=" << r.production_result.diagnostic.derivative_attach_failed_after_convergence
                          << " test_only_valid=" << r.test_only_valid
                          << " test_only_final_k=" << (r.test_only_valid ? r.test_only_result.branch.transfer_revolution : -1)
                          << " test_only_final_q=" << (r.test_only_valid ? r.test_only_result.branch.target_revolution : -1)
                          << " test_only_final_time=" << (r.test_only_valid ? r.test_only_result.branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN())
                          << " test_only_final_alpha=" << (r.test_only_valid ? r.test_only_result.branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN())
                          << " test_only_final_residual=" << r.test_only_result.final_residual_seconds
                          << " test_only_stop_reason=" << stop_reason_name(r.test_only_result.stop_reason)
                          << "\n";
            }

            for (const auto& [k, exact_group] : exact_by_k) {
                const int prod_before_count = prod_before_by_k.count(k) ? static_cast<int>(prod_before_by_k.at(k).size()) : 0;
                const int test_before_count = test_before_by_k.count(k) ? static_cast<int>(test_before_by_k.at(k).size()) : 0;
                const int prod_after_count = prod_after_by_k.count(k) ? static_cast<int>(prod_after_by_k.at(k).size()) : 0;
                const int test_after_count = test_after_by_k.count(k) ? static_cast<int>(test_after_by_k.at(k).size()) : 0;
                const int prod_cov = engineering_coverage_for_k(exact_group, prod_after_by_k.count(k) ? prod_after_by_k.at(k) : std::vector<IndexedBranch>{});
                const int test_cov = engineering_coverage_for_k(exact_group, test_after_by_k.count(k) ? test_after_by_k.at(k) : std::vector<IndexedBranch>{});
                std::cout << "k_summary k=" << k
                          << " exact_count_by_k=" << exact_group.size()
                          << " production_before_dedup_count_by_k=" << prod_before_count
                          << " test_only_before_dedup_count_by_k=" << test_before_count
                          << " production_after_dedup_count_by_k=" << prod_after_count
                          << " test_only_after_dedup_count_by_k=" << test_after_count
                          << " production_engineering_coverage_by_k=" << prod_cov
                          << " test_only_engineering_coverage_by_k=" << test_cov
                          << "\n";
                if (prod_before_count != test_before_count || prod_after_count != test_after_count || prod_cov != test_cov) {
                    print_candidate_table("production_candidates_before_dedup", prod_before_by_k.count(k) ? prod_before_by_k.at(k) : std::vector<IndexedBranch>{}, records, true);
                    print_candidate_table("test_only_candidates_before_dedup", test_before_by_k.count(k) ? test_before_by_k.at(k) : std::vector<IndexedBranch>{}, records, false);
                    print_candidate_table("production_candidates_after_dedup", prod_after_by_k.count(k) ? prod_after_by_k.at(k) : std::vector<IndexedBranch>{}, records, true);
                    print_candidate_table("test_only_candidates_after_dedup", test_after_by_k.count(k) ? test_after_by_k.at(k) : std::vector<IndexedBranch>{}, records, false);
                    for (const auto& exact_entry : exact_group) {
                        bool prod_match = false;
                        if (prod_after_by_k.count(k)) {
                            for (const auto& cand : prod_after_by_k.at(k)) {
                                if (is_engineering_match(cand.branch, exact_entry.branch)) {
                                    prod_match = true;
                                    break;
                                }
                            }
                        }
                        bool test_match = false;
                        if (test_after_by_k.count(k)) {
                            for (const auto& cand : test_after_by_k.at(k)) {
                                if (is_engineering_match(cand.branch, exact_entry.branch)) {
                                    test_match = true;
                                    break;
                                }
                            }
                        }
                        if (test_match && !prod_match) {
                            double best_prod_time = std::numeric_limits<double>::infinity();
                            double best_prod_alpha = std::numeric_limits<double>::infinity();
                            double best_prod_residual = std::numeric_limits<double>::quiet_NaN();
                            if (prod_after_by_k.count(k)) {
                                for (const auto& cand : prod_after_by_k.at(k)) {
                                    const double terr = std::abs(cand.branch.time_of_flight_seconds - exact_entry.branch.time_of_flight_seconds);
                                    const double aerr = wrapped_alpha_distance(cand.branch.encounter_global_angle, exact_entry.branch.encounter_global_angle);
                                    if (terr < best_prod_time) {
                                        best_prod_time = terr;
                                        best_prod_alpha = aerr;
                                        best_prod_residual = cand.branch.residual_seconds;
                                    }
                                }
                            }
                            double best_test_time = std::numeric_limits<double>::infinity();
                            double best_test_alpha = std::numeric_limits<double>::infinity();
                            double best_test_residual = std::numeric_limits<double>::quiet_NaN();
                            if (test_after_by_k.count(k)) {
                                for (const auto& cand : test_after_by_k.at(k)) {
                                    const double terr = std::abs(cand.branch.time_of_flight_seconds - exact_entry.branch.time_of_flight_seconds);
                                    const double aerr = wrapped_alpha_distance(cand.branch.encounter_global_angle, exact_entry.branch.encounter_global_angle);
                                    if (terr < best_test_time) {
                                        best_test_time = terr;
                                        best_test_alpha = aerr;
                                        best_test_residual = cand.branch.residual_seconds;
                                    }
                                }
                            }
                            std::cout << "missing_in_production_but_covered_in_test_only"
                                      << " k=" << k
                                      << " exact_time=" << exact_entry.branch.time_of_flight_seconds
                                      << " exact_alpha=" << exact_entry.branch.encounter_global_angle
                                      << " nearest_production_time_error=" << best_prod_time
                                      << " nearest_production_alpha_error=" << best_prod_alpha
                                      << " nearest_production_residual=" << best_prod_residual
                                      << " nearest_test_time_error=" << best_test_time
                                      << " nearest_test_alpha_error=" << best_test_alpha
                                      << " nearest_test_residual=" << best_test_residual
                                      << "\n";
                        }
                    }
                }
            }
        }
    }

    std::cout << "midcell_total_exact=" << total_exact << "\n";
    std::cout << "midcell_total_production_engineering_coverage=" << total_prod_coverage << "\n";
    std::cout << "midcell_total_testonly_engineering_coverage=" << total_test_coverage << "\n";

    assert(!samples.empty());
    return 0;
}
