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
#include <set>
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

struct PipelineRecord {
    int source_node_branch_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch source_node_branch;
    bool attach_derivative_ok = false;
    bool raw_prediction_finite = false;
    double raw_predicted_alpha = std::numeric_limits<double>::quiet_NaN();
    bool raw_residual_valid = false;
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1RootRefinementResult refinement_result;
};

double wrapped_alpha_distance(double a, double b) {
    return std::abs(normalize_angle_minus_pi_pi(a - b));
}

struct QSheetSelectionResult {
    int selected_q = -1;
    bool selection_failed = false;
    double selected_continuity_error = std::numeric_limits<double>::quiet_NaN();
    double source_q_continuity_error = std::numeric_limits<double>::quiet_NaN();
};

QSheetSelectionResult select_q_by_target_time_sheet_continuity_test_only(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int transfer_revolution,
    double alpha_linear,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source_branch,
    int max_target_revolution
) {
    QSheetSelectionResult result{};
    const double source_time_reference =
        std::isfinite(source_branch.target_time_seconds)
            ? source_branch.target_time_seconds
            : source_branch.time_of_flight_seconds;
    double best_error = std::numeric_limits<double>::infinity();
    int best_q = source_branch.target_revolution;
    bool any_valid = false;
    for (int q = 0; q <= max_target_revolution; ++q) {
        const auto residual = spaceship_cpp::problem1::evaluate_problem1_root_residual(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            alpha_linear,
            transfer_revolution,
            q);
        if (!residual.valid) {
            continue;
        }
        const double continuity_error = std::abs(residual.target_time_seconds - source_time_reference);
        if (q == source_branch.target_revolution) {
            result.source_q_continuity_error = continuity_error;
        }
        if (!any_valid || continuity_error < best_error) {
            any_valid = true;
            best_error = continuity_error;
            best_q = q;
        }
    }
    if (!any_valid) {
        result.selection_failed = true;
        result.selected_q = source_branch.target_revolution;
        result.selected_continuity_error = std::numeric_limits<double>::infinity();
        return result;
    }
    result.selected_q = best_q;
    result.selected_continuity_error = best_error;
    return result;
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

std::vector<QuerySample> build_midcell_samples(int samples_per_group) {
    std::vector<QuerySample> samples;
    const double mid_offset = kMidCellOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < samples_per_group; ++i) {
        const double base2_nu_A = normalize_angle_0_2pi(0.91 + static_cast<double>(i) * 1.8123456789);
        const double base2_nu_B = normalize_angle_0_2pi(1.41 + static_cast<double>(i) * 2.2718281828);
        const double base2_theta_A = normalize_angle_0_2pi(0.63 + static_cast<double>(i) * 1.1447298860);
        const auto mid_node = find_nearest_virtual_root_table_node(base2_nu_A, base2_nu_B, base2_theta_A);
        samples.push_back({"mid_cell",
                           normalize_angle_0_2pi(mid_node.nu_A_node + kMidCellOffsetDegrees * kPi / 180.0),
                           normalize_angle_0_2pi(mid_node.nu_B_node + kMidCellOffsetDegrees * kPi / 180.0),
                           normalize_angle_0_2pi(mid_node.theta_A_node + kMidCellOffsetDegrees * kPi / 180.0)});
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

std::vector<PipelineRecord> build_records_for_node(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    const Problem1RootVirtualGridNode& node,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& node_branches,
    int max_target_revolution = 1
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
            departure_planet,
            target_planet,
            node.nu_A_node,
            node.nu_B_node,
            node.theta_A_node,
            source,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        record.attach_derivative_ok = attached.valid && attached.derivatives_available;
        if (!record.attach_derivative_ok) {
            records.push_back(record);
            continue;
        }
        const double delta_alpha =
            attached.d_encounter_global_angle_d_nu_A * node.dnu_A +
            attached.d_encounter_global_angle_d_nu_B * node.dnu_B +
            attached.d_encounter_global_angle_d_theta_A * node.dtheta_A;
        record.raw_predicted_alpha = normalize_angle_0_2pi(attached.encounter_global_angle + delta_alpha);
        const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            source.transfer_revolution,
            record.raw_predicted_alpha,
            source,
            max_target_revolution);
        record.raw_prediction_finite = std::isfinite(record.raw_predicted_alpha);
        if (!record.raw_prediction_finite) {
            records.push_back(record);
            continue;
        }
        const auto raw_residual = problem1::evaluate_problem1_root_residual(
            departure_planet, target_planet,
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            record.raw_predicted_alpha,
            source.transfer_revolution, q_selection.selected_q);
        record.raw_residual_valid = raw_residual.valid;
        record.raw_residual_seconds = raw_residual.residual_seconds;
        record.refinement_result = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            departure_planet, target_planet,
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            source.transfer_revolution, q_selection.selected_q,
            record.raw_predicted_alpha,
            80, 1e-2, 1e-12,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        records.push_back(record);
    }
    return records;
}

std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> collect_valid_refined(
    const std::vector<PipelineRecord>& records
) {
    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> out;
    for (const auto& record : records) {
        if (record.refinement_result.valid) {
            out.push_back(record.refinement_result.branch);
        }
    }
    return out;
}

int engineering_coverage_for_k(
    const std::vector<IndexedBranch>& exact_group,
    const std::vector<IndexedBranch>& refined_group
) {
    const std::size_t n = std::min(exact_group.size(), refined_group.size());
    int covered = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (is_engineering_match(refined_group[i].branch, exact_group[i].branch)) {
            covered += 1;
        }
    }
    return covered;
}

std::string classify_failure_type(
    int exact_count,
    int refined_count,
    const std::vector<IndexedBranch>& exact_group,
    const std::vector<IndexedBranch>& refined_group
) {
    if (refined_count < exact_count) {
        return "missing";
    }
    if (refined_count > exact_count) {
        return "extra";
    }
    const std::size_t n = std::min(exact_group.size(), refined_group.size());
    for (std::size_t i = 0; i < n; ++i) {
        if (!is_engineering_match(refined_group[i].branch, exact_group[i].branch)) {
            return "wrong_root";
        }
    }
    return "unknown";
}

void print_sorted_table(
    const char* title,
    const std::vector<IndexedBranch>& group
) {
    std::cout << title << "\n";
    for (const auto& entry : group) {
        std::cout << "  original_index=" << entry.original_index
                  << " k=" << entry.branch.transfer_revolution
                  << " q=" << entry.branch.target_revolution
                  << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
                  << " alpha=" << entry.branch.encounter_global_angle
                  << " residual_seconds=" << entry.branch.residual_seconds
                  << "\n";
    }
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const auto samples = build_midcell_samples(12);

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "midcell_engineering_fail_debug\n";

    int fail_case_count = 0;

    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const auto& sample = samples[sample_index];
        const auto node = find_nearest_virtual_root_table_node(
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);

        auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet, target_planet,
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            max_transfer_revolution, max_target_revolution);
        for (auto& branch : exact_branches) {
            if (!branch.valid) {
                continue;
            }
            const auto polished = problem1::refine_problem1_root_branch_newton_seconds(
                departure_planet, target_planet,
                sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                branch.transfer_revolution, branch.target_revolution,
                branch.encounter_global_angle,
                30, 1e-6, 1e-14);
            if (polished.valid) {
                branch = polished;
            }
        }

        const auto node_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet, target_planet,
            node.nu_A_node, node.nu_B_node, node.theta_A_node,
            max_transfer_revolution, max_target_revolution);
        const auto records = build_records_for_node(
            departure_planet, target_planet, sample, node, node_branches);
        const auto before_dedup = collect_valid_refined(records);
        auto after_dedup = before_dedup;
        std::vector<problem1::Problem1SolutionBranch> deduped;
        for (const auto& b : after_dedup) {
            add_if_not_duplicate(&deduped, b);
        }
        after_dedup = deduped;

        const auto exact_by_k = group_by_k_sorted(exact_branches);
        const auto refined_by_k = group_by_k_sorted(after_dedup);

        bool sample_has_fail = false;
        for (const auto& [k, exact_group] : exact_by_k) {
            const auto refined_it = refined_by_k.find(k);
            const std::vector<IndexedBranch> empty;
            const auto& refined_group = refined_it == refined_by_k.end() ? empty : refined_it->second;
            const int exact_count = static_cast<int>(exact_group.size());
            const int refined_count = static_cast<int>(refined_group.size());
            if (engineering_coverage_for_k(exact_group, refined_group) < exact_count) {
                sample_has_fail = true;
            }
        }
        if (!sample_has_fail) {
            continue;
        }

        for (const auto& [k, exact_group] : exact_by_k) {
            const auto refined_it = refined_by_k.find(k);
            const std::vector<IndexedBranch> empty;
            const auto& refined_group = refined_it == refined_by_k.end() ? empty : refined_it->second;
            const int exact_count = static_cast<int>(exact_group.size());
            const int refined_count = static_cast<int>(refined_group.size());
            const int coverage = engineering_coverage_for_k(exact_group, refined_group);
            if (coverage == exact_count) {
                continue;
            }

            fail_case_count += 1;
            std::cout << "engineering_fail_case\n";
            std::cout << "sample_index=" << sample_index << "\n";
            std::cout << "query_nu_A=" << sample.query_nu_A << "\n";
            std::cout << "query_nu_B=" << sample.query_nu_B << "\n";
            std::cout << "query_theta_A=" << sample.query_theta_A << "\n";
            std::cout << "nearest_node_index=(" << node.nu_A_index << "," << node.nu_B_index << "," << node.theta_A_index << ")\n";
            std::cout << "node_offset_degrees=("
                      << (node.dnu_A * 180.0 / kPi) << ","
                      << (node.dnu_B * 180.0 / kPi) << ","
                      << (node.dtheta_A * 180.0 / kPi) << ")\n";
            std::cout << "k=" << k << "\n";
            std::cout << "exact_count_by_k=" << exact_count << "\n";
            std::cout << "refined_count_by_k=" << refined_count << "\n";
            std::cout << "failure_type=" << classify_failure_type(exact_count, refined_count, exact_group, refined_group) << "\n";

            print_sorted_table("exact_branches_sorted_by_time", exact_group);
            print_sorted_table("refined_branches_sorted_by_time", refined_group);

            std::cout << "paired_comparison\n";
            const std::size_t pair_count = std::min(exact_group.size(), refined_group.size());
            for (std::size_t i = 0; i < pair_count; ++i) {
                const auto& exact = exact_group[i].branch;
                const auto& refined = refined_group[i].branch;
                std::cout << "  pair_index=" << i
                          << " exact_time=" << exact.time_of_flight_seconds
                          << " refined_time=" << refined.time_of_flight_seconds
                          << " time_error=" << std::abs(refined.time_of_flight_seconds - exact.time_of_flight_seconds)
                          << " exact_alpha=" << exact.encounter_global_angle
                          << " refined_alpha=" << refined.encounter_global_angle
                          << " alpha_wrapped_error=" << wrapped_alpha_distance(refined.encounter_global_angle, exact.encounter_global_angle)
                          << " residual_seconds=" << refined.residual_seconds
                          << " q_exact=" << exact.target_revolution
                          << " q_refined=" << refined.target_revolution
                          << "\n";
            }

            std::vector<problem1::Problem1SolutionBranch> neighbor_pool;
            const int axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));
            for (int dA = -1; dA <= 1; ++dA) {
                for (int dB = -1; dB <= 1; ++dB) {
                    for (int dT = -1; dT <= 1; ++dT) {
                        auto wrap_index = [&](int base, int delta) {
                            int idx = (base + delta) % axis_count;
                            if (idx < 0) {
                                idx += axis_count;
                            }
                            return idx;
                        };
                        Problem1RootVirtualGridNode neighbor = node;
                        neighbor.nu_A_index = wrap_index(node.nu_A_index, dA);
                        neighbor.nu_B_index = wrap_index(node.nu_B_index, dB);
                        neighbor.theta_A_index = wrap_index(node.theta_A_index, dT);
                        neighbor.nu_A_node = normalize_angle_0_2pi(neighbor.nu_A_index * kVirtualGridStepRadians);
                        neighbor.nu_B_node = normalize_angle_0_2pi(neighbor.nu_B_index * kVirtualGridStepRadians);
                        neighbor.theta_A_node = normalize_angle_0_2pi(neighbor.theta_A_index * kVirtualGridStepRadians);
                        neighbor.dnu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(sample.query_nu_A) - neighbor.nu_A_node);
                        neighbor.dnu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(sample.query_nu_B) - neighbor.nu_B_node);
                        neighbor.dtheta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(sample.query_theta_A) - neighbor.theta_A_node);

                        const auto neighbor_branches = problem1::solve_problem1_from_departure_anomalies(
                            departure_planet, target_planet,
                            neighbor.nu_A_node, neighbor.nu_B_node, neighbor.theta_A_node,
                            max_transfer_revolution, max_target_revolution);
                        const auto neighbor_records = build_records_for_node(
                            departure_planet, target_planet, sample, neighbor, neighbor_branches);
                        const auto neighbor_valid = collect_valid_refined(neighbor_records);
                        for (const auto& branch : neighbor_valid) {
                            add_if_not_duplicate(&neighbor_pool, branch);
                        }
                    }
                }
            }

            const auto neighbor_by_k = group_by_k_sorted(neighbor_pool);
            const auto neighbor_it = neighbor_by_k.find(k);
            const auto& neighbor_group = neighbor_it == neighbor_by_k.end() ? empty : neighbor_it->second;
            const int single_cov = engineering_coverage_for_k(exact_group, refined_group);
            const int neighbor_cov = engineering_coverage_for_k(exact_group, neighbor_group);
            int recovered = 0;
            int neighbor_wrong = 0;
            const std::size_t neighbor_pair_count = std::min(exact_group.size(), neighbor_group.size());
            for (std::size_t i = 0; i < pair_count; ++i) {
                const bool single_match = is_engineering_match(refined_group[i].branch, exact_group[i].branch);
                const bool neighbor_match =
                    i < neighbor_pair_count && is_engineering_match(neighbor_group[i].branch, exact_group[i].branch);
                if (!single_match && neighbor_match) {
                    recovered += 1;
                }
            }
            for (std::size_t i = 0; i < neighbor_pair_count; ++i) {
                if (!is_engineering_match(neighbor_group[i].branch, exact_group[i].branch)) {
                    neighbor_wrong += 1;
                }
            }
            std::cout << "neighbor_pool_diagnostic\n";
            std::cout << "  single_node_coverage_by_k=" << single_cov << "/" << exact_count << "\n";
            std::cout << "  neighbor_pool_coverage_by_k=" << neighbor_cov << "/" << exact_count << "\n";
            std::cout << "  neighbor_pool_candidate_count_by_k=" << neighbor_group.size() << "\n";
            std::cout << "  neighbor_pool_recovered_missing_branch_count=" << recovered << "\n";
            std::cout << "  neighbor_pool_physical_wrong_root_count=" << neighbor_wrong << "\n";
        }
    }

    std::cout << "engineering_fail_case_count=" << fail_case_count << "\n";
    assert(fail_case_count > 0);
    return 0;
}
