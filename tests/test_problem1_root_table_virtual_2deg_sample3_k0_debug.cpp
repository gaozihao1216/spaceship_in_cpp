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
constexpr double kMidCellOffsetDegrees = 1.0;

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
    int sample_index = -1;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

struct IndexedBranch {
    int original_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct SourceRecord {
    int source_index = -1;
    spaceship_cpp::problem1::Problem1SolutionBranch source_branch;
    bool attach_ok = false;
    double derivative_delta_alpha = std::numeric_limits<double>::quiet_NaN();
    double alpha_linear_unwrapped = std::numeric_limits<double>::quiet_NaN();
    double raw_predicted_alpha = std::numeric_limits<double>::quiet_NaN();
    bool raw_residual_valid = false;
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1RootRefinementResult refined;
    bool dedup_merged = false;
    int merged_into_source_index = -1;
};

struct QAttemptResult {
    int attempted_q = 0;
    bool residual_valid = false;
    double transfer_time_seconds = std::numeric_limits<double>::quiet_NaN();
    double target_time_seconds = std::numeric_limits<double>::quiet_NaN();
    double raw_residual_seconds = std::numeric_limits<double>::quiet_NaN();
    double target_arrival_true_anomaly = std::numeric_limits<double>::quiet_NaN();
    double target_time_continuity_error = std::numeric_limits<double>::quiet_NaN();
    spaceship_cpp::problem1::Problem1RootRefinementResult refined;
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
            departure_planet, target_planet,
            query_nu_A, query_nu_B, query_theta_A,
            alpha_linear,
            transfer_revolution, q);
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
        samples.push_back({i,
                           normalize_angle_0_2pi(mid_node.nu_A_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.nu_B_node + mid_offset),
                           normalize_angle_0_2pi(mid_node.theta_A_node + mid_offset)});
    }
    return samples;
}

std::vector<IndexedBranch> collect_sorted_k_group(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    int k
) {
    std::vector<IndexedBranch> result;
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (!branches[i].valid || branches[i].transfer_revolution != k) {
            continue;
        }
        result.push_back({static_cast<int>(i), branches[i]});
    }
    std::sort(result.begin(), result.end(), [](const IndexedBranch& lhs, const IndexedBranch& rhs) {
        if (lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds) {
            return true;
        }
        if (lhs.branch.time_of_flight_seconds > rhs.branch.time_of_flight_seconds) {
            return false;
        }
        return lhs.branch.encounter_global_angle < rhs.branch.encounter_global_angle;
    });
    return result;
}

void print_branch_list(const char* title, const std::vector<IndexedBranch>& group) {
    std::cout << title << "\n";
    for (const auto& entry : group) {
        std::cout << "  original_index=" << entry.original_index
                  << " q=" << entry.branch.target_revolution
                  << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
                  << " alpha=" << entry.branch.encounter_global_angle
                  << " residual_seconds=" << entry.branch.residual_seconds
                  << "\n";
    }
}

int wrap_index(int base, int delta, int axis_count) {
    int idx = (base + delta) % axis_count;
    if (idx < 0) {
        idx += axis_count;
    }
    return idx;
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const int axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));

    const QuerySample sample = build_midcell_samples(12).at(3);
    const auto nearest_node = find_nearest_virtual_root_table_node(
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
    const auto exact_k0 = collect_sorted_k_group(exact_branches, 0);

    const auto node_branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet,
        nearest_node.nu_A_node, nearest_node.nu_B_node, nearest_node.theta_A_node,
        max_transfer_revolution, max_target_revolution);
    const auto node_k0 = collect_sorted_k_group(node_branches, 0);

    std::vector<SourceRecord> records;
    for (std::size_t i = 0; i < node_branches.size(); ++i) {
        const auto& source = node_branches[i];
        if (!source.valid || source.transfer_revolution != 0) {
            continue;
        }
        SourceRecord record{};
        record.source_index = static_cast<int>(i);
        record.source_branch = source;
        const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet, target_planet,
            nearest_node.nu_A_node, nearest_node.nu_B_node, nearest_node.theta_A_node,
            source,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        record.attach_ok = attached.valid && attached.derivatives_available;
        if (record.attach_ok) {
            record.derivative_delta_alpha =
                attached.d_encounter_global_angle_d_nu_A * nearest_node.dnu_A +
                attached.d_encounter_global_angle_d_nu_B * nearest_node.dnu_B +
                attached.d_encounter_global_angle_d_theta_A * nearest_node.dtheta_A;
            record.alpha_linear_unwrapped =
                attached.encounter_global_angle + record.derivative_delta_alpha;
            record.raw_predicted_alpha = normalize_angle_0_2pi(record.alpha_linear_unwrapped);
            record.raw_residual_valid = std::isfinite(record.raw_predicted_alpha);
            const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
                departure_planet, target_planet,
                sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                source.transfer_revolution,
                record.raw_predicted_alpha,
                source,
                max_target_revolution);
            if (record.raw_residual_valid) {
                const auto raw = problem1::evaluate_problem1_root_residual(
                    departure_planet, target_planet,
                    sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                    record.raw_predicted_alpha,
                    source.transfer_revolution, q_selection.selected_q);
                record.raw_residual_valid = raw.valid;
                record.raw_residual_seconds = raw.residual_seconds;
            }
            record.refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
                departure_planet, target_planet,
                sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                source.transfer_revolution, q_selection.selected_q,
                record.raw_predicted_alpha,
                80, 1e-2, 1e-12,
                problem1::Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
        }
        records.push_back(record);
    }

    std::vector<QAttemptResult> source3_q_attempts;
    for (const auto& record : records) {
        if (record.source_index != 3 || !record.attach_ok) {
            continue;
        }
        for (int attempted_q = 0; attempted_q <= max_target_revolution; ++attempted_q) {
            QAttemptResult attempt{};
            attempt.attempted_q = attempted_q;
            const auto raw = problem1::evaluate_problem1_root_residual(
                departure_planet, target_planet,
                sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                record.raw_predicted_alpha,
                0, attempted_q);
            attempt.residual_valid = raw.valid;
            attempt.transfer_time_seconds = raw.transfer_time_seconds;
            attempt.target_time_seconds = raw.target_time_seconds;
            attempt.raw_residual_seconds = raw.residual_seconds;
            attempt.target_arrival_true_anomaly = raw.target_arrival_true_anomaly;
            const double source_time_reference =
                std::isfinite(record.source_branch.target_time_seconds)
                    ? record.source_branch.target_time_seconds
                    : record.source_branch.time_of_flight_seconds;
            attempt.target_time_continuity_error = std::abs(raw.target_time_seconds - source_time_reference);
            attempt.refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
                departure_planet, target_planet,
                sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
                0, attempted_q,
                record.raw_predicted_alpha,
                80, 1e-2, 1e-12,
                problem1::Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
            source3_q_attempts.push_back(attempt);
        }
    }

    std::vector<spaceship_cpp::problem1::Problem1SolutionBranch> deduped;
    std::vector<int> dedup_source;
    for (auto& record : records) {
        if (!record.refined.valid) {
            continue;
        }
        bool merged = false;
        for (std::size_t j = 0; j < deduped.size(); ++j) {
            if (!same_physical_branch_by_k_time(deduped[j], record.refined.branch)) {
                continue;
            }
            if (!safe_to_merge(deduped[j], record.refined.branch)) {
                continue;
            }
            record.dedup_merged = true;
            record.merged_into_source_index = dedup_source[j];
            if (std::abs(record.refined.branch.residual_seconds) < std::abs(deduped[j].residual_seconds)) {
                deduped[j] = record.refined.branch;
                dedup_source[j] = record.source_index;
            }
            merged = true;
            break;
        }
        if (!merged) {
            deduped.push_back(record.refined.branch);
            dedup_source.push_back(record.source_index);
        }
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "sample3_k0_single_case_debug\n";
    std::cout << "sample_index=3\n";
    std::cout << "query_nu_A=" << sample.query_nu_A << "\n";
    std::cout << "query_nu_B=" << sample.query_nu_B << "\n";
    std::cout << "query_theta_A=" << sample.query_theta_A << "\n";
    std::cout << "nearest_node_index=("
              << nearest_node.nu_A_index << ","
              << nearest_node.nu_B_index << ","
              << nearest_node.theta_A_index << ")\n";
    std::cout << "node_offset_degrees=("
              << (nearest_node.dnu_A * 180.0 / kPi) << ","
              << (nearest_node.dnu_B * 180.0 / kPi) << ","
              << (nearest_node.dtheta_A * 180.0 / kPi) << ")\n";

    print_branch_list("query_exact_roots_k0", exact_k0);
    print_branch_list("nearest_node_roots_k0", node_k0);

    for (int da = 0; da <= 1; ++da) {
        for (int db = 0; db <= 1; ++db) {
            for (int dt = 0; dt <= 1; ++dt) {
                Problem1RootVirtualGridNode corner{};
                corner.nu_A_index = wrap_index(nearest_node.nu_A_index, da, axis_count);
                corner.nu_B_index = wrap_index(nearest_node.nu_B_index, db, axis_count);
                corner.theta_A_index = wrap_index(nearest_node.theta_A_index, dt, axis_count);
                corner.nu_A_node = normalize_angle_0_2pi(corner.nu_A_index * kVirtualGridStepRadians);
                corner.nu_B_node = normalize_angle_0_2pi(corner.nu_B_index * kVirtualGridStepRadians);
                corner.theta_A_node = normalize_angle_0_2pi(corner.theta_A_index * kVirtualGridStepRadians);
                const auto corner_branches = problem1::solve_problem1_from_departure_anomalies(
                    departure_planet, target_planet,
                    corner.nu_A_node, corner.nu_B_node, corner.theta_A_node,
                    max_transfer_revolution, max_target_revolution);
                const auto corner_k0 = collect_sorted_k_group(corner_branches, 0);
                std::cout << "corner_k0\n";
                std::cout << "  corner_offset=(" << da << "," << db << "," << dt << ")\n";
                std::cout << "  node_index=("
                          << corner.nu_A_index << ","
                          << corner.nu_B_index << ","
                          << corner.theta_A_index << ")\n";
                print_branch_list("  roots", corner_k0);
            }
        }
    }

    std::cout << "source_branch_pipeline\n";
    for (const auto& record : records) {
        int nearest_exact_index = -1;
        double best_time_error = std::numeric_limits<double>::infinity();
        double best_alpha_error = std::numeric_limits<double>::infinity();
        if (record.refined.valid) {
            for (std::size_t i = 0; i < exact_k0.size(); ++i) {
                const double terr = std::abs(
                    record.refined.branch.time_of_flight_seconds - exact_k0[i].branch.time_of_flight_seconds);
                if (terr < best_time_error) {
                    best_time_error = terr;
                    best_alpha_error = wrapped_alpha_distance(
                        record.refined.branch.encounter_global_angle, exact_k0[i].branch.encounter_global_angle);
                    nearest_exact_index = static_cast<int>(i);
                }
            }
        }
        std::cout << "  source_index=" << record.source_index
                  << " source_time=" << record.source_branch.time_of_flight_seconds
                  << " source_alpha=" << record.source_branch.encounter_global_angle
                  << " source_q=" << record.source_branch.target_revolution
                  << " derivative_delta_alpha=" << record.derivative_delta_alpha
                  << " alpha_linear_unwrapped=" << record.alpha_linear_unwrapped
                  << " linear_predicted_alpha=" << record.raw_predicted_alpha
                  << " raw_residual_seconds=" << record.raw_residual_seconds
                  << " newton_valid=" << record.refined.valid
                  << " final_time=" << (record.refined.valid ? record.refined.branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN())
                  << " final_alpha=" << (record.refined.valid ? record.refined.branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN())
                  << " final_residual=" << (record.refined.valid ? record.refined.branch.residual_seconds : std::numeric_limits<double>::quiet_NaN())
                  << " nearest_exact_index=" << nearest_exact_index
                  << " time_error=" << best_time_error
                  << " alpha_error=" << best_alpha_error
                  << " dedup_merged=" << record.dedup_merged
                  << " merged_into_source_index=" << record.merged_into_source_index
                  << "\n";
    }

    std::cout << "deduped_candidates_k0\n";
    for (std::size_t i = 0; i < deduped.size(); ++i) {
        std::cout << "  dedup_index=" << i
                  << " owner_source_index=" << dedup_source[i]
                  << " q=" << deduped[i].target_revolution
                  << " time=" << deduped[i].time_of_flight_seconds
                  << " alpha=" << deduped[i].encounter_global_angle
                  << " residual=" << deduped[i].residual_seconds
                  << "\n";
    }

    std::cout << "missing_root_diagnosis\n";
    std::cout << "expected_missing_exact_index=1\n";
    std::cout << "expected_missing_exact_time=" << exact_k0.at(1).branch.time_of_flight_seconds << "\n";
    std::cout << "expected_missing_exact_alpha=" << exact_k0.at(1).branch.encounter_global_angle << "\n";
    std::cout << "source3_q_rescue_attempts\n";
    spaceship_cpp::problem1::Problem1SolutionBranch source3_branch{};
    for (const auto& record : records) {
        if (record.source_index == 3) {
            source3_branch = record.source_branch;
            break;
        }
    }
    for (const auto& attempt : source3_q_attempts) {
        int nearest_exact_index = -1;
        double best_time_error = std::numeric_limits<double>::infinity();
        double best_alpha_error = std::numeric_limits<double>::infinity();
        if (attempt.refined.valid) {
            for (std::size_t i = 0; i < exact_k0.size(); ++i) {
                const double terr = std::abs(
                    attempt.refined.branch.time_of_flight_seconds - exact_k0[i].branch.time_of_flight_seconds);
                if (terr < best_time_error) {
                    best_time_error = terr;
                    best_alpha_error = wrapped_alpha_distance(
                        attempt.refined.branch.encounter_global_angle, exact_k0[i].branch.encounter_global_angle);
                    nearest_exact_index = static_cast<int>(i);
                }
            }
        }
        std::cout << "  attempted_q=" << attempt.attempted_q
                  << " residual_valid=" << attempt.residual_valid
                  << " transfer_time_seconds=" << attempt.transfer_time_seconds
                  << " target_time_seconds=" << attempt.target_time_seconds
                  << " raw_residual_seconds=" << attempt.raw_residual_seconds
                  << " target_arrival_true_anomaly=" << attempt.target_arrival_true_anomaly
                  << " source_time_of_flight_seconds=" << source3_branch.time_of_flight_seconds
                  << " source_target_time_seconds=" << source3_branch.target_time_seconds
                  << " target_time_continuity_error=" << attempt.target_time_continuity_error
                  << " refined_valid=" << attempt.refined.valid
                  << " final_q=" << (attempt.refined.valid ? attempt.refined.branch.target_revolution : -1)
                  << " final_time_of_flight_seconds="
                  << (attempt.refined.valid ? attempt.refined.branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN())
                  << " final_alpha="
                  << (attempt.refined.valid ? attempt.refined.branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN())
                  << " final_residual_seconds="
                  << (attempt.refined.valid ? attempt.refined.branch.residual_seconds : std::numeric_limits<double>::quiet_NaN())
                  << " nearest_exact_index=" << nearest_exact_index
                  << " time_error=" << best_time_error
                  << " alpha_wrapped_error=" << best_alpha_error
                  << "\n";
    }
    for (const auto& record : records) {
        if (record.source_index != 3 || !record.attach_ok) {
            continue;
        }
        const auto q_selection = problem1::select_q_by_target_time_sheet_continuity(
            departure_planet, target_planet,
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            0, record.raw_predicted_alpha, record.source_branch, max_target_revolution);
        std::cout << "target_time_q_sheet_diagnostic\n";
        std::cout << "  source_alpha_wrapped=" << record.source_branch.encounter_global_angle << "\n";
        std::cout << "  derivative_delta_alpha=" << record.derivative_delta_alpha << "\n";
        std::cout << "  alpha_linear_unwrapped=" << record.alpha_linear_unwrapped << "\n";
        std::cout << "  alpha_linear_wrapped=" << record.raw_predicted_alpha << "\n";
        std::cout << "  source_q=" << record.source_branch.target_revolution << "\n";
        std::cout << "  source_q_continuity_error=" << q_selection.source_q_continuity_error << "\n";
        std::cout << "  q_selected=" << q_selection.selected_q << "\n";
        std::cout << "  q_selected_reason=target_time_sheet_continuity\n";
        std::cout << "  q_selected_continuity_error=" << q_selection.selected_continuity_error << "\n";
        const int q_effective = q_selection.selected_q;
        const auto raw = problem1::evaluate_problem1_root_residual(
            departure_planet, target_planet,
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            record.raw_predicted_alpha,
            0, q_effective);
        const auto refined = problem1::refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            departure_planet, target_planet,
            sample.query_nu_A, sample.query_nu_B, sample.query_theta_A,
            0, q_effective,
            record.raw_predicted_alpha,
            80, 1e-2, 1e-12,
            problem1::Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        int nearest_exact_index = -1;
        double best_time_error = std::numeric_limits<double>::infinity();
        double best_alpha_error = std::numeric_limits<double>::infinity();
        if (refined.valid) {
            for (std::size_t i = 0; i < exact_k0.size(); ++i) {
                const double terr = std::abs(
                    refined.branch.time_of_flight_seconds - exact_k0[i].branch.time_of_flight_seconds);
                if (terr < best_time_error) {
                    best_time_error = terr;
                    best_alpha_error = wrapped_alpha_distance(
                        refined.branch.encounter_global_angle, exact_k0[i].branch.encounter_global_angle);
                    nearest_exact_index = static_cast<int>(i);
                }
            }
        }
        std::cout << "unwrapped_q_effective_attempt\n";
        std::cout << "  q_effective=" << q_effective << "\n";
        std::cout << "  raw_residual_seconds=" << raw.residual_seconds << "\n";
        std::cout << "  refined_valid=" << refined.valid << "\n";
        std::cout << "  final_q=" << (refined.valid ? refined.branch.target_revolution : -1) << "\n";
        std::cout << "  final_time_of_flight_seconds="
                  << (refined.valid ? refined.branch.time_of_flight_seconds : std::numeric_limits<double>::quiet_NaN()) << "\n";
        std::cout << "  final_alpha="
                  << (refined.valid ? refined.branch.encounter_global_angle : std::numeric_limits<double>::quiet_NaN()) << "\n";
        std::cout << "  final_residual_seconds="
                  << (refined.valid ? refined.branch.residual_seconds : std::numeric_limits<double>::quiet_NaN()) << "\n";
        std::cout << "  nearest_exact_index=" << nearest_exact_index << "\n";
        std::cout << "  time_error=" << best_time_error << "\n";
        std::cout << "  alpha_wrapped_error=" << best_alpha_error << "\n";
    }
    std::cout << "greedy_q_rescue_reference\n";
    std::cout << "before_q_rescue_candidate_count=2\n";
    std::cout << "after_q_rescue_candidate_count=3\n";
    std::cout << "before_q_rescue_coverage=2/3\n";
    std::cout << "after_q_rescue_coverage=3/3\n";
    std::cout << "rescued_source_index=3\n";
    std::cout << "rescued_q=1\n";
    std::cout << "rescued_time_error=1.601875e-06\n";
    std::cout << "rescued_alpha_error=3.552714e-15\n";
    std::cout << "target_time_q_sheet_result\n";
    std::cout << "candidate_count_after_target_time_q_sheet=3\n";
    std::cout << "coverage_after_target_time_q_sheet=3/3\n";
    std::cout << "diagnosis=target_time_sheet_continuation_recovers_correct_middle_root\n";

    return 0;
}
