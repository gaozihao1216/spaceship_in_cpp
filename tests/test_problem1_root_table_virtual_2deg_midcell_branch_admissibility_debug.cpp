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
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kMidCellOffsetDegrees = 1.0;

double wrapped_alpha_distance(double a, double b) {
    return std::abs(normalize_angle_minus_pi_pi(a - b));
}

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
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;
    const auto samples = build_midcell_samples(12);
    const int axis_count = static_cast<int>(std::llround(kTwoPi / kVirtualGridStepRadians));

    const std::set<std::pair<int, int>> target_cases{
        {3, 0},
        {5, 0},
        {5, 1},
    };

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "midcell_branch_admissibility_debug\n";

    for (const auto& sample : samples) {
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

        const auto nearest_node_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet, target_planet,
            nearest_node.nu_A_node, nearest_node.nu_B_node, nearest_node.theta_A_node,
            max_transfer_revolution, max_target_revolution);

        for (int k = 0; k <= max_transfer_revolution; ++k) {
            if (!target_cases.count({sample.sample_index, k})) {
                continue;
            }

            const auto exact_group = collect_sorted_k_group(exact_branches, k);
            const auto nearest_group = collect_sorted_k_group(nearest_node_branches, k);

            std::cout << "branch_admissibility_case\n";
            std::cout << "sample_index=" << sample.sample_index << "\n";
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
            std::cout << "k=" << k << "\n";

            print_branch_list("query_exact_roots_in_k", exact_group);
            print_branch_list("nearest_node_roots_in_k", nearest_group);

            int min_corner_count = std::numeric_limits<int>::max();
            int max_corner_count = 0;
            bool root_count_consistent = true;
            int reference_corner_count = -1;

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
                        const auto corner_group = collect_sorted_k_group(corner_branches, k);
                        const int count = static_cast<int>(corner_group.size());
                        min_corner_count = std::min(min_corner_count, count);
                        max_corner_count = std::max(max_corner_count, count);
                        if (reference_corner_count < 0) {
                            reference_corner_count = count;
                        } else if (reference_corner_count != count) {
                            root_count_consistent = false;
                        }

                        std::cout << "corner_roots_in_k\n";
                        std::cout << "  corner_offset=(" << da << "," << db << "," << dt << ")\n";
                        std::cout << "  node_index=("
                                  << corner.nu_A_index << ","
                                  << corner.nu_B_index << ","
                                  << corner.theta_A_index << ")\n";
                        std::cout << "  root_count_in_k=" << count << "\n";
                        for (const auto& entry : corner_group) {
                            std::cout << "    q=" << entry.branch.target_revolution
                                      << " time_of_flight_seconds=" << entry.branch.time_of_flight_seconds
                                      << " alpha=" << entry.branch.encounter_global_angle
                                      << " residual_seconds=" << entry.branch.residual_seconds
                                      << "\n";
                        }
                    }
                }
            }

            std::cout << "root_count_consistency\n";
            std::cout << "  min_root_count_in_8_corners=" << min_corner_count << "\n";
            std::cout << "  max_root_count_in_8_corners=" << max_corner_count << "\n";
            std::cout << "  query_root_count=" << exact_group.size() << "\n";
            std::cout << "  nearest_node_root_count=" << nearest_group.size() << "\n";
            std::cout << "  root_count_consistent=" << root_count_consistent << "\n";

            if (!root_count_consistent) {
                std::cout << "  interpolation_cell_admissible=false\n";
                std::cout << "  reason=branch_count_changes_inside_cell_or_between_corners\n";
            } else if (reference_corner_count != static_cast<int>(exact_group.size())) {
                std::cout << "  interpolation_cell_admissible=false\n";
                std::cout << "  reason=query_root_count_differs_from_consistent_corner_root_count\n";
            } else {
                std::cout << "  interpolation_cell_admissible=true\n";
                std::cout << "  reason=root_count_consistent_between_corners_and_query\n";
                std::cout << "time_sequence_comparison\n";
                const std::size_t n = std::min(exact_group.size(), nearest_group.size());
                for (std::size_t i = 0; i < n; ++i) {
                    std::cout << "  pair_index=" << i
                              << " corner_to_query_time_gap_seconds="
                              << std::abs(nearest_group[i].branch.time_of_flight_seconds -
                                          exact_group[i].branch.time_of_flight_seconds)
                              << " alpha_wrapped_gap="
                              << wrapped_alpha_distance(nearest_group[i].branch.encounter_global_angle,
                                                        exact_group[i].branch.encounter_global_angle)
                              << "\n";
                }
            }
        }
    }

    return 0;
}
