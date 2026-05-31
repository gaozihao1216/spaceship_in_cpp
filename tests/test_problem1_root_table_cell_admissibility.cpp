#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr double kVirtualGridStepRadians = kPi / 90.0;
constexpr double kMidCellOffsetDegrees = 1.0;

struct QuerySample {
    int sample_index = -1;
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
    node.dnu_A = spaceship_cpp::common::normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - node.nu_A_node);
    node.dnu_B = spaceship_cpp::common::normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - node.nu_B_node);
    node.dtheta_A = spaceship_cpp::common::normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - node.theta_A_node);
    return node;
}

std::vector<QuerySample> build_midcell_samples(int samples_per_group) {
    std::vector<QuerySample> samples;
    const double mid_offset = kMidCellOffsetDegrees * kPi / 180.0;
    for (int i = 0; i < samples_per_group; ++i) {
        const double base_nu_A = normalize_angle_0_2pi(0.91 + static_cast<double>(i) * 1.8123456789);
        const double base_nu_B = normalize_angle_0_2pi(1.41 + static_cast<double>(i) * 2.2718281828);
        const double base_theta_A = normalize_angle_0_2pi(0.63 + static_cast<double>(i) * 1.1447298860);
        const auto nearest = find_nearest_virtual_root_table_node(base_nu_A, base_nu_B, base_theta_A);
        samples.push_back({
            i,
            normalize_angle_0_2pi(nearest.nu_A_node + mid_offset),
            normalize_angle_0_2pi(nearest.nu_B_node + mid_offset),
            normalize_angle_0_2pi(nearest.theta_A_node + mid_offset),
        });
    }
    return samples;
}

int lower_corner_index(int nearest_index, double offset_radians) {
    return offset_radians >= 0.0 ? nearest_index : nearest_index - 1;
}

spaceship_cpp::problem1::Problem1RootTable build_local_cell_table_for_query(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const QuerySample& sample,
    int max_transfer_revolution,
    int max_target_revolution,
    int* lower_i,
    int* lower_j,
    int* lower_k
) {
    const auto nearest = find_nearest_virtual_root_table_node(
        sample.query_nu_A, sample.query_nu_B, sample.query_theta_A);
    *lower_i = lower_corner_index(nearest.nu_A_index, nearest.dnu_A);
    *lower_j = lower_corner_index(nearest.nu_B_index, nearest.dnu_B);
    *lower_k = lower_corner_index(nearest.theta_A_index, nearest.dtheta_A);

    spaceship_cpp::problem1::Problem1RootTableConfig config{};
    config.departure_planet = departure_planet;
    config.target_planet = target_planet;
    config.nu_A_start = static_cast<double>(*lower_i) * kVirtualGridStepRadians;
    config.nu_A_step = kVirtualGridStepRadians;
    config.nu_A_count = 2;
    config.nu_B_depart_start = static_cast<double>(*lower_j) * kVirtualGridStepRadians;
    config.nu_B_depart_step = kVirtualGridStepRadians;
    config.nu_B_depart_count = 2;
    config.theta_A_start = static_cast<double>(*lower_k) * kVirtualGridStepRadians;
    config.theta_A_step = kVirtualGridStepRadians;
    config.theta_A_count = 2;
    config.max_transfer_revolution = max_transfer_revolution;
    config.max_target_revolution = max_target_revolution;
    return spaceship_cpp::problem1::build_problem1_root_table(config);
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    namespace planet_params = spaceship_cpp::planet_params;

    std::cout << std::setprecision(6) << std::scientific;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const int max_transfer_revolution = 1;
    const int max_target_revolution = 1;

    problem1::Problem1RootTableConfig simple_config{};
    simple_config.departure_planet = departure_planet;
    simple_config.target_planet = target_planet;
    simple_config.nu_A_start = 0.0;
    simple_config.nu_A_step = kVirtualGridStepRadians;
    simple_config.nu_A_count = 2;
    simple_config.nu_B_depart_start = 0.0;
    simple_config.nu_B_depart_step = kVirtualGridStepRadians;
    simple_config.nu_B_depart_count = 2;
    simple_config.theta_A_start = 0.0;
    simple_config.theta_A_step = kVirtualGridStepRadians;
    simple_config.theta_A_count = 2;
    simple_config.max_transfer_revolution = 0;
    simple_config.max_target_revolution = 1;
    problem1::Problem1RootTable synthetic(simple_config);
    for (auto& cell : synthetic.mutable_cells()) {
        cell.solved = true;
        problem1::Problem1SolutionBranch a{};
        a.valid = true;
        a.transfer_revolution = 0;
        a.target_revolution = 0;
        a.time_of_flight_seconds = 10.0;
        problem1::Problem1SolutionBranch b = a;
        b.time_of_flight_seconds = 20.0;
        cell.solutions_sorted_by_time_of_flight = {a, b};
    }
    synthetic.mutable_cells().back().solutions_sorted_by_time_of_flight.pop_back();
    const auto synthetic_result = problem1::evaluate_problem1_root_cell_admissibility(
        synthetic, 0, 0, 0, 0, 1);
    std::cout << "synthetic_nonadmissible=" << synthetic_result.admissible << '\n';
    std::cout << "synthetic_reason=" << synthetic_result.reason << '\n';
    assert(!synthetic_result.admissible);
    assert(synthetic_result.reason == "branch_count_inconsistent_across_corners");

    const auto midcell_samples = build_midcell_samples(12);

    int lower_i = 0;
    int lower_j = 0;
    int lower_k = 0;
    const auto sample3_table = build_local_cell_table_for_query(
        departure_planet,
        target_planet,
        midcell_samples.at(3),
        max_transfer_revolution,
        max_target_revolution,
        &lower_i,
        &lower_j,
        &lower_k);
    const auto sample3_result = problem1::evaluate_problem1_root_cell_admissibility(
        sample3_table, 0, 0, 0, 0, max_target_revolution);
    std::cout << "sample3_lower_indices=(" << lower_i << "," << lower_j << "," << lower_k << ")\n";
    std::cout << "sample3_admissible=" << sample3_result.admissible << '\n';
    std::cout << "sample3_reason=" << sample3_result.reason << '\n';
    assert(sample3_result.admissible);
    assert(sample3_result.reason == "corner_branch_count_consistent");

    const auto sample5_table = build_local_cell_table_for_query(
        departure_planet,
        target_planet,
        midcell_samples.at(5),
        max_transfer_revolution,
        max_target_revolution,
        &lower_i,
        &lower_j,
        &lower_k);
    const auto sample5_result = problem1::evaluate_problem1_root_cell_admissibility(
        sample5_table, 0, 0, 0, 0, max_target_revolution);
    std::cout << "sample5_lower_indices=(" << lower_i << "," << lower_j << "," << lower_k << ")\n";
    std::cout << "sample5_admissible=" << sample5_result.admissible << '\n';
    std::cout << "sample5_reason=" << sample5_result.reason << '\n';
    assert(!sample5_result.admissible);
    assert(sample5_result.reason == "branch_count_inconsistent_across_corners");

    std::cout << "cell_admissibility_checker_ok\n";
    return 0;
}
