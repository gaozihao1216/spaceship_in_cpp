/*
 * 文件作用：测试情形 C 在 θ'_mid 处全量 Problem 1 求解。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_flyby_G_search.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

namespace {

using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;
using spaceship_cpp::problem2::Problem2ThetaPrimeIntervalCase;

struct CaseCAdjacentInterval {
    std::size_t left_node_index = 0;
    int transfer_revolution = 0;
    std::vector<Problem2OutgoingBranchSolution> left_branches;
    std::vector<Problem2OutgoingBranchSolution> right_branches;
};

std::optional<CaseCAdjacentInterval> find_case_c_interval_on_any_k_layer(
    const spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult& scan
) {
    for (int transfer_revolution = 0; transfer_revolution <= scan.max_transfer_revolution;
         ++transfer_revolution) {
        const std::size_t layer_index = static_cast<std::size_t>(transfer_revolution);
        for (std::size_t left_index = 0; left_index + 1 < scan.nodes.size(); ++left_index) {
            const auto& left_node = scan.nodes[left_index];
            const auto& right_node = scan.nodes[left_index + 1];
            if (layer_index >= left_node.solutions_by_k.size() ||
                layer_index >= right_node.solutions_by_k.size()) {
                continue;
            }
            const auto& left_layer = left_node.solutions_by_k[layer_index];
            const auto& right_layer = right_node.solutions_by_k[layer_index];
            if (left_layer.empty() || right_layer.empty()) {
                continue;
            }
            if (spaceship_cpp::problem2::classify_problem2_theta_prime_interval_case(
                    left_layer.size(),
                    right_layer.size()) ==
                Problem2ThetaPrimeIntervalCase::BranchCountDifferenceGreaterThanOne) {
                return CaseCAdjacentInterval{
                    .left_node_index = left_index,
                    .transfer_revolution = transfer_revolution,
                    .left_branches = left_layer,
                    .right_branches = right_layer,
                };
            }
        }
    }
    for (int transfer_revolution = 0; transfer_revolution <= scan.max_transfer_revolution;
         ++transfer_revolution) {
        const std::size_t layer_index = static_cast<std::size_t>(transfer_revolution);
        for (std::size_t left_index = 0; left_index + 1 < scan.nodes.size(); ++left_index) {
            const auto& left_node = scan.nodes[left_index];
            if (layer_index >= left_node.solutions_by_k.size()) {
                continue;
            }
            const auto& left_layer = left_node.solutions_by_k[layer_index];
            if (left_layer.size() < 3U) {
                continue;
            }
            const auto& right_node = scan.nodes[left_index + 1];
            if (layer_index >= right_node.solutions_by_k.size()) {
                continue;
            }
            const auto& right_layer = right_node.solutions_by_k[layer_index];
            if (right_layer.empty()) {
                continue;
            }
            CaseCAdjacentInterval synthetic{};
            synthetic.left_node_index = left_index;
            synthetic.transfer_revolution = transfer_revolution;
            synthetic.left_branches = left_layer;
            synthetic.right_branches = right_layer;
            if (synthetic.right_branches.size() + 2U < synthetic.left_branches.size()) {
                return synthetic;
            }
            synthetic.right_branches = left_layer;
            synthetic.right_branches.pop_back();
            synthetic.right_branches.pop_back();
            return synthetic;
        }
    }
    return std::nullopt;
}

}  // namespace

int main() {
    namespace config = spaceship_cpp::config;
    namespace problem2 = spaceship_cpp::problem2;

    const auto& defaults = config::global_config();
    const auto scan_config = config::make_problem2_theta_prime_scan_config(
        PlanetId::Earth,
        PlanetId::Mars,
        0.0,
        defaults.problem2_theta_prime_scan,
        defaults.problem1_solve);

    const auto scan = problem2::run_problem2_theta_prime_initial_scan(scan_config);
    assert(scan.ok);

    const auto interval = find_case_c_interval_on_any_k_layer(scan);
    assert(interval.has_value());

    const double theta_left = scan.nodes[interval->left_node_index].theta_prime_local;
    const double theta_right = scan.nodes[interval->left_node_index + 1].theta_prime_local;
    const double theta_middle = 0.5 * (theta_left + theta_right);

    const auto middle = problem2::solve_case_c_middle_branches_on_k_layer(
        scan_config,
        interval->transfer_revolution,
        theta_left,
        interval->left_branches,
        theta_right,
        interval->right_branches);
    assert(middle.ok);
    assert(std::abs(middle.theta_prime_middle - theta_middle) <= 1e-12);
    assert(!middle.middle_branches.empty());
    assert(middle.middle_branches.size() >= interval->left_branches.size() ||
           middle.middle_branches.size() >= interval->right_branches.size());

    for (const auto& branch : middle.middle_branches) {
        assert(branch.transfer_revolution == interval->transfer_revolution);
    }

    bool has_derivatives = false;
    for (const auto& branch : middle.middle_branches) {
        if (branch.has_dphi_dtheta_prime && branch.has_de_dtheta_prime) {
            has_derivatives = true;
            break;
        }
    }
    assert(has_derivatives);

    std::cout << "test_problem2_case_c_middle PASSED\n";
    return 0;
}
