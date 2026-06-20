/*
 * 文件作用：验证 θ' 初扫在离散点直接按 k 分桶，且各 k 层导数与扁平列表一致。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem2/problem2_theta_prime_scan.hpp"

#include <cassert>
#include <iostream>
#include <unordered_map>
#include <utility>

namespace {

using spaceship_cpp::config::Problem2ThetaPrimeScanDefaults;
using spaceship_cpp::planet_params::PlanetId;
using spaceship_cpp::problem2::Problem2OutgoingBranchSolution;
using spaceship_cpp::problem2::Problem2ThetaPrimeScanConfig;

constexpr double kDayInSeconds = 86400.0;

struct BranchKey {
    int transfer_revolution = 0;
    int target_revolution = 0;
    double encounter_global_angle = 0.0;

    bool operator==(const BranchKey& other) const {
        return transfer_revolution == other.transfer_revolution &&
            target_revolution == other.target_revolution &&
            encounter_global_angle == other.encounter_global_angle;
    }
};

struct BranchKeyHash {
    std::size_t operator()(const BranchKey& key) const {
        const auto angle_bits = static_cast<std::size_t>(std::hash<double>{}(key.encounter_global_angle));
        return static_cast<std::size_t>(key.transfer_revolution) ^
            (static_cast<std::size_t>(key.target_revolution) << 8U) ^
            (angle_bits << 16U);
    }
};

BranchKey make_branch_key(const Problem2OutgoingBranchSolution& solution) {
    return BranchKey{
        .transfer_revolution = solution.transfer_revolution,
        .target_revolution = solution.target_revolution,
        .encounter_global_angle = solution.encounter_global_angle,
    };
}

Problem2ThetaPrimeScanConfig make_test_scan_config(
    PlanetId flyby_planet,
    PlanetId target_planet,
    double flyby_time_seconds,
    int theta_prime_count
) {
    const auto& defaults = spaceship_cpp::config::global_config();
    const Problem2ThetaPrimeScanDefaults scan_defaults{
        theta_prime_count,
        defaults.problem1_solve.phi_scan_count,
        defaults.problem2_theta_prime_scan.branch_phi_pairing_max_gap,
    };
    return spaceship_cpp::config::make_problem2_theta_prime_scan_config(
        flyby_planet,
        target_planet,
        flyby_time_seconds,
        scan_defaults,
        defaults.problem1_solve);
}

std::unordered_map<BranchKey, int, BranchKeyHash> count_solutions(
    const std::vector<Problem2OutgoingBranchSolution>& solutions
) {
    std::unordered_map<BranchKey, int, BranchKeyHash> counts;
    for (const auto& solution : solutions) {
        ++counts[make_branch_key(solution)];
    }
    return counts;
}

std::vector<Problem2OutgoingBranchSolution> merge_solutions_from_k_layers(
    const std::vector<std::vector<Problem2OutgoingBranchSolution>>& solutions_by_k
) {
    std::vector<Problem2OutgoingBranchSolution> merged;
    for (const auto& layer : solutions_by_k) {
        merged.insert(merged.end(), layer.begin(), layer.end());
    }
    return merged;
}

bool verify_k_grouping_for_scan(
    const spaceship_cpp::problem2::Problem2ThetaPrimeInitialScanResult& scan
) {
    if (!scan.ok || scan.nodes.empty()) {
        return false;
    }

    bool found_nonempty_k_layer = false;
    bool found_k_layer_derivative = false;

    for (const auto& node : scan.nodes) {
        if (static_cast<int>(node.solutions_by_k.size()) != scan.max_transfer_revolution + 1) {
            return false;
        }

        std::size_t bucket_count = 0;
        for (const auto& layer : node.solutions_by_k) {
            bucket_count += layer.size();
        }
        if (bucket_count != node.solutions.size()) {
            return false;
        }

        for (std::size_t k = 0; k < node.solutions_by_k.size(); ++k) {
            for (const auto& solution : node.solutions_by_k[k]) {
                if (solution.transfer_revolution != static_cast<int>(k)) {
                    return false;
                }
                if (solution.has_dphi_dtheta_prime && solution.has_de_dtheta_prime) {
                    found_k_layer_derivative = true;
                }
            }
            if (!node.solutions_by_k[k].empty()) {
                found_nonempty_k_layer = true;
            }
        }

        const auto flat_counts = count_solutions(node.solutions);
        const auto merged_counts = count_solutions(merge_solutions_from_k_layers(node.solutions_by_k));
        if (flat_counts != merged_counts) {
            return false;
        }

        for (const auto& solution : node.solutions) {
            const auto& layer = node.solutions_by_k[static_cast<std::size_t>(solution.transfer_revolution)];
            bool found_in_layer = false;
            for (const auto& bucketed_solution : layer) {
                if (make_branch_key(solution) == make_branch_key(bucketed_solution)) {
                    found_in_layer = true;
                    if (bucketed_solution.has_dphi_dtheta_prime != solution.has_dphi_dtheta_prime ||
                        bucketed_solution.has_de_dtheta_prime != solution.has_de_dtheta_prime ||
                        bucketed_solution.dphi_dtheta_prime != solution.dphi_dtheta_prime ||
                        bucketed_solution.de_dtheta_prime != solution.de_dtheta_prime) {
                        return false;
                    }
                    break;
                }
            }
            if (!found_in_layer) {
                return false;
            }
        }
    }

    return found_nonempty_k_layer && found_k_layer_derivative;
}

bool try_verify_k_grouping_for_scenario(
    PlanetId flyby_planet,
    PlanetId target_planet,
    double flyby_time_seconds
) {
    const auto config = make_test_scan_config(flyby_planet, target_planet, flyby_time_seconds, 32);
    const auto scan = spaceship_cpp::problem2::run_problem2_theta_prime_initial_scan(config);
    if (!verify_k_grouping_for_scan(scan)) {
        return false;
    }
    std::cout << "k_grouping_verified_for_flyby_scenario\n";
    return true;
}

}  // namespace

int main() {
    namespace problem2 = spaceship_cpp::problem2;

    {
        const auto config = make_test_scan_config(PlanetId::Earth, PlanetId::Mars, 0.0, 32);
        const auto scan = problem2::run_problem2_theta_prime_initial_scan(config);
        assert(scan.ok);
        assert(scan.max_transfer_revolution == config.problem1_solve.max_transfer_revolution);
        assert(verify_k_grouping_for_scan(scan));
    }

    {
        bool verified = false;
        const std::pair<PlanetId, PlanetId> planet_pairs[] = {
            {PlanetId::Earth, PlanetId::Mars},
            {PlanetId::Mars, PlanetId::Earth},
            {PlanetId::Venus, PlanetId::Earth},
        };
        const double flyby_times[] = {0.0, 100.0 * kDayInSeconds};
        for (const auto& [flyby_planet, target_planet] : planet_pairs) {
            for (const double flyby_time : flyby_times) {
                if (try_verify_k_grouping_for_scenario(flyby_planet, target_planet, flyby_time)) {
                    verified = true;
                    break;
                }
            }
            if (verified) {
                break;
            }
        }
        assert(verified);
    }

    std::cout << "test_problem2_theta_prime_k_grouping PASSED\n";
    return 0;
}
