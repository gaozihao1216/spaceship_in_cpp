#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

struct PairCase {
    spaceship_cpp::planet_params::PlanetId from;
    spaceship_cpp::planet_params::PlanetId to;
};

int env_int(const char* name, int fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    const int value = std::atoi(raw);
    return value > 0 ? value : fallback;
}

unsigned env_uint(const char* name, unsigned fallback) {
    const char* raw = std::getenv(name);
    if (raw == nullptr || *raw == '\0') return fallback;
    const unsigned long value = std::strtoul(raw, nullptr, 10);
    return value > 0 ? static_cast<unsigned>(value) : fallback;
}

double angle_diff(double lhs, double rhs) {
    double diff = std::fmod(lhs - rhs + spaceship_cpp::common::kPi, spaceship_cpp::common::kTwoPi);
    if (diff < 0.0) diff += spaceship_cpp::common::kTwoPi;
    return std::abs(diff - spaceship_cpp::common::kPi);
}

bool branches_equal(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& lhs,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& rhs
) {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const auto& a = lhs[i];
        const auto& b = rhs[i];
        if (a.valid != b.valid ||
            a.transfer_revolution != b.transfer_revolution ||
            a.target_revolution != b.target_revolution ||
            angle_diff(a.target_arrival_true_anomaly, b.target_arrival_true_anomaly) > 1e-12 ||
            angle_diff(a.encounter_global_angle, b.encounter_global_angle) > 1e-12 ||
            std::abs(a.target_time_seconds - b.target_time_seconds) > 1e-9 ||
            std::abs(a.residual_seconds - b.residual_seconds) > 1e-12) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const int sample_count = env_int("ROOT_TABLE_5DEG_TEST_SAMPLE_COUNT", 20);
    const unsigned seed = env_uint("ROOT_TABLE_5DEG_TEST_SEED", 12345);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> angle_dist(0.0, common::kTwoPi);

    {
        std::vector<problem1::Problem1SolutionBranch> baseline(2);
        baseline[0].valid = true;
        baseline[0].encounter_global_angle = 0.25;
        baseline[0].target_arrival_true_anomaly = 0.75;
        baseline[0].transfer_revolution = 0;
        baseline[0].target_revolution = 0;
        baseline[0].time_of_flight_seconds = 1000.0;
        baseline[0].target_time_seconds = 1000.0;
        baseline[0].residual_seconds = 1e-9;

        baseline[1].valid = true;
        baseline[1].encounter_global_angle = 1.25;
        baseline[1].target_arrival_true_anomaly = 1.75;
        baseline[1].transfer_revolution = 1;
        baseline[1].target_revolution = 0;
        baseline[1].time_of_flight_seconds = 2000.0;
        baseline[1].target_time_seconds = 2000.0;
        baseline[1].residual_seconds = -2e-9;

        std::vector<problem1::Problem1SolutionBranch> shuffled{baseline[1], baseline[0]};
        const auto report = problem1::compare_problem1_solution_branches_by_identity(
            shuffled, baseline, 0.1, 0.2, 0.3);
        assert(report.matched_count == 2);
        assert(report.missing_count == 0);
        assert(report.extra_count == 0);
        assert(report.missing_branches.empty());
        assert(report.extra_branches.empty());

        const std::vector<problem1::Problem1SolutionBranch> missing_candidate{baseline[1]};
        const auto missing_report = problem1::compare_problem1_solution_branches_by_identity(
            missing_candidate, baseline, 0.1, 0.2, 0.3);
        assert(missing_report.matched_count == 1);
        assert(missing_report.missing_count == 1);
        assert(missing_report.extra_count == 0);
        assert(missing_report.missing_branches.size() == 1);
        assert(missing_report.missing_branches[0].transfer_revolution == baseline[0].transfer_revolution);
        assert(missing_report.missing_branches[0].target_revolution == baseline[0].target_revolution);
        assert(missing_report.missing_branches[0].nu_A_depart == 0.1);
        assert(missing_report.missing_branches[0].nu_B_depart == 0.2);
        assert(missing_report.missing_branches[0].theta_A == 0.3);

        std::vector<problem1::Problem1SolutionBranch> extra_candidate{baseline[0], baseline[1]};
        problem1::Problem1SolutionBranch extra = baseline[0];
        extra.encounter_global_angle = 2.25;
        extra.target_arrival_true_anomaly = 2.75;
        extra.time_of_flight_seconds = 3000.0;
        extra.target_time_seconds = 3000.0;
        extra_candidate.push_back(extra);
        const auto extra_report = problem1::compare_problem1_solution_branches_by_identity(
            extra_candidate, baseline, 0.1, 0.2, 0.3);
        assert(extra_report.matched_count == 2);
        assert(extra_report.missing_count == 0);
        assert(extra_report.extra_count == 1);
        assert(extra_report.extra_branches.size() == 1);
        assert(angle_diff(extra_report.extra_branches[0].encounter_global_angle, extra.encounter_global_angle) <= 1e-12);
    }

    const std::vector<PairCase> pairs{
        {planet::PlanetId::Earth, planet::PlanetId::Mars},
        {planet::PlanetId::Mars, planet::PlanetId::Earth},
        {planet::PlanetId::Venus, planet::PlanetId::Mercury},
    };

    int query_count = 0;
    int mismatch_count = 0;
    int missing_from_adaptive_total = 0;
    int new_in_adaptive_total = 0;
    double max_alpha_diff = 0.0;
    double max_theta_prime_diff = 0.0;
    double max_target_time_diff = 0.0;
    double max_residual_diff = 0.0;
    long long adaptive_coarse_scan_count_total = 0;
    long long adaptive_interval_total = 0;
    long long adaptive_sign_change_interval_total = 0;
    long long adaptive_near_zero_interval_total = 0;
    long long adaptive_local_min_interval_total = 0;
    long long adaptive_local_max_interval_total = 0;
    long long adaptive_valid_boundary_interval_total = 0;
    long long adaptive_wrap_interval_total = 0;
    long long adaptive_rapid_change_interval_total = 0;
    long long adaptive_residual_eval_count_total = 0;
    long long adaptive_refined_interval_total = 0;
    long long adaptive_bisection_interval_total = 0;
    long long adaptive_ternary_interval_total = 0;
    long long adaptive_local_fine_scan_interval_total = 0;
    long long adaptive_candidate_count_before_dedup_total = 0;
    long long adaptive_candidate_count_after_dedup_total = 0;
    long long adaptive_ternary_accept_total = 0;
    long long adaptive_ternary_reject_total = 0;
    long long adaptive_local_fine_scan_root_total = 0;
    long long fullscan_residual_eval_count_total = 0;
    int adaptive_fallback_to_fullscan_count = 0;

    for (const auto& pair : pairs) {
        int pair_mismatch_count = 0;
        int pair_branch_count = 0;
        for (int i = 0; i < sample_count; ++i) {
            const double nu_A = angle_dist(rng);
            const double nu_B = angle_dist(rng);
            const double theta_A = angle_dist(rng);

            const auto old_branches = problem1::solve_problem1_from_departure_anomalies(
                pair.from, pair.to, nu_A, nu_B, theta_A, 1, 1);

            problem1::Problem1SolveOptions options{};
            options.mode = problem1::Problem1SolveMode::AdaptiveScanDebugCompare;
            options.max_transfer_revolution = 1;
            options.max_target_revolution = 1;
            const auto mode_result = problem1::solve_problem1_from_departure_anomalies_with_mode(
                pair.from, pair.to, nu_A, nu_B, theta_A, options);

            const bool ok = branches_equal(old_branches, mode_result.branches) &&
                mode_result.mode_profile.mode_requested == problem1::Problem1SolveMode::AdaptiveScanDebugCompare &&
                mode_result.mode_profile.mode_used == problem1::Problem1SolveMode::FullScan2880 &&
                mode_result.mode_profile.adaptive_attempted &&
                !mode_result.mode_profile.fallback_used &&
                mode_result.mode_profile.debug_compare_enabled &&
                mode_result.mode_profile.adaptive_coarse_scan_count > 0 &&
                mode_result.mode_profile.adaptive_residual_eval_count ==
                    mode_result.mode_profile.adaptive_coarse_scan_count &&
                mode_result.mode_profile.adaptive_interval_total >= 0 &&
                mode_result.mode_profile.adaptive_refined_interval_count >=
                    mode_result.mode_profile.adaptive_interval_total &&
                mode_result.mode_profile.adaptive_candidate_count_after_dedup >= 0;
            if (!ok) {
                mismatch_count += 1;
                pair_mismatch_count += 1;
                std::cout << "Problem1SolveModeDebugComparePlaceholderMismatch\n";
                std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
                std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
                std::cout << "sample_index=" << i << '\n';
                std::cout << "old_branch_count=" << old_branches.size() << '\n';
                std::cout << "new_branch_count=" << mode_result.branches.size() << '\n';
                std::cout << "missing_from_adaptive=" << mode_result.mode_profile.debug_compare_missing_from_adaptive << '\n';
                std::cout << "new_in_adaptive=" << mode_result.mode_profile.debug_compare_new_in_adaptive << '\n';
            }
            if (!mode_result.mode_profile.debug_compare_missing_branches.empty()) {
                const auto& missing = mode_result.mode_profile.debug_compare_missing_branches.front();
                std::cout << "Problem1SolveModeDebugCompareFirstMissing\n";
                std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
                std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
                std::cout << "sample_index=" << i << '\n';
                std::cout << "nu_A=" << missing.nu_A_depart << '\n';
                std::cout << "nu_B=" << missing.nu_B_depart << '\n';
                std::cout << "theta_A=" << missing.theta_A << '\n';
                std::cout << "transfer_revolution=" << missing.transfer_revolution << '\n';
                std::cout << "target_revolution=" << missing.target_revolution << '\n';
                std::cout << "encounter_global_angle=" << missing.encounter_global_angle << '\n';
                std::cout << "time_of_flight_seconds=" << missing.time_of_flight_seconds << '\n';
                std::cout << "residual_seconds=" << missing.residual_seconds << '\n';
                std::cout << "relative_residual=" << missing.relative_residual << '\n';
            }
            if (!mode_result.mode_profile.debug_compare_extra_diagnostics.empty()) {
                const auto& extra = mode_result.mode_profile.debug_compare_extra_diagnostics.front();
                std::cout << "Problem1SolveModeDebugCompareFirstExtra\n";
                std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
                std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
                std::cout << "sample_index=" << i << '\n';
                std::cout << "nu_A=" << extra.branch.nu_A_depart << '\n';
                std::cout << "nu_B=" << extra.branch.nu_B_depart << '\n';
                std::cout << "theta_A=" << extra.branch.theta_A << '\n';
                std::cout << "transfer_revolution=" << extra.branch.transfer_revolution << '\n';
                std::cout << "target_revolution=" << extra.branch.target_revolution << '\n';
                std::cout << "encounter_global_angle=" << extra.branch.encounter_global_angle << '\n';
                std::cout << "target_arrival_true_anomaly=" << extra.branch.target_arrival_true_anomaly << '\n';
                std::cout << "time_of_flight_seconds=" << extra.branch.time_of_flight_seconds << '\n';
                std::cout << "residual_seconds=" << extra.branch.residual_seconds << '\n';
                std::cout << "relative_residual=" << extra.branch.relative_residual << '\n';
                std::cout << "root_bracket_width=" << extra.root_bracket_width << '\n';
                std::cout << "bisection_iterations=" << extra.bisection_iterations << '\n';
                std::cout << "adaptive_source_reason=" << extra.adaptive_source_reason << '\n';
                std::cout << "ternary_phi_star=" << extra.ternary_phi_star << '\n';
                std::cout << "ternary_residual=" << extra.ternary_residual << '\n';
                std::cout << "ternary_residual_sq=" << extra.ternary_residual_sq << '\n';
                std::cout << "local_interval_left_phi=" << extra.local_interval_left_phi << '\n';
                std::cout << "local_interval_right_phi=" << extra.local_interval_right_phi << '\n';
                std::cout << "local_subdivision_index=" << extra.local_subdivision_index << '\n';
                std::cout << "sweep_min_abs_residual=" << extra.sweep_min_abs_residual << '\n';
                std::cout << "sweep_phi_at_min_abs_residual=" << extra.sweep_phi_at_min_abs_residual << '\n';
                std::cout << "sweep_sign_change_exists=" << (extra.sweep_sign_change_exists ? 1 : 0) << '\n';
                std::cout << "sweep_valid_sample_count=" << extra.sweep_valid_sample_count << '\n';
                std::cout << "sweep_invalid_sample_count=" << extra.sweep_invalid_sample_count << '\n';
                std::cout << "sweep_residual_at_adaptive_phi=" << extra.sweep_residual_at_adaptive_phi << '\n';
                std::cout << "sweep_left_residual=" << extra.sweep_left_residual << '\n';
                std::cout << "sweep_right_residual=" << extra.sweep_right_residual << '\n';
                std::cout << "nearest_baseline_transfer_revolution=" << extra.nearest_baseline_transfer_revolution << '\n';
                std::cout << "nearest_baseline_target_revolution=" << extra.nearest_baseline_target_revolution << '\n';
                std::cout << "nearest_baseline_angle_diff=" << extra.nearest_baseline_angle_diff << '\n';
                std::cout << "nearest_baseline_target_anomaly_diff=" << extra.nearest_baseline_target_anomaly_diff << '\n';
                std::cout << "nearest_baseline_time_diff=" << extra.nearest_baseline_time_diff << '\n';
                std::cout << "nearest_baseline_residual_diff=" << extra.nearest_baseline_residual_diff << '\n';
                std::cout << "nearest_baseline_match_failure_reason=" << extra.nearest_baseline_match_failure_reason << '\n';
                std::cout << "high_res_oracle_enabled=" << (extra.high_res_oracle_enabled ? 1 : 0) << '\n';
                std::cout << "high_res_scan_count=" << extra.high_res_scan_count << '\n';
                std::cout << "high_res_oracle_matched=" << (extra.high_res_oracle_matched ? 1 : 0) << '\n';
                std::cout << "high_res_oracle_branch_count=" << extra.high_res_oracle_branch_count << '\n';
            }

            query_count += 1;
            pair_branch_count += static_cast<int>(old_branches.size());
            missing_from_adaptive_total += mode_result.mode_profile.debug_compare_missing_from_adaptive;
            new_in_adaptive_total += mode_result.mode_profile.debug_compare_new_in_adaptive;
            max_alpha_diff = std::max(max_alpha_diff, mode_result.mode_profile.debug_compare_max_alpha_diff);
            max_theta_prime_diff = std::max(max_theta_prime_diff, mode_result.mode_profile.debug_compare_max_theta_prime_diff);
            max_target_time_diff = std::max(max_target_time_diff, mode_result.mode_profile.debug_compare_max_target_time_diff);
            max_residual_diff = std::max(max_residual_diff, mode_result.mode_profile.debug_compare_max_residual_diff);
            adaptive_coarse_scan_count_total += mode_result.mode_profile.adaptive_coarse_scan_count;
            adaptive_interval_total += mode_result.mode_profile.adaptive_interval_total;
            adaptive_sign_change_interval_total += mode_result.mode_profile.adaptive_sign_change_interval_count;
            adaptive_near_zero_interval_total += mode_result.mode_profile.adaptive_near_zero_interval_count;
            adaptive_local_min_interval_total += mode_result.mode_profile.adaptive_local_min_interval_count;
            adaptive_local_max_interval_total += mode_result.mode_profile.adaptive_local_max_interval_count;
            adaptive_valid_boundary_interval_total += mode_result.mode_profile.adaptive_valid_boundary_interval_count;
            adaptive_wrap_interval_total += mode_result.mode_profile.adaptive_wrap_interval_count;
            adaptive_rapid_change_interval_total += mode_result.mode_profile.adaptive_rapid_change_interval_count;
            adaptive_residual_eval_count_total += mode_result.mode_profile.adaptive_residual_eval_count;
            adaptive_refined_interval_total += mode_result.mode_profile.adaptive_refined_interval_count;
            adaptive_bisection_interval_total += mode_result.mode_profile.adaptive_bisection_interval_count;
            adaptive_ternary_interval_total += mode_result.mode_profile.adaptive_ternary_interval_count;
            adaptive_local_fine_scan_interval_total += mode_result.mode_profile.adaptive_local_fine_scan_interval_count;
            adaptive_candidate_count_before_dedup_total += mode_result.mode_profile.adaptive_candidate_count_before_dedup;
            adaptive_candidate_count_after_dedup_total += mode_result.mode_profile.adaptive_candidate_count_after_dedup;
            adaptive_ternary_accept_total += mode_result.mode_profile.adaptive_ternary_accept_count;
            adaptive_ternary_reject_total += mode_result.mode_profile.adaptive_ternary_reject_count;
            adaptive_local_fine_scan_root_total += mode_result.mode_profile.adaptive_local_fine_scan_root_count;
            fullscan_residual_eval_count_total += mode_result.mode_profile.residual_evaluation_count;
            if (mode_result.mode_profile.adaptive_fallback_to_fullscan) {
                adaptive_fallback_to_fullscan_count += 1;
            }
        }

        std::cout << "Problem1SolveModeDebugComparePlaceholderPairSummary\n";
        std::cout << "from_planet=" << planet::planet_name(pair.from) << '\n';
        std::cout << "to_planet=" << planet::planet_name(pair.to) << '\n';
        std::cout << "query_count=" << sample_count << '\n';
        std::cout << "branch_count=" << pair_branch_count << '\n';
        std::cout << "mismatch_count=" << pair_mismatch_count << '\n';
    }

    std::cout << "Problem1SolveModeDebugComparePlaceholderSummary\n";
    std::cout << "query_count=" << query_count << '\n';
    std::cout << "mismatch_count=" << mismatch_count << '\n';
    std::cout << "missing_from_adaptive_total=" << missing_from_adaptive_total << '\n';
    std::cout << "new_in_adaptive_total=" << new_in_adaptive_total << '\n';
    std::cout << "max_alpha_diff=" << max_alpha_diff << '\n';
    std::cout << "max_theta_prime_diff=" << max_theta_prime_diff << '\n';
    std::cout << "max_target_time_diff=" << max_target_time_diff << '\n';
    std::cout << "max_residual_diff=" << max_residual_diff << '\n';
    std::cout << "adaptive_coarse_scan_count_total=" << adaptive_coarse_scan_count_total << '\n';
    std::cout << "adaptive_interval_total=" << adaptive_interval_total << '\n';
    std::cout << "adaptive_sign_change_interval_total=" << adaptive_sign_change_interval_total << '\n';
    std::cout << "adaptive_near_zero_interval_total=" << adaptive_near_zero_interval_total << '\n';
    std::cout << "adaptive_local_min_interval_total=" << adaptive_local_min_interval_total << '\n';
    std::cout << "adaptive_local_max_interval_total=" << adaptive_local_max_interval_total << '\n';
    std::cout << "adaptive_valid_boundary_interval_total=" << adaptive_valid_boundary_interval_total << '\n';
    std::cout << "adaptive_wrap_interval_total=" << adaptive_wrap_interval_total << '\n';
    std::cout << "adaptive_rapid_change_interval_total=" << adaptive_rapid_change_interval_total << '\n';
    std::cout << "adaptive_residual_eval_count_total=" << adaptive_residual_eval_count_total << '\n';
    std::cout << "adaptive_refined_interval_total=" << adaptive_refined_interval_total << '\n';
    std::cout << "adaptive_bisection_interval_total=" << adaptive_bisection_interval_total << '\n';
    std::cout << "adaptive_ternary_interval_total=" << adaptive_ternary_interval_total << '\n';
    std::cout << "adaptive_local_fine_scan_interval_total=" << adaptive_local_fine_scan_interval_total << '\n';
    std::cout << "adaptive_candidate_count_before_dedup_total=" << adaptive_candidate_count_before_dedup_total << '\n';
    std::cout << "adaptive_candidate_count_after_dedup_total=" << adaptive_candidate_count_after_dedup_total << '\n';
    std::cout << "adaptive_ternary_accept_total=" << adaptive_ternary_accept_total << '\n';
    std::cout << "adaptive_ternary_reject_total=" << adaptive_ternary_reject_total << '\n';
    std::cout << "adaptive_local_fine_scan_root_total=" << adaptive_local_fine_scan_root_total << '\n';
    std::cout << "adaptive_fallback_to_fullscan_count=" << adaptive_fallback_to_fullscan_count << '\n';
    std::cout << "fullscan_residual_eval_count_total=" << fullscan_residual_eval_count_total << '\n';
    std::cout << "placeholder_ok=" << (mismatch_count == 0 ? 1 : 0) << '\n';
    return mismatch_count == 0 ? 0 : 1;
}
