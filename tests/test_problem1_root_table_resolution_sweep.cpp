#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace {

bool approx_equal(double lhs, double rhs, double abs_tol = 1e-5, double rel_tol = 1e-5) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= abs_tol + rel_tol * scale;
}

std::optional<spaceship_cpp::problem1::Problem1SolutionBranch> find_same_kq_nearest_branch(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    int transfer_revolution,
    int target_revolution,
    double encounter_global_angle
) {
    namespace common = spaceship_cpp::common;
    const spaceship_cpp::problem1::Problem1SolutionBranch* best = nullptr;
    double best_distance = std::numeric_limits<double>::infinity();
    for (const auto& branch : branches) {
        if (!branch.valid ||
            branch.transfer_revolution != transfer_revolution ||
            branch.target_revolution != target_revolution) {
            continue;
        }
        const double distance = std::abs(
            common::normalize_angle_minus_pi_pi(branch.encounter_global_angle - encounter_global_angle));
        if (distance < best_distance) {
            best_distance = distance;
            best = &branch;
        }
    }
    if (best == nullptr) {
        return std::nullopt;
    }
    return *best;
}

double average_or_zero(double sum, int count) {
    if (count <= 0) {
        return 0.0;
    }
    return sum / static_cast<double>(count);
}

void update_max(double value, double* current_max) {
    assert(current_max != nullptr);
    if (!std::isfinite(value)) {
        return;
    }
    if (value > *current_max) {
        *current_max = value;
    }
}

struct ResolutionCase {
    std::string name;
    int nu_A_count = 0;
    int nu_B_depart_count = 0;
    int theta_A_count = 0;
};

struct QuerySample {
    spaceship_cpp::planet_params::PlanetId departure_planet = spaceship_cpp::planet_params::PlanetId::Earth;
    spaceship_cpp::planet_params::PlanetId target_planet = spaceship_cpp::planet_params::PlanetId::Mars;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
};

struct SweepSummary {
    spaceship_cpp::planet_params::PlanetId departure_planet = spaceship_cpp::planet_params::PlanetId::Earth;
    spaceship_cpp::planet_params::PlanetId target_planet = spaceship_cpp::planet_params::PlanetId::Mars;
    ResolutionCase resolution{};
    int total_nodes = 0;
    double nu_A_step = 0.0;
    double nu_B_depart_step = 0.0;
    double theta_A_step = 0.0;
    int max_transfer_revolution = 0;
    int max_target_revolution = 0;

    double build_time_seconds = 0.0;
    double query_time_seconds_total = 0.0;
    double exact_solve_time_seconds_total = 0.0;

    int table_nonempty_cell_count = 0;
    int table_empty_cell_count = 0;
    int table_total_solution_count = 0;
    double average_solution_count_per_nonempty_cell = 0.0;

    int sample_count = 0;
    int exact_solve_branch_count = 0;
    int route_a_candidate_count = 0;
    int route_a_matched_branch_count = 0;
    int route_a_missing_branch_count = 0;
    double route_a_coverage_ratio = 0.0;
    int matched_error_count = 0;
    double route_a_max_alpha_error = 0.0;
    double route_a_avg_alpha_error = 0.0;
    double route_a_max_time_error = 0.0;
    double route_a_avg_time_error = 0.0;
    double route_a_max_residual = 0.0;
    double route_a_avg_residual = 0.0;
};

spaceship_cpp::problem1::Problem1RootTableConfig make_root_table_config(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const ResolutionCase& resolution
);

std::vector<QuerySample> build_deterministic_query_samples() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::vector<QuerySample> samples;
    const std::vector<std::pair<planet_params::PlanetId, planet_params::PlanetId>> planet_pairs{
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Venus},
    };
    for (const auto& [departure_planet, target_planet] : planet_pairs) {
        const ResolutionCase reference_resolution{"8x8x12_reference", 8, 8, 12};
        const problem1::Problem1RootTableConfig config =
            make_root_table_config(departure_planet, target_planet, reference_resolution);
        const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
        const double offset_fraction = 0.005;
        int pair_sample_count = 0;
        for (const auto& cell : table.cells()) {
            if (!cell.solved || cell.solutions_sorted_by_time_of_flight.empty()) {
                continue;
            }
            QuerySample sample{};
            sample.departure_planet = departure_planet;
            sample.target_planet = target_planet;
            sample.query_nu_A = common::normalize_angle_0_2pi(
                cell.nu_A_depart + offset_fraction * config.nu_A_step);
            sample.query_nu_B = common::normalize_angle_0_2pi(
                cell.nu_B_depart + offset_fraction * config.nu_B_depart_step);
            sample.query_theta_A = common::normalize_angle_0_2pi(
                cell.theta_A + offset_fraction * config.theta_A_step);
            samples.push_back(sample);
            pair_sample_count += 1;
            if (pair_sample_count >= 6) {
                break;
            }
        }
    }
    return samples;
}

spaceship_cpp::problem1::Problem1RootTableConfig make_root_table_config(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const ResolutionCase& resolution
) {
    namespace common = spaceship_cpp::common;
    return spaceship_cpp::problem1::Problem1RootTableConfig{
        departure_planet,
        target_planet,
        0.0,
        common::kTwoPi / static_cast<double>(resolution.nu_A_count),
        resolution.nu_A_count,
        0.0,
        common::kTwoPi / static_cast<double>(resolution.nu_B_depart_count),
        resolution.nu_B_depart_count,
        0.0,
        common::kTwoPi / static_cast<double>(resolution.theta_A_count),
        resolution.theta_A_count,
        1,
        1,
        "problem1_root_table_draft_v0",
    };
}

SweepSummary run_resolution_case(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const ResolutionCase& resolution,
    const std::vector<QuerySample>& all_samples
) {
    namespace problem1 = spaceship_cpp::problem1;
    using clock = std::chrono::steady_clock;
    const double newton_residual_tolerance_seconds = 1e-6;
    const double newton_residual_tolerance_scale_free =
        problem1::problem1_residual_seconds_to_scale_free(newton_residual_tolerance_seconds);

    const problem1::Problem1RootTableConfig config =
        make_root_table_config(departure_planet, target_planet, resolution);
    SweepSummary summary{};
    summary.resolution = resolution;
    summary.departure_planet = departure_planet;
    summary.target_planet = target_planet;
    summary.nu_A_step = config.nu_A_step;
    summary.nu_B_depart_step = config.nu_B_depart_step;
    summary.theta_A_step = config.theta_A_step;
    summary.total_nodes = config.nu_A_count * config.nu_B_depart_count * config.theta_A_count;
    summary.max_transfer_revolution = config.max_transfer_revolution;
    summary.max_target_revolution = config.max_target_revolution;

    const auto build_start = clock::now();
    const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
    const auto build_end = clock::now();
    summary.build_time_seconds = std::chrono::duration<double>(build_end - build_start).count();

    for (const auto& cell : table.cells()) {
        const int branch_count = static_cast<int>(cell.solutions_sorted_by_time_of_flight.size());
        if (branch_count > 0) {
            summary.table_nonempty_cell_count += 1;
            summary.table_total_solution_count += branch_count;
        } else {
            summary.table_empty_cell_count += 1;
        }
    }
    summary.average_solution_count_per_nonempty_cell = average_or_zero(
        static_cast<double>(summary.table_total_solution_count),
        summary.table_nonempty_cell_count);

    double alpha_error_sum = 0.0;
    double time_error_sum = 0.0;
    double residual_sum = 0.0;

    for (const QuerySample& sample : all_samples) {
        if (sample.departure_planet != departure_planet || sample.target_planet != target_planet) {
            continue;
        }
        summary.sample_count += 1;

        const auto exact_start = clock::now();
        const auto exact_branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            config.max_transfer_revolution,
            config.max_target_revolution);
        const auto exact_end = clock::now();
        summary.exact_solve_time_seconds_total += std::chrono::duration<double>(exact_end - exact_start).count();

        summary.exact_solve_branch_count += static_cast<int>(exact_branches.size());
        if (exact_branches.empty()) {
            continue;
        }

        const auto query_start = clock::now();
        const auto route_a = problem1::query_problem1_root_table_linear_newton(
            table,
            sample.query_nu_A,
            sample.query_nu_B,
            sample.query_theta_A,
            30,
            newton_residual_tolerance_scale_free,
            1e-12);
        const auto query_end = clock::now();
        summary.query_time_seconds_total += std::chrono::duration<double>(query_end - query_start).count();

        assert(route_a.valid);
        summary.route_a_candidate_count += static_cast<int>(route_a.branches.size());

        for (const auto& exact_branch : exact_branches) {
            const auto candidate = find_same_kq_nearest_branch(
                route_a.branches,
                exact_branch.transfer_revolution,
                exact_branch.target_revolution,
                exact_branch.encounter_global_angle);
            if (!candidate.has_value()) {
                summary.route_a_missing_branch_count += 1;
                continue;
            }
            const double alpha_error = std::abs(
                spaceship_cpp::common::normalize_angle_minus_pi_pi(
                    candidate->encounter_global_angle - exact_branch.encounter_global_angle));
            const double time_error =
                std::abs(candidate->time_of_flight_seconds - exact_branch.time_of_flight_seconds);
            const double residual_abs = std::abs(candidate->residual_seconds);

            summary.route_a_matched_branch_count += 1;
            summary.matched_error_count += 1;
            alpha_error_sum += alpha_error;
            time_error_sum += time_error;
            residual_sum += residual_abs;
            update_max(alpha_error, &summary.route_a_max_alpha_error);
            update_max(time_error, &summary.route_a_max_time_error);
            update_max(residual_abs, &summary.route_a_max_residual);
        }
    }

    summary.route_a_coverage_ratio = average_or_zero(
        static_cast<double>(summary.route_a_matched_branch_count),
        summary.exact_solve_branch_count);
    summary.route_a_avg_alpha_error = average_or_zero(alpha_error_sum, summary.matched_error_count);
    summary.route_a_avg_time_error = average_or_zero(time_error_sum, summary.matched_error_count);
    summary.route_a_avg_residual = average_or_zero(residual_sum, summary.matched_error_count);
    return summary;
}

void print_summary(const SweepSummary& summary) {
    std::cout << "planet_pair=" <<
        spaceship_cpp::planet_params::planet_name(summary.departure_planet) << "->" <<
        spaceship_cpp::planet_params::planet_name(summary.target_planet) << '\n';
    std::cout << "resolution_name=" << summary.resolution.name << '\n';
    std::cout << "nu_A_count=" << summary.resolution.nu_A_count << '\n';
    std::cout << "nu_A_step=" << summary.nu_A_step << '\n';
    std::cout << "nu_B_depart_count=" << summary.resolution.nu_B_depart_count << '\n';
    std::cout << "nu_B_depart_step=" << summary.nu_B_depart_step << '\n';
    std::cout << "theta_A_count=" << summary.resolution.theta_A_count << '\n';
    std::cout << "theta_A_step=" << summary.theta_A_step << '\n';
    std::cout << "total_nodes=" << summary.total_nodes << '\n';
    std::cout << "max_transfer_revolution=" << summary.max_transfer_revolution << '\n';
    std::cout << "max_target_revolution=" << summary.max_target_revolution << '\n';
    std::cout << "build_time_seconds=" << summary.build_time_seconds << '\n';
    std::cout << "query_time_seconds_total=" << summary.query_time_seconds_total << '\n';
    std::cout << "exact_solve_time_seconds_total=" << summary.exact_solve_time_seconds_total << '\n';
    std::cout << "table_nonempty_cell_count=" << summary.table_nonempty_cell_count << '\n';
    std::cout << "table_empty_cell_count=" << summary.table_empty_cell_count << '\n';
    std::cout << "table_total_solution_count=" << summary.table_total_solution_count << '\n';
    std::cout << "average_solution_count_per_nonempty_cell=" <<
        summary.average_solution_count_per_nonempty_cell << '\n';
    std::cout << "sample_count=" << summary.sample_count << '\n';
    std::cout << "exact_solve_branch_count=" << summary.exact_solve_branch_count << '\n';
    std::cout << "route_a_candidate_count=" << summary.route_a_candidate_count << '\n';
    std::cout << "route_a_matched_branch_count=" << summary.route_a_matched_branch_count << '\n';
    std::cout << "route_a_missing_branch_count=" << summary.route_a_missing_branch_count << '\n';
    std::cout << "route_a_coverage_ratio=" << summary.route_a_coverage_ratio << '\n';
    std::cout << "route_a_max_alpha_error=" << summary.route_a_max_alpha_error << '\n';
    std::cout << "route_a_avg_alpha_error=" << summary.route_a_avg_alpha_error << '\n';
    std::cout << "route_a_max_time_error=" << summary.route_a_max_time_error << '\n';
    std::cout << "route_a_avg_time_error=" << summary.route_a_avg_time_error << '\n';
    std::cout << "route_a_max_residual=" << summary.route_a_max_residual << '\n';
    std::cout << "route_a_avg_residual=" << summary.route_a_avg_residual << '\n';
    std::cout << "average_newton_iterations=unavailable\n";
    std::cout << "max_newton_iterations=unavailable\n";
    std::cout << std::flush;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;

    const std::vector<QuerySample> all_samples = build_deterministic_query_samples();
    assert(!all_samples.empty());

    const std::vector<ResolutionCase> lightweight_cases{
        {"8x8x12_smoke", 8, 8, 12},
        {"12x12x24_light", 12, 12, 24},
    };
    const std::vector<ResolutionCase> expensive_cases{
        {"16x16x32_expensive", 16, 16, 32},
        {"24x24x48_expensive", 24, 24, 48},
    };
    const bool run_expensive = [] {
        const char* env = std::getenv("RUN_EXPENSIVE_ROOT_TABLE_SWEEP");
        return env != nullptr && std::string(env) == "1";
    }();

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1RootTable resolution sweep\n";
    std::cout << "Current root-table query smoke/functionality tests use 8x8x12.\n";
    std::cout << "This is a smoke-test resolution only, not a production-quality root table.\n";
    if (!run_expensive) {
        std::cout <<
            "RUN_EXPENSIVE_ROOT_TABLE_SWEEP is not set to 1; running lightweight Earth->Mars sweep only.\n";
    }

    std::vector<SweepSummary> summaries;
    for (const auto& resolution : lightweight_cases) {
        std::cout << "\n=== Earth->Mars ===\n";
        SweepSummary mars_summary = run_resolution_case(
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            resolution,
            all_samples);
        print_summary(mars_summary);
        summaries.push_back(mars_summary);
    }

    if (run_expensive) {
        for (const auto& resolution : lightweight_cases) {
            std::cout << "\n=== Earth->Venus ===\n";
            summaries.push_back(run_resolution_case(
                planet_params::PlanetId::Earth,
                planet_params::PlanetId::Venus,
                resolution,
                all_samples));
            print_summary(summaries.back());
        }
        for (const auto& resolution : expensive_cases) {
            std::cout << "\n=== Earth->Mars ===\n";
            summaries.push_back(run_resolution_case(
                planet_params::PlanetId::Earth,
                planet_params::PlanetId::Mars,
                resolution,
                all_samples));
            print_summary(summaries.back());

            std::cout << "\n=== Earth->Venus ===\n";
            summaries.push_back(run_resolution_case(
                planet_params::PlanetId::Earth,
                planet_params::PlanetId::Venus,
                resolution,
                all_samples));
            print_summary(summaries.back());
        }
        std::cout << "TODO: 32x32x64 sweep is intentionally skipped in default expensive mode to avoid excessive runtime.\n";
    }

    bool any_matched = false;
    for (const auto& summary : summaries) {
        any_matched = any_matched || (summary.route_a_matched_branch_count > 0);
    }
    assert(any_matched);

    const auto check_monotonicish = [&](planet_params::PlanetId target_planet) {
        const SweepSummary* coarse = nullptr;
        const SweepSummary* finer = nullptr;
        for (const auto& summary : summaries) {
            if (summary.resolution.name == "8x8x12_smoke" &&
                summary.target_planet == target_planet) {
                coarse = &summary;
                break;
            }
        }
        for (const auto& summary : summaries) {
            if (summary.resolution.name == "12x12x24_light" &&
                summary.target_planet == target_planet) {
                finer = &summary;
                break;
            }
        }
        assert(coarse != nullptr);
        assert(finer != nullptr);
        if (coarse->route_a_coverage_ratio > 0.0) {
            assert(finer->route_a_coverage_ratio >= coarse->route_a_coverage_ratio * 0.8);
        }
        if (finer->route_a_coverage_ratio > coarse->route_a_coverage_ratio + 0.05 &&
            finer->route_a_missing_branch_count + 2 < coarse->route_a_missing_branch_count) {
            std::cout << "Likely table resolution is a major factor.\n";
        } else {
            std::cout <<
                "Coverage is still low; likely due to nearest-node seed limitation or branch matching ambiguity.\n";
        }
    };

    check_monotonicish(planet_params::PlanetId::Mars);
    if (run_expensive) {
        check_monotonicish(planet_params::PlanetId::Venus);
    }
    std::cout << std::flush;
    for (const auto& summary : summaries) {
        if (summary.route_a_matched_branch_count > 0) {
            assert(summary.route_a_max_residual <= 1e-6);
            assert(std::isfinite(summary.route_a_max_alpha_error));
            assert(std::isfinite(summary.route_a_max_time_error));
        }
    }
    std::cout << std::flush;
    return 0;
}
