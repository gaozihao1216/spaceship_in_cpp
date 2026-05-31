#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"

#include "spaceship_cpp/common/common.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

namespace spaceship_cpp::problem1 {
namespace {

using common::kPi;
using common::normalize_angle_0_2pi;
using common::normalize_angle_minus_pi_pi;

using Clock = std::chrono::steady_clock;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = kPi / 90.0;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

int nearest_grid_index(double angle) {
    const double wrapped = normalize_angle_0_2pi(angle);
    int index = static_cast<int>(std::llround(wrapped / kGridStepRadians)) % kSamplesPerDimension;
    if (index < 0) {
        index += kSamplesPerDimension;
    }
    return index;
}

double grid_angle(int index) {
    return normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(index));
}

bool duplicate_branch(const Problem1SolutionBranch& lhs, const Problem1SolutionBranch& rhs) {
    return lhs.transfer_revolution == rhs.transfer_revolution &&
           std::abs(lhs.time_of_flight_seconds - rhs.time_of_flight_seconds) < 1e-3 &&
           std::abs(normalize_angle_minus_pi_pi(lhs.encounter_global_angle - rhs.encounter_global_angle)) < 1e-9;
}

int lower_grid_index(double angle) {
    const double wrapped = normalize_angle_0_2pi(angle);
    int index = static_cast<int>(std::floor(wrapped / kGridStepRadians)) % kSamplesPerDimension;
    if (index < 0) {
        index += kSamplesPerDimension;
    }
    return index;
}

double grid_distance_squared(
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int nu_A_index,
    int nu_B_index,
    int theta_A_index
) {
    const double d_a = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - grid_angle(nu_A_index)) /
        kGridStepRadians;
    const double d_b = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - grid_angle(nu_B_index)) /
        kGridStepRadians;
    const double d_t = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - grid_angle(theta_A_index)) /
        kGridStepRadians;
    return d_a * d_a + d_b * d_b + d_t * d_t;
}

void dedup_refined_branches(std::vector<Problem1SolutionBranch>* branches) {
    std::vector<Problem1SolutionBranch> deduped;
    for (const auto& branch : *branches) {
        auto it = std::find_if(deduped.begin(), deduped.end(), [&](const auto& existing) {
            return duplicate_branch(branch, existing);
        });
        if (it == deduped.end()) {
            deduped.push_back(branch);
        } else if (std::abs(branch.residual_seconds) < std::abs(it->residual_seconds)) {
            *it = branch;
        }
    }
    std::sort(deduped.begin(), deduped.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
    });
    *branches = std::move(deduped);
}

void fallback_direct_solve(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_transfer_revolution,
    int max_target_revolution,
    Problem1NearestNodeQueryResult* result
) {
    result->profile.fallback_direct_solve_count += 1;
    const auto direct_start = Clock::now();
    result->branches = solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        query_nu_A,
        query_nu_B,
        query_theta_A,
        max_transfer_revolution,
        max_target_revolution);
    result->profile.fallback_direct_solve_ms += elapsed_ms(direct_start, Clock::now());
    result->profile.fallback_direct_branch_count += static_cast<int>(result->branches.size());
    for (auto& branch : result->branches) {
        if (!branch.valid) {
            continue;
        }
        const auto attached = attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            branch,
            Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        if (attached.valid && attached.derivatives_available) {
            branch = attached;
        }
    }
    result->used_direct_solve_fallback = true;
    result->valid = !result->branches.empty();
    if (!result->valid) {
        result->invalid_reason = "direct_solve_returned_no_branches";
    } else {
        result->invalid_reason.clear();
    }
}

struct CellVertexSeedCandidate {
    bool valid = false;
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double distance_squared = 0.0;
    Problem1RootTable2DegLoadedNode node;
    int valid_branch_count = 0;
};

std::vector<CellVertexSeedCandidate> select_cell_vertex_seed_candidates(
    const Problem1RootTable2DegLoader& loader,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_transfer_revolution,
    int max_target_revolution,
    const Problem1NearestNodeQueryOptions& options,
    Problem1NearestNodeQueryProfile* profile
) {
    const auto start = Clock::now();
    const int lower_a = lower_grid_index(query_nu_A);
    const int lower_b = lower_grid_index(query_nu_B);
    const int lower_t = lower_grid_index(query_theta_A);
    std::vector<CellVertexSeedCandidate> candidates;
    candidates.reserve(8);
    for (int da = 0; da <= 1; ++da) {
        for (int db = 0; db <= 1; ++db) {
            for (int dt = 0; dt <= 1; ++dt) {
                CellVertexSeedCandidate candidate{};
                candidate.nu_A_index = (lower_a + da) % kSamplesPerDimension;
                candidate.nu_B_index = (lower_b + db) % kSamplesPerDimension;
                candidate.theta_A_index = (lower_t + dt) % kSamplesPerDimension;
                candidate.distance_squared = grid_distance_squared(
                    query_nu_A, query_nu_B, query_theta_A,
                    candidate.nu_A_index, candidate.nu_B_index, candidate.theta_A_index);
                const auto linear = Problem1RootTable2DegLoader::linear_index_from_indices(
                    candidate.nu_A_index, candidate.nu_B_index, candidate.theta_A_index);
                candidate.node = loader.load_node_by_linear_index(linear);
                profile->cell_vertex_loaded_vertex_count += 1;
                profile->table_node_read_ms += candidate.node.profile.node_read_ms;
                profile->table_offset_build_ms += candidate.node.profile.offset_build_ms;
                if (!candidate.node.valid) {
                    continue;
                }
                for (const auto& branch : candidate.node.branches) {
                    if (branch.valid &&
                        branch.derivatives_available &&
                        branch.transfer_revolution <= max_transfer_revolution &&
                        branch.target_revolution <= max_target_revolution) {
                        candidate.valid_branch_count += 1;
                    }
                }
                candidate.valid = candidate.valid_branch_count > 0;
                if (candidate.valid) {
                    candidates.push_back(std::move(candidate));
                }
            }
        }
    }
    if (candidates.empty()) {
        profile->cell_vertex_routeA_no_valid_vertex_count += 1;
        profile->cell_vertex_selection_ms += elapsed_ms(start, Clock::now());
        return candidates;
    }
    std::sort(candidates.begin(), candidates.end(), [&](const auto& lhs, const auto& rhs) {
        if (options.cell_vertex_seed_policy == 2) {
            return lhs.distance_squared < rhs.distance_squared;
        }
        if (options.cell_vertex_seed_policy == 3) {
            if (lhs.valid_branch_count != rhs.valid_branch_count) {
                return lhs.valid_branch_count > rhs.valid_branch_count;
            }
            return lhs.distance_squared < rhs.distance_squared;
        }
        if (lhs.valid_branch_count != rhs.valid_branch_count) {
            return lhs.valid_branch_count > rhs.valid_branch_count;
        }
        return lhs.distance_squared < rhs.distance_squared;
    });
    const int keep = std::max(0, std::min<int>(options.cell_vertex_top_k_vertices, candidates.size()));
    candidates.resize(static_cast<std::size_t>(keep));
    profile->cell_vertex_selected_vertex_count += keep;
    profile->cell_vertex_selection_ms += elapsed_ms(start, Clock::now());
    return candidates;
}

bool try_cell_vertex_routeA(
    const Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_transfer_revolution,
    int max_target_revolution,
    const Problem1NearestNodeQueryOptions& options,
    Problem1NearestNodeQueryResult* result
) {
    result->profile.cell_vertex_routeA_attempt_count += 1;
    const auto selected = select_cell_vertex_seed_candidates(
        loader, query_nu_A, query_nu_B, query_theta_A, max_transfer_revolution, max_target_revolution,
        options, &result->profile);
    if (selected.empty()) {
        result->profile.cell_vertex_routeA_failure_count += 1;
        return false;
    }

    int attempted = 0;
    const auto route_start = Clock::now();
    for (const auto& vertex : selected) {
        const double vertex_nu_A = grid_angle(vertex.nu_A_index);
        const double vertex_nu_B = grid_angle(vertex.nu_B_index);
        const double vertex_theta_A = grid_angle(vertex.theta_A_index);
        const double delta_nu_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_A) - vertex_nu_A);
        const double delta_nu_B = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_nu_B) - vertex_nu_B);
        const double delta_theta_A = normalize_angle_minus_pi_pi(normalize_angle_0_2pi(query_theta_A) - vertex_theta_A);
        for (const auto& branch : vertex.node.branches) {
            if (attempted >= options.cell_vertex_max_branch_attempts) {
                break;
            }
            if (!branch.valid ||
                !branch.derivatives_available ||
                branch.transfer_revolution > max_transfer_revolution ||
                branch.target_revolution > max_target_revolution) {
                continue;
            }
            attempted += 1;
            result->profile.cell_vertex_attempted_branch_count += 1;
            const double alpha_seed = normalize_angle_0_2pi(
                branch.encounter_global_angle +
                branch.d_encounter_global_angle_d_nu_A * delta_nu_A +
                branch.d_encounter_global_angle_d_nu_B * delta_nu_B +
                branch.d_encounter_global_angle_d_theta_A * delta_theta_A);
            const auto refined = refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
                departure_planet,
                target_planet,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                branch.transfer_revolution,
                branch.target_revolution,
                alpha_seed,
                options.max_newton_iterations,
                options.residual_tolerance_seconds,
                1e-12,
                Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
            const auto validation_start = Clock::now();
            const bool valid = refined.valid &&
                refined.branch.valid &&
                std::abs(refined.branch.residual_seconds) <= options.residual_tolerance_seconds &&
                std::isfinite(refined.branch.time_of_flight_seconds) &&
                refined.branch.time_of_flight_seconds > 0.0 &&
                std::isfinite(refined.branch.encounter_global_angle);
            result->profile.cell_vertex_validation_ms += elapsed_ms(validation_start, Clock::now());
            if (!valid) {
                result->profile.cell_vertex_routeA_validation_failure_count += 1;
                continue;
            }
            auto attached = attach_problem1_root_derivatives_with_mode(
                departure_planet,
                target_planet,
                query_nu_A,
                query_nu_B,
                query_theta_A,
                refined.branch,
                Problem1RootDerivativeMode::AnalyticOnly,
                1e-6);
            if (!(attached.valid && attached.derivatives_available)) {
                attached = refined.branch;
            }
            result->branches.push_back(attached);
        }
        if (attempted >= options.cell_vertex_max_branch_attempts) {
            break;
        }
    }
    result->profile.cell_vertex_routeA_ms += elapsed_ms(route_start, Clock::now());
    if (attempted == 0) {
        result->profile.cell_vertex_routeA_no_valid_branch_count += 1;
    }
    dedup_refined_branches(&result->branches);
    if (result->branches.empty()) {
        result->profile.cell_vertex_routeA_failure_count += 1;
        return false;
    }
    result->profile.cell_vertex_routeA_success_count += 1;
    result->profile.direct_fallback_avoided_count += 1;
    result->valid = true;
    result->invalid_reason.clear();
    return true;
}

void handle_query_fallback(
    const Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_transfer_revolution,
    int max_target_revolution,
    const Problem1NearestNodeQueryOptions& options,
    Problem1NearestNodeQueryResult* result
) {
    const bool try_cell = options.cell_vertex_routeA_enabled ||
        options.fallback_mode == Problem1FallbackMode::CellVertexRouteA ||
        options.fallback_mode == Problem1FallbackMode::ScoutNoDirectFallback ||
        options.scout_no_direct_fallback;
    if (try_cell &&
        options.cell_vertex_seed_policy != 0 &&
        try_cell_vertex_routeA(
            loader, departure_planet, target_planet, query_nu_A, query_nu_B, query_theta_A,
            max_transfer_revolution, max_target_revolution, options, result)) {
        return;
    }

    const bool scout_no_direct =
        options.scout_no_direct_fallback || options.fallback_mode == Problem1FallbackMode::ScoutNoDirectFallback;
    const bool allow_direct = options.fallback_direct_solve &&
        !scout_no_direct &&
        (options.fallback_mode == Problem1FallbackMode::DirectSolve ||
         options.fallback_mode == Problem1FallbackMode::CellVertexRouteA);
    if (allow_direct) {
        fallback_direct_solve(
            departure_planet,
            target_planet,
            query_nu_A,
            query_nu_B,
            query_theta_A,
            max_transfer_revolution,
            max_target_revolution,
            result);
    } else {
        result->profile.direct_fallback_skipped_count += 1;
        result->valid = false;
    }
}

}  // namespace

Problem1NearestNodeQueryResult query_problem1_from_2deg_nearest_node(
    const Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double query_nu_A,
    double query_nu_B,
    double query_theta_A,
    int max_transfer_revolution,
    int max_target_revolution,
    const Problem1NearestNodeQueryOptions& options
) {
    const auto total_start = Clock::now();
    Problem1NearestNodeQueryResult result{};
    result.nearest_nu_A_index = nearest_grid_index(query_nu_A);
    result.nearest_nu_B_index = nearest_grid_index(query_nu_B);
    result.nearest_theta_A_index = nearest_grid_index(query_theta_A);
    result.nearest_linear_index = Problem1RootTable2DegLoader::linear_index_from_indices(
        result.nearest_nu_A_index,
        result.nearest_nu_B_index,
        result.nearest_theta_A_index);

    const double node_nu_A = grid_angle(result.nearest_nu_A_index);
    const double node_nu_B = grid_angle(result.nearest_nu_B_index);
    const double node_theta_A = grid_angle(result.nearest_theta_A_index);
    const double wrapped_query_nu_A = normalize_angle_0_2pi(query_nu_A);
    const double wrapped_query_nu_B = normalize_angle_0_2pi(query_nu_B);
    const double wrapped_query_theta_A = normalize_angle_0_2pi(query_theta_A);
    result.delta_nu_A = normalize_angle_minus_pi_pi(wrapped_query_nu_A - node_nu_A);
    result.delta_nu_B = normalize_angle_minus_pi_pi(wrapped_query_nu_B - node_nu_B);
    result.delta_theta_A = normalize_angle_minus_pi_pi(wrapped_query_theta_A - node_theta_A);

    result.profile.table_load_node_count += 1;
    const auto table_load_start = Clock::now();
    const auto loaded = loader.load_node_by_linear_index(result.nearest_linear_index);
    result.profile.table_load_ms += elapsed_ms(table_load_start, Clock::now());
    if (loaded.profile.offset_built_this_call) {
        result.profile.table_offset_build_count += 1;
    }
    result.profile.table_offset_build_ms += loaded.profile.offset_build_ms;
    result.profile.table_node_read_ms += loaded.profile.node_read_ms;
    if (!loaded.valid) {
        result.invalid_reason = "nearest_node_load_failed:" + loaded.invalid_reason;
        handle_query_fallback(
            loader,
            departure_planet,
            target_planet,
            wrapped_query_nu_A,
            wrapped_query_nu_B,
            wrapped_query_theta_A,
            max_transfer_revolution,
            max_target_revolution,
            options,
            &result);
        result.profile.total_ms = elapsed_ms(total_start, Clock::now());
        return result;
    }
    result.profile.table_loaded_branch_count += static_cast<int>(loaded.branches.size());

    for (const auto& branch : loaded.branches) {
        if (!branch.valid ||
            !branch.derivatives_available ||
            branch.transfer_revolution > max_transfer_revolution ||
            branch.target_revolution > max_target_revolution) {
            continue;
        }
        result.table_seed_count += 1;
        result.profile.seed_branch_count += 1;
        const auto seed_start = Clock::now();
        const double alpha_seed = normalize_angle_0_2pi(
            branch.encounter_global_angle +
            branch.d_encounter_global_angle_d_nu_A * result.delta_nu_A +
            branch.d_encounter_global_angle_d_nu_B * result.delta_nu_B +
            branch.d_encounter_global_angle_d_theta_A * result.delta_theta_A);
        result.profile.seed_generation_ms += elapsed_ms(seed_start, Clock::now());
        result.profile.refine_attempt_count += 1;
        const auto refine_start = Clock::now();
        const auto refined = refine_problem1_root_branch_newton_residual_first_diagnostic_seconds(
            departure_planet,
            target_planet,
            wrapped_query_nu_A,
            wrapped_query_nu_B,
            wrapped_query_theta_A,
            branch.transfer_revolution,
            branch.target_revolution,
            alpha_seed,
            options.max_newton_iterations,
            options.residual_tolerance_seconds,
            1e-12,
            Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        result.profile.refine_ms += elapsed_ms(refine_start, Clock::now());
        if (!refined.valid ||
            !refined.branch.valid ||
            std::abs(refined.branch.residual_seconds) > options.residual_tolerance_seconds) {
            result.refined_fail_count += 1;
            result.profile.refine_fail_count += 1;
            continue;
        }
        result.profile.refine_success_count += 1;
        result.profile.derivative_attach_attempt_count += 1;
        const auto derivative_attach_start = Clock::now();
        auto attached = attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            wrapped_query_nu_A,
            wrapped_query_nu_B,
            wrapped_query_theta_A,
            refined.branch,
            Problem1RootDerivativeMode::AnalyticOnly,
            1e-6);
        result.profile.derivative_attach_ms += elapsed_ms(derivative_attach_start, Clock::now());
        if (!(attached.valid && attached.derivatives_available)) {
            result.profile.derivative_attach_fail_count += 1;
            attached = refined.branch;
        } else {
            result.profile.derivative_attach_success_count += 1;
        }
        result.branches.push_back(attached);
        result.refined_success_count += 1;
    }

    result.profile.dedup_input_count += static_cast<int>(result.branches.size());
    const auto dedup_start = Clock::now();
    dedup_refined_branches(&result.branches);
    result.profile.dedup_ms += elapsed_ms(dedup_start, Clock::now());
    result.profile.dedup_output_count += static_cast<int>(result.branches.size());
    if (result.branches.empty()) {
        result.invalid_reason = "nearest_node_refine_produced_no_branches";
        handle_query_fallback(
            loader,
            departure_planet,
            target_planet,
            wrapped_query_nu_A,
            wrapped_query_nu_B,
            wrapped_query_theta_A,
            max_transfer_revolution,
            max_target_revolution,
            options,
            &result);
        result.profile.total_ms = elapsed_ms(total_start, Clock::now());
        return result;
    }

    result.valid = true;
    result.invalid_reason.clear();
    result.profile.total_ms = elapsed_ms(total_start, Clock::now());
    return result;
}

}  // namespace spaceship_cpp::problem1
