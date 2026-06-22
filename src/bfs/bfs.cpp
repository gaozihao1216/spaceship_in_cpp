/*
 * 文件作用：BFS 轨迹搜索公共 API 实现。
 * 主要工作：封装 Steps 1–4，按给定发射时刻与起终点行星返回最佳转移轨道描述。
 */
#include "spaceship_cpp/bfs/bfs.hpp"

#include "spaceship_cpp/bfs/step3_top_k_sequences.hpp"
#include "spaceship_cpp/bfs/step4_fixed_sequence_fine_search.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace spaceship_cpp::bfs {
namespace {

using planet_params::PlanetId;

TransferLegDescriptor leg_from_edge(const TrajectorySearchEdge& edge) {
    TransferLegDescriptor leg{};
    if (!edge.valid) {
        return leg;
    }
    leg.from_planet = edge.from_planet;
    leg.to_planet = edge.to_planet;
    leg.departure_time_seconds_since_j2000 = edge.departure_time;
    leg.arrival_time_seconds_since_j2000 = edge.arrival_time;
    leg.time_of_flight_seconds = edge.transfer_time_seconds;
    leg.eccentricity = edge.outgoing_e;
    leg.semi_latus_rectum_au = edge.outgoing_p;
    leg.perihelion_angle_global_rad = edge.outgoing_theta;
    leg.transfer_revolution = edge.transfer_revolution;
    leg.target_revolution = edge.target_revolution;
    leg.flyby_theta_prime_local_rad = edge.theta_prime;
    leg.encounter_angle_rad = edge.alpha;
    leg.flyby_constraint_G = edge.slingshot_residual;
    return leg;
}

std::vector<TransferLegDescriptor> legs_from_edges(
    const std::vector<TrajectorySearchEdge>& edges
) {
    std::vector<TransferLegDescriptor> legs;
    legs.reserve(edges.size());
    for (const TrajectorySearchEdge& edge : edges) {
        if (edge.valid) {
            legs.push_back(leg_from_edge(edge));
        }
    }
    return legs;
}

bool sequence_ends_at(
    const std::vector<PlanetId>& sequence,
    PlanetId destination
) {
    return !sequence.empty() && sequence.back() == destination;
}

TrajectorySearchOutput make_output_from_step4_entry(
    const Step4SequenceFineSearchResult& refined,
    PlanetId departure,
    PlanetId destination
) {
    TrajectorySearchOutput output{};
    output.ok = true;
    output.departure_planet = departure;
    output.destination_planet = destination;
    output.visit_sequence = refined.planet_sequence;
    output.leg0_theta_global_rad = refined.best_theta_global;

    if (!refined.found_solution) {
        output.found_solution = false;
        output.error_message = refined.error_message.empty()
            ? "no_feasible_trajectory_for_destination"
            : refined.error_message;
        return output;
    }

    const FixedSequenceBfsResult& bfs = refined.best_bfs;
    output.found_solution = true;
    output.score = refined.best_score;
    output.launch_v_inf_mps = bfs.best_launch_v_inf;
    output.arrival_v_inf_mps = bfs.best_arrival_v_inf;
    output.total_time_seconds = bfs.best_total_time_seconds;
    output.legs = legs_from_edges(bfs.best_edges);
    return output;
}

TrajectorySearchOutput make_output_from_step3_entry(
    const RankedPlanetSequence& ranked,
    PlanetId departure,
    PlanetId destination
) {
    TrajectorySearchOutput output{};
    output.ok = true;
    output.found_solution = true;
    output.departure_planet = departure;
    output.destination_planet = destination;
    output.visit_sequence = ranked.planet_sequence;
    output.score = ranked.best_score;
    output.launch_v_inf_mps = ranked.best_solution.launch_v_inf;
    output.arrival_v_inf_mps = ranked.best_solution.arrival_v_inf;
    output.total_time_seconds = ranked.best_solution.total_time_seconds;
    output.leg0_theta_global_rad = ranked.best_seed.transfer_theta_global;
    output.legs = legs_from_edges(ranked.best_solution.edges);
    return output;
}

TrajectorySearchGlobalConfig config_for_mission(
    const TrajectorySearchInput& input,
    const TrajectorySearchGlobalConfig& base
) {
    TrajectorySearchGlobalConfig config = base;
    config.mission.launch_time_seconds_since_j2000 = input.launch_time_seconds_since_j2000;
    return config;
}

bool validate_input(
    const TrajectorySearchInput& input,
    std::string& error_message
) {
    if (input.departure_planet != PlanetId::Earth) {
        error_message = "departure_planet_must_be_earth";
        return false;
    }
    if (input.destination_planet == PlanetId::Earth) {
        error_message = "destination_planet_cannot_be_earth";
        return false;
    }
    return true;
}

}  // namespace

TrajectorySearchOutput search_best_trajectory(
    const TrajectorySearchInput& input,
    const TrajectorySearchGlobalConfig& config
) {
    TrajectorySearchOutput output{};
    output.departure_planet = input.departure_planet;
    output.destination_planet = input.destination_planet;

    std::string error_message;
    if (!validate_input(input, error_message)) {
        output.error_message = std::move(error_message);
        return output;
    }

    const TrajectorySearchGlobalConfig mission_config = config_for_mission(input, config);

    const Leg0MultiTargetThetaResult leg0 =
        find_leg0_feasible_theta_for_first_leg_targets(mission_config);
    if (!leg0.ok) {
        output.error_message = leg0.error_message.empty()
            ? "step1_leg0_scan_failed"
            : leg0.error_message;
        return output;
    }

    const std::vector<Leg0ThetaSeed> seeds =
        discretize_leg0_theta_seeds(mission_config, leg0);

    const Step2FreePathSearchResult step2 =
        run_step2_free_path_search(mission_config, seeds);
    if (!step2.ok) {
        output.error_message = step2.error_message.empty()
            ? "step2_free_path_bfs_failed"
            : step2.error_message;
        return output;
    }

    const Step3TopKSequencesResult step3 = select_top_k_planet_sequences_reaching(
        step2,
        mission_config.discretization.top_k_sequences,
        input.destination_planet);
    if (!step3.ok) {
        output.error_message = step3.error_message.empty()
            ? "step3_top_k_failed"
            : step3.error_message;
        return output;
    }

    output.ok = true;

    if (step3.sequences.empty()) {
        output.found_solution = false;
        output.error_message = "no_candidate_sequences_from_step2";
        return output;
    }

    const Step4FineSearchResult step4 =
        run_step4_fixed_sequence_fine_search(mission_config, step3, leg0);

    const PlanetId destination = input.destination_planet;

    const Step4SequenceFineSearchResult* best_refined = nullptr;
    double best_score = std::numeric_limits<double>::infinity();
    for (const Step4SequenceFineSearchResult& refined : step4.by_sequence) {
        if (!sequence_ends_at(refined.planet_sequence, destination)) {
            continue;
        }
        if (!refined.found_solution) {
            continue;
        }
        if (refined.best_score < best_score) {
            best_score = refined.best_score;
            best_refined = &refined;
        }
    }

    if (best_refined != nullptr) {
        return make_output_from_step4_entry(*best_refined, input.departure_planet, destination);
    }

    const RankedPlanetSequence* best_coarse = nullptr;
    best_score = std::numeric_limits<double>::infinity();
    for (const RankedPlanetSequence& ranked : step3.sequences) {
        if (!sequence_ends_at(ranked.planet_sequence, destination)) {
            continue;
        }
        if (ranked.best_score < best_score) {
            best_score = ranked.best_score;
            best_coarse = &ranked;
        }
    }

    if (best_coarse != nullptr) {
        TrajectorySearchOutput coarse =
            make_output_from_step3_entry(*best_coarse, input.departure_planet, destination);
        coarse.error_message =
            "fine_search_found_no_solution; returning_step3_coarse_best";
        return coarse;
    }

    output.found_solution = false;
    output.error_message = "no_trajectory_reaching_destination_planet";
    return output;
}

}  // namespace spaceship_cpp::bfs
