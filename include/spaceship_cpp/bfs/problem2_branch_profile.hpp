#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <limits>

namespace spaceship_cpp::bfs {

enum class Problem2BranchRejectReason {
    None,
    InvalidEdge,
    InvalidNextState,
    NonFiniteArrivalTime,
    NonFiniteIncomingOrbit,
    NonIncreasingTime,
    NonPositiveTransferTime,
    NonPositiveOutgoingP,
    SlingshotResidualTooLarge,
    Problem1ResidualTooLarge,
    ArrivalTimeTooLarge,
    FlybyEvaluatorInvalid,
    FlybyPhysicalInfeasible
};

struct Problem2SolveBranchStats {
    bool valid = false;

    int depth = 0;
    int parent_node_index = -1;

    planet_params::PlanetId from_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId to_planet = planet_params::PlanetId::Mercury;

    double parent_time = 0.0;
    double parent_incoming_e = 0.0;
    double parent_incoming_theta = 0.0;

    int raw_edge_count = 0;
    int accepted_edge_count = 0;

    int reject_invalid_edge_count = 0;
    int reject_invalid_next_state_count = 0;
    int reject_nonfinite_arrival_time_count = 0;
    int reject_nonfinite_incoming_orbit_count = 0;
    int reject_non_increasing_time_count = 0;
    int reject_non_positive_transfer_time_count = 0;
    int reject_non_positive_outgoing_p_count = 0;
    int reject_slingshot_residual_count = 0;
    int reject_problem1_residual_count = 0;
    int reject_arrival_time_count = 0;
    int reject_flyby_evaluator_invalid_count = 0;
    int reject_flyby_physical_infeasible_count = 0;

    int flyby_evaluated_count = 0;
    int flyby_valid_count = 0;
    int flyby_feasible_count = 0;
    int flyby_infeasible_count = 0;
    int flyby_invalid_count = 0;

    int would_reject_flyby_evaluator_invalid_count = 0;
    int would_reject_flyby_physical_infeasible_count = 0;

    double min_v_inf_in = std::numeric_limits<double>::infinity();
    double max_v_inf_in = 0.0;
    double min_v_inf_out = std::numeric_limits<double>::infinity();
    double max_v_inf_out = 0.0;
    double min_v_inf_mismatch = std::numeric_limits<double>::infinity();
    double max_v_inf_mismatch = 0.0;
    double min_v_inf_mismatch_tolerance = std::numeric_limits<double>::infinity();
    double max_v_inf_mismatch_tolerance = 0.0;

    double min_arrival_time = std::numeric_limits<double>::infinity();
    double max_arrival_time = -std::numeric_limits<double>::infinity();
    double min_transfer_time_seconds = std::numeric_limits<double>::infinity();
    double max_transfer_time_seconds = -std::numeric_limits<double>::infinity();

    double min_abs_slingshot_residual = std::numeric_limits<double>::infinity();
    double max_abs_slingshot_residual = 0.0;
    double min_abs_problem1_residual_seconds = std::numeric_limits<double>::infinity();
    double max_abs_problem1_residual_seconds = 0.0;

    double min_flyby_required_periapsis_radius_m = std::numeric_limits<double>::infinity();
    double max_flyby_required_periapsis_radius_m = 0.0;
    double min_flyby_allowed_periapsis_radius_m = std::numeric_limits<double>::infinity();
    double max_flyby_allowed_periapsis_radius_m = 0.0;
    double min_flyby_turn_angle_rad = std::numeric_limits<double>::infinity();
    double max_flyby_turn_angle_rad = 0.0;
    double min_flyby_max_turn_angle_rad = std::numeric_limits<double>::infinity();
    double max_flyby_max_turn_angle_rad = 0.0;

    int rejected_by_vinf_mismatch_count = 0;
    int rejected_by_turn_angle_count = 0;
    int rejected_by_periapsis_radius_count = 0;

    double solve_wall_time_ms = 0.0;
};

struct Problem2LayerBranchStats {
    int depth = 0;

    planet_params::PlanetId from_planet = planet_params::PlanetId::Mercury;
    planet_params::PlanetId to_planet = planet_params::PlanetId::Mercury;

    int input_state_count = 0;
    int problem2_solve_count = 0;

    int raw_edge_count_total = 0;
    int accepted_edge_count_total = 0;
    int output_state_count_after_beam = 0;
    int beam_pruned_count = 0;

    int raw_edge_count_min = 0;
    int raw_edge_count_max = 0;
    double raw_edge_count_mean = 0.0;

    int accepted_edge_count_min = 0;
    int accepted_edge_count_max = 0;
    double accepted_edge_count_mean = 0.0;

    double layer_wall_time_ms = 0.0;
};

}  // namespace spaceship_cpp::bfs
