/*
 * 文件作用：声明 Problem 1 的直接残差评估和求解接口。
 * 主要工作：定义输入、残差结果、候选解，并暴露扫描加二分细化的求解函数。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <vector>

namespace spaceship_cpp::problem1 {

enum class Problem1ResidualStatus {
    Success,
    InvalidInput,
    SingularGeometry,
    InvalidTransferOrbit,
    InvalidBranch,
    InvalidTimeOfFlight,
};

struct Problem1ResidualInput {
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;
    double launch_time_seconds_since_j2000;
    double transfer_perihelion_angle;
    double encounter_global_angle;
    int transfer_revolution;
    int target_revolution;
};

struct Problem1ResidualResult {
    Problem1ResidualStatus status;
    double residual;

    double r1;
    double r2;

    double xi1;
    double xi2;
    double target_theta_start;
    double target_theta_end;

    double transfer_e_raw;
    double transfer_e;
    double transfer_perihelion_angle_used;
    double transfer_p;

    double deltaF_transfer;
    double deltaF_target;

    double transfer_time_scale_free;
    double target_time_scale_free;
};

struct Problem1SolveInput {
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;

    double launch_time_seconds_since_j2000;
    double transfer_perihelion_angle;

    int max_transfer_revolution;
    int max_target_revolution;

    int phi_scan_count;

    double phi_tolerance;
    double residual_tolerance;
    int max_bisection_iterations;
    double max_candidate_relative_residual = 1e-6;
};

struct Problem1Candidate {
    Problem1ResidualResult residual_result;

    double encounter_global_angle;
    double launch_time_seconds_since_j2000;
    double time_of_flight_seconds;
    double arrival_time_seconds_since_j2000;

    double residual_scale;
    double relative_residual;

    double root_bracket_width;
    int bisection_iterations;
    bool refined_by_bisection;

    int transfer_revolution;
    int target_revolution;
};

Problem1ResidualResult evaluate_problem1_residual(const Problem1ResidualInput& input);

double compute_transfer_e_from_two_points(
    double r1,
    double xi1,
    double r2,
    double xi2
);

double compute_transfer_p_from_departure(
    double r1,
    double e_transfer,
    double xi1
);

std::vector<Problem1Candidate> solve_problem1(const Problem1SolveInput& input);

}  // namespace spaceship_cpp::problem1
