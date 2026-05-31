#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct ThetaLaunchFeasibilityScoutOptions {
    int theta_scout_count = 720;

    double max_launch_v_inf = 7000.0;
    double near_v_inf_buffer = 1000.0;

    int max_transfer_revolution = 1;
    int max_target_revolution = 1;

    double problem1_residual_tolerance_seconds = 1e-2;
    int problem1_max_newton_iterations = 80;

    bool refine_interval_boundaries = true;
    int max_boundary_refine_iterations = 30;
    double theta_boundary_tolerance = 1e-5;
};

struct ThetaLaunchFeasibilitySample {
    bool valid = false;

    double theta = 0.0;

    int raw_branch_count = 0;
    int accepted_candidate_count = 0;

    double min_launch_v_inf = std::numeric_limits<double>::infinity();
    double second_min_launch_v_inf = std::numeric_limits<double>::infinity();

    planet_params::PlanetId best_target_planet = planet_params::PlanetId::Mercury;
    int best_transfer_revolution = 0;
    int best_target_revolution = 0;

    double best_arrival_time = 0.0;
    double best_transfer_time_seconds = 0.0;
};

struct ThetaLaunchFeasibilityInterval {
    bool valid = false;

    double theta_left = 0.0;
    double theta_right = 0.0;

    double min_v_inf_inside = std::numeric_limits<double>::infinity();
    double theta_at_min = 0.0;

    int sample_count_inside = 0;
};

struct ThetaLaunchFeasibilityScoutResult {
    bool ok = false;
    std::string error_message;

    std::vector<ThetaLaunchFeasibilitySample> samples;

    std::vector<ThetaLaunchFeasibilityInterval> valid_intervals;
    std::vector<ThetaLaunchFeasibilityInterval> near_valid_intervals;

    ThetaLaunchFeasibilitySample best_sample;

    int theta_scout_count = 0;
    int valid_sample_count = 0;
    int near_valid_sample_count = 0;
    int no_candidate_sample_count = 0;

    double global_min_launch_v_inf = std::numeric_limits<double>::infinity();
    double theta_at_global_min = 0.0;
};

ThetaLaunchFeasibilityScoutResult scout_theta_launch_feasibility_with_problem1_table(
    const problem1::Problem1RootTable2DegLoader& loader,
    planet_params::PlanetId launch_planet,
    double launch_time,
    const std::vector<planet_params::PlanetId>& allowed_first_targets,
    const ThetaLaunchFeasibilityScoutOptions& options
);

}  // namespace spaceship_cpp::bfs
