#pragma once

#include "spaceship_cpp/bfs/theta_launch_feasibility_scout.hpp"

#include <limits>
#include <string>
#include <vector>

namespace spaceship_cpp::bfs {

struct ThetaCandidateSelectorOptions {
    int max_theta_candidates = 50;

    int min_samples_per_valid_interval = 5;
    int max_samples_per_valid_interval = 30;

    bool include_interval_boundaries = true;
    bool include_interval_midpoint = true;
    bool include_interval_minimum = true;

    bool include_near_valid_intervals = false;
    int max_near_valid_extra_candidates = 10;

    double duplicate_theta_tolerance = 1e-8;
};

struct ThetaCandidate {
    bool valid = false;

    double theta = 0.0;
    double estimated_min_launch_v_inf = std::numeric_limits<double>::infinity();

    bool from_valid_interval = false;
    bool from_near_valid_interval = false;
    bool is_interval_boundary = false;
    bool is_interval_midpoint = false;
    bool is_interval_minimum = false;

    int interval_index = -1;
};

struct ThetaCandidateSelectionResult {
    bool ok = false;
    std::string error_message;

    std::vector<ThetaCandidate> candidates;

    int valid_interval_count = 0;
    int near_valid_interval_count = 0;
    int candidate_count = 0;
};

ThetaCandidateSelectionResult select_theta_candidates_from_launch_scout(
    const ThetaLaunchFeasibilityScoutResult& scout,
    const ThetaCandidateSelectorOptions& options
);

}  // namespace spaceship_cpp::bfs
