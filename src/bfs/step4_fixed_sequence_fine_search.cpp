/*
 * 文件作用：实现 Step 4 固定序列 θ 细扫 + BFS 精搜。
 */
#include "spaceship_cpp/bfs/step4_fixed_sequence_fine_search.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace spaceship_cpp::bfs {
namespace {

using planet_params::PlanetId;

bool near_equal_theta(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-12;
}

bool validate_planet_sequence(
    const std::vector<PlanetId>& planet_sequence,
    std::string& error_message
) {
    if (planet_sequence.size() < 2U) {
        error_message = "planet_sequence_must_have_at_least_two_planets";
        return false;
    }
    if (planet_sequence.front() != PlanetId::Earth) {
        error_message = "planet_sequence_must_start_at_earth";
        return false;
    }
    return true;
}

bool validate_step3_input(
    const Step3TopKSequencesResult& step3,
    std::string& error_message
) {
    if (!step3.ok) {
        error_message = step3.error_message.empty()
            ? "step3_result_not_ok"
            : step3.error_message;
        return false;
    }
    if (step3.sequences.empty()) {
        error_message = "step3_has_no_sequences";
        return false;
    }
    return true;
}

bool validate_leg0_input(
    const Leg0MultiTargetThetaResult& leg0,
    std::string& error_message
) {
    if (!leg0.ok) {
        error_message = leg0.error_message.empty()
            ? "leg0_result_not_ok"
            : leg0.error_message;
        return false;
    }
    return true;
}

void append_unique_theta(std::vector<double>& samples, double theta) {
    for (double existing : samples) {
        if (near_equal_theta(existing, theta)) {
            return;
        }
    }
    samples.push_back(theta);
}

std::vector<double> build_theta_samples_for_sequence(
    const TrajectorySearchGlobalConfig& config,
    const Leg0MultiTargetThetaResult& leg0,
    const RankedPlanetSequence& ranked_sequence
) {
    if (ranked_sequence.planet_sequence.size() < 2U) {
        return {};
    }

    const PlanetId first_leg_target = ranked_sequence.planet_sequence[1];
    std::vector<double> samples = discretize_leg0_theta_samples_for_target(
        leg0,
        first_leg_target,
        config.discretization.leg0_theta_fine_scan_count);
    append_unique_theta(samples, ranked_sequence.best_seed.transfer_theta_global);
    return samples;
}

void update_global_best(Step4FineSearchResult& result) {
    result.best_sequence_index = -1;
    result.global_best_score = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < result.by_sequence.size(); ++index) {
        const Step4SequenceFineSearchResult& entry = result.by_sequence[index];
        if (!entry.found_solution) {
            continue;
        }
        if (entry.best_score < result.global_best_score) {
            result.global_best_score = entry.best_score;
            result.best_sequence_index = static_cast<int>(index);
        }
    }
}

}  // namespace

Step4SequenceFineSearchResult search_fixed_sequence_over_theta_samples(
    const TrajectorySearchGlobalConfig& config,
    const std::vector<PlanetId>& planet_sequence,
    const std::vector<double>& theta_samples
) {
    Step4SequenceFineSearchResult result{};
    result.planet_sequence = planet_sequence;

    std::string error_message;
    if (!validate_planet_sequence(planet_sequence, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }
    if (theta_samples.empty()) {
        result.error_message = "theta_samples_must_not_be_empty";
        return result;
    }
    if (!(config.discretization.leg0_theta_fine_scan_count > 0)) {
        result.error_message = "leg0_theta_fine_scan_count_must_be_positive";
        return result;
    }

    for (double theta : theta_samples) {
        ++result.theta_samples_tried;
        const FixedSequenceBfsConfig bfs_config =
            make_fixed_sequence_bfs_config(config, planet_sequence, theta);
        const FixedSequenceBfsResult bfs_result = search_fixed_sequence_bfs(bfs_config);
        if (!bfs_result.ok) {
            result.error_message = bfs_result.error_message;
            return result;
        }
        if (!bfs_result.found_solution) {
            continue;
        }

        ++result.theta_samples_with_solution;
        if (bfs_result.best_score < result.best_score) {
            result.found_solution = true;
            result.best_score = bfs_result.best_score;
            result.best_theta_global = theta;
            result.best_bfs = bfs_result;
        }
    }

    if (!result.found_solution) {
        result.error_message = "no_feasible_fixed_sequence_over_theta_samples";
    }
    return result;
}

Step4FineSearchResult run_step4_fixed_sequence_fine_search(
    const TrajectorySearchGlobalConfig& config,
    const Step3TopKSequencesResult& step3,
    const Leg0MultiTargetThetaResult& leg0
) {
    Step4FineSearchResult result{};

    std::string error_message;
    if (!validate_step3_input(step3, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }
    if (!validate_leg0_input(leg0, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }
    if (!(config.discretization.leg0_theta_fine_scan_count > 0)) {
        result.error_message = "leg0_theta_fine_scan_count_must_be_positive";
        return result;
    }

    result.by_sequence.reserve(step3.sequences.size());
    for (const RankedPlanetSequence& ranked_sequence : step3.sequences) {
        ++result.stats.sequences_processed;
        Step4SequenceFineSearchResult sequence_result{};
        sequence_result.planet_sequence = ranked_sequence.planet_sequence;

        if (!validate_planet_sequence(ranked_sequence.planet_sequence, error_message)) {
            sequence_result.error_message = error_message;
            result.by_sequence.push_back(std::move(sequence_result));
            continue;
        }

        const std::vector<double> theta_samples =
            build_theta_samples_for_sequence(config, leg0, ranked_sequence);
        if (theta_samples.empty()) {
            sequence_result.error_message = "no_theta_samples_for_first_leg_target";
            result.by_sequence.push_back(std::move(sequence_result));
            continue;
        }

        result.stats.total_theta_samples += theta_samples.size();
        sequence_result = search_fixed_sequence_over_theta_samples(
            config,
            ranked_sequence.planet_sequence,
            theta_samples);
        result.stats.total_bfs_calls += sequence_result.theta_samples_tried;

        if (sequence_result.found_solution) {
            ++result.stats.sequences_with_solution;
        }
        result.by_sequence.push_back(std::move(sequence_result));
    }

    update_global_best(result);
    result.ok = true;
    if (result.stats.sequences_with_solution == 0U) {
        result.error_message = "step4_found_no_feasible_sequence";
    }
    return result;
}

Step4FineSearchResult run_step4_fixed_sequence_fine_search(
    const TrajectorySearchGlobalConfig& config,
    const Step3TopKSequencesResult& step3
) {
    const Leg0MultiTargetThetaResult leg0 =
        find_leg0_feasible_theta_for_first_leg_targets(config);
    Step4FineSearchResult result{};
    if (!leg0.ok) {
        result.error_message = leg0.error_message;
        return result;
    }
    return run_step4_fixed_sequence_fine_search(config, step3, leg0);
}

}  // namespace spaceship_cpp::bfs
