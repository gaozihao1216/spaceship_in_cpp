/*
 * 文件作用：实现 Step 3 行星序列 Top-K 筛选。
 */
#include "spaceship_cpp/bfs/step3_top_k_sequences.hpp"

#include <algorithm>
#include <map>
#include <utility>

namespace spaceship_cpp::bfs {
namespace {

using planet_params::PlanetId;

bool validate_top_k(int top_k, std::string& error_message) {
    if (!(top_k > 0)) {
        error_message = "top_k_must_be_positive";
        return false;
    }
    return true;
}

bool validate_step2_input(
    const Step2FreePathSearchResult& step2,
    std::string& error_message
) {
    if (!step2.ok) {
        error_message = step2.error_message.empty()
            ? "step2_result_not_ok"
            : step2.error_message;
        return false;
    }
    if (step2.by_seed.empty()) {
        error_message = "step2_has_no_seed_results";
        return false;
    }
    return true;
}

struct SequenceBestEntry {
    double best_score = 0.0;
    Leg0ThetaSeed best_seed{};
    FreePathBfsSolution best_solution{};
};

void consider_solution(
    const Leg0ThetaSeed& seed,
    const FreePathBfsSolution& solution,
    std::map<std::vector<PlanetId>, SequenceBestEntry>& best_by_sequence
) {
    if (solution.planet_sequence.empty()) {
        return;
    }

    const auto found = best_by_sequence.find(solution.planet_sequence);
    if (found == best_by_sequence.end() || solution.score < found->second.best_score) {
        best_by_sequence[solution.planet_sequence] = SequenceBestEntry{
            .best_score = solution.score,
            .best_seed = seed,
            .best_solution = solution,
        };
    }
}

std::vector<RankedPlanetSequence> build_ranked_sequences(
    std::map<std::vector<PlanetId>, SequenceBestEntry>&& best_by_sequence
) {
    std::vector<RankedPlanetSequence> ranked;
    ranked.reserve(best_by_sequence.size());
    for (auto& [planet_sequence, entry] : best_by_sequence) {
        ranked.push_back(RankedPlanetSequence{
            .planet_sequence = std::move(planet_sequence),
            .best_score = entry.best_score,
            .best_seed = entry.best_seed,
            .best_solution = std::move(entry.best_solution),
        });
    }

    std::sort(
        ranked.begin(),
        ranked.end(),
        [](const RankedPlanetSequence& lhs, const RankedPlanetSequence& rhs) {
            if (lhs.best_score != rhs.best_score) {
                return lhs.best_score < rhs.best_score;
            }
            return lhs.planet_sequence < rhs.planet_sequence;
        });
    return ranked;
}

}  // namespace

Step3TopKSequencesResult select_top_k_planet_sequences(
    const Step2FreePathSearchResult& step2,
    int top_k
) {
    Step3TopKSequencesResult result{};
    result.stats.top_k_requested = top_k;

    std::string error_message;
    if (!validate_top_k(top_k, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }
    if (!validate_step2_input(step2, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }

    std::map<std::vector<PlanetId>, SequenceBestEntry> best_by_sequence;
    for (const FreePathBfsResult& seed_result : step2.by_seed) {
        for (const FreePathBfsSolution& solution : seed_result.solutions) {
            ++result.stats.solutions_pooled;
            consider_solution(seed_result.seed, solution, best_by_sequence);
        }
    }

    result.stats.unique_sequences = best_by_sequence.size();
    std::vector<RankedPlanetSequence> ranked =
        build_ranked_sequences(std::move(best_by_sequence));

    const std::size_t keep_count = std::min(
        ranked.size(),
        static_cast<std::size_t>(top_k));
    result.sequences.assign(ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(keep_count));
    result.stats.sequences_returned = result.sequences.size();
    result.ok = true;
    if (result.sequences.empty()) {
        result.error_message = "no_planet_sequences_from_step2";
    }
    return result;
}

Step3TopKSequencesResult select_top_k_planet_sequences_reaching(
    const Step2FreePathSearchResult& step2,
    int top_k,
    PlanetId destination_planet
) {
    Step3TopKSequencesResult result{};
    result.stats.top_k_requested = top_k;

    std::string error_message;
    if (!validate_top_k(top_k, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }
    if (!validate_step2_input(step2, error_message)) {
        result.error_message = std::move(error_message);
        return result;
    }

    std::map<std::vector<PlanetId>, SequenceBestEntry> best_by_sequence;
    for (const FreePathBfsResult& seed_result : step2.by_seed) {
        for (const FreePathBfsSolution& solution : seed_result.solutions) {
            if (solution.planet_sequence.empty() ||
                solution.planet_sequence.back() != destination_planet) {
                continue;
            }
            ++result.stats.solutions_pooled;
            consider_solution(seed_result.seed, solution, best_by_sequence);
        }
    }

    result.stats.unique_sequences = best_by_sequence.size();
    std::vector<RankedPlanetSequence> ranked =
        build_ranked_sequences(std::move(best_by_sequence));

    const std::size_t keep_count = std::min(
        ranked.size(),
        static_cast<std::size_t>(top_k));
    result.sequences.assign(ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(keep_count));
    result.stats.sequences_returned = result.sequences.size();
    result.ok = true;
    if (result.sequences.empty()) {
        result.error_message = "no_planet_sequences_reaching_destination";
    }
    return result;
}

Step3TopKSequencesResult run_step3_select_top_k_sequences(
    const TrajectorySearchGlobalConfig& config,
    const Step2FreePathSearchResult& step2
) {
    return select_top_k_planet_sequences(step2, config.discretization.top_k_sequences);
}

}  // namespace spaceship_cpp::bfs
