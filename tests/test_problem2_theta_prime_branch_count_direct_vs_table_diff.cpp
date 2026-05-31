#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kThetaSampleCount = 32;
constexpr int kPreviewLimit = 10;

struct ResidualBranch {
    bool problem1_valid = false;
    bool slingshot_valid = false;
    std::string problem1_invalid_reason;
    std::string slingshot_invalid_reason;
    double theta_prime = 0.0;
    double alpha = 0.0;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double problem1_residual_seconds = 0.0;
    double slingshot_residual = 0.0;
};

struct SampleSummary {
    double theta_prime = 0.0;
    std::vector<ResidualBranch> branches;
    std::map<int, int> raw_count_by_k;
    std::map<int, int> slingshot_valid_count_by_k;
    bool table_fallback_used = false;
};

struct PairDiff {
    double max_alpha_diff = 0.0;
    double max_time_diff_seconds = 0.0;
    double max_slingshot_residual_diff = 0.0;
};

struct SignChangeStats {
    int sign_change_pair_count = 0;
    int residual_invalid_pair_count = 0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

bool residual_branch_less(const ResidualBranch& lhs, const ResidualBranch& rhs) {
    return std::tie(lhs.transfer_revolution, lhs.target_revolution, lhs.time_of_flight_seconds) <
           std::tie(rhs.transfer_revolution, rhs.target_revolution, rhs.time_of_flight_seconds);
}

std::string format_count_by_k(const std::map<int, int>& count_by_k) {
    std::ostringstream os;
    os << "{";
    bool first = true;
    for (const auto& [k, count] : count_by_k) {
        if (!first) {
            os << ",";
        }
        first = false;
        os << k << ":" << count;
    }
    os << "}";
    return os.str();
}

int problem1_valid_count(const std::vector<ResidualBranch>& branches) {
    int count = 0;
    for (const auto& branch : branches) {
        if (branch.problem1_valid) {
            count += 1;
        }
    }
    return count;
}

int slingshot_valid_count(const std::vector<ResidualBranch>& branches) {
    int count = 0;
    for (const auto& branch : branches) {
        if (branch.slingshot_valid) {
            count += 1;
        }
    }
    return count;
}

std::vector<ResidualBranch> convert_problem1_to_problem2_residuals(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_anomaly_phi,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem2 = spaceship_cpp::problem2;
    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);

    std::vector<ResidualBranch> results;
    results.reserve(branches.size());
    for (const auto& branch : branches) {
        ResidualBranch residual_branch{};
        residual_branch.theta_prime = theta_prime;
        residual_branch.alpha = branch.target_arrival_true_anomaly;
        residual_branch.transfer_revolution = branch.transfer_revolution;
        residual_branch.target_revolution = branch.target_revolution;
        residual_branch.time_of_flight_seconds = branch.time_of_flight_seconds;
        residual_branch.target_time_seconds = branch.target_time_seconds;
        residual_branch.problem1_residual_seconds = branch.residual_seconds;
        if (!branch.valid) {
            residual_branch.problem1_invalid_reason =
                branch.invalid_reason.empty() ? "problem1_branch_invalid" : branch.invalid_reason;
            results.push_back(residual_branch);
            continue;
        }
        residual_branch.problem1_valid = true;
        const auto theta_alpha_residual = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            encounter_anomaly_phi,
            branch.target_arrival_true_anomaly,
            incoming_e,
            incoming_theta,
            theta_prime);
        if (!theta_alpha_residual.valid) {
            residual_branch.slingshot_invalid_reason =
                theta_alpha_residual.invalid_reason.empty() ? "slingshot_invalid" : theta_alpha_residual.invalid_reason;
            results.push_back(residual_branch);
            continue;
        }
        residual_branch.slingshot_valid = true;
        residual_branch.slingshot_residual = theta_alpha_residual.slingshot_residual;
        results.push_back(residual_branch);
    }
    std::sort(results.begin(), results.end(), residual_branch_less);
    return results;
}

SampleSummary make_sample_summary(double theta_prime, std::vector<ResidualBranch> branches, bool table_fallback_used) {
    SampleSummary summary{};
    summary.theta_prime = theta_prime;
    summary.branches = std::move(branches);
    summary.table_fallback_used = table_fallback_used;
    for (const auto& branch : summary.branches) {
        if (branch.problem1_valid) {
            summary.raw_count_by_k[branch.transfer_revolution] += 1;
        }
        if (branch.slingshot_valid) {
            summary.slingshot_valid_count_by_k[branch.transfer_revolution] += 1;
        }
    }
    return summary;
}

SampleSummary build_direct_summary(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet, phi, beta, theta_A, 1, 1);
    return make_sample_summary(
        theta_prime,
        convert_problem1_to_problem2_residuals(
            departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime, branches),
        false);
}

SampleSummary build_table_summary(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto query = problem1::query_problem1_from_2deg_nearest_node(
        loader, departure_planet, target_planet, phi, beta, theta_A, 1, 1, options);
    return make_sample_summary(
        theta_prime,
        convert_problem1_to_problem2_residuals(
            departure_planet, target_planet, phi, incoming_e, incoming_theta, theta_prime, query.branches),
        query.used_direct_solve_fallback);
}

PairDiff compute_pair_diff(const SampleSummary& direct, const SampleSummary& table) {
    PairDiff diff{};
    const int pair_count = std::min<int>(direct.branches.size(), table.branches.size());
    for (int i = 0; i < pair_count; ++i) {
        const auto& d = direct.branches[static_cast<std::size_t>(i)];
        const auto& t = table.branches[static_cast<std::size_t>(i)];
        diff.max_alpha_diff = std::max(diff.max_alpha_diff, std::abs(normalize_angle_minus_pi_pi(d.alpha - t.alpha)));
        diff.max_time_diff_seconds =
            std::max(diff.max_time_diff_seconds, std::abs(d.time_of_flight_seconds - t.time_of_flight_seconds));
        if (d.slingshot_valid && t.slingshot_valid) {
            diff.max_slingshot_residual_diff =
                std::max(diff.max_slingshot_residual_diff, std::abs(d.slingshot_residual - t.slingshot_residual));
        }
    }
    return diff;
}

std::map<int, std::vector<ResidualBranch>> raw_problem1_branches_by_k(const SampleSummary& sample) {
    std::map<int, std::vector<ResidualBranch>> by_k;
    for (const auto& branch : sample.branches) {
        if (branch.problem1_valid) {
            by_k[branch.transfer_revolution].push_back(branch);
        }
    }
    for (auto& [k, branches] : by_k) {
        (void)k;
        std::sort(branches.begin(), branches.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.time_of_flight_seconds < rhs.time_of_flight_seconds;
        });
    }
    return by_k;
}

SignChangeStats count_sign_change_pairs(const SampleSummary& left, const SampleSummary& right) {
    const auto left_by_k = raw_problem1_branches_by_k(left);
    const auto right_by_k = raw_problem1_branches_by_k(right);
    SignChangeStats stats{};
    for (const auto& [k, left_branches] : left_by_k) {
        const auto right_it = right_by_k.find(k);
        if (right_it == right_by_k.end()) {
            continue;
        }
        const auto& right_branches = right_it->second;
        const int pair_count = std::min<int>(left_branches.size(), right_branches.size());
        for (int rank = 0; rank < pair_count; ++rank) {
            const auto& a = left_branches[static_cast<std::size_t>(rank)];
            const auto& b = right_branches[static_cast<std::size_t>(rank)];
            if (!a.slingshot_valid || !b.slingshot_valid) {
                stats.residual_invalid_pair_count += 1;
                continue;
            }
            if ((a.slingshot_residual < 0.0 && b.slingshot_residual > 0.0) ||
                (a.slingshot_residual > 0.0 && b.slingshot_residual < 0.0)) {
                stats.sign_change_pair_count += 1;
            }
        }
    }
    return stats;
}

void print_sample_summary(int i, const SampleSummary& direct, const SampleSummary& table, const PairDiff& diff) {
    std::cout << "Problem2ThetaPrimeDirectVsTableSample"
              << " i=" << i
              << " theta_prime=" << direct.theta_prime
              << " direct_raw_branch_count=" << direct.branches.size()
              << " table_raw_branch_count=" << table.branches.size()
              << " direct_problem1_valid_count=" << problem1_valid_count(direct.branches)
              << " table_problem1_valid_count=" << problem1_valid_count(table.branches)
              << " direct_slingshot_valid_count=" << slingshot_valid_count(direct.branches)
              << " table_slingshot_valid_count=" << slingshot_valid_count(table.branches)
              << " direct_raw_count_by_k=" << format_count_by_k(direct.raw_count_by_k)
              << " table_raw_count_by_k=" << format_count_by_k(table.raw_count_by_k)
              << " direct_slingshot_valid_count_by_k=" << format_count_by_k(direct.slingshot_valid_count_by_k)
              << " table_slingshot_valid_count_by_k=" << format_count_by_k(table.slingshot_valid_count_by_k)
              << " table_fallback_used=" << (table.table_fallback_used ? 1 : 0)
              << " max_alpha_diff_sorted=" << diff.max_alpha_diff
              << " max_time_diff_sorted=" << diff.max_time_diff_seconds
              << " max_slingshot_residual_diff_sorted=" << diff.max_slingshot_residual_diff
              << '\n';
}

void print_branch_preview(double theta_prime, const char* source, const SampleSummary& sample) {
    const int preview_count = std::min<int>(kPreviewLimit, sample.branches.size());
    for (int rank = 0; rank < preview_count; ++rank) {
        const auto& branch = sample.branches[static_cast<std::size_t>(rank)];
        std::cout << "Problem2ThetaPrimeDirectVsTableBranchPreview"
                  << " theta_prime=" << theta_prime
                  << " source=" << source
                  << " rank=" << rank
                  << " k=" << branch.transfer_revolution
                  << " q=" << branch.target_revolution
                  << " alpha=" << branch.alpha
                  << " time_of_flight_seconds=" << branch.time_of_flight_seconds
                  << " problem1_residual_seconds=" << branch.problem1_residual_seconds
                  << " slingshot_residual=" << branch.slingshot_residual
                  << " problem1_valid=" << (branch.problem1_valid ? 1 : 0)
                  << " slingshot_valid=" << (branch.slingshot_valid ? 1 : 0)
                  << " problem1_invalid_reason=" << branch.problem1_invalid_reason
                  << " slingshot_invalid_reason=" << branch.slingshot_invalid_reason
                  << '\n';
    }
}

std::string mismatch_reason(
    bool raw_count_by_k_mismatch,
    bool stable_flag_mismatch,
    bool sign_change_count_mismatch,
    bool residual_invalid_pair_count_mismatch,
    bool slingshot_valid_count_by_k_mismatch
) {
    if (raw_count_by_k_mismatch) {
        return "raw_count_by_k_mismatch";
    }
    if (stable_flag_mismatch) {
        return "stable_flag_mismatch";
    }
    if (sign_change_count_mismatch) {
        return "sign_change_count_mismatch";
    }
    if (residual_invalid_pair_count_mismatch) {
        return "residual_invalid_pair_count_mismatch";
    }
    if (slingshot_valid_count_by_k_mismatch) {
        return "slingshot_valid_count_by_k_mismatch";
    }
    return "none";
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_direct_vs_table_diff_skipped_missing_table\n";
        return 0;
    }

    try {
        const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
        problem1::Problem1NearestNodeQueryOptions options{};
        options.residual_tolerance_seconds = 1e-2;
        options.max_newton_iterations = 80;
        options.fallback_direct_solve = true;

        const auto departure_planet = planet_params::PlanetId::Earth;
        const auto target_planet = planet_params::PlanetId::Mars;
        const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
        const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
        const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
        const double phi = departure_state.varphi;
        const double beta = target_state.varphi;
        const double incoming_e = 0.3;
        const double incoming_theta = 0.4;

        std::vector<SampleSummary> direct_samples;
        std::vector<SampleSummary> table_samples;
        std::vector<PairDiff> sample_diffs;
        direct_samples.reserve(kThetaSampleCount);
        table_samples.reserve(kThetaSampleCount);
        sample_diffs.reserve(kThetaSampleCount);

        int raw_branch_count_mismatch_count = 0;
        int raw_count_by_k_mismatch_count = 0;
        int slingshot_valid_count_by_k_mismatch_count = 0;
        int table_fallback_count = 0;
        int table_non_empty_sample_count = 0;
        std::map<std::string, int> direct_slingshot_invalid_reason_count;
        std::map<std::string, int> table_slingshot_invalid_reason_count;
        double max_alpha_diff = 0.0;
        double max_time_diff_seconds = 0.0;
        double max_slingshot_residual_diff = 0.0;

        for (int i = 0; i < kThetaSampleCount; ++i) {
            const double theta_prime = kTwoPi * static_cast<double>(i) / static_cast<double>(kThetaSampleCount);
            direct_samples.push_back(build_direct_summary(
                departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta, theta_prime));
            table_samples.push_back(build_table_summary(
                loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta,
                theta_prime));
            const auto& direct = direct_samples.back();
            const auto& table = table_samples.back();
            const PairDiff diff = compute_pair_diff(direct, table);
            sample_diffs.push_back(diff);
            max_alpha_diff = std::max(max_alpha_diff, diff.max_alpha_diff);
            max_time_diff_seconds = std::max(max_time_diff_seconds, diff.max_time_diff_seconds);
            max_slingshot_residual_diff =
                std::max(max_slingshot_residual_diff, diff.max_slingshot_residual_diff);
            if (direct.branches.size() != table.branches.size()) {
                raw_branch_count_mismatch_count += 1;
            }
            if (direct.raw_count_by_k != table.raw_count_by_k) {
                raw_count_by_k_mismatch_count += 1;
            }
            if (direct.slingshot_valid_count_by_k != table.slingshot_valid_count_by_k) {
                slingshot_valid_count_by_k_mismatch_count += 1;
            }
            if (table.table_fallback_used) {
                table_fallback_count += 1;
            }
            for (const auto& branch : direct.branches) {
                if (branch.problem1_valid && !branch.slingshot_valid) {
                    direct_slingshot_invalid_reason_count[
                        branch.slingshot_invalid_reason.empty() ? "unknown" : branch.slingshot_invalid_reason] += 1;
                }
            }
            for (const auto& branch : table.branches) {
                if (branch.problem1_valid && !branch.slingshot_valid) {
                    table_slingshot_invalid_reason_count[
                        branch.slingshot_invalid_reason.empty() ? "unknown" : branch.slingshot_invalid_reason] += 1;
                }
            }
            if (!table.branches.empty()) {
                table_non_empty_sample_count += 1;
            }
            print_sample_summary(i, direct, table, diff);
        }

        int interval_mismatch_count = 0;
        int stable_flag_mismatch_count = 0;
        int sign_change_count_mismatch_count = 0;
        int residual_invalid_pair_count_mismatch_count = 0;

        for (int i = 0; i + 1 < kThetaSampleCount; ++i) {
            const auto& direct_left = direct_samples[static_cast<std::size_t>(i)];
            const auto& direct_right = direct_samples[static_cast<std::size_t>(i + 1)];
            const auto& table_left = table_samples[static_cast<std::size_t>(i)];
            const auto& table_right = table_samples[static_cast<std::size_t>(i + 1)];

            const bool direct_stable = direct_left.raw_count_by_k == direct_right.raw_count_by_k;
            const bool table_stable = table_left.raw_count_by_k == table_right.raw_count_by_k;
            const SignChangeStats direct_sign_change = count_sign_change_pairs(direct_left, direct_right);
            const SignChangeStats table_sign_change = count_sign_change_pairs(table_left, table_right);
            const bool raw_count_by_k_mismatch =
                direct_left.raw_count_by_k != table_left.raw_count_by_k ||
                direct_right.raw_count_by_k != table_right.raw_count_by_k;
            const bool slingshot_valid_count_by_k_mismatch =
                direct_left.slingshot_valid_count_by_k != table_left.slingshot_valid_count_by_k ||
                direct_right.slingshot_valid_count_by_k != table_right.slingshot_valid_count_by_k;
            const bool stable_flag_mismatch = direct_stable != table_stable;
            const bool sign_change_count_mismatch =
                direct_sign_change.sign_change_pair_count != table_sign_change.sign_change_pair_count;
            const bool residual_invalid_pair_count_mismatch =
                direct_sign_change.residual_invalid_pair_count != table_sign_change.residual_invalid_pair_count;
            if (!raw_count_by_k_mismatch && !stable_flag_mismatch && !sign_change_count_mismatch &&
                !residual_invalid_pair_count_mismatch && !slingshot_valid_count_by_k_mismatch) {
                continue;
            }

            interval_mismatch_count += 1;
            if (stable_flag_mismatch) {
                stable_flag_mismatch_count += 1;
            }
            if (sign_change_count_mismatch) {
                sign_change_count_mismatch_count += 1;
            }
            if (residual_invalid_pair_count_mismatch) {
                residual_invalid_pair_count_mismatch_count += 1;
            }
            std::cout << "Problem2ThetaPrimeDirectVsTableIntervalMismatch"
                      << " i_left=" << i
                      << " i_right=" << (i + 1)
                      << " theta_left=" << direct_left.theta_prime
                      << " theta_right=" << direct_right.theta_prime
                      << " direct_left_raw_count_by_k=" << format_count_by_k(direct_left.raw_count_by_k)
                      << " direct_right_raw_count_by_k=" << format_count_by_k(direct_right.raw_count_by_k)
                      << " table_left_raw_count_by_k=" << format_count_by_k(table_left.raw_count_by_k)
                      << " table_right_raw_count_by_k=" << format_count_by_k(table_right.raw_count_by_k)
                      << " direct_left_slingshot_valid_count_by_k="
                      << format_count_by_k(direct_left.slingshot_valid_count_by_k)
                      << " direct_right_slingshot_valid_count_by_k="
                      << format_count_by_k(direct_right.slingshot_valid_count_by_k)
                      << " table_left_slingshot_valid_count_by_k="
                      << format_count_by_k(table_left.slingshot_valid_count_by_k)
                      << " table_right_slingshot_valid_count_by_k="
                      << format_count_by_k(table_right.slingshot_valid_count_by_k)
                      << " direct_stable=" << (direct_stable ? 1 : 0)
                      << " table_stable=" << (table_stable ? 1 : 0)
                      << " direct_sign_change_pair_count=" << direct_sign_change.sign_change_pair_count
                      << " table_sign_change_pair_count=" << table_sign_change.sign_change_pair_count
                      << " direct_residual_invalid_pair_count=" << direct_sign_change.residual_invalid_pair_count
                      << " table_residual_invalid_pair_count=" << table_sign_change.residual_invalid_pair_count
                      << " reason="
                      << mismatch_reason(
                             raw_count_by_k_mismatch,
                             stable_flag_mismatch,
                             sign_change_count_mismatch,
                             residual_invalid_pair_count_mismatch,
                             slingshot_valid_count_by_k_mismatch)
                      << '\n';
            print_branch_preview(direct_left.theta_prime, "direct", direct_left);
            print_branch_preview(direct_left.theta_prime, "table", table_left);
            print_branch_preview(direct_right.theta_prime, "direct", direct_right);
            print_branch_preview(direct_right.theta_prime, "table", table_right);
        }

        for (const auto& [reason, count] : direct_slingshot_invalid_reason_count) {
            std::cout << "Problem2ThetaPrimeSlingshotInvalidReasonSummary"
                      << " source=direct"
                      << " reason=" << reason
                      << " count=" << count
                      << '\n';
        }
        for (const auto& [reason, count] : table_slingshot_invalid_reason_count) {
            std::cout << "Problem2ThetaPrimeSlingshotInvalidReasonSummary"
                      << " source=table"
                      << " reason=" << reason
                      << " count=" << count
                      << '\n';
        }

        const bool ok = table_non_empty_sample_count > 0;
        std::cout << "Problem2ThetaPrimeDirectVsTableDiffSummary\n";
        std::cout << "theta_sample_count=" << kThetaSampleCount << '\n';
        std::cout << "raw_branch_count_mismatch_count=" << raw_branch_count_mismatch_count << '\n';
        std::cout << "raw_count_by_k_mismatch_count=" << raw_count_by_k_mismatch_count << '\n';
        std::cout << "slingshot_valid_count_by_k_mismatch_count="
                  << slingshot_valid_count_by_k_mismatch_count << '\n';
        std::cout << "interval_mismatch_count=" << interval_mismatch_count << '\n';
        std::cout << "stable_flag_mismatch_count=" << stable_flag_mismatch_count << '\n';
        std::cout << "sign_change_count_mismatch_count=" << sign_change_count_mismatch_count << '\n';
        std::cout << "residual_invalid_pair_count_mismatch_count="
                  << residual_invalid_pair_count_mismatch_count << '\n';
        std::cout << "table_fallback_count=" << table_fallback_count << '\n';
        std::cout << "max_alpha_diff=" << max_alpha_diff << '\n';
        std::cout << "max_time_diff_seconds=" << max_time_diff_seconds << '\n';
        std::cout << "max_slingshot_residual_diff=" << max_slingshot_residual_diff << '\n';
        std::cout << "diagnostic_ok=" << (ok ? 1 : 0) << '\n';
        return ok ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cout << "problem2_theta_prime_direct_vs_table_diff_error=" << ex.what() << '\n';
        return 1;
    }
}
