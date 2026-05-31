#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = kPi / 90.0;

struct GridAngles {
    std::array<double, kSamplesPerDimension> angles{};
};

struct NodeCoordinates {
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
};

struct BranchBinaryRecord {
    std::uint8_t valid = 0;
    std::int32_t transfer_revolution = 0;
    std::int32_t target_revolution = 0;
    double encounter_global_angle = 0.0;
    double target_arrival_true_anomaly = 0.0;
    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double residual_seconds = 0.0;
    double transfer_e = 0.0;
    double transfer_p = 0.0;
    double transfer_a = 0.0;
    double theta_B = 0.0;
    std::uint8_t derivatives_available = 0;
    double d_encounter_global_angle_d_nu_A = 0.0;
    double d_encounter_global_angle_d_nu_B = 0.0;
    double d_encounter_global_angle_d_theta_A = 0.0;
};

struct ChunkHeader {
    std::uint32_t version = 0;
    double grid_step_deg = 0.0;
    std::int32_t samples_per_dimension = 0;
    long long start_node = 0;
    long long end_node = 0;
    std::int32_t max_transfer_revolution = 0;
    std::int32_t max_target_revolution = 0;
    std::uint64_t node_count = 0;
};

struct FailureRecord {
    std::uint64_t linear_index = 0;
    NodeCoordinates node{};
    BranchBinaryRecord branch{};
};

std::string env_string(const char* name, const std::string& default_value) {
    if (const char* raw = std::getenv(name)) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return default_value;
}

template <class T>
bool read_scalar(std::ifstream& in, T* value) {
    in.read(reinterpret_cast<char*>(value), sizeof(*value));
    return static_cast<bool>(in);
}

GridAngles make_grid_angles() {
    GridAngles grid{};
    for (int i = 0; i < kSamplesPerDimension; ++i) {
        grid.angles[static_cast<std::size_t>(i)] =
            normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(i));
    }
    return grid;
}

NodeCoordinates decode_linear_index(long long linear_index, const GridAngles& grid) {
    NodeCoordinates node{};
    node.theta_A_index = static_cast<int>(linear_index % kSamplesPerDimension);
    linear_index /= kSamplesPerDimension;
    node.nu_B_index = static_cast<int>(linear_index % kSamplesPerDimension);
    linear_index /= kSamplesPerDimension;
    node.nu_A_index = static_cast<int>(linear_index % kSamplesPerDimension);
    node.nu_A = grid.angles[static_cast<std::size_t>(node.nu_A_index)];
    node.nu_B = grid.angles[static_cast<std::size_t>(node.nu_B_index)];
    node.theta_A = grid.angles[static_cast<std::size_t>(node.theta_A_index)];
    return node;
}

spaceship_cpp::problem1::Problem1SolutionBranch to_solution_branch(const BranchBinaryRecord& record) {
    spaceship_cpp::problem1::Problem1SolutionBranch branch{};
    branch.valid = record.valid != 0;
    branch.transfer_revolution = record.transfer_revolution;
    branch.target_revolution = record.target_revolution;
    branch.encounter_global_angle = record.encounter_global_angle;
    branch.target_arrival_true_anomaly = record.target_arrival_true_anomaly;
    branch.time_of_flight_seconds = record.time_of_flight_seconds;
    branch.target_time_seconds = record.target_time_seconds;
    branch.residual_seconds = record.residual_seconds;
    branch.transfer_e = record.transfer_e;
    branch.transfer_p = record.transfer_p;
    branch.transfer_a = record.transfer_a;
    branch.theta_B = record.theta_B;
    branch.derivatives_available = record.derivatives_available != 0;
    branch.d_encounter_global_angle_d_nu_A = record.d_encounter_global_angle_d_nu_A;
    branch.d_encounter_global_angle_d_nu_B = record.d_encounter_global_angle_d_nu_B;
    branch.d_encounter_global_angle_d_theta_A = record.d_encounter_global_angle_d_theta_A;
    return branch;
}

double angle_distance(double lhs, double rhs) {
    double delta = std::abs(normalize_angle_0_2pi(lhs) - normalize_angle_0_2pi(rhs));
    return std::min(delta, 2.0 * kPi - delta);
}

const spaceship_cpp::problem1::Problem1SolutionBranch* find_matching_rerun_branch(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    const BranchBinaryRecord& original
) {
    const spaceship_cpp::problem1::Problem1SolutionBranch* best = nullptr;
    double best_score = std::numeric_limits<double>::infinity();
    for (const auto& branch : branches) {
        if (!branch.valid ||
            branch.transfer_revolution != original.transfer_revolution ||
            branch.target_revolution != original.target_revolution) {
            continue;
        }
        const double score = angle_distance(branch.encounter_global_angle, original.encounter_global_angle) +
                             1e-6 * std::abs(branch.residual_seconds - original.residual_seconds) +
                             1e-9 * std::abs(branch.time_of_flight_seconds - original.time_of_flight_seconds);
        if (score < best_score) {
            best_score = score;
            best = &branch;
        }
    }
    return best;
}

bool read_header(std::ifstream& in, ChunkHeader* header, std::string* error) {
    char magic[8] = {};
    in.read(magic, sizeof(magic));
    if (!in || std::string(magic, sizeof(magic)) != "P1RT2DEG") {
        *error = "invalid_magic";
        return false;
    }
    if (!read_scalar(in, &header->version) ||
        !read_scalar(in, &header->grid_step_deg) ||
        !read_scalar(in, &header->samples_per_dimension) ||
        !read_scalar(in, &header->start_node) ||
        !read_scalar(in, &header->end_node) ||
        !read_scalar(in, &header->max_transfer_revolution) ||
        !read_scalar(in, &header->max_target_revolution) ||
        !read_scalar(in, &header->node_count)) {
        *error = "invalid_header";
        return false;
    }
    if (header->samples_per_dimension != kSamplesPerDimension) {
        *error = "unexpected_samples_per_dimension";
        return false;
    }
    return true;
}

std::vector<FailureRecord> scan_failures(
    const std::filesystem::path& bin_path,
    const GridAngles& grid,
    ChunkHeader* header,
    std::string* error
) {
    std::vector<FailureRecord> failures;
    std::ifstream in(bin_path, std::ios::binary);
    if (!in) {
        *error = "missing_or_unreadable_bin";
        return failures;
    }
    if (!read_header(in, header, error)) {
        return failures;
    }

    for (std::uint64_t node_index = 0; node_index < header->node_count; ++node_index) {
        std::uint64_t linear_index = 0;
        std::uint32_t branch_count = 0;
        if (!read_scalar(in, &linear_index) || !read_scalar(in, &branch_count)) {
            *error = "truncated_node_record";
            return failures;
        }
        const NodeCoordinates node = decode_linear_index(static_cast<long long>(linear_index), grid);
        for (std::uint32_t branch_index = 0; branch_index < branch_count; ++branch_index) {
            BranchBinaryRecord branch{};
            in.read(reinterpret_cast<char*>(&branch), sizeof(branch));
            if (!in) {
                *error = "truncated_branch_record";
                return failures;
            }
            if (branch.valid != 0 && branch.derivatives_available == 0) {
                failures.push_back(FailureRecord{linear_index, node, branch});
            }
        }
    }
    return failures;
}

void print_failure(const ChunkHeader& header, const FailureRecord& failure) {
    const BranchBinaryRecord& branch = failure.branch;
    std::cout << "Problem1RootTable2DegDerivativeFailure\n";
    std::cout << "chunk_start=" << header.start_node << '\n';
    std::cout << "chunk_end=" << header.end_node << '\n';
    std::cout << "linear_index=" << failure.linear_index << '\n';
    std::cout << "nu_A_index=" << failure.node.nu_A_index << '\n';
    std::cout << "nu_B_index=" << failure.node.nu_B_index << '\n';
    std::cout << "theta_A_index=" << failure.node.theta_A_index << '\n';
    std::cout << "transfer_revolution=" << branch.transfer_revolution << '\n';
    std::cout << "target_revolution=" << branch.target_revolution << '\n';
    std::cout << "encounter_global_angle=" << branch.encounter_global_angle << '\n';
    std::cout << "target_arrival_true_anomaly=" << branch.target_arrival_true_anomaly << '\n';
    std::cout << "time_of_flight_seconds=" << branch.time_of_flight_seconds << '\n';
    std::cout << "target_time_seconds=" << branch.target_time_seconds << '\n';
    std::cout << "residual_seconds=" << branch.residual_seconds << '\n';
    std::cout << "transfer_e=" << branch.transfer_e << '\n';
    std::cout << "transfer_p=" << branch.transfer_p << '\n';
    std::cout << "transfer_a=" << branch.transfer_a << '\n';
    std::cout << "theta_B=" << branch.theta_B << '\n';
}

void print_reproduction(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const ChunkHeader& header,
    const FailureRecord& failure
) {
    namespace problem1 = spaceship_cpp::problem1;

    const auto branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet,
        target_planet,
        failure.node.nu_A,
        failure.node.nu_B,
        failure.node.theta_A,
        header.max_transfer_revolution,
        header.max_target_revolution);
    const problem1::Problem1SolutionBranch* matched =
        find_matching_rerun_branch(branches, failure.branch);

    problem1::Problem1SolutionBranch original_branch = to_solution_branch(failure.branch);
    problem1::Problem1SolutionBranch branch_for_attach = matched != nullptr ? *matched : original_branch;
    const auto analytic = problem1::attach_problem1_root_derivatives_with_mode(
        departure_planet,
        target_planet,
        failure.node.nu_A,
        failure.node.nu_B,
        failure.node.theta_A,
        branch_for_attach,
        problem1::Problem1RootDerivativeMode::AnalyticOnly,
        1e-6);
    const auto finite_difference = problem1::attach_problem1_root_derivatives_with_mode(
        departure_planet,
        target_planet,
        failure.node.nu_A,
        failure.node.nu_B,
        failure.node.theta_A,
        branch_for_attach,
        problem1::Problem1RootDerivativeMode::FiniteDifferenceOnly,
        1e-6);

    std::cout << "Problem1DerivativeFailureReproduction\n";
    std::cout << "linear_index=" << failure.linear_index << '\n';
    std::cout << "k=" << failure.branch.transfer_revolution << '\n';
    std::cout << "q=" << failure.branch.target_revolution << '\n';
    std::cout << "original_residual_seconds=" << failure.branch.residual_seconds << '\n';
    std::cout << "rerun_branch_found=" << (matched != nullptr ? 1 : 0) << '\n';
    std::cout << "rerun_branch_residual_seconds="
              << (matched != nullptr ? matched->residual_seconds : std::numeric_limits<double>::quiet_NaN()) << '\n';
    std::cout << "analytic_derivative_success="
              << (analytic.valid && analytic.derivatives_available ? 1 : 0) << '\n';
    std::cout << "analytic_failure_reason=" << analytic.invalid_reason << '\n';
    std::cout << "finite_difference_derivative_success="
              << (finite_difference.valid && finite_difference.derivatives_available ? 1 : 0) << '\n';
    std::cout << "finite_difference_failure_reason=" << finite_difference.invalid_reason << '\n';
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;

    const auto bin_path = std::filesystem::path(env_string("ROOT_TABLE_CHUNK_BIN", ""));
    const auto meta_path = std::filesystem::path(env_string("ROOT_TABLE_CHUNK_META", ""));
    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;

    std::cout << std::setprecision(12) << std::scientific;
    std::cout << "Problem1RootTable2DegDerivativeFailureScanConfig\n";
    std::cout << "bin_path=" << bin_path.string() << '\n';
    std::cout << "meta_path=" << meta_path.string() << '\n';

    if (bin_path.empty()) {
        std::cout << "scan_error=missing_ROOT_TABLE_CHUNK_BIN\n";
        return 1;
    }

    const GridAngles grid = make_grid_angles();
    ChunkHeader header{};
    std::string error;
    const std::vector<FailureRecord> failures = scan_failures(bin_path, grid, &header, &error);
    if (!error.empty()) {
        std::cout << "scan_error=" << error << '\n';
        return 1;
    }

    std::cout << "Problem1RootTable2DegDerivativeFailureScanSummary\n";
    std::cout << "chunk_start=" << header.start_node << '\n';
    std::cout << "chunk_end=" << header.end_node << '\n';
    std::cout << "node_count=" << header.node_count << '\n';
    std::cout << "failure_count=" << failures.size() << '\n';

    for (const FailureRecord& failure : failures) {
        print_failure(header, failure);
        print_reproduction(departure_planet, target_planet, header, failure);
    }

    return 0;
}
