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

struct MetaInfo {
    long long start_node = 0;
    long long end_node = 0;
    long long completed_node_count = 0;
    double wall_time_seconds = 0.0;
    double mean_node_ms = 0.0;
    long long total_branch_count = 0;
    long long derivative_success_count = 0;
    long long derivative_fail_count = 0;
    long long invalid_node_count = 0;
};

struct FailureRecord {
    std::uint64_t linear_index = 0;
    std::streamoff branch_offset = 0;
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
bool read_scalar(std::istream& in, T* value) {
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

bool read_header(std::istream& in, ChunkHeader* header, std::string* error) {
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

bool read_meta(const std::filesystem::path& meta_path, MetaInfo* meta, std::string* error) {
    std::ifstream in(meta_path);
    if (!in) {
        *error = "missing_or_unreadable_meta";
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        auto read_ll = [&](const std::string& prefix, long long* out) {
            if (line.rfind(prefix, 0) == 0) {
                *out = std::atoll(line.substr(prefix.size()).c_str());
            }
        };
        auto read_double = [&](const std::string& prefix, double* out) {
            if (line.rfind(prefix, 0) == 0) {
                *out = std::atof(line.substr(prefix.size()).c_str());
            }
        };
        read_ll("start_node=", &meta->start_node);
        read_ll("end_node=", &meta->end_node);
        read_ll("completed_node_count=", &meta->completed_node_count);
        read_double("wall_time_seconds=", &meta->wall_time_seconds);
        read_double("mean_node_ms=", &meta->mean_node_ms);
        read_ll("total_branch_count=", &meta->total_branch_count);
        read_ll("derivative_success_count=", &meta->derivative_success_count);
        read_ll("derivative_fail_count=", &meta->derivative_fail_count);
        read_ll("invalid_node_count=", &meta->invalid_node_count);
    }
    return true;
}

bool write_meta(const std::filesystem::path& meta_path, const MetaInfo& meta, std::string* error) {
    std::ofstream out(meta_path);
    if (!out) {
        *error = "meta_write_failed";
        return false;
    }
    out << std::setprecision(6) << std::scientific;
    out << "start_node=" << meta.start_node << '\n';
    out << "end_node=" << meta.end_node << '\n';
    out << "completed_node_count=" << meta.completed_node_count << '\n';
    out << "wall_time_seconds=" << meta.wall_time_seconds << '\n';
    out << "mean_node_ms=" << meta.mean_node_ms << '\n';
    out << "total_branch_count=" << meta.total_branch_count << '\n';
    out << "derivative_success_count=" << meta.derivative_success_count << '\n';
    out << "derivative_fail_count=" << meta.derivative_fail_count << '\n';
    out << "invalid_node_count=" << meta.invalid_node_count << '\n';
    return true;
}

std::vector<FailureRecord> scan_failures(
    std::fstream& io,
    const GridAngles& grid,
    ChunkHeader* header,
    std::string* error
) {
    std::vector<FailureRecord> failures;
    io.seekg(0, std::ios::beg);
    if (!read_header(io, header, error)) {
        return failures;
    }

    for (std::uint64_t node_index = 0; node_index < header->node_count; ++node_index) {
        std::uint64_t linear_index = 0;
        std::uint32_t branch_count = 0;
        if (!read_scalar(io, &linear_index) || !read_scalar(io, &branch_count)) {
            *error = "truncated_node_record";
            return failures;
        }
        const NodeCoordinates node = decode_linear_index(static_cast<long long>(linear_index), grid);
        for (std::uint32_t branch_index = 0; branch_index < branch_count; ++branch_index) {
            const std::streamoff branch_offset = static_cast<std::streamoff>(io.tellg());
            BranchBinaryRecord branch{};
            io.read(reinterpret_cast<char*>(&branch), sizeof(branch));
            if (!io) {
                *error = "truncated_branch_record";
                return failures;
            }
            if (branch.valid != 0 && branch.derivatives_available == 0) {
                failures.push_back(FailureRecord{linear_index, branch_offset, node, branch});
            }
        }
    }
    return failures;
}

bool copy_backup_if_needed(const std::filesystem::path& path) {
    const auto backup_path = path.string() + ".bak";
    if (std::filesystem::exists(backup_path)) {
        std::cout << "backup_warning=backup_exists_not_overwritten:" << backup_path << '\n';
        return true;
    }
    std::error_code ec;
    std::filesystem::copy_file(path, backup_path, std::filesystem::copy_options::none, ec);
    if (ec) {
        std::cout << "backup_error=" << ec.message() << '\n';
        return false;
    }
    std::cout << "backup_created=" << backup_path << '\n';
    return true;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const auto bin_path = std::filesystem::path(env_string("ROOT_TABLE_CHUNK_BIN", ""));
    const auto meta_path = std::filesystem::path(env_string("ROOT_TABLE_CHUNK_META", ""));
    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;

    std::cout << std::setprecision(12) << std::scientific;
    if (bin_path.empty() || meta_path.empty()) {
        std::cout << "patch_error=missing_ROOT_TABLE_CHUNK_BIN_or_META\n";
        return 1;
    }

    MetaInfo meta{};
    std::string error;
    if (!read_meta(meta_path, &meta, &error)) {
        std::cout << "patch_error=" << error << '\n';
        return 1;
    }
    if (!copy_backup_if_needed(bin_path) || !copy_backup_if_needed(meta_path)) {
        return 1;
    }

    std::fstream io(bin_path, std::ios::binary | std::ios::in | std::ios::out);
    if (!io) {
        std::cout << "patch_error=missing_or_unwritable_bin\n";
        return 1;
    }

    const GridAngles grid = make_grid_angles();
    ChunkHeader header{};
    const std::vector<FailureRecord> failures = scan_failures(io, grid, &header, &error);
    if (!error.empty()) {
        std::cout << "patch_error=" << error << '\n';
        return 1;
    }
    if (meta.derivative_fail_count - static_cast<long long>(failures.size()) < 0) {
        std::cout << "patch_error=meta_derivative_fail_count_would_be_negative\n";
        return 1;
    }

    const long long success_before = meta.derivative_success_count;
    const long long fail_before = meta.derivative_fail_count;
    long long patched_count = 0;
    long long patch_failed_count = 0;

    for (const FailureRecord& failure : failures) {
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
        if (matched == nullptr) {
            patch_failed_count += 1;
            continue;
        }
        const auto finite_difference = problem1::attach_problem1_root_derivatives_with_mode(
            departure_planet,
            target_planet,
            failure.node.nu_A,
            failure.node.nu_B,
            failure.node.theta_A,
            *matched,
            problem1::Problem1RootDerivativeMode::FiniteDifferenceOnly,
            1e-6);
        if (!finite_difference.valid || !finite_difference.derivatives_available) {
            patch_failed_count += 1;
            continue;
        }

        BranchBinaryRecord patched = failure.branch;
        patched.derivatives_available = 1;
        patched.d_encounter_global_angle_d_nu_A = finite_difference.d_encounter_global_angle_d_nu_A;
        patched.d_encounter_global_angle_d_nu_B = finite_difference.d_encounter_global_angle_d_nu_B;
        patched.d_encounter_global_angle_d_theta_A = finite_difference.d_encounter_global_angle_d_theta_A;

        io.clear();
        io.seekp(failure.branch_offset, std::ios::beg);
        io.write(reinterpret_cast<const char*>(&patched), sizeof(patched));
        if (!io) {
            std::cout << "patch_error=bin_write_failed\n";
            return 1;
        }
        patched_count += 1;

        std::cout << "Problem1RootTable2DegDerivativePatchedBranch\n";
        std::cout << "linear_index=" << failure.linear_index << '\n';
        std::cout << "nu_A_index=" << failure.node.nu_A_index << '\n';
        std::cout << "nu_B_index=" << failure.node.nu_B_index << '\n';
        std::cout << "theta_A_index=" << failure.node.theta_A_index << '\n';
        std::cout << "k=" << failure.branch.transfer_revolution << '\n';
        std::cout << "q=" << failure.branch.target_revolution << '\n';
        std::cout << "original_residual_seconds=" << failure.branch.residual_seconds << '\n';
        std::cout << "finite_difference_success=1\n";
        std::cout << "d_encounter_global_angle_d_nu_A=" << patched.d_encounter_global_angle_d_nu_A << '\n';
        std::cout << "d_encounter_global_angle_d_nu_B=" << patched.d_encounter_global_angle_d_nu_B << '\n';
        std::cout << "d_encounter_global_angle_d_theta_A=" << patched.d_encounter_global_angle_d_theta_A << '\n';
    }
    io.flush();

    meta.derivative_success_count += patched_count;
    meta.derivative_fail_count -= patched_count;
    if (meta.derivative_fail_count < 0) {
        std::cout << "patch_error=meta_derivative_fail_count_would_be_negative\n";
        return 1;
    }
    if (patched_count > 0 && !write_meta(meta_path, meta, &error)) {
        std::cout << "patch_error=" << error << '\n';
        return 1;
    }

    std::cout << "Problem1RootTable2DegDerivativePatchSummary\n";
    std::cout << "bin_path=" << bin_path.string() << '\n';
    std::cout << "meta_path=" << meta_path.string() << '\n';
    std::cout << "failure_count_before=" << failures.size() << '\n';
    std::cout << "patched_count=" << patched_count << '\n';
    std::cout << "patch_failed_count=" << patch_failed_count << '\n';
    std::cout << "derivative_success_count_before=" << success_before << '\n';
    std::cout << "derivative_success_count_after=" << meta.derivative_success_count << '\n';
    std::cout << "derivative_fail_count_before=" << fail_before << '\n';
    std::cout << "derivative_fail_count_after=" << meta.derivative_fail_count << '\n';

    return patch_failed_count == 0 ? 0 : 1;
}
