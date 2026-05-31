#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using spaceship_cpp::common::kPi;
using spaceship_cpp::common::normalize_angle_0_2pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepDegrees = 2.0;
constexpr double kGridStepRadians = kPi / 90.0;
constexpr long long kTotalNodes =
    1LL * kSamplesPerDimension * kSamplesPerDimension * kSamplesPerDimension;
constexpr std::uint32_t kFormatVersion = 1;

struct BuildConfig {
    long long start_node = 0;
    long long end_node = kTotalNodes;
    long long chunk_size = 100000;
    int thread_count = 1;
    int max_transfer_revolution = 1;
    int max_target_revolution = 1;
    std::filesystem::path output_dir;
};

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

struct NodeBinaryRecord {
    std::uint64_t linear_index = 0;
    std::vector<BranchBinaryRecord> branches;
};

struct ChunkStats {
    long long completed_node_count = 0;
    long long total_branch_count = 0;
    long long derivative_success_count = 0;
    long long derivative_fail_count = 0;
    long long invalid_node_count = 0;
    double wall_time_seconds = 0.0;
    double mean_node_ms = 0.0;
};

struct WorkerResult {
    std::vector<NodeBinaryRecord> nodes;
    long long total_branch_count = 0;
    long long derivative_success_count = 0;
    long long derivative_fail_count = 0;
    long long invalid_node_count = 0;
};

struct ChunkHeaderValidation {
    bool valid = false;
    std::string invalid_reason;
    std::uint32_t version = 0;
    double grid_step_deg = 0.0;
    std::int32_t samples_per_dimension = 0;
    long long start_node = 0;
    long long end_node = 0;
    std::int32_t max_transfer_revolution = 0;
    std::int32_t max_target_revolution = 0;
    std::uint64_t node_count = 0;
};

int env_int(const char* name, int default_value) {
    if (const char* raw = std::getenv(name)) {
        const int value = std::atoi(raw);
        if (value > 0) {
            return value;
        }
    }
    return default_value;
}

long long env_long_long(const char* name, long long default_value) {
    if (const char* raw = std::getenv(name)) {
        const long long value = std::atoll(raw);
        if (value >= 0) {
            return value;
        }
    }
    return default_value;
}

std::string env_string(const char* name, const std::string& default_value) {
    if (const char* raw = std::getenv(name)) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return default_value;
}

BuildConfig read_config() {
    BuildConfig config{};
    config.start_node = env_long_long("ROOT_TABLE_START_NODE", 0);
    config.end_node = env_long_long("ROOT_TABLE_END_NODE", kTotalNodes);
    config.chunk_size = env_long_long("ROOT_TABLE_CHUNK_SIZE", 100000);
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    config.thread_count = env_int("ROOT_TABLE_THREAD_COUNT", static_cast<int>(std::min(16u, hw)));
    config.max_transfer_revolution = env_int("ROOT_TABLE_MAX_TRANSFER_REVOLUTION", 1);
    config.max_target_revolution = env_int("ROOT_TABLE_MAX_TARGET_REVOLUTION", 1);
    config.output_dir = env_string("ROOT_TABLE_OUTPUT_DIR", "./results/problem1_root_table_2deg");
    config.start_node = std::clamp<long long>(config.start_node, 0, kTotalNodes);
    config.end_node = std::clamp<long long>(config.end_node, config.start_node, kTotalNodes);
    config.chunk_size = std::max<long long>(1, config.chunk_size);
    config.thread_count = std::max(1, config.thread_count);
    return config;
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

BranchBinaryRecord to_branch_record(const spaceship_cpp::problem1::Problem1SolutionBranch& branch) {
    BranchBinaryRecord record{};
    record.valid = static_cast<std::uint8_t>(branch.valid ? 1 : 0);
    record.transfer_revolution = branch.transfer_revolution;
    record.target_revolution = branch.target_revolution;
    record.encounter_global_angle = branch.encounter_global_angle;
    record.target_arrival_true_anomaly = branch.target_arrival_true_anomaly;
    record.time_of_flight_seconds = branch.time_of_flight_seconds;
    record.target_time_seconds = branch.target_time_seconds;
    record.residual_seconds = branch.residual_seconds;
    record.transfer_e = branch.transfer_e;
    record.transfer_p = branch.transfer_p;
    record.transfer_a = branch.transfer_a;
    record.theta_B = branch.theta_B;
    record.derivatives_available = static_cast<std::uint8_t>(branch.derivatives_available ? 1 : 0);
    record.d_encounter_global_angle_d_nu_A = branch.d_encounter_global_angle_d_nu_A;
    record.d_encounter_global_angle_d_nu_B = branch.d_encounter_global_angle_d_nu_B;
    record.d_encounter_global_angle_d_theta_A = branch.d_encounter_global_angle_d_theta_A;
    return record;
}

WorkerResult process_range(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const GridAngles& grid,
    long long begin_node,
    long long end_node,
    int max_transfer_revolution,
    int max_target_revolution
) {
    namespace problem1 = spaceship_cpp::problem1;

    WorkerResult result{};
    result.nodes.reserve(static_cast<std::size_t>(std::max<long long>(0, end_node - begin_node)));
    for (long long linear_index = begin_node; linear_index < end_node; ++linear_index) {
        const NodeCoordinates node = decode_linear_index(linear_index, grid);
        NodeBinaryRecord node_record{};
        node_record.linear_index = static_cast<std::uint64_t>(linear_index);

        auto branches = problem1::solve_problem1_from_departure_anomalies(
            departure_planet,
            target_planet,
            node.nu_A,
            node.nu_B,
            node.theta_A,
            max_transfer_revolution,
            max_target_revolution);
        result.total_branch_count += static_cast<long long>(branches.size());
        bool has_valid_branch = false;
        for (auto& branch : branches) {
            if (branch.valid) {
                has_valid_branch = true;
                const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
                    departure_planet,
                    target_planet,
                    node.nu_A,
                    node.nu_B,
                    node.theta_A,
                    branch,
                    problem1::Problem1RootDerivativeMode::AnalyticOnly,
                    1e-6);
                if (attached.valid && attached.derivatives_available) {
                    branch = attached;
                    result.derivative_success_count += 1;
                } else {
                    result.derivative_fail_count += 1;
                }
            }
            node_record.branches.push_back(to_branch_record(branch));
        }
        if (!has_valid_branch) {
            result.invalid_node_count += 1;
        }
        result.nodes.push_back(std::move(node_record));
    }
    return result;
}

void write_scalar(std::ofstream& out, const auto& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template <class T>
bool read_scalar(std::ifstream& in, T* value) {
    in.read(reinterpret_cast<char*>(value), sizeof(*value));
    return static_cast<bool>(in);
}

ChunkHeaderValidation validate_chunk_header(
    const std::filesystem::path& bin_path,
    long long expected_start,
    long long expected_end,
    long long expected_count
) {
    ChunkHeaderValidation validation{};
    if (!std::filesystem::exists(bin_path)) {
        validation.invalid_reason = "missing_bin";
        return validation;
    }

    std::ifstream in(bin_path, std::ios::binary);
    char magic[8] = {};
    in.read(magic, sizeof(magic));
    if (!in || std::string(magic, sizeof(magic)) != "P1RT2DEG") {
        validation.invalid_reason = "invalid_bin_header";
        return validation;
    }
    if (!read_scalar(in, &validation.version) ||
        !read_scalar(in, &validation.grid_step_deg) ||
        !read_scalar(in, &validation.samples_per_dimension) ||
        !read_scalar(in, &validation.start_node) ||
        !read_scalar(in, &validation.end_node) ||
        !read_scalar(in, &validation.max_transfer_revolution) ||
        !read_scalar(in, &validation.max_target_revolution) ||
        !read_scalar(in, &validation.node_count)) {
        validation.invalid_reason = "invalid_bin_header";
        return validation;
    }
    if (validation.start_node != expected_start || validation.end_node != expected_end) {
        validation.invalid_reason = "range_mismatch";
        return validation;
    }
    if (validation.node_count != static_cast<std::uint64_t>(expected_count)) {
        validation.invalid_reason = "node_count_mismatch";
        return validation;
    }
    validation.valid = true;
    return validation;
}

void write_chunk_binary(
    const std::filesystem::path& output_path,
    const BuildConfig& config,
    long long chunk_start,
    long long chunk_end,
    const std::vector<NodeBinaryRecord>& nodes
) {
    std::ofstream out(output_path, std::ios::binary);
    const char magic[8] = {'P', '1', 'R', 'T', '2', 'D', 'E', 'G'};
    out.write(magic, sizeof(magic));
    write_scalar(out, kFormatVersion);
    write_scalar(out, kGridStepDegrees);
    write_scalar(out, kSamplesPerDimension);
    write_scalar(out, chunk_start);
    write_scalar(out, chunk_end);
    write_scalar(out, config.max_transfer_revolution);
    write_scalar(out, config.max_target_revolution);
    const std::uint64_t node_count = static_cast<std::uint64_t>(nodes.size());
    write_scalar(out, node_count);
    for (const auto& node : nodes) {
        write_scalar(out, node.linear_index);
        const std::uint32_t branch_count = static_cast<std::uint32_t>(node.branches.size());
        write_scalar(out, branch_count);
        for (const auto& branch : node.branches) {
            out.write(reinterpret_cast<const char*>(&branch), sizeof(branch));
        }
    }
}

void write_chunk_meta(
    const std::filesystem::path& meta_path,
    long long chunk_start,
    long long chunk_end,
    const ChunkStats& stats
) {
    std::ofstream out(meta_path);
    out << std::setprecision(6) << std::scientific;
    out << "start_node=" << chunk_start << '\n';
    out << "end_node=" << chunk_end << '\n';
    out << "completed_node_count=" << stats.completed_node_count << '\n';
    out << "wall_time_seconds=" << stats.wall_time_seconds << '\n';
    out << "mean_node_ms=" << stats.mean_node_ms << '\n';
    out << "total_branch_count=" << stats.total_branch_count << '\n';
    out << "derivative_success_count=" << stats.derivative_success_count << '\n';
    out << "derivative_fail_count=" << stats.derivative_fail_count << '\n';
    out << "invalid_node_count=" << stats.invalid_node_count << '\n';
}

bool chunk_already_complete(
    const std::filesystem::path& bin_path,
    const std::filesystem::path& meta_path,
    long long expected_start,
    long long expected_end,
    long long expected_count,
    std::string* invalid_reason
) {
    if (!std::filesystem::exists(meta_path)) {
        return false;
    }
    std::ifstream in(meta_path);
    std::string line;
    bool count_ok = false;
    while (std::getline(in, line)) {
        const std::string prefix = "completed_node_count=";
        if (line.rfind(prefix, 0) == 0) {
            const long long completed = std::atoll(line.substr(prefix.size()).c_str());
            count_ok = completed == expected_count;
            break;
        }
    }
    if (!count_ok) {
        if (invalid_reason != nullptr) {
            *invalid_reason = "meta_count_mismatch";
        }
        return false;
    }
    const auto header_validation = validate_chunk_header(bin_path, expected_start, expected_end, expected_count);
    if (!header_validation.valid) {
        if (invalid_reason != nullptr) {
            *invalid_reason = header_validation.invalid_reason;
        }
        return false;
    }
    return true;
}

ChunkStats build_chunk(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    const GridAngles& grid,
    const BuildConfig& config,
    long long chunk_start,
    long long chunk_end
) {
    using clock = std::chrono::steady_clock;

    const long long node_count = chunk_end - chunk_start;
    const int thread_count = std::max(1, std::min<int>(config.thread_count, static_cast<int>(node_count)));
    std::vector<WorkerResult> worker_results(static_cast<std::size_t>(thread_count));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    const auto begin = clock::now();
    long long cursor = chunk_start;
    const long long base = node_count / thread_count;
    const long long extra = node_count % thread_count;
    for (int worker_index = 0; worker_index < thread_count; ++worker_index) {
        const long long count = base + (worker_index < extra ? 1 : 0);
        const long long worker_begin = cursor;
        const long long worker_end = cursor + count;
        cursor = worker_end;
        workers.emplace_back([&, worker_index, worker_begin, worker_end]() {
            worker_results[static_cast<std::size_t>(worker_index)] = process_range(
                departure_planet,
                target_planet,
                grid,
                worker_begin,
                worker_end,
                config.max_transfer_revolution,
                config.max_target_revolution);
        });
    }
    for (auto& worker : workers) {
        worker.join();
    }
    const auto end = clock::now();

    std::vector<NodeBinaryRecord> all_nodes;
    all_nodes.reserve(static_cast<std::size_t>(node_count));
    ChunkStats stats{};
    stats.completed_node_count = node_count;
    for (auto& worker_result : worker_results) {
        stats.total_branch_count += worker_result.total_branch_count;
        stats.derivative_success_count += worker_result.derivative_success_count;
        stats.derivative_fail_count += worker_result.derivative_fail_count;
        stats.invalid_node_count += worker_result.invalid_node_count;
        std::move(
            worker_result.nodes.begin(),
            worker_result.nodes.end(),
            std::back_inserter(all_nodes));
    }
    std::sort(
        all_nodes.begin(),
        all_nodes.end(),
        [](const NodeBinaryRecord& lhs, const NodeBinaryRecord& rhs) { return lhs.linear_index < rhs.linear_index; });

    const auto chunk_stem = "problem1_root_table_2deg_chunk_" + std::to_string(chunk_start) + "_" +
                            std::to_string(chunk_end);
    const auto bin_path = config.output_dir / (chunk_stem + ".bin");
    const auto meta_path = config.output_dir / (chunk_stem + ".meta.txt");
    write_chunk_binary(bin_path, config, chunk_start, chunk_end, all_nodes);
    stats.wall_time_seconds = std::chrono::duration<double>(end - begin).count();
    stats.mean_node_ms =
        node_count > 0 ? 1000.0 * stats.wall_time_seconds / static_cast<double>(node_count) : 0.0;
    write_chunk_meta(meta_path, chunk_start, chunk_end, stats);
    return stats;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;

    const BuildConfig config = read_config();
    const GridAngles grid = make_grid_angles();
    std::filesystem::create_directories(config.output_dir);

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1RootTable2DegBuildConfig\n";
    std::cout << "grid_step_deg=" << kGridStepDegrees << '\n';
    std::cout << "samples_per_dimension=" << kSamplesPerDimension << '\n';
    std::cout << "total_nodes=" << kTotalNodes << '\n';
    std::cout << "start_node=" << config.start_node << '\n';
    std::cout << "end_node=" << config.end_node << '\n';
    std::cout << "thread_count=" << config.thread_count << '\n';
    std::cout << "max_transfer_revolution=" << config.max_transfer_revolution << '\n';
    std::cout << "max_target_revolution=" << config.max_target_revolution << '\n';
    std::cout << "output_dir=" << config.output_dir.string() << '\n';

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    for (long long chunk_start = config.start_node; chunk_start < config.end_node; chunk_start += config.chunk_size) {
        const long long chunk_end = std::min(config.end_node, chunk_start + config.chunk_size);
        const long long expected_count = chunk_end - chunk_start;
        const auto chunk_stem = "problem1_root_table_2deg_chunk_" + std::to_string(chunk_start) + "_" +
                                std::to_string(chunk_end);
        const auto bin_path = config.output_dir / (chunk_stem + ".bin");
        const auto meta_path = config.output_dir / (chunk_stem + ".meta.txt");
        std::string checkpoint_invalid_reason;
        if (chunk_already_complete(
                bin_path,
                meta_path,
                chunk_start,
                chunk_end,
                expected_count,
                &checkpoint_invalid_reason)) {
            std::cout << "Problem1RootTable2DegChunkSkipped\n";
            std::cout << "start_node=" << chunk_start << '\n';
            std::cout << "end_node=" << chunk_end << '\n';
            std::cout << "reason=checkpoint_complete\n";
            continue;
        }
        if (std::filesystem::exists(meta_path) || std::filesystem::exists(bin_path)) {
            std::cout << "Problem1RootTable2DegCheckpointInvalid\n";
            std::cout << "start_node=" << chunk_start << '\n';
            std::cout << "end_node=" << chunk_end << '\n';
            std::cout << "reason=" << checkpoint_invalid_reason << '\n';
            std::cout << "action=rebuild\n";
        }

        const ChunkStats stats = build_chunk(
            departure_planet, target_planet, grid, config, chunk_start, chunk_end);
        std::cout << "Problem1RootTable2DegChunkSummary\n";
        std::cout << "start_node=" << chunk_start << '\n';
        std::cout << "end_node=" << chunk_end << '\n';
        std::cout << "completed_node_count=" << stats.completed_node_count << '\n';
        std::cout << "wall_time_seconds=" << stats.wall_time_seconds << '\n';
        std::cout << "mean_node_ms=" << stats.mean_node_ms << '\n';
        std::cout << "total_branch_count=" << stats.total_branch_count << '\n';
        std::cout << "derivative_success_count=" << stats.derivative_success_count << '\n';
        std::cout << "derivative_fail_count=" << stats.derivative_fail_count << '\n';
        std::cout << "invalid_node_count=" << stats.invalid_node_count << '\n';
    }

    std::cout << "build_problem1_root_table_2deg_ok\n";
    return 0;
}
