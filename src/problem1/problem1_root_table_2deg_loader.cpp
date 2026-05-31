#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace spaceship_cpp::problem1 {
namespace {

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepDegrees = 2.0;
constexpr std::uint32_t kFormatVersion = 1;
constexpr const char* kMagic = "P1RT2DEG";

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

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

template <class T>
bool read_scalar(std::istream& in, T* value) {
    in.read(reinterpret_cast<char*>(value), sizeof(*value));
    return static_cast<bool>(in);
}

ChunkHeader read_chunk_header(const std::filesystem::path& bin_path) {
    std::ifstream in(bin_path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("failed to open chunk bin: " + bin_path.string());
    }

    char magic[8] = {};
    in.read(magic, sizeof(magic));
    if (!in || std::string(magic, sizeof(magic)) != kMagic) {
        throw std::runtime_error("invalid chunk magic: " + bin_path.string());
    }

    ChunkHeader header{};
    if (!read_scalar(in, &header.version) ||
        !read_scalar(in, &header.grid_step_deg) ||
        !read_scalar(in, &header.samples_per_dimension) ||
        !read_scalar(in, &header.start_node) ||
        !read_scalar(in, &header.end_node) ||
        !read_scalar(in, &header.max_transfer_revolution) ||
        !read_scalar(in, &header.max_target_revolution) ||
        !read_scalar(in, &header.node_count)) {
        throw std::runtime_error("invalid chunk header: " + bin_path.string());
    }
    if (header.version != kFormatVersion) {
        throw std::runtime_error("unsupported chunk version: " + bin_path.string());
    }
    if (header.samples_per_dimension != kSamplesPerDimension) {
        throw std::runtime_error("unexpected samples_per_dimension: " + bin_path.string());
    }
    if (std::abs(header.grid_step_deg - kGridStepDegrees) > 1e-12) {
        throw std::runtime_error("unexpected grid_step_deg: " + bin_path.string());
    }
    if (header.end_node < header.start_node ||
        header.node_count != static_cast<std::uint64_t>(header.end_node - header.start_node)) {
        throw std::runtime_error("inconsistent chunk range/count: " + bin_path.string());
    }
    return header;
}

Problem1SolutionBranch to_solution_branch(const BranchBinaryRecord& record) {
    Problem1SolutionBranch branch{};
    branch.valid = record.valid != 0;
    branch.encounter_global_angle = record.encounter_global_angle;
    branch.target_arrival_true_anomaly = record.target_arrival_true_anomaly;
    branch.transfer_revolution = record.transfer_revolution;
    branch.target_revolution = record.target_revolution;
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

Problem1RootTable2DegLoadedNode invalid_node(long long linear_index, std::string reason) {
    Problem1RootTable2DegLoadedNode node{};
    node.linear_index = linear_index;
    node.invalid_reason = std::move(reason);
    Problem1RootTable2DegLoader::indices_from_linear_index(
        linear_index, &node.nu_A_index, &node.nu_B_index, &node.theta_A_index);
    return node;
}

}  // namespace

Problem1RootTable2DegLoader::Problem1RootTable2DegLoader(
    std::vector<Problem1RootTable2DegChunkInfo> chunks
) : chunks_(std::move(chunks)), runtime_indices_(chunks_.size()) {
    for (const auto& chunk : chunks_) {
        total_nodes_ += static_cast<long long>(chunk.node_count);
    }
}

bool Problem1RootTable2DegLoader::ensure_chunk_offsets_built(
    std::size_t chunk_index,
    std::string* invalid_reason,
    Problem1RootTable2DegLoadProfile* profile
) const {
    if (chunk_index >= chunks_.size()) {
        if (invalid_reason != nullptr) {
            *invalid_reason = "chunk_index_out_of_range";
        }
        return false;
    }
    auto& runtime = runtime_indices_[chunk_index];
    if (runtime.offsets_built) {
        return true;
    }

    const auto offset_start = Clock::now();
    const auto& chunk = chunks_[chunk_index];
    std::ifstream in(chunk.bin_path, std::ios::binary);
    if (!in) {
        if (invalid_reason != nullptr) {
            *invalid_reason = "failed_to_open_chunk";
        }
        return false;
    }

    char magic[8] = {};
    in.read(magic, sizeof(magic));
    std::uint32_t version = 0;
    double grid_step_deg = 0.0;
    std::int32_t samples_per_dimension = 0;
    long long start_node = 0;
    long long end_node = 0;
    std::int32_t max_transfer_revolution = 0;
    std::int32_t max_target_revolution = 0;
    std::uint64_t node_count = 0;
    if (!in || std::string(magic, sizeof(magic)) != kMagic ||
        !read_scalar(in, &version) ||
        !read_scalar(in, &grid_step_deg) ||
        !read_scalar(in, &samples_per_dimension) ||
        !read_scalar(in, &start_node) ||
        !read_scalar(in, &end_node) ||
        !read_scalar(in, &max_transfer_revolution) ||
        !read_scalar(in, &max_target_revolution) ||
        !read_scalar(in, &node_count)) {
        if (invalid_reason != nullptr) {
            *invalid_reason = "invalid_chunk_header";
        }
        return false;
    }
    if (version != kFormatVersion ||
        samples_per_dimension != kSamplesPerDimension ||
        std::abs(grid_step_deg - kGridStepDegrees) > 1e-12 ||
        start_node != chunk.start_node ||
        end_node != chunk.end_node ||
        node_count != chunk.node_count) {
        if (invalid_reason != nullptr) {
            *invalid_reason = "chunk_header_mismatch";
        }
        return false;
    }

    std::vector<std::uint64_t> offsets;
    offsets.reserve(static_cast<std::size_t>(node_count));
    for (std::uint64_t node_i = 0; node_i < node_count; ++node_i) {
        const auto position = in.tellg();
        if (position < 0) {
            if (invalid_reason != nullptr) {
                *invalid_reason = "invalid_node_stream_position";
            }
            return false;
        }
        offsets.push_back(static_cast<std::uint64_t>(position));

        std::uint64_t node_linear_index = 0;
        std::uint32_t branch_count = 0;
        if (!read_scalar(in, &node_linear_index) || !read_scalar(in, &branch_count)) {
            if (invalid_reason != nullptr) {
                *invalid_reason = "truncated_node_record";
            }
            return false;
        }
        const long long expected_linear_index = chunk.start_node + static_cast<long long>(node_i);
        if (static_cast<long long>(node_linear_index) != expected_linear_index) {
            if (invalid_reason != nullptr) {
                std::ostringstream os;
                os << "unexpected_node_linear_index:" << node_linear_index;
                *invalid_reason = os.str();
            }
            return false;
        }
        in.seekg(static_cast<std::streamoff>(branch_count * sizeof(BranchBinaryRecord)), std::ios::cur);
        if (!in) {
            if (invalid_reason != nullptr) {
                *invalid_reason = "truncated_branch_records";
            }
            return false;
        }
    }

    runtime.node_offsets = std::move(offsets);
    runtime.offsets_built = true;
    if (profile != nullptr) {
        profile->offset_built_this_call = true;
        profile->offset_build_ms += elapsed_ms(offset_start, Clock::now());
    }
    return true;
}

Problem1RootTable2DegLoader Problem1RootTable2DegLoader::open(const std::filesystem::path& table_dir) {
    if (!std::filesystem::exists(table_dir)) {
        throw std::runtime_error("root table directory does not exist: " + table_dir.string());
    }

    const std::regex chunk_name_re(R"(problem1_root_table_2deg_chunk_([0-9]+)_([0-9]+)\.bin)");
    std::vector<Problem1RootTable2DegChunkInfo> chunks;
    for (const auto& entry : std::filesystem::directory_iterator(table_dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const std::string filename = entry.path().filename().string();
        std::smatch match;
        if (!std::regex_match(filename, match, chunk_name_re)) {
            continue;
        }

        const long long name_start = std::stoll(match[1].str());
        const long long name_end = std::stoll(match[2].str());
        const ChunkHeader header = read_chunk_header(entry.path());
        if (header.start_node != name_start || header.end_node != name_end) {
            throw std::runtime_error("chunk filename/header range mismatch: " + entry.path().string());
        }

        Problem1RootTable2DegChunkInfo info{};
        info.bin_path = entry.path();
        info.meta_path = entry.path().parent_path() /
            ("problem1_root_table_2deg_chunk_" + std::to_string(header.start_node) + "_" +
             std::to_string(header.end_node) + ".meta.txt");
        info.start_node = header.start_node;
        info.end_node = header.end_node;
        info.node_count = header.node_count;
        info.file_size_bytes = std::filesystem::file_size(entry.path());
        chunks.push_back(std::move(info));
    }

    std::sort(chunks.begin(), chunks.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.start_node < rhs.start_node;
    });
    for (std::size_t i = 0; i < chunks.size(); ++i) {
        if (!std::filesystem::exists(chunks[i].meta_path)) {
            throw std::runtime_error("missing chunk meta: " + chunks[i].meta_path.string());
        }
        if (i > 0 && chunks[i - 1].end_node > chunks[i].start_node) {
            throw std::runtime_error("overlapping chunk ranges");
        }
        if (i > 0 && chunks[i - 1].end_node != chunks[i].start_node) {
            throw std::runtime_error("non-contiguous chunk ranges");
        }
    }

    return Problem1RootTable2DegLoader(std::move(chunks));
}

Problem1RootTable2DegLoadedNode Problem1RootTable2DegLoader::load_node_by_indices(
    int nu_A_index,
    int nu_B_index,
    int theta_A_index
) const {
    if (nu_A_index < 0 || nu_A_index >= kSamplesPerDimension ||
        nu_B_index < 0 || nu_B_index >= kSamplesPerDimension ||
        theta_A_index < 0 || theta_A_index >= kSamplesPerDimension) {
        Problem1RootTable2DegLoadedNode node{};
        node.invalid_reason = "index_out_of_range";
        return node;
    }
    return load_node_by_linear_index(
        linear_index_from_indices(nu_A_index, nu_B_index, theta_A_index));
}

Problem1RootTable2DegLoadedNode Problem1RootTable2DegLoader::load_node_by_linear_index(
    long long linear_index
) const {
    if (linear_index < 0 || linear_index >= total_nodes_) {
        return invalid_node(linear_index, "linear_index_out_of_range");
    }

    const auto it = std::upper_bound(
        chunks_.begin(),
        chunks_.end(),
        linear_index,
        [](long long value, const Problem1RootTable2DegChunkInfo& chunk) {
            return value < chunk.start_node;
        });
    if (it == chunks_.begin()) {
        return invalid_node(linear_index, "chunk_not_found");
    }
    const std::size_t chunk_index = static_cast<std::size_t>(std::distance(chunks_.begin(), it - 1));
    const auto& chunk = chunks_[chunk_index];
    if (linear_index < chunk.start_node || linear_index >= chunk.end_node) {
        return invalid_node(linear_index, "chunk_not_found");
    }

    Problem1RootTable2DegLoadProfile load_profile{};
    std::string offset_error;
    if (!ensure_chunk_offsets_built(chunk_index, &offset_error, &load_profile)) {
        auto node = invalid_node(linear_index, offset_error.empty() ? "failed_to_build_chunk_offsets" : offset_error);
        node.profile = load_profile;
        return node;
    }
    const auto& runtime = runtime_indices_[chunk_index];
    const std::uint64_t local_index = static_cast<std::uint64_t>(linear_index - chunk.start_node);
    if (local_index >= runtime.node_offsets.size()) {
        auto node = invalid_node(linear_index, "node_offset_not_found");
        node.profile = load_profile;
        return node;
    }

    const auto node_read_start = Clock::now();
    std::ifstream in(chunk.bin_path, std::ios::binary);
    if (!in) {
        load_profile.node_read_ms += elapsed_ms(node_read_start, Clock::now());
        auto node = invalid_node(linear_index, "failed_to_open_chunk");
        node.profile = load_profile;
        return node;
    }

    in.seekg(static_cast<std::streamoff>(runtime.node_offsets[local_index]), std::ios::beg);
    if (!in) {
        load_profile.node_read_ms += elapsed_ms(node_read_start, Clock::now());
        auto node = invalid_node(linear_index, "failed_to_seek_node_offset");
        node.profile = load_profile;
        return node;
    }

    std::uint64_t node_linear_index = 0;
    std::uint32_t branch_count = 0;
    if (!read_scalar(in, &node_linear_index) || !read_scalar(in, &branch_count)) {
        load_profile.node_read_ms += elapsed_ms(node_read_start, Clock::now());
        auto node = invalid_node(linear_index, "truncated_node_record");
        node.profile = load_profile;
        return node;
    }
    if (static_cast<long long>(node_linear_index) != linear_index) {
        load_profile.node_read_ms += elapsed_ms(node_read_start, Clock::now());
        auto node = invalid_node(linear_index, "node_linear_index_mismatch");
        node.profile = load_profile;
        return node;
    }

    Problem1RootTable2DegLoadedNode node{};
    node.valid = true;
    node.linear_index = linear_index;
    indices_from_linear_index(linear_index, &node.nu_A_index, &node.nu_B_index, &node.theta_A_index);
    node.branches.reserve(branch_count);
    for (std::uint32_t branch_i = 0; branch_i < branch_count; ++branch_i) {
        BranchBinaryRecord record{};
        in.read(reinterpret_cast<char*>(&record), sizeof(record));
        if (!in) {
            load_profile.node_read_ms += elapsed_ms(node_read_start, Clock::now());
            auto invalid = invalid_node(linear_index, "truncated_branch_record");
            invalid.profile = load_profile;
            return invalid;
        }
        node.branches.push_back(to_solution_branch(record));
    }
    load_profile.node_read_ms += elapsed_ms(node_read_start, Clock::now());
    node.profile = load_profile;
    return node;
}

std::size_t Problem1RootTable2DegLoader::chunk_count() const {
    return chunks_.size();
}

long long Problem1RootTable2DegLoader::total_nodes() const {
    return total_nodes_;
}

const std::vector<Problem1RootTable2DegChunkInfo>& Problem1RootTable2DegLoader::chunks() const {
    return chunks_;
}

long long Problem1RootTable2DegLoader::linear_index_from_indices(
    int nu_A_index,
    int nu_B_index,
    int theta_A_index
) {
    return ((static_cast<long long>(nu_A_index) * kSamplesPerDimension) + nu_B_index) *
               kSamplesPerDimension +
           theta_A_index;
}

void Problem1RootTable2DegLoader::indices_from_linear_index(
    long long linear_index,
    int* nu_A_index,
    int* nu_B_index,
    int* theta_A_index
) {
    if (theta_A_index != nullptr) {
        *theta_A_index = static_cast<int>(linear_index % kSamplesPerDimension);
    }
    linear_index /= kSamplesPerDimension;
    if (nu_B_index != nullptr) {
        *nu_B_index = static_cast<int>(linear_index % kSamplesPerDimension);
    }
    linear_index /= kSamplesPerDimension;
    if (nu_A_index != nullptr) {
        *nu_A_index = static_cast<int>(linear_index);
    }
}

}  // namespace spaceship_cpp::problem1
