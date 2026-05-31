#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

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

struct MetaInfo {
    long long completed_node_count = -1;
    long long total_branch_count = -1;
    long long derivative_success_count = -1;
};

struct ValidationSummary {
    bool valid = false;
    std::string validation_error;
    std::uint32_t version = 0;
    double grid_step_deg = 0.0;
    std::int32_t samples_per_dimension = 0;
    long long start_node = 0;
    long long end_node = 0;
    std::uint64_t node_count_header = 0;
    std::uint64_t node_count_read = 0;
    std::uint64_t total_branch_count_read = 0;
    std::uint64_t valid_branch_count_read = 0;
    std::uint64_t derivatives_available_count_read = 0;
    long long meta_completed_node_count = -1;
    long long meta_total_branch_count = -1;
    long long meta_derivative_success_count = -1;
    std::uint64_t min_linear_index = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t max_linear_index = 0;
    double max_abs_residual_seconds = 0.0;
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

MetaInfo read_meta(const std::filesystem::path& meta_path) {
    MetaInfo meta{};
    std::ifstream in(meta_path);
    std::string line;
    while (std::getline(in, line)) {
        auto read_ll = [&](const std::string& prefix, long long* out) {
            if (line.rfind(prefix, 0) == 0) {
                *out = std::atoll(line.substr(prefix.size()).c_str());
            }
        };
        read_ll("completed_node_count=", &meta.completed_node_count);
        read_ll("total_branch_count=", &meta.total_branch_count);
        read_ll("derivative_success_count=", &meta.derivative_success_count);
    }
    return meta;
}

ValidationSummary validate_chunk(
    const std::filesystem::path& bin_path,
    const std::filesystem::path& meta_path
) {
    ValidationSummary summary{};
    if (!std::filesystem::exists(bin_path)) {
        summary.validation_error = "missing_bin";
        return summary;
    }
    if (!std::filesystem::exists(meta_path)) {
        summary.validation_error = "missing_meta";
        return summary;
    }

    const MetaInfo meta = read_meta(meta_path);
    summary.meta_completed_node_count = meta.completed_node_count;
    summary.meta_total_branch_count = meta.total_branch_count;
    summary.meta_derivative_success_count = meta.derivative_success_count;

    std::ifstream in(bin_path, std::ios::binary);
    char magic[8] = {};
    in.read(magic, sizeof(magic));
    const bool magic_ok = in && std::string(magic, sizeof(magic)) == "P1RT2DEG";
    if (!magic_ok) {
        summary.validation_error = "invalid_magic";
        return summary;
    }
    if (!read_scalar(in, &summary.version) ||
        !read_scalar(in, &summary.grid_step_deg) ||
        !read_scalar(in, &summary.samples_per_dimension) ||
        !read_scalar(in, &summary.start_node) ||
        !read_scalar(in, &summary.end_node) ||
        !read_scalar(in, &summary.node_count_header)) {
        summary.validation_error = "invalid_header";
        return summary;
    }

    std::int32_t max_transfer_revolution = 0;
    std::int32_t max_target_revolution = 0;
    summary.node_count_header = 0;
    in.clear();
    in.seekg(0, std::ios::beg);
    in.read(magic, sizeof(magic));
    read_scalar(in, &summary.version);
    read_scalar(in, &summary.grid_step_deg);
    read_scalar(in, &summary.samples_per_dimension);
    read_scalar(in, &summary.start_node);
    read_scalar(in, &summary.end_node);
    read_scalar(in, &max_transfer_revolution);
    read_scalar(in, &max_target_revolution);
    read_scalar(in, &summary.node_count_header);

    for (std::uint64_t node_index = 0; node_index < summary.node_count_header; ++node_index) {
        std::uint64_t linear_index = 0;
        std::uint32_t branch_count = 0;
        if (!read_scalar(in, &linear_index) || !read_scalar(in, &branch_count)) {
            summary.validation_error = "truncated_node_record";
            return summary;
        }
        summary.node_count_read += 1;
        summary.min_linear_index = std::min(summary.min_linear_index, linear_index);
        summary.max_linear_index = std::max(summary.max_linear_index, linear_index);
        for (std::uint32_t branch_index = 0; branch_index < branch_count; ++branch_index) {
            BranchBinaryRecord branch{};
            in.read(reinterpret_cast<char*>(&branch), sizeof(branch));
            if (!in) {
                summary.validation_error = "truncated_branch_record";
                return summary;
            }
            summary.total_branch_count_read += 1;
            if (branch.valid != 0) {
                summary.valid_branch_count_read += 1;
            }
            if (branch.derivatives_available != 0) {
                summary.derivatives_available_count_read += 1;
            }
            summary.max_abs_residual_seconds =
                std::max(summary.max_abs_residual_seconds, std::abs(branch.residual_seconds));
        }
    }

    if (summary.meta_completed_node_count != static_cast<long long>(summary.node_count_read)) {
        summary.validation_error = "meta_completed_node_count_mismatch";
        return summary;
    }
    if (summary.meta_total_branch_count != static_cast<long long>(summary.total_branch_count_read)) {
        summary.validation_error = "meta_total_branch_count_mismatch";
        return summary;
    }
    if (summary.meta_derivative_success_count != static_cast<long long>(summary.derivatives_available_count_read)) {
        summary.validation_error = "meta_derivative_success_count_mismatch";
        return summary;
    }
    summary.valid = true;
    return summary;
}

}  // namespace

int main() {
    const auto bin_path = std::filesystem::path(
        env_string("ROOT_TABLE_CHUNK_BIN", "/tmp/problem1_root_table_2deg_smoke/problem1_root_table_2deg_chunk_0_1000.bin"));
    const auto meta_path = std::filesystem::path(
        env_string("ROOT_TABLE_CHUNK_META", "/tmp/problem1_root_table_2deg_smoke/problem1_root_table_2deg_chunk_0_1000.meta.txt"));

    const auto summary = validate_chunk(bin_path, meta_path);

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "Problem1RootTable2DegChunkValidationSummary\n";
    std::cout << "valid=" << (summary.valid ? 1 : 0) << '\n';
    std::cout << "bin_path=" << bin_path.string() << '\n';
    std::cout << "meta_path=" << meta_path.string() << '\n';
    std::cout << "magic_ok=" << (summary.validation_error != "invalid_magic" ? 1 : 0) << '\n';
    std::cout << "version=" << summary.version << '\n';
    std::cout << "grid_step_deg=" << summary.grid_step_deg << '\n';
    std::cout << "samples_per_dimension=" << summary.samples_per_dimension << '\n';
    std::cout << "start_node=" << summary.start_node << '\n';
    std::cout << "end_node=" << summary.end_node << '\n';
    std::cout << "node_count_header=" << summary.node_count_header << '\n';
    std::cout << "node_count_read=" << summary.node_count_read << '\n';
    std::cout << "total_branch_count_read=" << summary.total_branch_count_read << '\n';
    std::cout << "valid_branch_count_read=" << summary.valid_branch_count_read << '\n';
    std::cout << "derivatives_available_count_read=" << summary.derivatives_available_count_read << '\n';
    std::cout << "meta_completed_node_count=" << summary.meta_completed_node_count << '\n';
    std::cout << "meta_total_branch_count=" << summary.meta_total_branch_count << '\n';
    std::cout << "meta_derivative_success_count=" << summary.meta_derivative_success_count << '\n';
    std::cout << "max_abs_residual_seconds=" << summary.max_abs_residual_seconds << '\n';
    std::cout << "validation_error=" << summary.validation_error << '\n';

    return summary.valid ? 0 : 1;
}
