#pragma once

#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

struct Problem1RootTable2DegChunkInfo {
    std::filesystem::path bin_path;
    std::filesystem::path meta_path;
    long long start_node = 0;
    long long end_node = 0;
    std::uint64_t node_count = 0;
    std::uint64_t file_size_bytes = 0;
};

struct Problem1RootTable2DegLoadProfile {
    bool offset_built_this_call = false;
    double offset_build_ms = 0.0;
    double node_read_ms = 0.0;
};

struct Problem1RootTable2DegLoadedNode {
    bool valid = false;
    std::string invalid_reason;
    long long linear_index = 0;
    int nu_A_index = 0;
    int nu_B_index = 0;
    int theta_A_index = 0;
    std::vector<Problem1SolutionBranch> branches;
    Problem1RootTable2DegLoadProfile profile;
};

struct Problem1RootTable2DegChunkRuntimeIndex {
    bool offsets_built = false;
    std::vector<std::uint64_t> node_offsets;
};

class Problem1RootTable2DegLoader {
public:
    static Problem1RootTable2DegLoader open(const std::filesystem::path& table_dir);

    Problem1RootTable2DegLoadedNode load_node_by_indices(
        int nu_A_index,
        int nu_B_index,
        int theta_A_index
    ) const;

    Problem1RootTable2DegLoadedNode load_node_by_linear_index(long long linear_index) const;

    std::size_t chunk_count() const;
    long long total_nodes() const;
    const std::vector<Problem1RootTable2DegChunkInfo>& chunks() const;

    static long long linear_index_from_indices(int nu_A_index, int nu_B_index, int theta_A_index);
    static void indices_from_linear_index(
        long long linear_index,
        int* nu_A_index,
        int* nu_B_index,
        int* theta_A_index
    );

private:
    explicit Problem1RootTable2DegLoader(std::vector<Problem1RootTable2DegChunkInfo> chunks);

    bool ensure_chunk_offsets_built(
        std::size_t chunk_index,
        std::string* invalid_reason,
        Problem1RootTable2DegLoadProfile* profile
    ) const;

    std::vector<Problem1RootTable2DegChunkInfo> chunks_;
    mutable std::vector<Problem1RootTable2DegChunkRuntimeIndex> runtime_indices_;
    long long total_nodes_ = 0;
};

}  // namespace spaceship_cpp::problem1
