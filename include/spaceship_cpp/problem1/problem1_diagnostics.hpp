#pragma once

#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <map>
#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

struct Problem1SolveSummary {
    int candidate_count;
    double min_time_of_flight_seconds;
    double max_time_of_flight_seconds;
    double min_relative_residual;
    double max_relative_residual;
};

Problem1SolveSummary summarize_problem1_candidates(
    const std::vector<Problem1Candidate>& candidates
);

void write_problem1_candidates_csv(
    const std::vector<Problem1Candidate>& candidates,
    const std::string& output_path
);

struct Problem1TableSummary {
    int cell_count;
    int valid_geometry_cell_count;
    int total_branch_count;
    int valid_branch_count;

    double min_departure_radius;
    double max_departure_radius;
    double min_target_radius;
    double max_target_radius;

    double min_time_of_flight_seconds;
    double max_time_of_flight_seconds;
    double min_time_of_flight_days;
    double max_time_of_flight_days;

    double min_transfer_e;
    double max_transfer_e;
    double min_transfer_p;
    double max_transfer_p;
};

struct Problem1TableBranchDiagnostics {
    int cell_count;
    int valid_cell_count;
    int invalid_cell_count;
    int min_transfer_revolution;
    int max_transfer_revolution;
    int min_target_revolution;
    int max_target_revolution;
    std::map<std::string, int> invalid_reason_counts;
    std::map<std::string, int> conic_type_counts;
    std::map<int, int> branch_count_distribution;
    std::map<int, int> branch_count_by_transfer_revolution;
    std::map<int, int> branch_count_by_target_revolution;
    std::map<std::string, int> valid_branch_count_by_pair;
    std::map<std::string, int> invalid_branch_count_by_pair;
    std::map<std::string, int> invalid_reason_counts_by_pair;
};

Problem1TableSummary summarize_problem1_table(const Problem1Table& table);

Problem1TableBranchDiagnostics summarize_problem1_table_branches(const Problem1Table& table);

void write_problem1_table_csv(
    const Problem1Table& table,
    const std::string& output_path
);

void write_problem1_table_branch_summary_csv(
    const Problem1Table& table,
    const std::string& output_path
);

}  // namespace spaceship_cpp::problem1
