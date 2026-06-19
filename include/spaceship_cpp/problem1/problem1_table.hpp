#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

enum class Problem1TransferConicType {
    Invalid,
    Elliptic,
    Parabolic,
    Hyperbolic,
};

struct Problem1TimeOfFlightBranch {
    bool valid;
    int transfer_revolution;
    int target_revolution;
    double theta_arrival_branch;
    double target_true_anomaly_start;
    double target_true_anomaly_end_branch;
    double deltaF_transfer;
    double deltaF_target;
    double time_of_flight_scale_free;
    double target_time_of_flight_scale_free;
    double time_of_flight_seconds;
    double target_time_of_flight_seconds;
    double residual_scale_free;
    double residual_seconds;
    std::string invalid_reason;
};

struct Problem1TableMetadata {
    std::string schema_version = "planet_angle_pair_table_v2";
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;
    std::string departure_planet_name;
    std::string target_planet_name;

    double departure_orbit_p;
    double departure_orbit_e;
    double departure_orbit_theta_0;
    double target_orbit_p;
    double target_orbit_e;
    double target_orbit_theta_0;

    std::string axis1_definition;
    std::string axis2_definition;
    std::string axis3_definition;
    std::string angle_convention;
    std::string derivative_method;
    std::string target_branch_validity_note;

    int departure_true_anomaly_count;
    int target_true_anomaly_count;
    int transfer_theta_departure_count;
    int max_transfer_revolution;
    int max_target_revolution;
};

struct Problem1TableConfig {
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;

    double departure_true_anomaly_start;
    double departure_true_anomaly_step;
    int departure_true_anomaly_count;

    double target_true_anomaly_start;
    double target_true_anomaly_step;
    int target_true_anomaly_count;

    double transfer_theta_departure_start;
    double transfer_theta_departure_step;
    int transfer_theta_departure_count;

    int max_transfer_revolution;
    int max_target_revolution;
};

struct Problem1TableCell {
    int departure_true_anomaly_index;
    int target_true_anomaly_index;
    int transfer_theta_departure_index;

    double departure_true_anomaly_input;
    double target_true_anomaly_input;
    double transfer_theta_departure_input;

    double departure_true_anomaly;
    double target_true_anomaly;
    double transfer_theta_departure;
    double transfer_theta_arrival;
    double launch_time_phase_seconds_since_j2000;
    double target_true_anomaly_at_departure;
    bool target_branches_are_launch_phase_specific;

    double departure_radius;
    double target_radius;
    double departure_global_angle;
    double target_global_angle;
    double delta_global_angle;

    double transfer_perihelion_angle_global_raw;
    double transfer_perihelion_angle_global_used;

    double transfer_e_raw;
    double transfer_e;
    double transfer_p;
    double transfer_a;
    bool normalized_negative_e;

    Problem1TransferConicType conic_type;
    bool valid;
    std::string invalid_reason;

    bool derivatives_available;
    double dtransfer_e_dnu_departure;
    double dtransfer_e_dnu_target;
    double dtransfer_e_dtheta_departure;

    std::vector<Problem1TimeOfFlightBranch> time_of_flight_branches;
};

struct Problem1TableQueryResult {
    Problem1TableCell cell;
    bool sourced_from_grid = false;
};

struct Problem1TableBranchQueryResult {
    Problem1TableCell cell;
    Problem1TimeOfFlightBranch branch;
    bool branch_found = false;
    bool sourced_from_grid = false;
};

struct Problem1TableVertexIndex {
    int departure_true_anomaly_index;
    int target_true_anomaly_index;
    int transfer_theta_departure_index;
};

struct Problem1TableInterpolationAdmissibility {
    bool admissible = false;
    std::string reason;
    std::array<Problem1TableVertexIndex, 8> vertex_indices{};
    std::array<Problem1TableCell, 8> vertices{};
    int transfer_revolution = 0;
    double local_departure_true_anomaly = 0.0;
    double local_target_true_anomaly = 0.0;
    double local_transfer_theta_departure = 0.0;
};

struct Problem1TransferBranchView {
    bool valid = false;
    int transfer_revolution = 0;
    double theta_arrival_branch = 0.0;
    double deltaF_transfer = 0.0;
    double time_of_flight_scale_free = 0.0;
    double time_of_flight_seconds = 0.0;
    std::string invalid_reason;
};

struct Problem1TransferTimeQueryResult {
    bool ok = false;
    bool interpolation_admissible = false;
    std::string method;
    std::string reason;
    int transfer_revolution = 0;
    double theta_arrival_branch = 0.0;
    double deltaF_transfer = 0.0;
    double time_of_flight_scale_free = 0.0;
    double time_of_flight_seconds = 0.0;
    Problem1TableInterpolationAdmissibility admissibility{};
};

void validate_problem1_table_metadata(const Problem1TableMetadata& metadata);

Problem1TableCell evaluate_problem1_table_cell_geometry(
    double departure_radius,
    double departure_global_angle,
    double target_radius,
    double target_global_angle,
    double transfer_theta_departure_input,
    int max_transfer_revolution
);

Problem1TableCell evaluate_problem1_table_cell_for_planets(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double departure_true_anomaly_input,
    double target_true_anomaly_input,
    double transfer_theta_departure_input,
    int max_transfer_revolution,
    int max_target_revolution
);

class Problem1Table {
public:
    // 当前实现是 endpoint-table / transfer-time-table 方向的实验性中间产物。
    // 它只保留几何和飞行时间分支查询，不再承担预计算根表职责。
    explicit Problem1Table(Problem1TableConfig config);

    const Problem1TableConfig& config() const;
    const Problem1TableMetadata& metadata() const;

    int departure_true_anomaly_count() const;
    int target_true_anomaly_count() const;
    int transfer_theta_departure_count() const;

    std::size_t flat_index(
        int departure_true_anomaly_index,
        int target_true_anomaly_index,
        int transfer_theta_departure_index
    ) const;

    Problem1TableCell& at(
        int departure_true_anomaly_index,
        int target_true_anomaly_index,
        int transfer_theta_departure_index
    );
    const Problem1TableCell& at(
        int departure_true_anomaly_index,
        int target_true_anomaly_index,
        int transfer_theta_departure_index
    ) const;

    const std::vector<Problem1TableCell>& cells() const;

private:
    Problem1TableConfig config_;
    Problem1TableMetadata metadata_;
    std::vector<Problem1TableCell> cells_;
};

Problem1Table build_problem1_table(const Problem1TableConfig& config);

Problem1TableQueryResult query_problem1_table_exact(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure
);

Problem1TableBranchQueryResult query_problem1_table_exact_branch(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution,
    int target_revolution
);  // 仅返回 representative / launch-phase-specific q branch；不应用于未来插值准入检查。

Problem1TableBranchQueryResult query_problem1_table_exact_branch_at_departure_time(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    double departure_time_seconds_since_j2000,
    int transfer_revolution,
    int target_revolution
);

Problem1TimeOfFlightBranch evaluate_problem1_table_branch_with_target_departure_true_anomaly(
    const Problem1Table& table,
    const Problem1TableCell& cell,
    double target_true_anomaly_at_departure,
    int transfer_revolution,
    int target_revolution
);

const Problem1TimeOfFlightBranch* find_problem1_table_branch(
    const Problem1TableCell& cell,
    int transfer_revolution,
    int target_revolution
);

std::string problem1_table_branch_signature(const Problem1TableCell& cell);

bool is_problem1_transfer_branch_valid_for_interpolation(
    const Problem1TimeOfFlightBranch& branch,
    int transfer_revolution
);

Problem1TransferBranchView get_problem1_transfer_branch_view(
    const Problem1TableCell& cell,
    int transfer_revolution
);

Problem1TransferTimeQueryResult query_problem1_transfer_time_with_interpolation_stub(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution
);  // 基于旧 endpoint-table 的 future interpolator stub；不代表最终的 root-solution-table 方向。

Problem1TableInterpolationAdmissibility check_problem1_table_transfer_branch_interpolation_admissibility(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution
);  // 仅适用于当前 endpoint-table / transfer-time-table 设计；未来 root-solution-table 需要新的 branch matching 准入规则。

}  // namespace spaceship_cpp::problem1
