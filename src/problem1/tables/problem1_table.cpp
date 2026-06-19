/*
 * 文件作用：实现 Problem 1 endpoint transfer-time table。
 * 主要工作：生成表格单元几何、计算飞行时间分支，并提供最近网格和插值查询。
 */
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/common/orbit_math.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <limits>
#include <stdexcept>

namespace spaceship_cpp::problem1 {

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kDefaultEpsilon;
using spaceship_cpp::common::kPi;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;
using spaceship_cpp::common::orbit_F;

// 返回 quiet_NaN 哨兵值。
double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

// 构造一个标记为 invalid 的空单元，所有数值字段填 NaN。
Problem1TableCell make_invalid_cell(const std::string& reason) {
    const double nan = nan_value();
    Problem1TableCell cell{};
    cell.departure_true_anomaly_index = 0;
    cell.target_true_anomaly_index = 0;
    cell.transfer_theta_departure_index = 0;
    cell.departure_true_anomaly_input = nan;
    cell.target_true_anomaly_input = nan;
    cell.transfer_theta_departure_input = nan;
    cell.departure_true_anomaly = nan;
    cell.target_true_anomaly = nan;
    cell.transfer_theta_departure = nan;
    cell.transfer_theta_arrival = nan;
    cell.launch_time_phase_seconds_since_j2000 = nan;
    cell.target_true_anomaly_at_departure = nan;
    cell.target_branches_are_launch_phase_specific = true;
    cell.departure_radius = nan;
    cell.target_radius = nan;
    cell.departure_global_angle = nan;
    cell.target_global_angle = nan;
    cell.delta_global_angle = nan;
    cell.transfer_perihelion_angle_global_raw = nan;
    cell.transfer_perihelion_angle_global_used = nan;
    cell.transfer_e_raw = nan;
    cell.transfer_e = nan;
    cell.transfer_p = nan;
    cell.transfer_a = nan;
    cell.normalized_negative_e = false;
    cell.conic_type = Problem1TransferConicType::Invalid;
    cell.valid = false;
    cell.invalid_reason = reason;
    cell.derivatives_available = false;
    cell.dtransfer_e_dnu_departure = nan;
    cell.dtransfer_e_dnu_target = nan;
    cell.dtransfer_e_dtheta_departure = nan;
    return cell;
}

// 检查 1+e·cos(angle) > 0，确保轨道半径公式分母为正。
bool positive_orbit_denominator(double e, double angle) {
    const double denominator = 1.0 + e * std::cos(angle);
    return is_finite(denominator) && denominator > kDefaultEpsilon;
}

// 校验周期轴步长为正；步长非法时抛异常。
double positive_mod_step(double start, double step, int count) {
    if (!(count > 0) || !is_finite(step) || !(step > 0.0)) {
        throw std::invalid_argument("invalid periodic axis configuration");
    }
    return step;
}

// 校验 Problem1TableConfig 三维轴参数和圈数上限的合法性。
void validate_problem1_table_config(const Problem1TableConfig& config) {
    if (config.departure_true_anomaly_count <= 0 ||
        config.target_true_anomaly_count <= 0 ||
        config.transfer_theta_departure_count <= 0 ||
        !is_finite(config.departure_true_anomaly_start) ||
        !is_finite(config.target_true_anomaly_start) ||
        !is_finite(config.transfer_theta_departure_start) ||
        !is_finite(config.departure_true_anomaly_step) ||
        !is_finite(config.target_true_anomaly_step) ||
        !is_finite(config.transfer_theta_departure_step) ||
        !(config.departure_true_anomaly_step > 0.0) ||
        !(config.target_true_anomaly_step > 0.0) ||
        !(config.transfer_theta_departure_step > 0.0) ||
        config.max_transfer_revolution < 0 ||
        config.max_target_revolution < 0) {
        throw std::invalid_argument("invalid Problem1TableConfig");
    }
}

// 校验后返回配置副本（用于构造时确保参数合法）。
Problem1TableConfig validated_problem1_table_config(Problem1TableConfig config) {
    validate_problem1_table_config(config);
    return config;
}

// 计算表格总单元数 = nu_A × nu_B × theta_A。
std::size_t problem1_table_cell_count(const Problem1TableConfig& config) {
    return static_cast<std::size_t>(config.departure_true_anomaly_count) *
        static_cast<std::size_t>(config.target_true_anomaly_count) *
        static_cast<std::size_t>(config.transfer_theta_departure_count);
}

// 由轴起点、步长和索引计算该网格点的角度值（归一化到 [0,2π)）。
double axis_value(double start, double step, int index) {
    return normalize_angle_0_2pi(start + static_cast<double>(index) * step);
}

struct PeriodicAxisLocation {
    int lower_index;
    int upper_index;
    double lower_unwrapped;
    double upper_unwrapped;
    double query_unwrapped;
};

// 在周期轴上定位查询角所在的 [lower, upper] 网格区间；
// 解决周期边界（0↔2π）处插值需要正确包裹的问题。
PeriodicAxisLocation locate_periodic_axis(
    double start,
    double step,
    int count,
    double query_angle
) {
    if (count <= 0 || !is_finite(start) || !is_finite(step) || !(step > 0.0) || !is_finite(query_angle)) {
        throw std::invalid_argument("invalid periodic axis locator inputs");
    }

    const double query_normalized = normalize_angle_0_2pi(query_angle);
    const double start_normalized = normalize_angle_0_2pi(start);
    const double offset = spaceship_cpp::common::positive_mod(query_normalized - start_normalized, kTwoPi);
    const double position = offset / step;
    const int lower_index = static_cast<int>(std::floor(position));
    const int wrapped_lower_index = ((lower_index % count) + count) % count;
    const int wrapped_upper_index = (wrapped_lower_index + 1) % count;
    const double lower_unwrapped = start + static_cast<double>(lower_index) * step;
    const double upper_unwrapped = lower_unwrapped + step;

    return PeriodicAxisLocation{
        wrapped_lower_index,
        wrapped_upper_index,
        lower_unwrapped,
        upper_unwrapped,
        start + position * step,
    };
}

// 根据偏心率分类转移圆锥曲线类型。
Problem1TransferConicType classify_transfer_conic(double e) {
    if (!is_finite(e) || e < 0.0) {
        return Problem1TransferConicType::Invalid;
    }
    if (std::abs(e - 1.0) <= 1e-12) {
        return Problem1TransferConicType::Parabolic;
    }
    if (e < 1.0) {
        return Problem1TransferConicType::Elliptic;
    }
    return Problem1TransferConicType::Hyperbolic;
}

double reduced_time_from_true_anomaly(
    planet_params::PlanetId planet_id,
    double true_anomaly
) {
    // 中文注释：由出发行星真近点角反推一个“模本轨道周期”的发射相位时间，用于恢复目标行星在同一时刻的起始相位。
    const planet_params::PlanetParams& planet = planet_params::get_planet_params(planet_id);
    const double mu = planet_params::get_solar_system_physical_params().GM_sun;
    const double orbit_scale = std::sqrt(mu / std::pow(planet.orbit.p, 3.0));
    const double raw_time =
        (orbit_F(planet.orbit.e, true_anomaly) - orbit_F(planet.orbit.e, planet.orbit.varphi_0)) / orbit_scale;
    return std::remainder(raw_time, planet_params::planet_orbital_period(planet_id));
}

// 评估单个 (k,q) 分支的飞行时间和残差；是表格分支枚举的核心单元。
Problem1TimeOfFlightBranch evaluate_problem1_table_branch(
    double transfer_e,
    double transfer_theta_departure,
    double transfer_theta_arrival,
    double transfer_p,
    planet_params::PlanetId target_planet,
    double target_true_anomaly_start,
    double target_true_anomaly_end_base,
    int transfer_revolution,
    int target_revolution
) {
    const double mu = planet_params::get_solar_system_physical_params().GM_sun;
    const planet_params::PlanetParams& target = planet_params::get_planet_params(target_planet);
    Problem1TimeOfFlightBranch branch{};
    branch.valid = false;
    branch.transfer_revolution = transfer_revolution;
    branch.target_revolution = target_revolution;
    branch.theta_arrival_branch = nan_value();
    branch.target_true_anomaly_start = target_true_anomaly_start;
    branch.target_true_anomaly_end_branch = nan_value();
    branch.deltaF_transfer = nan_value();
    branch.deltaF_target = nan_value();
    branch.time_of_flight_scale_free = nan_value();
    branch.target_time_of_flight_scale_free = nan_value();
    branch.time_of_flight_seconds = nan_value();
    branch.target_time_of_flight_seconds = nan_value();
    branch.residual_scale_free = nan_value();
    branch.residual_seconds = nan_value();

    const Problem1TransferConicType conic_type = classify_transfer_conic(transfer_e);
    if (conic_type == Problem1TransferConicType::Invalid ||
        conic_type == Problem1TransferConicType::Parabolic ||
        !is_finite(transfer_p) || !(transfer_p > 0.0) ||
        !is_finite(target_true_anomaly_start) ||
        !is_finite(target_true_anomaly_end_base)) {
        branch.invalid_reason = "invalid branch inputs";
        return branch;
    }

    double theta_departure_forward = transfer_theta_departure;
    double theta_arrival_forward = transfer_theta_arrival;
    if (conic_type == Problem1TransferConicType::Elliptic) {
        while (theta_arrival_forward <= theta_departure_forward) {
            theta_arrival_forward += kTwoPi;
        }
    } else {
        theta_departure_forward = normalize_angle_minus_pi_pi(transfer_theta_departure);
        theta_arrival_forward = normalize_angle_minus_pi_pi(transfer_theta_arrival);
    }

    double target_theta_end_base = target_true_anomaly_end_base;
    while (target_theta_end_base <= target_true_anomaly_start) {
        target_theta_end_base += kTwoPi;
    }

    if (conic_type == Problem1TransferConicType::Hyperbolic && transfer_revolution != 0) {
        // 中文注释：双曲线没有额外绕日整圈的物理分支，但仍保留这个 (k,q) 记录，便于完整诊断覆盖率。
        branch.invalid_reason = "hyperbolic transfer does not support k > 0";
        return branch;
    }

    if (conic_type == Problem1TransferConicType::Hyperbolic &&
        theta_arrival_forward <= theta_departure_forward) {
        branch.invalid_reason = "hyperbolic arrival angle is not forward of departure angle";
        branch.theta_arrival_branch = theta_arrival_forward;
        return branch;
    }

    branch.theta_arrival_branch =
        theta_arrival_forward + static_cast<double>(transfer_revolution) * kTwoPi;
    branch.target_true_anomaly_end_branch =
        target_theta_end_base + static_cast<double>(target_revolution) * kTwoPi;

    try {
        branch.deltaF_transfer =
            orbit_F(transfer_e, branch.theta_arrival_branch) - orbit_F(transfer_e, theta_departure_forward);
        branch.deltaF_target =
            orbit_F(target.orbit.e, branch.target_true_anomaly_end_branch) -
            orbit_F(target.orbit.e, target_true_anomaly_start);
        branch.time_of_flight_scale_free = std::pow(transfer_p, 1.5) * branch.deltaF_transfer;
        branch.target_time_of_flight_scale_free =
            std::pow(target.orbit.p, 1.5) * branch.deltaF_target;
        branch.time_of_flight_seconds = branch.time_of_flight_scale_free / std::sqrt(mu);
        branch.target_time_of_flight_seconds =
            branch.target_time_of_flight_scale_free / std::sqrt(mu);
        branch.residual_scale_free =
            branch.time_of_flight_scale_free - branch.target_time_of_flight_scale_free;
        branch.residual_seconds =
            branch.time_of_flight_seconds - branch.target_time_of_flight_seconds;
        branch.valid =
            is_finite(branch.deltaF_transfer) && branch.deltaF_transfer > 0.0 &&
            is_finite(branch.deltaF_target) && branch.deltaF_target > 0.0 &&
            is_finite(branch.time_of_flight_seconds) && branch.time_of_flight_seconds > 0.0 &&
            is_finite(branch.target_time_of_flight_seconds) && branch.target_time_of_flight_seconds > 0.0 &&
            is_finite(branch.residual_scale_free) && is_finite(branch.residual_seconds);
        if (!branch.valid) {
            branch.invalid_reason = "non-positive or non-finite branch time of flight";
        }
    } catch (const std::domain_error& ex) {
        branch.invalid_reason = ex.what();
    }

    return branch;
}

// 枚举所有 (k,q) 组合并计算各自的时间分支列表。
std::vector<Problem1TimeOfFlightBranch> build_problem1_table_branches(
    double transfer_e,
    double transfer_theta_departure,
    double transfer_theta_arrival,
    double transfer_p,
    planet_params::PlanetId target_planet,
    double target_true_anomaly_start,
    double target_true_anomaly_end_base,
    int max_transfer_revolution,
    int max_target_revolution
) {
    std::vector<Problem1TimeOfFlightBranch> branches;
    for (int transfer_revolution = 0; transfer_revolution <= max_transfer_revolution; ++transfer_revolution) {
        for (int target_revolution = 0; target_revolution <= max_target_revolution; ++target_revolution) {
            branches.push_back(evaluate_problem1_table_branch(
                transfer_e,
                transfer_theta_departure,
                transfer_theta_arrival,
                transfer_p,
                target_planet,
                target_true_anomaly_start,
                target_true_anomaly_end_base,
                transfer_revolution,
                target_revolution));
        }
    }
    return branches;
}

// 在单元分支列表中查找指定 k 的插值有效分支指针。
const Problem1TimeOfFlightBranch* find_valid_transfer_branch_by_k(
    const Problem1TableCell& cell,
    int transfer_revolution
) {
    const Problem1TransferBranchView view = get_problem1_transfer_branch_view(cell, transfer_revolution);
    if (!view.valid) {
        return nullptr;
    }
    for (const Problem1TimeOfFlightBranch& branch : cell.time_of_flight_branches) {
        if (is_problem1_transfer_branch_valid_for_interpolation(branch, transfer_revolution)) {
            return &branch;
        }
    }
    return nullptr;
}

// 由表格配置生成元数据（轴定义、行星参数快照、schema 版本等）。
Problem1TableMetadata make_problem1_table_metadata(const Problem1TableConfig& config) {
    const planet_params::PlanetParams& departure =
        planet_params::get_planet_params(config.departure_planet);
    const planet_params::PlanetParams& target =
        planet_params::get_planet_params(config.target_planet);
    return Problem1TableMetadata{
        "planet_angle_pair_table_v2",
        config.departure_planet,
        config.target_planet,
        departure.name,
        target.name,
        departure.orbit.p,
        departure.orbit.e,
        departure.orbit.theta_0,
        target.orbit.p,
        target.orbit.e,
        target.orbit.theta_0,
        "nu_A: departure planet true anomaly in [0, 2pi)",
        "nu_B: target planet true anomaly in [0, 2pi)",
        "theta_A: transfer true anomaly at departure relative to transfer perihelion",
        "theta_B = theta_A + wrap_to_2pi(lambda_B - lambda_A)",
        "TODO: finite-difference / analytic derivatives not implemented in C++ table v2",
        "Stored q branches are launch-phase-specific: target_true_anomaly_at_departure is reconstructed from "
        "nu_A modulo T_A, so q is not universally valid without a real departure time or nu_B_depart input.",
        config.departure_true_anomaly_count,
        config.target_true_anomaly_count,
        config.transfer_theta_departure_count,
        config.max_transfer_revolution,
        config.max_target_revolution,
    };
}

}  // namespace

// 校验表格元数据 schema 版本是否为当前代码支持的 v2。
void validate_problem1_table_metadata(const Problem1TableMetadata& metadata) {
    if (metadata.schema_version != "planet_angle_pair_table_v2") {
        throw std::invalid_argument(
            "unsupported Problem1Table metadata schema; expected planet_angle_pair_table_v2");
    }
}

// 仅从两端点太阳中心几何构造表格单元（不含行星相位和时间分支）；
// 隔离纯几何求值，便于单元测试。
Problem1TableCell evaluate_problem1_table_cell_geometry(
    double departure_radius,
    double departure_global_angle,
    double target_radius,
    double target_global_angle,
    double transfer_theta_departure_input,
    int max_transfer_revolution
) {
    Problem1TableCell cell = make_invalid_cell("");
    cell.departure_radius = departure_radius;
    cell.target_radius = target_radius;
    cell.departure_global_angle = normalize_angle_0_2pi(departure_global_angle);
    cell.target_global_angle = normalize_angle_0_2pi(target_global_angle);
    cell.departure_true_anomaly_input = nan_value();
    cell.target_true_anomaly_input = nan_value();
    cell.transfer_theta_departure_input = normalize_angle_0_2pi(transfer_theta_departure_input);
    cell.launch_time_phase_seconds_since_j2000 = nan_value();
    cell.target_true_anomaly_at_departure = nan_value();
    cell.target_branches_are_launch_phase_specific = true;
    cell.derivatives_available = false;
    cell.dtransfer_e_dnu_departure = nan_value();
    cell.dtransfer_e_dnu_target = nan_value();
    cell.dtransfer_e_dtheta_departure = nan_value();

    if (!is_finite(departure_radius) || !(departure_radius > 0.0) ||
        !is_finite(target_radius) || !(target_radius > 0.0) ||
        !is_finite(departure_global_angle) ||
        !is_finite(target_global_angle) ||
        !is_finite(transfer_theta_departure_input) ||
        max_transfer_revolution < 0) {
        cell.invalid_reason = "invalid geometry inputs";
        return cell;
    }

    cell.delta_global_angle = normalize_angle_0_2pi(cell.target_global_angle - cell.departure_global_angle);
    cell.transfer_theta_arrival =
        normalize_angle_0_2pi(cell.transfer_theta_departure_input + cell.delta_global_angle);

    const double raw_perihelion_global =
        normalize_angle_0_2pi(cell.departure_global_angle - cell.transfer_theta_departure_input);
    cell.transfer_perihelion_angle_global_raw = raw_perihelion_global;

    const double denominator =
        departure_radius * std::cos(cell.transfer_theta_departure_input) -
        target_radius * std::cos(cell.transfer_theta_arrival);
    if (!is_finite(denominator) || std::abs(denominator) <= kDefaultEpsilon) {
        cell.invalid_reason = "singular geometry denominator";
        return cell;
    }

    const double transfer_e_raw = (target_radius - departure_radius) / denominator;
    double transfer_e = transfer_e_raw;
    cell.normalized_negative_e = false;
    cell.transfer_perihelion_angle_global_used = raw_perihelion_global;
    cell.transfer_theta_departure = cell.transfer_theta_departure_input;
    double transfer_theta_arrival = cell.transfer_theta_arrival;
    if (transfer_e < 0.0) {
        transfer_e = -transfer_e;
        cell.normalized_negative_e = true;
        cell.transfer_perihelion_angle_global_used = normalize_angle_0_2pi(raw_perihelion_global + kPi);
        cell.transfer_theta_departure = normalize_angle_0_2pi(cell.transfer_theta_departure_input - kPi);
        transfer_theta_arrival = normalize_angle_0_2pi(cell.transfer_theta_arrival - kPi);
    }
    cell.transfer_theta_arrival = transfer_theta_arrival;
    cell.transfer_e_raw = transfer_e_raw;
    cell.transfer_e = transfer_e;

    const double orbit_factor = 1.0 + transfer_e * std::cos(cell.transfer_theta_departure);
    if (!(is_finite(orbit_factor) && orbit_factor > kDefaultEpsilon)) {
        cell.invalid_reason = "non-positive departure orbit factor";
        return cell;
    }
    cell.transfer_p = departure_radius * orbit_factor;
    if (!(is_finite(cell.transfer_p) && cell.transfer_p > 0.0) ||
        !positive_orbit_denominator(transfer_e, cell.transfer_theta_departure) ||
        !positive_orbit_denominator(transfer_e, cell.transfer_theta_arrival)) {
        cell.invalid_reason = "invalid transfer conic geometry";
        return cell;
    }

    cell.conic_type = classify_transfer_conic(transfer_e);
    if (cell.conic_type == Problem1TransferConicType::Parabolic) {
        cell.transfer_a = nan_value();
        cell.invalid_reason = "parabolic transfer is unsupported in table v2";
        cell.time_of_flight_branches.clear();
        cell.valid = false;
        return cell;
    }
    if (cell.conic_type == Problem1TransferConicType::Invalid) {
        cell.transfer_a = nan_value();
        cell.invalid_reason = "invalid transfer eccentricity";
        return cell;
    }

    const double one_minus_e2 = 1.0 - transfer_e * transfer_e;
    if (std::abs(one_minus_e2) <= kDefaultEpsilon) {
        cell.transfer_a = nan_value();
    } else {
        cell.transfer_a = cell.transfer_p / one_minus_e2;
    }

    // 中文注释：纯几何求值阶段不展开 (k,q) 时间分支；完整 branch 列表在带行星上下文的 evaluate_problem1_table_cell_for_planets 中生成。
    cell.time_of_flight_branches.clear();
    cell.valid = true;
    cell.invalid_reason.clear();
    return cell;
}

// 从行星局部真近点角 (nu_A, nu_B, theta_A) 生成完整表格单元（含时间分支）；
// 解决"该行星相位下表格值是什么"的问题。
Problem1TableCell evaluate_problem1_table_cell_for_planets(
    planet_params::PlanetId departure_planet,
    planet_params::PlanetId target_planet,
    double departure_true_anomaly_input,
    double target_true_anomaly_input,
    double transfer_theta_departure_input,
    int max_transfer_revolution,
    int max_target_revolution
) {
    const planet_params::PlanetParams& departure = planet_params::get_planet_params(departure_planet);
    const planet_params::PlanetParams& target = planet_params::get_planet_params(target_planet);

    const double departure_true_anomaly = normalize_angle_0_2pi(departure_true_anomaly_input);
    const double target_true_anomaly = normalize_angle_0_2pi(target_true_anomaly_input);
    const double departure_radius =
        planet_params::planet_radius_at_true_anomaly(departure_planet, departure_true_anomaly);
    const double target_radius =
        planet_params::planet_radius_at_true_anomaly(target_planet, target_true_anomaly);
    const double departure_global_angle =
        normalize_angle_0_2pi(departure.orbit.theta_0 + departure_true_anomaly);
    const double target_global_angle =
        normalize_angle_0_2pi(target.orbit.theta_0 + target_true_anomaly);

    Problem1TableCell cell = evaluate_problem1_table_cell_geometry(
        departure_radius,
        departure_global_angle,
        target_radius,
        target_global_angle,
        transfer_theta_departure_input,
        max_transfer_revolution);
    cell.departure_true_anomaly_input = departure_true_anomaly_input;
    cell.target_true_anomaly_input = target_true_anomaly_input;
    cell.transfer_theta_departure_input = normalize_angle_0_2pi(transfer_theta_departure_input);
    cell.departure_true_anomaly = departure_true_anomaly;
    cell.target_true_anomaly = target_true_anomaly;
    cell.launch_time_phase_seconds_since_j2000 =
        reduced_time_from_true_anomaly(departure_planet, departure_true_anomaly);
    cell.target_true_anomaly_at_departure =
        planet_params::planet_true_anomaly_at_time(target_planet, cell.launch_time_phase_seconds_since_j2000);
    cell.target_branches_are_launch_phase_specific = true;
    if (cell.valid) {
        // 中文注释：几何由 (nu_A, nu_B, theta_A) 唯一决定；这里存下来的 q 分支只对“由 nu_A 模 T_A 反推出的代表性发射相位”有效。
        cell.time_of_flight_branches = build_problem1_table_branches(
            cell.transfer_e,
            cell.transfer_theta_departure,
            cell.transfer_theta_arrival,
            cell.transfer_p,
            target_planet,
            cell.target_true_anomaly_at_departure,
            cell.target_true_anomaly,
            max_transfer_revolution,
            max_target_revolution);
    } else {
        cell.time_of_flight_branches.clear();
    }
    return cell;
}

Problem1Table::Problem1Table(Problem1TableConfig config)
    : config_(validated_problem1_table_config(config)),
      metadata_(make_problem1_table_metadata(config_)),
      cells_(problem1_table_cell_count(config_)) {
    validate_problem1_table_metadata(metadata_);
}

const Problem1TableConfig& Problem1Table::config() const {
    return config_;
}

const Problem1TableMetadata& Problem1Table::metadata() const {
    return metadata_;
}

int Problem1Table::departure_true_anomaly_count() const {
    return config_.departure_true_anomaly_count;
}

int Problem1Table::target_true_anomaly_count() const {
    return config_.target_true_anomaly_count;
}

int Problem1Table::transfer_theta_departure_count() const {
    return config_.transfer_theta_departure_count;
}

std::size_t Problem1Table::flat_index(
    int departure_true_anomaly_index,
    int target_true_anomaly_index,
    int transfer_theta_departure_index
) const {
    if (departure_true_anomaly_index < 0 ||
        departure_true_anomaly_index >= config_.departure_true_anomaly_count ||
        target_true_anomaly_index < 0 ||
        target_true_anomaly_index >= config_.target_true_anomaly_count ||
        transfer_theta_departure_index < 0 ||
        transfer_theta_departure_index >= config_.transfer_theta_departure_count) {
        throw std::out_of_range("Problem1Table index out of range");
    }
    return (
        static_cast<std::size_t>(departure_true_anomaly_index) *
            static_cast<std::size_t>(config_.target_true_anomaly_count) +
        static_cast<std::size_t>(target_true_anomaly_index)
    ) * static_cast<std::size_t>(config_.transfer_theta_departure_count) +
        static_cast<std::size_t>(transfer_theta_departure_index);
}

Problem1TableCell& Problem1Table::at(
    int departure_true_anomaly_index,
    int target_true_anomaly_index,
    int transfer_theta_departure_index
) {
    return cells_.at(flat_index(
        departure_true_anomaly_index,
        target_true_anomaly_index,
        transfer_theta_departure_index));
}

const Problem1TableCell& Problem1Table::at(
    int departure_true_anomaly_index,
    int target_true_anomaly_index,
    int transfer_theta_departure_index
) const {
    return cells_.at(flat_index(
        departure_true_anomaly_index,
        target_true_anomaly_index,
        transfer_theta_departure_index));
}

const std::vector<Problem1TableCell>& Problem1Table::cells() const {
    return cells_;
}

// 遍历三维网格填充所有单元，构建完整 endpoint transfer-time 表格。
Problem1Table build_problem1_table(const Problem1TableConfig& config) {
    validate_problem1_table_config(config);

    Problem1Table table(config);
    for (int departure_index = 0; departure_index < config.departure_true_anomaly_count; ++departure_index) {
        const double departure_true_anomaly = axis_value(
            config.departure_true_anomaly_start,
            config.departure_true_anomaly_step,
            departure_index);
        for (int target_index = 0; target_index < config.target_true_anomaly_count; ++target_index) {
            const double target_true_anomaly = axis_value(
                config.target_true_anomaly_start,
                config.target_true_anomaly_step,
                target_index);
            for (int theta_index = 0; theta_index < config.transfer_theta_departure_count; ++theta_index) {
                const double transfer_theta_departure = axis_value(
                    config.transfer_theta_departure_start,
                    config.transfer_theta_departure_step,
                    theta_index);
                Problem1TableCell cell = evaluate_problem1_table_cell_for_planets(
                    config.departure_planet,
                    config.target_planet,
                    departure_true_anomaly,
                    target_true_anomaly,
                    transfer_theta_departure,
                    config.max_transfer_revolution,
                    config.max_target_revolution);
                cell.departure_true_anomaly_index = departure_index;
                cell.target_true_anomaly_index = target_index;
                cell.transfer_theta_departure_index = theta_index;
                table.at(departure_index, target_index, theta_index) = std::move(cell);
            }
        }
    }
    return table;
}

// 即时计算（非查表）指定角度的表格单元；sourced_from_grid=false。
Problem1TableQueryResult query_problem1_table_exact(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure
) {
    const Problem1TableConfig& config = table.config();
    Problem1TableCell cell = evaluate_problem1_table_cell_for_planets(
        config.departure_planet,
        config.target_planet,
        normalize_angle_0_2pi(departure_true_anomaly),
        normalize_angle_0_2pi(target_true_anomaly),
        normalize_angle_0_2pi(transfer_theta_departure),
        config.max_transfer_revolution,
        config.max_target_revolution);
    cell.departure_true_anomaly_index = -1;
    cell.target_true_anomaly_index = -1;
    cell.transfer_theta_departure_index = -1;
    return Problem1TableQueryResult{cell, false};
}

Problem1TimeOfFlightBranch evaluate_problem1_table_branch_with_target_departure_true_anomaly(
    const Problem1Table& table,
    const Problem1TableCell& cell,
    double target_true_anomaly_at_departure,
    int transfer_revolution,
    int target_revolution
) {
    if (!cell.valid) {
        Problem1TimeOfFlightBranch branch{};
        branch.valid = false;
        branch.transfer_revolution = transfer_revolution;
        branch.target_revolution = target_revolution;
        branch.theta_arrival_branch = nan_value();
        branch.target_true_anomaly_start = target_true_anomaly_at_departure;
        branch.target_true_anomaly_end_branch = nan_value();
        branch.deltaF_transfer = nan_value();
        branch.deltaF_target = nan_value();
        branch.time_of_flight_scale_free = nan_value();
        branch.target_time_of_flight_scale_free = nan_value();
        branch.time_of_flight_seconds = nan_value();
        branch.target_time_of_flight_seconds = nan_value();
        branch.residual_scale_free = nan_value();
        branch.residual_seconds = nan_value();
        branch.invalid_reason = "invalid geometry cell";
        return branch;
    }

    // 中文注释：这个接口显式要求真实的目标行星发射时刻 anomaly，因此得到的 q 分支才具有对应 launch time 的物理意义。
    return evaluate_problem1_table_branch(
        cell.transfer_e,
        cell.transfer_theta_departure,
        cell.transfer_theta_arrival,
        cell.transfer_p,
        table.config().target_planet,
        normalize_angle_0_2pi(target_true_anomaly_at_departure),
        cell.target_true_anomaly,
        transfer_revolution,
        target_revolution);
}

Problem1TableBranchQueryResult query_problem1_table_exact_branch(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution,
    int target_revolution
) {
    // 中文注释：这个接口返回的是 representative / launch-phase-specific q branch，
    // 因为它内部仍使用了由 nu_A 模 T_A 反推出的目标行星出发 anomaly；
    // future interpolation admissibility 不应依赖这个接口。
    Problem1TableQueryResult query = query_problem1_table_exact(
        table,
        departure_true_anomaly,
        target_true_anomaly,
        transfer_theta_departure);
    const Problem1TimeOfFlightBranch* branch = find_problem1_table_branch(
        query.cell,
        transfer_revolution,
        target_revolution);
    if (branch == nullptr) {
        if (query.cell.valid) {
            const Problem1TableConfig& config = table.config();
            const Problem1TimeOfFlightBranch direct_branch = evaluate_problem1_table_branch(
                query.cell.transfer_e,
                query.cell.transfer_theta_departure,
                query.cell.transfer_theta_arrival,
                query.cell.transfer_p,
                config.target_planet,
                query.cell.target_true_anomaly_at_departure,
                query.cell.target_true_anomaly,
                transfer_revolution,
                target_revolution);
            return Problem1TableBranchQueryResult{std::move(query.cell), direct_branch, true, false};
        }
        Problem1TimeOfFlightBranch invalid_branch{};
        invalid_branch.valid = false;
        invalid_branch.transfer_revolution = transfer_revolution;
        invalid_branch.target_revolution = target_revolution;
        invalid_branch.theta_arrival_branch = nan_value();
        invalid_branch.target_true_anomaly_start = nan_value();
        invalid_branch.target_true_anomaly_end_branch = nan_value();
        invalid_branch.deltaF_transfer = nan_value();
        invalid_branch.deltaF_target = nan_value();
        invalid_branch.time_of_flight_scale_free = nan_value();
        invalid_branch.target_time_of_flight_scale_free = nan_value();
        invalid_branch.time_of_flight_seconds = nan_value();
        invalid_branch.target_time_of_flight_seconds = nan_value();
        invalid_branch.residual_scale_free = nan_value();
        invalid_branch.residual_seconds = nan_value();
        invalid_branch.invalid_reason = "requested (k,q) branch is not present in this cell";
        return Problem1TableBranchQueryResult{std::move(query.cell), invalid_branch, false, false};
    }
    return Problem1TableBranchQueryResult{std::move(query.cell), *branch, true, false};
}

Problem1TableBranchQueryResult query_problem1_table_exact_branch_at_departure_time(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    double departure_time_seconds_since_j2000,
    int transfer_revolution,
    int target_revolution
) {
    Problem1TableQueryResult query = query_problem1_table_exact(
        table,
        departure_true_anomaly,
        target_true_anomaly,
        transfer_theta_departure);
    const double target_true_anomaly_at_departure = planet_params::planet_true_anomaly_at_time(
        table.config().target_planet,
        departure_time_seconds_since_j2000);
    Problem1TimeOfFlightBranch branch = evaluate_problem1_table_branch_with_target_departure_true_anomaly(
        table,
        query.cell,
        target_true_anomaly_at_departure,
        transfer_revolution,
        target_revolution);
    return Problem1TableBranchQueryResult{std::move(query.cell), branch, true, false};
}

// 在单元已有分支列表中按 (k,q) 精确查找；未找到返回 nullptr。
const Problem1TimeOfFlightBranch* find_problem1_table_branch(
    const Problem1TableCell& cell,
    int transfer_revolution,
    int target_revolution
) {
    for (const Problem1TimeOfFlightBranch& branch : cell.time_of_flight_branches) {
        if (branch.transfer_revolution == transfer_revolution &&
            branch.target_revolution == target_revolution) {
            return &branch;
        }
    }
    return nullptr;
}

std::string problem1_table_branch_signature(const Problem1TableCell& cell) {
    // 中文注释：branch signature 显式编码 (k,q,valid)，避免未来插值阶段误把“branch_count 一样”当成“分支结构一样”。
    std::vector<std::string> parts;
    parts.reserve(cell.time_of_flight_branches.size());
    for (const Problem1TimeOfFlightBranch& branch : cell.time_of_flight_branches) {
        std::ostringstream builder;
        builder
            << "k=" << branch.transfer_revolution
            << ",q=" << branch.target_revolution
            << ",valid=" << (branch.valid ? 1 : 0);
        parts.push_back(builder.str());
    }
    std::sort(parts.begin(), parts.end());

    std::ostringstream joined;
    for (std::size_t index = 0; index < parts.size(); ++index) {
        if (index > 0) {
            joined << '|';
        }
        joined << parts[index];
    }
    return joined.str();
}

bool is_problem1_transfer_branch_valid_for_interpolation(
    const Problem1TimeOfFlightBranch& branch,
    int transfer_revolution
) {
    // 中文注释：future interpolation admissibility 只检查 transfer-side k 有效性；
    // 不检查 q，不检查 target-side TOF，也不检查 representative residual。
    return
        branch.transfer_revolution == transfer_revolution &&
        is_finite(branch.theta_arrival_branch) &&
        is_finite(branch.deltaF_transfer) &&
        branch.deltaF_transfer > 0.0 &&
        is_finite(branch.time_of_flight_scale_free) &&
        is_finite(branch.time_of_flight_seconds) &&
        branch.time_of_flight_seconds >= 0.0;
}

Problem1TransferBranchView get_problem1_transfer_branch_view(
    const Problem1TableCell& cell,
    int transfer_revolution
) {
    Problem1TransferBranchView view{};
    view.transfer_revolution = transfer_revolution;

    bool found_transfer_side_branch = false;
    for (const Problem1TimeOfFlightBranch& branch : cell.time_of_flight_branches) {
        if (!is_problem1_transfer_branch_valid_for_interpolation(branch, transfer_revolution)) {
            continue;
        }

        if (!found_transfer_side_branch) {
            // 中文注释：对同一个 geometry cell 和同一个 k，只提取纯 transfer-side 视图；
            // q、target_time、residual 都不进入这个 view。
            view.valid = true;
            view.theta_arrival_branch = branch.theta_arrival_branch;
            view.deltaF_transfer = branch.deltaF_transfer;
            view.time_of_flight_scale_free = branch.time_of_flight_scale_free;
            view.time_of_flight_seconds = branch.time_of_flight_seconds;
            found_transfer_side_branch = true;
            continue;
        }

        const bool consistent =
            std::abs(normalize_angle_minus_pi_pi(branch.theta_arrival_branch - view.theta_arrival_branch)) <= 1e-12 &&
            std::abs(branch.deltaF_transfer - view.deltaF_transfer) <= 1e-12 &&
            std::abs(branch.time_of_flight_scale_free - view.time_of_flight_scale_free) <= 1e-12 &&
            std::abs(branch.time_of_flight_seconds - view.time_of_flight_seconds) <= 1e-12;
        if (!consistent) {
            view.valid = false;
            view.invalid_reason = "same-k branches disagree on transfer-side invariants across q";
            return view;
        }
    }

    if (!found_transfer_side_branch) {
        view.valid = false;
        view.invalid_reason = "no transfer-side-valid branch found for requested k";
    } else {
        view.invalid_reason.clear();
    }
    return view;
}

Problem1TransferTimeQueryResult query_problem1_transfer_time_with_interpolation_stub(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution
) {
    Problem1TransferTimeQueryResult result{};
    result.transfer_revolution = transfer_revolution;
    result.method = "exact_fallback";

    // 中文注释：当前没有真正插值；这里只是固定未来 transfer-time interpolator 的外部接口。
    result.admissibility = check_problem1_table_transfer_branch_interpolation_admissibility(
        table,
        departure_true_anomaly,
        target_true_anomaly,
        transfer_theta_departure,
        transfer_revolution);
    result.interpolation_admissible = result.admissibility.admissible;

    if (!result.admissibility.admissible) {
        result.ok = false;
        result.reason = result.admissibility.reason.empty()
            ? "interpolation admissibility failed"
            : result.admissibility.reason;
        return result;
    }

    const Problem1TableQueryResult query = query_problem1_table_exact(
        table,
        departure_true_anomaly,
        target_true_anomaly,
        transfer_theta_departure);
    const Problem1TransferBranchView view =
        get_problem1_transfer_branch_view(query.cell, transfer_revolution);
    if (!view.valid) {
        result.ok = false;
        result.reason = view.invalid_reason.empty()
            ? "transfer branch view is unavailable"
            : view.invalid_reason;
        return result;
    }

    result.ok = true;
    result.method = "exact_stub";
    result.reason.clear();
    result.theta_arrival_branch = view.theta_arrival_branch;
    result.deltaF_transfer = view.deltaF_transfer;
    result.time_of_flight_scale_free = view.time_of_flight_scale_free;
    result.time_of_flight_seconds = view.time_of_flight_seconds;
    return result;
}

Problem1TableInterpolationAdmissibility check_problem1_table_transfer_branch_interpolation_admissibility(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution
) {
    Problem1TableInterpolationAdmissibility result{};
    result.transfer_revolution = transfer_revolution;

    if (transfer_revolution < 0) {
        result.reason = "transfer_revolution must be non-negative";
        return result;
    }

    const Problem1TableConfig& config = table.config();
    if (transfer_revolution > config.max_transfer_revolution) {
        result.reason = "transfer_revolution exceeds table max_transfer_revolution";
        return result;
    }

    PeriodicAxisLocation departure_axis{};
    PeriodicAxisLocation target_axis{};
    PeriodicAxisLocation theta_axis{};
    try {
        // 中文注释：这里显式只定位连续几何轴，不接收 q；q 不参与 future interpolation admissibility。
        departure_axis = locate_periodic_axis(
            config.departure_true_anomaly_start,
            config.departure_true_anomaly_step,
            config.departure_true_anomaly_count,
            departure_true_anomaly);
        target_axis = locate_periodic_axis(
            config.target_true_anomaly_start,
            config.target_true_anomaly_step,
            config.target_true_anomaly_count,
            target_true_anomaly);
        theta_axis = locate_periodic_axis(
            config.transfer_theta_departure_start,
            config.transfer_theta_departure_step,
            config.transfer_theta_departure_count,
            transfer_theta_departure);
    } catch (const std::invalid_argument& ex) {
        result.reason = ex.what();
        return result;
    }

    result.local_departure_true_anomaly = departure_axis.query_unwrapped;
    result.local_target_true_anomaly = target_axis.query_unwrapped;
    result.local_transfer_theta_departure = theta_axis.query_unwrapped;

    const int departure_indices[2] = {departure_axis.lower_index, departure_axis.upper_index};
    const int target_indices[2] = {target_axis.lower_index, target_axis.upper_index};
    const int theta_indices[2] = {theta_axis.lower_index, theta_axis.upper_index};

    Problem1TransferConicType reference_conic_type = Problem1TransferConicType::Invalid;
    bool reference_conic_type_set = false;

    int vertex_slot = 0;
    for (int departure_corner = 0; departure_corner < 2; ++departure_corner) {
        for (int target_corner = 0; target_corner < 2; ++target_corner) {
            for (int theta_corner = 0; theta_corner < 2; ++theta_corner) {
                const Problem1TableVertexIndex vertex_index{
                    departure_indices[departure_corner],
                    target_indices[target_corner],
                    theta_indices[theta_corner],
                };
                result.vertex_indices[vertex_slot] = vertex_index;

                const Problem1TableCell& vertex = table.at(
                    vertex_index.departure_true_anomaly_index,
                    vertex_index.target_true_anomaly_index,
                    vertex_index.transfer_theta_departure_index);
                result.vertices[vertex_slot] = vertex;

                if (!vertex.valid) {
                    result.reason = "one or more interpolation vertices have invalid geometry";
                    return result;
                }

                if (!reference_conic_type_set) {
                    reference_conic_type = vertex.conic_type;
                    reference_conic_type_set = true;
                } else if (vertex.conic_type != reference_conic_type) {
                    result.reason = "interpolation vertices do not share a common conic_type";
                    return result;
                }

                const Problem1TransferBranchView branch_view =
                    get_problem1_transfer_branch_view(vertex, transfer_revolution);
                if (!branch_view.valid) {
                    result.reason = "one or more interpolation vertices do not contain a valid transfer branch for k";
                    return result;
                }
                if (branch_view.transfer_revolution != transfer_revolution) {
                    result.reason = "interpolation vertex transfer branch k does not match requested k";
                    return result;
                }

                // 中文注释：这里故意不检查 q；q 不是 universal table branch，只能在真实 t_depart 已知后在线计算。
                ++vertex_slot;
            }
        }
    }

    result.admissible = true;
    result.reason = "admissible";
    return result;
}

}  // namespace spaceship_cpp::problem1
