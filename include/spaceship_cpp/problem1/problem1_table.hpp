/*
 * 文件作用：声明 Problem 1 endpoint transfer-time table 数据结构和查询接口。
 * 主要工作：描述表格元数据、单元几何、飞行时间分支和插值可行性检查。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace spaceship_cpp::problem1 {

enum class Problem1TransferConicType {
    // 当前几何无法构造有效圆锥曲线。
    Invalid,
    // 椭圆转移，是当前时间积分最常用的有效类型。
    Elliptic,
    // 抛物线边界分类，帮助诊断接近 e=1 的情况。
    Parabolic,
    // 双曲线分类；当前表格记录类型但主要测试仍围绕椭圆转移。
    Hyperbolic,
};

struct Problem1TimeOfFlightBranch {
    // 该多圈时间分支是否物理可用。
    bool valid;
    // k/q 分支编号：分别表示转移轨道和目标轨道额外绕行圈数。
    int transfer_revolution;
    int target_revolution;
    // 目标到达角所属的分支角，用于区分不同绕行路径。
    double theta_arrival_branch;
    // 目标行星在发射和到达时的局部真近点角。
    double target_true_anomaly_start;
    double target_true_anomaly_end_branch;
    // 转移轨道和目标轨道走过的异常角差。
    double deltaF_transfer;
    double deltaF_target;
    // 两套时间表示：scale_free 便于理论比较，seconds 便于物理输出。
    double time_of_flight_scale_free;
    double target_time_of_flight_scale_free;
    double time_of_flight_seconds;
    double target_time_of_flight_seconds;
    // 两段时间差，用于判断该分支是否满足 Problem 1 时间一致性。
    double residual_scale_free;
    double residual_seconds;
    std::string invalid_reason;
};

struct Problem1TableMetadata {
    // 表格语义版本，防止旧数据被误当成当前轴定义读取。
    std::string schema_version = "planet_angle_pair_table_v2";
    // 表格固定服务一对出发/目标行星。
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;
    std::string departure_planet_name;
    std::string target_planet_name;

    // 生成表格时使用的两条行星轨道参数快照。
    double departure_orbit_p;
    double departure_orbit_e;
    double departure_orbit_theta_0;
    double target_orbit_p;
    double target_orbit_e;
    double target_orbit_theta_0;

    // 三个表格轴和角度约定的文字说明，用于诊断输出和文档化。
    std::string axis1_definition;
    std::string axis2_definition;
    std::string axis3_definition;
    std::string angle_convention;
    std::string derivative_method;
    std::string target_branch_validity_note;

    // 表格维度和每个单元枚举的多圈上限。
    int departure_true_anomaly_count;
    int target_true_anomaly_count;
    int transfer_theta_departure_count;
    int max_transfer_revolution;
    int max_target_revolution;
};

struct Problem1TableConfig {
    // 指定本表格覆盖的行星对。
    planet_params::PlanetId departure_planet;
    planet_params::PlanetId target_planet;

    // 轴 1：出发行星真近点角采样。
    double departure_true_anomaly_start;
    double departure_true_anomaly_step;
    int departure_true_anomaly_count;

    // 轴 2：目标行星真近点角采样。
    double target_true_anomaly_start;
    double target_true_anomaly_step;
    int target_true_anomaly_count;

    // 轴 3：转移轨道在出发点的局部角采样。
    double transfer_theta_departure_start;
    double transfer_theta_departure_step;
    int transfer_theta_departure_count;

    // 每个单元中枚举的最大转移圈数和目标圈数。
    int max_transfer_revolution;
    int max_target_revolution;
};

struct Problem1TableCell {
    // 单元在三维表格中的离散索引。
    int departure_true_anomaly_index;
    int target_true_anomaly_index;
    int transfer_theta_departure_index;

    // 原始查询/采样输入，保留未归一化前的调用侧信息。
    double departure_true_anomaly_input;
    double target_true_anomaly_input;
    double transfer_theta_departure_input;

    // 归一化后的三轴角度，以及由相位反推的发射时间信息。
    double departure_true_anomaly;
    double target_true_anomaly;
    double transfer_theta_departure;
    double transfer_theta_arrival;
    double launch_time_phase_seconds_since_j2000;
    double target_true_anomaly_at_departure;
    bool target_branches_are_launch_phase_specific;

    // 两个端点的太阳中心几何，是解转移圆锥曲线的基础数据。
    double departure_radius;
    double target_radius;
    double departure_global_angle;
    double target_global_angle;
    double delta_global_angle;

    // 转移轨道近日点全局角；raw 是直接几何结果，used 是规范化后实际使用值。
    double transfer_perihelion_angle_global_raw;
    double transfer_perihelion_angle_global_used;

    // 转移圆锥曲线参数。
    double transfer_e_raw;
    double transfer_e;
    double transfer_p;
    double transfer_a;
    bool normalized_negative_e;

    // 单元整体几何分类和有效性。
    Problem1TransferConicType conic_type;
    bool valid;
    std::string invalid_reason;

    // 偏心率对三个输入轴的数值导数，用来辅助判断局部变化是否平滑。
    bool derivatives_available;
    double dtransfer_e_dnu_departure;
    double dtransfer_e_dnu_target;
    double dtransfer_e_dtheta_departure;

    // 该几何单元下所有枚举出的飞行时间分支。
    std::vector<Problem1TimeOfFlightBranch> time_of_flight_branches;
};

struct Problem1TableQueryResult {
    // 查询返回的单元；可能来自已构建网格，也可能由调用方后续即时计算。
    Problem1TableCell cell;
    // true 表示 cell 是从离散表格中取出的精确网格点。
    bool sourced_from_grid = false;
};

struct Problem1TableBranchQueryResult {
    // 单元和匹配到的时间分支一起返回，方便调用方同时检查几何与时间。
    Problem1TableCell cell;
    Problem1TimeOfFlightBranch branch;
    // false 表示单元存在但没有指定 k/q 分支。
    bool branch_found = false;
    bool sourced_from_grid = false;
};

struct Problem1TableVertexIndex {
    // 插值立方体一个顶点的三维索引。
    int departure_true_anomaly_index;
    int target_true_anomaly_index;
    int transfer_theta_departure_index;
};

struct Problem1TableInterpolationAdmissibility {
    // 说明当前查询点周围的八个顶点是否允许在同一分支上插值。
    bool admissible = false;
    // 不可插值时记录第一处失败原因，便于诊断表格拓扑问题。
    std::string reason;
    // 八个顶点的索引和完整单元数据。
    std::array<Problem1TableVertexIndex, 8> vertex_indices{};
    std::array<Problem1TableCell, 8> vertices{};
    int transfer_revolution = 0;
    // 查询点在所在网格小立方体中的局部坐标。
    double local_departure_true_anomaly = 0.0;
    double local_target_true_anomaly = 0.0;
    double local_transfer_theta_departure = 0.0;
};

struct Problem1TransferBranchView {
    // 从完整分支中抽取插值需要的最小信息。
    bool valid = false;
    int transfer_revolution = 0;
    double theta_arrival_branch = 0.0;
    double deltaF_transfer = 0.0;
    double time_of_flight_scale_free = 0.0;
    double time_of_flight_seconds = 0.0;
    std::string invalid_reason;
};

struct Problem1TransferTimeQueryResult {
    // ok 表示当前查询拿到了可用转移时间。
    bool ok = false;
    // 标记成功是否来自通过准入检查的插值路径。
    bool interpolation_admissible = false;
    // method/reason 让测试知道走的是插值、回退还是失败路径。
    std::string method;
    std::string reason;
    int transfer_revolution = 0;
    double theta_arrival_branch = 0.0;
    double deltaF_transfer = 0.0;
    double time_of_flight_scale_free = 0.0;
    double time_of_flight_seconds = 0.0;
    Problem1TableInterpolationAdmissibility admissibility{};
};

// 解决“表格元数据是否符合当前代码约定”的问题。
void validate_problem1_table_metadata(const Problem1TableMetadata& metadata);

// 只从两端点太阳中心几何构造单元，隔离行星相位相关逻辑。
Problem1TableCell evaluate_problem1_table_cell_geometry(
    double departure_radius,
    double departure_global_angle,
    double target_radius,
    double target_global_angle,
    double transfer_theta_departure_input,
    int max_transfer_revolution
);

// 从行星局部真近点角生成完整单元，解决“该行星相位下表格值是什么”的问题。
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

    // 暴露构建配置和元数据，供诊断程序确认表格语义。
    const Problem1TableConfig& config() const;
    const Problem1TableMetadata& metadata() const;

    // 返回三个轴的尺寸，调用方遍历表格前用它们确定边界。
    int departure_true_anomaly_count() const;
    int target_true_anomaly_count() const;
    int transfer_theta_departure_count() const;

    // 将三维索引映射到一维数组位置，统一表格内部存储布局。
    std::size_t flat_index(
        int departure_true_anomaly_index,
        int target_true_anomaly_index,
        int transfer_theta_departure_index
    ) const;

    // 读写指定网格单元；非 const 版本用于构建阶段填表。
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

    // 返回所有单元，供批量统计和诊断遍历使用。
    const std::vector<Problem1TableCell>& cells() const;

private:
    Problem1TableConfig config_;
    Problem1TableMetadata metadata_;
    std::vector<Problem1TableCell> cells_;
};

// 构建完整 endpoint transfer-time 表格。
Problem1Table build_problem1_table(const Problem1TableConfig& config);

// 精确查询离散网格点，解决“输入是否正好落在已建表节点”的问题。
Problem1TableQueryResult query_problem1_table_exact(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure
);

// 精确查询网格点上的指定 k/q 分支。
Problem1TableBranchQueryResult query_problem1_table_exact_branch(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution,
    int target_revolution
);  // 仅返回 representative / launch-phase-specific q branch；不用于跨相位插值准入判断。

// 在给定真实出发时间时查询分支，解决目标初始相位依赖的问题。
Problem1TableBranchQueryResult query_problem1_table_exact_branch_at_departure_time(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    double departure_time_seconds_since_j2000,
    int transfer_revolution,
    int target_revolution
);

// 对指定单元重新按目标出发真近点角评估分支时间。
Problem1TimeOfFlightBranch evaluate_problem1_table_branch_with_target_departure_true_anomaly(
    const Problem1Table& table,
    const Problem1TableCell& cell,
    double target_true_anomaly_at_departure,
    int transfer_revolution,
    int target_revolution
);

// 在单元已有分支列表中寻找指定 k/q 分支。
const Problem1TimeOfFlightBranch* find_problem1_table_branch(
    const Problem1TableCell& cell,
    int transfer_revolution,
    int target_revolution
);

// 生成分支有效性签名，帮助测试快速比较单元拓扑是否一致。
std::string problem1_table_branch_signature(const Problem1TableCell& cell);

// 判断单个时间分支是否适合作为插值顶点。
bool is_problem1_transfer_branch_valid_for_interpolation(
    const Problem1TimeOfFlightBranch& branch,
    int transfer_revolution
);

// 获取指定转移圈数的轻量分支视图。
Problem1TransferBranchView get_problem1_transfer_branch_view(
    const Problem1TableCell& cell,
    int transfer_revolution
);

// 当前仅做准入检查和保守回退，解决“插值路径是否可以安全启用”的实验接口。
Problem1TransferTimeQueryResult query_problem1_transfer_time_with_interpolation_stub(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution
);  // 基于当前 endpoint-table 的插值占位接口；失败时不会伪造插值结果。

// 检查八个邻近顶点是否能组成同一转移时间分支。
Problem1TableInterpolationAdmissibility check_problem1_table_transfer_branch_interpolation_admissibility(
    const Problem1Table& table,
    double departure_true_anomaly,
    double target_true_anomaly,
    double transfer_theta_departure,
    int transfer_revolution
);  // 仅适用于当前 endpoint-table / transfer-time-table 设计。

}  // namespace spaceship_cpp::problem1
