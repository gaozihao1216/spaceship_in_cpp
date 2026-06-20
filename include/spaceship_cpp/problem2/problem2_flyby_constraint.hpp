/*
 * 文件作用：声明 Problem 2 弹弓约束 F 函数与可缓存残差接口。
 * 主要工作：提供 G = F_in - F_out 的求值，以及 θ' 多分支搜索所需的数据结构。
 *
 * θ' 搜索算法设计（不可对 θ' 直接做单根二分）：
 *
 * 1. 固定飞掠时刻 t_J，缓存入射侧 F_in = F(e_in, theta_g_in)。
 * 2. 扫描 θ' 网格；对每个 θ' 调用 Problem 1 all-roots，得到若干 branch：
 *    (branch_id, phi_out, e_out, p_out, outgoing_residual)。
 * 3. 对每个 branch 计算 G(theta', branch) = F_in - F(e_out, theta_g_out(theta'))。
 * 4. 按 branch_id 分组（禁止跨 branch 合并样本），在组内识别：
 *    - near_zero_node：|G| <= threshold
 *    - sign_change_interval：相邻 θ' 节点 G 变号
 * 5. refine 只在同一 branch 的 sign_change_interval 上做；near_zero 节点可直接收录。
 *    每个区间可能仍有多个根（fold），需复用 Problem 1 的 fold 处理，而非假设单峰。
 */
#pragma once

#include <string>
#include <vector>

namespace spaceship_cpp::problem2 {

struct FlybyConstraintFResult {
    // valid 表示 F 在当前几何下可计算。
    bool valid = false;
    std::string invalid_reason;
    // F(e, theta_global) 在飞掠点的约束值。
    double value = 0.0;
    // 1 + e·cos(phi - theta_global)，接近 0 时几何退化。
    double denominator = 0.0;
};

// 入射侧 F 值缓存：t_J 与入射轨道固定后，搜索 θ' 时只需重复计算出射 F。
struct FlybyConstraintIncomingCache {
    bool valid = false;
    std::string invalid_reason;

    double flyby_true_anomaly_phi = 0.0;
    double flyby_planet_eccentricity = 0.0;
    double incoming_eccentricity = 0.0;
    double incoming_theta_global = 0.0;

    double incoming_F = 0.0;
};

struct FlybyConstraintResidualResult {
    bool valid = false;
    std::string invalid_reason;
    // G = F_in - F_out；G = 0 表示满足弹弓约束。
    double residual = 0.0;
    double incoming_F = 0.0;
    double outgoing_F = 0.0;
};

enum class FlybyThetaPrimeCandidateType {
    // 离散网格节点上 |G| 已足够小。
    NearZeroNode,
    // 同一 branch 在相邻 θ' 节点间 G 变号（含端点为零）。
    SignChangeInterval,
};

// 固定 θ' 下、某一 Problem 1 出射分支的完整样本。
struct FlybyThetaPrimeBranchSample {
    double theta_prime = 0.0;
    int branch_id = -1;

    double outgoing_phi = 0.0;
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
    double outgoing_problem1_residual = 0.0;

    double flyby_constraint_residual = 0.0;
    bool flyby_constraint_valid = false;
};

// 由离散扫描识别出的 θ' 候选；同一 branch 上才可做区间 refine。
struct FlybyThetaPrimeCandidate {
    FlybyThetaPrimeCandidateType type = FlybyThetaPrimeCandidateType::NearZeroNode;
    int branch_id = -1;

    double theta_prime_node = 0.0;
    double theta_prime_left = 0.0;
    double theta_prime_right = 0.0;

    double G_node = 0.0;
    double G_left = 0.0;
    double G_right = 0.0;

    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
};

// 计算 F(e, theta_global) = [1 + e_J cos(phi) + e(cos(phi-theta_g) + e_J cos(theta_g))] / sqrt(1 + e cos(phi-theta_g))。
FlybyConstraintFResult evaluate_flyby_constraint_F(
    double orbit_eccentricity,
    double orbit_perihelion_angle_global,
    double flyby_true_anomaly_phi,
    double flyby_planet_eccentricity
);

// 将局部角（以飞掠点方向为参考，Python 约定）或全局近日点角转换为全局角。
double flyby_orbit_theta_global_from_input(
    double orbit_theta,
    double flyby_true_anomaly_phi,
    bool input_theta_is_local
);

// 构造入射侧缓存；input_theta_is_local=true 时 orbit_theta 为局部角。
FlybyConstraintIncomingCache build_flyby_constraint_incoming_cache(
    double incoming_eccentricity,
    double incoming_theta,
    double flyby_true_anomaly_phi,
    double flyby_planet_eccentricity,
    bool input_theta_is_local
);

// 直接计算 G = F_in - F_out（无缓存）。
FlybyConstraintResidualResult evaluate_flyby_constraint_residual(
    double incoming_eccentricity,
    double incoming_theta,
    double outgoing_eccentricity,
    double outgoing_theta,
    double flyby_true_anomaly_phi,
    double flyby_planet_eccentricity,
    bool input_angles_are_local
);

// 使用入射缓存计算 G；outgoing_theta 的坐标系约定与 build 时相同。
FlybyConstraintResidualResult evaluate_flyby_constraint_residual_from_incoming_cache(
    const FlybyConstraintIncomingCache& incoming_cache,
    double outgoing_eccentricity,
    double outgoing_theta,
    bool outgoing_theta_is_local
);

// 在同一 branch 的相邻 θ' 样本间识别 near-zero 节点与变号区间。
// 不能跨 branch 合并：每个 θ' 可能对应多个 P1 根，必须按 branch_id 分组后再检测。
std::vector<FlybyThetaPrimeCandidate> detect_flyby_theta_prime_candidates_from_branch_samples(
    const std::vector<FlybyThetaPrimeBranchSample>& branch_samples,
    double near_zero_threshold
);

}  // namespace spaceship_cpp::problem2
