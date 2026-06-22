/*
 * 文件作用：实现全局默认配置。
 * 主要工作：集中返回 Problem 1、Problem 2 与轨迹搜索管线的默认参数。
 */
#include "spaceship_cpp/config/global_config.hpp"

namespace spaceship_cpp::config {

const GlobalConfig& global_config() {
    static const GlobalConfig kConfig{
        Problem1SolveDefaults{
            0,      // 默认不允许转移轨道额外绕日，先聚焦最直接的前向转移分支
            0,      // 默认不允许目标行星多绕额外圈数，先只看最近一次前向相遇
            120,    // 单次 solve 中 phi 扫描 120 个点，约等于每 3 度采样一次
            1e-10,  // 二分求根的角度容差
            0.0,    // residual 绝对容差当前保持为 0
            80,     // 每个变号区间最多二分 80 次
            1e-6,   // 候选解相对残差阈值
        },
        Problem2ThetaPrimeScanDefaults{
            64,     // Problem 2 出射段 θ' 初扫默认 64 个离散点
            64,     // 初扫专用 phi 扫描 64 点（全精度 solve 仍为 120）；初扫只需发现 branch 拓扑
            0.75,   // branch 配对默认允许最大约 43° 的 φ 跳变
        },
        Problem2RouteANewtonDefaults{
            20,     // Route A Newton 默认最多 20 次迭代
            1e-6,   // 相对残差收敛阈值，与 Problem 1 候选过滤一致
            1e-10,  // φ 更新步长容差
            1e-8,   // df/dφ 数值差分步长
            true,   // 残差增大即判定发散
        },
        TrajectorySearchDefaults{
            60,     // Step 1：θ 粗扫约 6° 一步，只划分可行区间
            180,    // Step 4：固定序列 θ 细扫
            3,      // Step 2：每个首段目标（Earth→P₁）的 θ 种子数上限（粗筛序列，细 θ 在 Step 4）
            6,      // Step 2：自由路径最多 6 段转移（Earth→P₁→…→P₆）
            3,      // Step 3：Top-K 行星序列
            1,      // 搜索管线 max k：至少 1，k=0 会漏掉多圈转移分支
            1,      // 搜索管线 max q：至少 1，q=0 会漏掉目标多圈相遇分支
        },
    };

    return kConfig;
}

problem1::Problem1SolveInput make_problem1_solve_input(
    planet_params::PlanetId departure,
    planet_params::PlanetId target,
    double launch_time_seconds_since_j2000,
    double transfer_perihelion_angle_global,
    const Problem1SolveDefaults& defaults
) {
    return problem1::Problem1SolveInput{
        departure,
        target,
        launch_time_seconds_since_j2000,
        transfer_perihelion_angle_global,
        defaults.max_transfer_revolution,
        defaults.max_target_revolution,
        defaults.phi_scan_count,
        defaults.phi_tolerance,
        defaults.residual_tolerance,
        defaults.max_bisection_iterations,
        defaults.max_candidate_relative_residual,
    };
}

problem2::Problem2ThetaPrimeScanConfig make_problem2_theta_prime_scan_config(
    planet_params::PlanetId flyby_planet,
    planet_params::PlanetId target_planet,
    double flyby_time_seconds_since_j2000,
    const Problem2ThetaPrimeScanDefaults& scan_defaults,
    const Problem1SolveDefaults& problem1_defaults
) {
    problem2::Problem2ThetaPrimeScanConfig config{};
    config.flyby_planet = flyby_planet;
    config.target_planet = target_planet;
    config.flyby_time_seconds_since_j2000 = flyby_time_seconds_since_j2000;
    config.theta_prime_count = scan_defaults.theta_prime_count;
    config.branch_phi_pairing_max_gap = scan_defaults.branch_phi_pairing_max_gap;
    config.problem1_solve = make_problem1_solve_input(
        flyby_planet,
        target_planet,
        flyby_time_seconds_since_j2000,
        0.0,
        problem1_defaults);
    if (scan_defaults.phi_scan_count >= 3) {
        config.problem1_solve.phi_scan_count = scan_defaults.phi_scan_count;
    }
    return config;
}

problem2::Problem2RouteANewtonOptions make_problem2_route_a_newton_options(
    const Problem2RouteANewtonDefaults& defaults,
    const Problem1SolveDefaults& problem1_defaults
) {
    problem2::Problem2RouteANewtonOptions options{};
    options.max_iterations = defaults.max_iterations;
    options.residual_tolerance = problem1_defaults.residual_tolerance;
    options.max_relative_residual = defaults.max_relative_residual;
    options.phi_tolerance = defaults.phi_tolerance;
    options.phi_derivative_step = defaults.phi_derivative_step;
    options.reject_on_residual_increase = defaults.reject_on_residual_increase;
    return options;
}

problem2::Problem2FlybySolveConfig make_problem2_flyby_solve_config(
    planet_params::PlanetId flyby_planet,
    planet_params::PlanetId target_planet,
    double flyby_time_seconds_since_j2000,
    const Problem2ThetaPrimeScanDefaults& scan_defaults,
    const Problem2RouteANewtonDefaults& route_a_defaults,
    const Problem1SolveDefaults& problem1_defaults
) {
    problem2::Problem2FlybySolveConfig config{};
    config.g_search_config.scan_config = make_problem2_theta_prime_scan_config(
        flyby_planet,
        target_planet,
        flyby_time_seconds_since_j2000,
        scan_defaults,
        problem1_defaults);
    config.g_search_config.route_a_newton_options = make_problem2_route_a_newton_options(
        route_a_defaults,
        problem1_defaults);
    config.g_search_config.g_newton_options.max_iterations = route_a_defaults.max_iterations;
    config.g_search_config.g_newton_options.G_tolerance = 1e-8;
    config.g_search_config.g_newton_options.theta_prime_tolerance = 1e-10;
    config.g_search_config.near_zero_G_threshold = 1e-4;
    config.g_search_config.solution_theta_prime_tolerance = 1e-8;
    config.g_search_config.theta_prime_tolerance = 1e-8;
    config.g_search_config.max_recursion_depth = 16;
    return config;
}

}  // namespace spaceship_cpp::config
