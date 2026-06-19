/*
 * 文件作用：测试全局默认配置。
 * 主要工作：检查默认行星、Problem 1 参数和表格配置是否满足基本有效性。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace {

// 中文说明：提供带绝对/相对容差的浮点近似相等判断，避免配置步长比较受舍入误差影响。
bool approx_equal(double a, double b, double abs_tol = 1e-12, double rel_tol = 1e-12) {
    const double diff = std::abs(a - b);
    if (diff <= abs_tol) {
        return true;
    }
    return diff <= rel_tol * std::max(std::abs(a), std::abs(b));
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace config = spaceship_cpp::config;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const config::GlobalConfig& cfg = config::global_config();

    // 中文说明：验证 Problem 1 求解器默认参数（圈数上限、扫描点数、容差等）均为正且满足最小约束。
    assert(cfg.problem1_solve.max_transfer_revolution >= 0);
    assert(cfg.problem1_solve.max_target_revolution >= 0);
    assert(cfg.problem1_solve.phi_scan_count >= 3);
    assert(cfg.problem1_solve.phi_tolerance > 0.0);
    assert(cfg.problem1_solve.max_bisection_iterations > 0);
    assert(cfg.problem1_solve.max_candidate_relative_residual > 0.0);

    // 中文说明：验证 Problem 1 表格 smoke 默认配置的采样计数与圈数上限合法。
    assert(cfg.problem1_table_smoke.departure_true_anomaly_count > 0);
    assert(cfg.problem1_table_smoke.target_true_anomaly_count > 0);
    assert(cfg.problem1_table_smoke.transfer_theta_departure_count > 0);
    assert(cfg.problem1_table_smoke.max_transfer_revolution >= 0);
    assert(cfg.problem1_table_smoke.max_target_revolution >= 0);

    {
        // 中文说明：验证 make_problem1_solve_input 将全局默认与行星/时间/转移角正确组装为求解输入。
        const problem1::Problem1SolveInput input = config::make_problem1_solve_input(
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            0.5,
            cfg.problem1_solve);
        assert(input.departure_planet == planet_params::PlanetId::Earth);
        assert(input.target_planet == planet_params::PlanetId::Mars);
        assert(input.launch_time_seconds_since_j2000 == 0.0);
        assert(input.transfer_perihelion_angle == 0.5);
        assert(input.phi_scan_count == cfg.problem1_solve.phi_scan_count);
        assert(input.max_candidate_relative_residual ==
            cfg.problem1_solve.max_candidate_relative_residual);
    }

    {
        // 中文说明：验证 make_problem1_table_config 按采样计数推导三轴步长，并传递圈数上限。
        const problem1::Problem1TableConfig table_config = config::make_problem1_table_config(
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            cfg.problem1_table_smoke);
        assert(table_config.departure_planet == planet_params::PlanetId::Earth);
        assert(table_config.target_planet == planet_params::PlanetId::Mars);
        assert(approx_equal(
            table_config.departure_true_anomaly_step,
            common::kTwoPi / static_cast<double>(cfg.problem1_table_smoke.departure_true_anomaly_count)));
        assert(approx_equal(
            table_config.target_true_anomaly_step,
            common::kTwoPi / static_cast<double>(cfg.problem1_table_smoke.target_true_anomaly_count)));
        assert(approx_equal(
            table_config.transfer_theta_departure_step,
            common::kTwoPi / static_cast<double>(cfg.problem1_table_smoke.transfer_theta_departure_count)));
        assert(table_config.max_transfer_revolution == cfg.problem1_table_smoke.max_transfer_revolution);
        assert(table_config.max_target_revolution == cfg.problem1_table_smoke.max_target_revolution);
    }

    {
        // 中文说明：验证缩小采样网格后仍能通过 build_problem1_table 成功生成非空表格。
        config::Problem1TableDefaults small_defaults = cfg.problem1_table_smoke;
        small_defaults.departure_true_anomaly_count = 2;
        small_defaults.target_true_anomaly_count = 2;
        small_defaults.transfer_theta_departure_count = 3;
        const problem1::Problem1TableConfig small_table_config = config::make_problem1_table_config(
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            small_defaults);
        const problem1::Problem1Table table = problem1::build_problem1_table(small_table_config);
        assert(!table.cells().empty());
    }

    return 0;
}
