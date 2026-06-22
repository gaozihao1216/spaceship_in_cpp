/*
 * 文件作用：测试全局默认配置。
 * 主要工作：检查 Problem 1、Problem 2 与轨迹搜索默认参数是否满足基本有效性。
 */
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"

#include <cassert>

int main() {
    namespace config = spaceship_cpp::config;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const config::GlobalConfig& cfg = config::global_config();

    assert(cfg.problem1_solve.max_transfer_revolution >= 0);
    assert(cfg.problem1_solve.max_target_revolution >= 0);
    assert(cfg.problem1_solve.phi_scan_count >= 3);
    assert(cfg.problem1_solve.phi_tolerance > 0.0);
    assert(cfg.problem1_solve.max_bisection_iterations > 0);
    assert(cfg.problem1_solve.max_candidate_relative_residual > 0.0);
    assert(cfg.problem2_theta_prime_scan.theta_prime_count >= 3);
    assert(cfg.problem2_theta_prime_scan.phi_scan_count >= 3);
    assert(cfg.problem2_theta_prime_scan.branch_phi_pairing_max_gap > 0.0);
    assert(cfg.trajectory_search.leg0_theta_coarse_scan_count >= 2);
    assert(cfg.trajectory_search.leg0_theta_fine_scan_count >= 2);
    assert(cfg.trajectory_search.leg0_theta_seeds_per_first_leg_target == 3);
    assert(cfg.trajectory_search.max_search_legs >= 1);
    assert(cfg.trajectory_search.top_k_sequences > 0);
    assert(cfg.trajectory_search.max_transfer_revolution == 1);
    assert(cfg.trajectory_search.max_target_revolution == 1);

    {
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
        const auto scan_config = config::make_problem2_theta_prime_scan_config(
            planet_params::PlanetId::Mars,
            planet_params::PlanetId::Earth,
            100.0 * 86400.0,
            cfg.problem2_theta_prime_scan,
            cfg.problem1_solve);
        assert(scan_config.flyby_planet == planet_params::PlanetId::Mars);
        assert(scan_config.target_planet == planet_params::PlanetId::Earth);
        assert(scan_config.theta_prime_count == cfg.problem2_theta_prime_scan.theta_prime_count);
        assert(scan_config.branch_phi_pairing_max_gap ==
            cfg.problem2_theta_prime_scan.branch_phi_pairing_max_gap);
        assert(scan_config.problem1_solve.phi_scan_count ==
            cfg.problem2_theta_prime_scan.phi_scan_count);
    }

    return 0;
}
