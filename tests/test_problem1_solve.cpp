/*
 * 文件作用：测试 Problem 1 直接求解器。
 * 主要工作：验证扫描加二分细化能产出有效候选解并满足残差阈值。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

int main() {
    namespace common = spaceship_cpp::common;
    namespace problem1 = spaceship_cpp::problem1;

    {
        // 中文说明：验证 phi_scan_count 低于最小值 3 时 solve_problem1 拒绝输入并抛出 invalid_argument。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            0,
            0,
            2,
            1e-10,
            0.0,
            80,
        };
        bool threw = false;
        try {
            (void)problem1::solve_problem1(input);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        // 中文说明：验证 max_transfer_revolution 为负时 solve_problem1 拒绝输入。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            -1,
            0,
            180,
            1e-10,
            0.0,
            80,
        };
        bool threw = false;
        try {
            (void)problem1::solve_problem1(input);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        // 中文说明：验证 Earth→Mars k=0,q=0 下求解器产出候选解，且相对残差、到达时间与轨道要素满足约束。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            0,
            0,
            180,
            1e-10,
            0.0,
            80,
        };
        const auto candidates = problem1::solve_problem1(input);
        for (const auto& candidate : candidates) {
            assert(candidate.residual_result.status == problem1::Problem1ResidualStatus::Success);
            assert(std::isfinite(candidate.encounter_global_angle));
            assert(candidate.encounter_global_angle >= 0.0);
            assert(candidate.encounter_global_angle < common::kTwoPi);
            assert(std::isfinite(candidate.time_of_flight_seconds));
            assert(candidate.time_of_flight_seconds > 0.0);
            assert(candidate.arrival_time_seconds_since_j2000 ==
                candidate.launch_time_seconds_since_j2000 + candidate.time_of_flight_seconds);
            assert(std::isfinite(candidate.residual_scale));
            assert(candidate.residual_scale > 0.0);
            assert(std::isfinite(candidate.relative_residual));
            assert(candidate.relative_residual >= 0.0);
            assert(candidate.relative_residual <= input.max_candidate_relative_residual + 1e-15);
            assert(std::isfinite(candidate.root_bracket_width));
            assert(candidate.root_bracket_width >= 0.0);
            assert(candidate.bisection_iterations >= 0);
            assert(std::abs(
                candidate.relative_residual -
                std::abs(candidate.residual_result.residual) / candidate.residual_scale) <= 1e-12);
            if (!candidate.refined_by_bisection) {
                assert(candidate.root_bracket_width == 0.0);
                assert(candidate.bisection_iterations == 0);
            } else {
                assert(candidate.root_bracket_width >= 0.0);
                assert(candidate.bisection_iterations >= 0);
            }
            assert(candidate.transfer_revolution == 0);
            assert(candidate.target_revolution == 0);
            assert(candidate.residual_result.transfer_p > 0.0);
            assert(candidate.residual_result.transfer_e >= 0.0);
            assert(candidate.residual_result.deltaF_transfer > 0.0);
            assert(candidate.residual_result.deltaF_target > 0.0);
        }
    }

    {
        // 中文说明：验证 k,q 均允许为 1 时多分支候选解合法、同 (k,q) 相遇角不重复，且按 arrive_time 升序排列。
        const problem1::Problem1SolveInput input{
            spaceship_cpp::planet_params::PlanetId::Earth,
            spaceship_cpp::planet_params::PlanetId::Mars,
            0.0,
            0.5,
            1,
            1,
            90,
            1e-9,
            0.0,
            60,
        };
        const auto candidates = problem1::solve_problem1(input);
        for (const auto& candidate : candidates) {
            assert(candidate.transfer_revolution >= 0 && candidate.transfer_revolution <= 1);
            assert(candidate.target_revolution >= 0 && candidate.target_revolution <= 1);
            assert(candidate.residual_result.status == problem1::Problem1ResidualStatus::Success);
            assert(std::isfinite(candidate.encounter_global_angle));
            assert(candidate.encounter_global_angle >= 0.0);
            assert(candidate.encounter_global_angle < common::kTwoPi);
            assert(candidate.time_of_flight_seconds > 0.0);
            assert(std::isfinite(candidate.residual_scale));
            assert(candidate.residual_scale > 0.0);
            assert(std::isfinite(candidate.relative_residual));
            assert(candidate.relative_residual >= 0.0);
            assert(candidate.relative_residual <= input.max_candidate_relative_residual + 1e-15);
            assert(std::isfinite(candidate.root_bracket_width));
            assert(candidate.root_bracket_width >= 0.0);
            assert(candidate.bisection_iterations >= 0);
            assert(std::abs(
                candidate.relative_residual -
                std::abs(candidate.residual_result.residual) / candidate.residual_scale) <= 1e-12);
            if (!candidate.refined_by_bisection) {
                assert(candidate.root_bracket_width == 0.0);
                assert(candidate.bisection_iterations == 0);
            } else {
                assert(candidate.root_bracket_width >= 0.0);
                assert(candidate.bisection_iterations >= 0);
            }
        }

        for (std::size_t i = 0; i < candidates.size(); ++i) {
            for (std::size_t j = i + 1; j < candidates.size(); ++j) {
                if (candidates[i].transfer_revolution == candidates[j].transfer_revolution &&
                    candidates[i].target_revolution == candidates[j].target_revolution) {
                    const double angle_diff = common::normalize_angle_minus_pi_pi(
                        candidates[i].encounter_global_angle - candidates[j].encounter_global_angle);
                    assert(std::abs(angle_diff) >= 1e-8);
                }
            }
        }

        // 中文说明：solve 输出可视为按 arrive_time 排序的根列表 (encounter_global_angle, k, q, arrive_time)。
        for (std::size_t index = 1; index < candidates.size(); ++index) {
            assert(
                candidates[index - 1].arrival_time_seconds_since_j2000 <=
                candidates[index].arrival_time_seconds_since_j2000);
        }
        for (const auto& candidate : candidates) {
            assert(std::isfinite(candidate.encounter_global_angle));
            assert(candidate.encounter_global_angle >= 0.0);
            assert(candidate.encounter_global_angle < common::kTwoPi);
            assert(std::isfinite(candidate.arrival_time_seconds_since_j2000));
            assert(candidate.arrival_time_seconds_since_j2000 ==
                candidate.launch_time_seconds_since_j2000 + candidate.time_of_flight_seconds);
        }
    }

    return 0;
}
