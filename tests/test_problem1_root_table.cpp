#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

namespace {

bool approx_equal(double lhs, double rhs, double abs_tol = 1e-6, double rel_tol = 1e-9) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= abs_tol + rel_tol * scale;
}

void compare_solution_lists(
    spaceship_cpp::planet_params::PlanetId target_planet,
    const std::vector<spaceship_cpp::problem1::Problem1Candidate>& candidates,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& direct
) {
    namespace common = spaceship_cpp::common;
    namespace problem1 = spaceship_cpp::problem1;

    const auto converted = problem1::convert_problem1_candidates_to_solution_branches(target_planet, candidates);
    assert(converted.size() == direct.size());
    for (std::size_t index = 0; index < direct.size(); ++index) {
        const auto& lhs = converted[index];
        const auto& rhs = direct[index];
        assert(lhs.transfer_revolution == rhs.transfer_revolution);
        assert(lhs.target_revolution == rhs.target_revolution);
        assert(std::abs(
            common::normalize_angle_minus_pi_pi(lhs.encounter_global_angle - rhs.encounter_global_angle)) <= 1e-8);
        assert(std::abs(
            common::normalize_angle_minus_pi_pi(
                lhs.target_arrival_true_anomaly - rhs.target_arrival_true_anomaly)) <= 1e-8);
        assert(approx_equal(lhs.time_of_flight_seconds, rhs.time_of_flight_seconds));
        assert(approx_equal(lhs.target_time_seconds, rhs.target_time_seconds));
        assert(approx_equal(lhs.residual_seconds, rhs.residual_seconds));
        assert(std::isnan(rhs.arrival_time_seconds_since_j2000));
    }
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    {
        // 中文注释：先验证 root-table 草案数据结构能被正常构建和索引。
        const problem1::Problem1RootTableConfig config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kHalfPi,
            2,
            0.25,
            common::kHalfPi,
            2,
            0.5,
            common::kPi,
            2,
            1,
            1,
            "problem1_root_table_draft_v0",
        };

        const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
        assert(table.config().nu_A_count == 2);
        assert(table.config().nu_B_depart_count == 2);
        assert(table.config().theta_A_count == 2);
        assert(table.cells().size() == 8);

        for (int i = 0; i < config.nu_A_count; ++i) {
            for (int j = 0; j < config.nu_B_depart_count; ++j) {
                for (int k = 0; k < config.theta_A_count; ++k) {
                    const auto& cell = table.at(i, j, k);
                    assert(std::isfinite(cell.nu_A_depart));
                    assert(std::isfinite(cell.nu_B_depart));
                    assert(std::isfinite(cell.theta_A));
                    assert(cell.nu_A_depart >= 0.0 && cell.nu_A_depart < common::kTwoPi);
                    assert(cell.nu_B_depart >= 0.0 && cell.nu_B_depart < common::kTwoPi);
                    assert(cell.theta_A >= 0.0 && cell.theta_A < common::kTwoPi);
                    assert(cell.solved);
                    // 中文注释：solved=true 只表示该 cell 已执行过 solve；即使无解，也允许 solutions 为空。
                    assert(cell.invalid_reason.empty());
                    for (std::size_t index = 1; index < cell.solutions_sorted_by_time_of_flight.size(); ++index) {
                        assert(
                            cell.solutions_sorted_by_time_of_flight[index - 1].time_of_flight_seconds <=
                            cell.solutions_sorted_by_time_of_flight[index].time_of_flight_seconds);
                    }
                }
            }
        }
    }

    {
        // 中文注释：当前真正可用的是把 solve_problem1 输出转换成 root-solution branch 列表。
        const problem1::Problem1SolveInput input{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
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
        const auto branches = problem1::convert_problem1_candidates_to_solution_branches(
            planet_params::PlanetId::Mars,
            candidates);
        assert(branches.size() == candidates.size());

        for (std::size_t index = 1; index < branches.size(); ++index) {
            assert(branches[index - 1].time_of_flight_seconds <= branches[index].time_of_flight_seconds);
        }
        for (const auto& branch : branches) {
            assert(branch.valid);
            assert(std::isfinite(branch.encounter_global_angle));
            assert(branch.encounter_global_angle >= 0.0);
            assert(branch.encounter_global_angle < common::kTwoPi);
            assert(std::isfinite(branch.target_arrival_true_anomaly));
            assert(branch.target_arrival_true_anomaly >= 0.0);
            assert(branch.target_arrival_true_anomaly < common::kTwoPi);
            assert(branch.transfer_revolution >= 0);
            assert(branch.target_revolution >= 0);
            assert(std::isfinite(branch.time_of_flight_seconds));
            assert(std::isfinite(branch.target_time_seconds));
            assert(std::isfinite(branch.residual_seconds));
            assert(std::isfinite(branch.arrival_time_seconds_since_j2000));
        }
    }

    {
        // 中文注释：扩大一致性测试，覆盖多组 launch_time、transfer_perihelion_angle、planet pair、(max_k,max_q)。
        const double earth_period = planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
        const std::vector<std::tuple<planet_params::PlanetId, planet_params::PlanetId>> planet_pairs{
            {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars},
            {planet_params::PlanetId::Earth, planet_params::PlanetId::Venus},
        };
        const std::vector<double> launch_times{
            0.0,
            0.25 * earth_period,
            0.5 * earth_period,
        };
        const std::vector<double> transfer_perihelion_angles{0.2, 0.5, 1.0};
        const std::vector<int> max_revolution_values{0, 1};

        for (const auto& [departure_planet, target_planet] : planet_pairs) {
            for (double launch_time_seconds_since_j2000 : launch_times) {
                const planet_params::PlanetState departure_state =
                    planet_params::planet_state_at_time(departure_planet, launch_time_seconds_since_j2000);
                const planet_params::PlanetState target_state =
                    planet_params::planet_state_at_time(target_planet, launch_time_seconds_since_j2000);

                for (double transfer_perihelion_angle : transfer_perihelion_angles) {
                    const double theta_A = common::normalize_angle_0_2pi(
                        departure_state.theta_global - transfer_perihelion_angle);

                    for (int max_revolutions : max_revolution_values) {
                        auto defaults = spaceship_cpp::config::global_config().problem1_solve;
                        defaults.max_transfer_revolution = max_revolutions;
                        defaults.max_target_revolution = max_revolutions;
                        const auto input = spaceship_cpp::config::make_problem1_solve_input(
                            departure_planet,
                            target_planet,
                            launch_time_seconds_since_j2000,
                            transfer_perihelion_angle,
                            defaults);
                        const auto candidates = problem1::solve_problem1(input);
                        const auto direct = problem1::solve_problem1_from_departure_anomalies(
                            departure_planet,
                            target_planet,
                            departure_state.varphi,
                            target_state.varphi,
                            theta_A,
                            defaults.max_transfer_revolution,
                            defaults.max_target_revolution);
                        compare_solution_lists(target_planet, candidates, direct);
                    }
                }
            }
        }
    }

    {
        // 中文注释：k,q 只存在于 solution branch 里，不是 root-table 连续轴。
        problem1::Problem1SolutionBranch branch{};
        branch.transfer_revolution = 3;
        branch.target_revolution = 4;
        assert(branch.transfer_revolution == 3);
        assert(branch.target_revolution == 4);
    }

    {
        // 中文注释：root table 里 k,q 仍然只是每个 solution 的 label，而不是表格轴。
        const problem1::Problem1RootTableConfig config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kTwoPi,
            1,
            0.0,
            common::kTwoPi,
            1,
            0.5,
            common::kTwoPi,
            1,
            1,
            1,
            "problem1_root_table_draft_v0",
        };
        const problem1::Problem1RootTable table = problem1::build_problem1_root_table(config);
        const auto& cell = table.at(0, 0, 0);
        for (const auto& branch : cell.solutions_sorted_by_time_of_flight) {
            assert(branch.transfer_revolution >= 0);
            assert(branch.target_revolution >= 0);
            assert(std::isfinite(branch.encounter_global_angle));
            assert(branch.encounter_global_angle >= 0.0);
            assert(branch.encounter_global_angle < common::kTwoPi);
            assert(std::isfinite(branch.target_arrival_true_anomaly));
            assert(branch.target_arrival_true_anomaly >= 0.0);
            assert(branch.target_arrival_true_anomaly < common::kTwoPi);
            assert(std::isfinite(branch.time_of_flight_seconds));
            assert(std::isfinite(branch.target_time_seconds));
            const double residual_scale_seconds = std::max(
                {1.0, std::abs(branch.time_of_flight_seconds), std::abs(branch.target_time_seconds)});
            assert(std::abs(branch.residual_seconds) <= 1e-6 * residual_scale_seconds + 1e-6);
        }
    }

    return 0;
}
