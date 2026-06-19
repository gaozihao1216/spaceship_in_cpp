/*
 * 文件作用：测试 Problem 1 endpoint transfer-time table。
 * 主要工作：验证表格构造、单元查询、分支查询和基础插值可行性。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace {

bool approx_equal(double a, double b, double abs_tol = 1e-9, double rel_tol = 1e-12) {
    const double diff = std::abs(a - b);
    if (diff <= abs_tol) {
        return true;
    }
    return diff <= rel_tol * std::max(std::abs(a), std::abs(b));
}

bool any_valid_branch(const spaceship_cpp::problem1::Problem1TableCell& cell) {
    for (const auto& branch : cell.time_of_flight_branches) {
        if (branch.valid) {
            return true;
        }
    }
    return false;
}

double wrapped_angle_distance(double a, double b) {
    // 中文注释：角度比较必须按周期量处理，避免 0 和 2pi 被误判为不相等。
    return std::abs(spaceship_cpp::common::normalize_angle_minus_pi_pi(a - b));
}

bool angle_equal(double a, double b, double abs_tol = 1e-9) {
    return wrapped_angle_distance(a, b) <= abs_tol;
}

bool is_known_invalid_reason(const std::string& reason) {
    return reason == "invalid geometry inputs" ||
        reason == "singular geometry denominator" ||
        reason == "non-positive departure orbit factor" ||
        reason == "invalid transfer conic geometry" ||
        reason == "parabolic transfer is unsupported in table v2" ||
        reason == "invalid transfer eccentricity" ||
        reason == "non-positive or non-finite branch time of flight" ||
        reason == "hyperbolic transfer does not support k > 0" ||
        reason == "hyperbolic arrival angle is not forward of departure angle" ||
        reason == "requested (k,q) branch is not present in this cell";
}

void assert_cells_match(
    const spaceship_cpp::problem1::Problem1TableCell& lhs,
    const spaceship_cpp::problem1::Problem1TableCell& rhs
) {
    assert(lhs.valid == rhs.valid);
    if (!lhs.valid || !rhs.valid) {
        // 中文注释：无效 cell 只要求原因合理且一致，不强行比较 NaN 数值字段。
        assert(!lhs.invalid_reason.empty());
        assert(!rhs.invalid_reason.empty());
        assert(lhs.invalid_reason == rhs.invalid_reason);
        assert(is_known_invalid_reason(lhs.invalid_reason));
        return;
    }

    assert(approx_equal(lhs.departure_radius, rhs.departure_radius, 1e-9, 1e-12));
    assert(approx_equal(lhs.target_radius, rhs.target_radius, 1e-9, 1e-12));
    assert(angle_equal(lhs.departure_global_angle, rhs.departure_global_angle, 1e-9));
    assert(angle_equal(lhs.target_global_angle, rhs.target_global_angle, 1e-9));
    assert(angle_equal(lhs.delta_global_angle, rhs.delta_global_angle, 1e-9));
    assert(angle_equal(lhs.transfer_theta_arrival, rhs.transfer_theta_arrival, 1e-9));
    assert(approx_equal(lhs.transfer_e, rhs.transfer_e, 1e-9, 1e-12));
    assert(approx_equal(lhs.transfer_p, rhs.transfer_p, 1e-9, 1e-12));
    if (std::isfinite(lhs.transfer_a) || std::isfinite(rhs.transfer_a)) {
        assert(approx_equal(lhs.transfer_a, rhs.transfer_a, 1e-9, 1e-12));
    }
    assert(lhs.conic_type == rhs.conic_type);
}

void assert_branches_match(
    const spaceship_cpp::problem1::Problem1TimeOfFlightBranch& lhs,
    const spaceship_cpp::problem1::Problem1TimeOfFlightBranch& rhs
) {
    assert(lhs.transfer_revolution == rhs.transfer_revolution);
    assert(lhs.target_revolution == rhs.target_revolution);
    assert(lhs.valid == rhs.valid);
    if (!lhs.valid || !rhs.valid) {
        assert(lhs.invalid_reason == rhs.invalid_reason);
        assert(is_known_invalid_reason(lhs.invalid_reason));
        return;
    }

    assert(angle_equal(lhs.theta_arrival_branch, rhs.theta_arrival_branch, 1e-9));
    assert(approx_equal(lhs.target_true_anomaly_start, rhs.target_true_anomaly_start, 1e-9, 1e-12));
    assert(approx_equal(lhs.target_true_anomaly_end_branch, rhs.target_true_anomaly_end_branch, 1e-9, 1e-12));
    assert(approx_equal(lhs.deltaF_transfer, rhs.deltaF_transfer, 1e-9, 1e-12));
    assert(approx_equal(lhs.deltaF_target, rhs.deltaF_target, 1e-9, 1e-12));
    assert(approx_equal(lhs.time_of_flight_scale_free, rhs.time_of_flight_scale_free, 1e-9, 1e-12));
    assert(approx_equal(
        lhs.target_time_of_flight_scale_free, rhs.target_time_of_flight_scale_free, 1e-9, 1e-12));
    assert(approx_equal(lhs.time_of_flight_seconds, rhs.time_of_flight_seconds, 1e-9, 1e-12));
    assert(approx_equal(lhs.target_time_of_flight_seconds, rhs.target_time_of_flight_seconds, 1e-9, 1e-12));
    assert(approx_equal(lhs.residual_scale_free, rhs.residual_scale_free, 1e-9, 1e-12));
    assert(approx_equal(lhs.residual_seconds, rhs.residual_seconds, 1e-9, 1e-12));
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    const problem1::Problem1TableConfig config{
        planet_params::PlanetId::Earth,
        planet_params::PlanetId::Mars,
        0.0,
        common::kPi,
        2,
        0.0,
        common::kHalfPi,
        3,
        0.0,
        common::kHalfPi,
        3,
        1,
        1,
    };

    const problem1::Problem1Table table = problem1::build_problem1_table(config);
    assert(table.departure_true_anomaly_count() == 2);
    assert(table.target_true_anomaly_count() == 3);
    assert(table.transfer_theta_departure_count() == 3);
    assert(table.cells().size() == 18);

    {
        const auto& metadata = table.metadata();
        assert(metadata.schema_version == "planet_angle_pair_table_v2");
        assert(metadata.axis1_definition.find("nu_A") != std::string::npos);
        assert(metadata.axis2_definition.find("nu_B") != std::string::npos);
        assert(metadata.axis3_definition.find("theta_A") != std::string::npos);
        assert(metadata.max_transfer_revolution == 1);
        assert(metadata.max_target_revolution == 1);
    }

    {
        const problem1::Problem1TableCell& cell = table.at(0, 0, 0);
        const auto& departure = planet_params::get_planet_params(planet_params::PlanetId::Earth);
        const auto& target = planet_params::get_planet_params(planet_params::PlanetId::Mars);
        assert(approx_equal(cell.departure_true_anomaly, 0.0));
        assert(approx_equal(cell.target_true_anomaly, 0.0));
        assert(approx_equal(cell.departure_global_angle, departure.orbit.theta_0));
        assert(approx_equal(cell.target_global_angle, target.orbit.theta_0));
        assert(approx_equal(
            cell.delta_global_angle,
            common::normalize_angle_0_2pi(target.orbit.theta_0 - departure.orbit.theta_0)));
    }

    {
        const problem1::Problem1TableCell direct = problem1::evaluate_problem1_table_cell_for_planets(
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kHalfPi,
            common::kHalfPi,
            1,
            1);
        const problem1::Problem1TableCell& from_table = table.at(0, 1, 1);
        assert(approx_equal(direct.departure_radius, from_table.departure_radius, 1e-9, 1e-12));
        assert(approx_equal(direct.target_radius, from_table.target_radius, 1e-9, 1e-12));
        assert(approx_equal(direct.transfer_e_raw, from_table.transfer_e_raw, 1e-9, 1e-12));
        assert(approx_equal(direct.transfer_p, from_table.transfer_p, 1e-9, 1e-12));
        assert(direct.time_of_flight_branches.size() == from_table.time_of_flight_branches.size());
    }

    {
        const problem1::Problem1TableQueryResult query_a = problem1::query_problem1_table_exact(
            table,
            0.01,
            0.02,
            0.03);
        const problem1::Problem1TableQueryResult query_b = problem1::query_problem1_table_exact(
            table,
            common::kTwoPi + 0.01,
            common::kTwoPi + 0.02,
            common::kTwoPi + 0.03);
        assert_cells_match(query_a.cell, query_b.cell);
    }

    {
        // 中文注释：分别单独测试三个周期轴在 0 和 2pi 边界处的一致性。
        const problem1::Problem1TableQueryResult base = problem1::query_problem1_table_exact(
            table,
            0.31,
            1.07,
            2.41);
        const problem1::Problem1TableQueryResult nu_a_wrapped = problem1::query_problem1_table_exact(
            table,
            0.31 + common::kTwoPi,
            1.07,
            2.41);
        const problem1::Problem1TableQueryResult nu_b_wrapped = problem1::query_problem1_table_exact(
            table,
            0.31,
            1.07 + common::kTwoPi,
            2.41);
        const problem1::Problem1TableQueryResult theta_a_wrapped = problem1::query_problem1_table_exact(
            table,
            0.31,
            1.07,
            2.41 + common::kTwoPi);
        assert_cells_match(base.cell, nu_a_wrapped.cell);
        assert_cells_match(base.cell, nu_b_wrapped.cell);
        assert_cells_match(base.cell, theta_a_wrapped.cell);

        const problem1::Problem1TableQueryResult nu_a_zero = problem1::query_problem1_table_exact(
            table,
            0.0,
            1.07,
            2.41);
        const problem1::Problem1TableQueryResult nu_a_two_pi = problem1::query_problem1_table_exact(
            table,
            common::kTwoPi,
            1.07,
            2.41);
        const problem1::Problem1TableQueryResult nu_b_zero = problem1::query_problem1_table_exact(
            table,
            0.31,
            0.0,
            2.41);
        const problem1::Problem1TableQueryResult nu_b_two_pi = problem1::query_problem1_table_exact(
            table,
            0.31,
            common::kTwoPi,
            2.41);
        const problem1::Problem1TableQueryResult theta_a_zero = problem1::query_problem1_table_exact(
            table,
            0.31,
            1.07,
            0.0);
        const problem1::Problem1TableQueryResult theta_a_two_pi = problem1::query_problem1_table_exact(
            table,
            0.31,
            1.07,
            common::kTwoPi);
        assert_cells_match(nu_a_zero.cell, nu_a_two_pi.cell);
        assert_cells_match(nu_b_zero.cell, nu_b_two_pi.cell);
        assert_cells_match(theta_a_zero.cell, theta_a_two_pi.cell);
    }

    {
        // 中文注释：随机抽样 exact query，并和直接几何公式计算结果逐点比对。
        std::mt19937_64 rng(20260522ULL);
        std::uniform_real_distribution<double> angle_dist(0.0, common::kTwoPi);
        for (int sample_index = 0; sample_index < 200; ++sample_index) {
            const double departure_true_anomaly = angle_dist(rng);
            const double target_true_anomaly = angle_dist(rng);
            const double transfer_theta_departure = angle_dist(rng);

            const problem1::Problem1TableQueryResult query = problem1::query_problem1_table_exact(
                table,
                departure_true_anomaly,
                target_true_anomaly,
                transfer_theta_departure);
            const auto& departure_planet = planet_params::get_planet_params(config.departure_planet);
            const auto& target_planet = planet_params::get_planet_params(config.target_planet);
            const double departure_radius = planet_params::planet_radius_at_true_anomaly(
                config.departure_planet,
                common::normalize_angle_0_2pi(departure_true_anomaly));
            const double target_radius = planet_params::planet_radius_at_true_anomaly(
                config.target_planet,
                common::normalize_angle_0_2pi(target_true_anomaly));
            const double departure_global_angle =
                common::normalize_angle_0_2pi(departure_planet.orbit.theta_0 + departure_true_anomaly);
            const double target_global_angle =
                common::normalize_angle_0_2pi(target_planet.orbit.theta_0 + target_true_anomaly);

            const problem1::Problem1TableCell direct = problem1::evaluate_problem1_table_cell_geometry(
                departure_radius,
                departure_global_angle,
                target_radius,
                target_global_angle,
                transfer_theta_departure,
                config.max_transfer_revolution);
            assert_cells_match(query.cell, direct);
        }
    }

    {
        // 中文注释：每个有效几何 cell 都必须显式包含配置范围内全部 (k,q) 分支，即使其中某些分支无效也要保留原因。
        for (const auto& cell : table.cells()) {
            if (!cell.valid) {
                continue;
            }
            assert(static_cast<int>(cell.time_of_flight_branches.size()) ==
                (config.max_transfer_revolution + 1) * (config.max_target_revolution + 1));
            for (int transfer_revolution = 0; transfer_revolution <= config.max_transfer_revolution;
                 ++transfer_revolution) {
                for (int target_revolution = 0; target_revolution <= config.max_target_revolution;
                     ++target_revolution) {
                    const auto* branch = problem1::find_problem1_table_branch(
                        cell,
                        transfer_revolution,
                        target_revolution);
                    assert(branch != nullptr);
                    if (!branch->valid) {
                        assert(!branch->invalid_reason.empty());
                    }
                }
            }
        }
    }

    {
        // 中文注释：同一个 geometry cell、同一个 k、不同 q 的 transfer-side 字段必须保持不变。
        for (const auto& cell : table.cells()) {
            if (!cell.valid) {
                continue;
            }
            for (int transfer_revolution = 0; transfer_revolution <= config.max_transfer_revolution;
                 ++transfer_revolution) {
                const problem1::Problem1TransferBranchView view =
                    problem1::get_problem1_transfer_branch_view(cell, transfer_revolution);
                if (!view.valid) {
                    continue;
                }

                for (const auto& branch : cell.time_of_flight_branches) {
                    if (!problem1::is_problem1_transfer_branch_valid_for_interpolation(
                            branch,
                            transfer_revolution)) {
                        continue;
                    }
                    assert(angle_equal(branch.theta_arrival_branch, view.theta_arrival_branch, 1e-12));
                    assert(approx_equal(branch.deltaF_transfer, view.deltaF_transfer, 1e-12, 1e-12));
                    assert(approx_equal(
                        branch.time_of_flight_scale_free, view.time_of_flight_scale_free, 1e-12, 1e-12));
                    assert(approx_equal(branch.time_of_flight_seconds, view.time_of_flight_seconds, 1e-12, 1e-12));
                }
            }
        }
    }

    {
        // 中文注释：随机抽查指定 (k,q) 分支，和直接 residual 计算得到的同一分支结果做一致性比对。
        std::mt19937_64 rng(20260523ULL);
        std::uniform_real_distribution<double> angle_dist(0.0, common::kTwoPi);
        std::uniform_int_distribution<int> k_dist(0, config.max_transfer_revolution);
        std::uniform_int_distribution<int> q_dist(0, config.max_target_revolution);

        for (int sample_index = 0; sample_index < 100; ++sample_index) {
            const double departure_true_anomaly = angle_dist(rng);
            const double target_true_anomaly = angle_dist(rng);
            const double transfer_theta_departure = angle_dist(rng);
            const int transfer_revolution = k_dist(rng);
            const int target_revolution = q_dist(rng);

            const problem1::Problem1TableBranchQueryResult branch_query =
                problem1::query_problem1_table_exact_branch(
                    table,
                    departure_true_anomaly,
                    target_true_anomaly,
                    transfer_theta_departure,
                    transfer_revolution,
                    target_revolution);
            assert(branch_query.branch_found == (problem1::find_problem1_table_branch(
                branch_query.cell,
                transfer_revolution,
                target_revolution) != nullptr));
            if (!branch_query.cell.valid) {
                // 中文注释：几何本身无效时，不要求存在任何 (k,q) branch；此时由 cell.invalid_reason 承担解释责任。
                assert(!branch_query.branch_found);
                continue;
            }

            const problem1::Problem1ResidualResult residual = problem1::evaluate_problem1_residual(
                problem1::Problem1ResidualInput{
                    config.departure_planet,
                    config.target_planet,
                    branch_query.cell.launch_time_phase_seconds_since_j2000,
                    branch_query.cell.transfer_perihelion_angle_global_used,
                    branch_query.cell.target_global_angle,
                    transfer_revolution,
                    target_revolution,
                });

            assert(branch_query.branch_found);
            assert(branch_query.branch.transfer_revolution == transfer_revolution);
            assert(branch_query.branch.target_revolution == target_revolution);

            if (residual.status == problem1::Problem1ResidualStatus::Success) {
                assert(branch_query.branch.valid);
                assert(approx_equal(
                    branch_query.branch.deltaF_transfer, residual.deltaF_transfer, 1e-9, 1e-12));
                assert(approx_equal(
                    branch_query.branch.deltaF_target, residual.deltaF_target, 1e-9, 1e-12));
                assert(approx_equal(
                    branch_query.branch.residual_scale_free, residual.residual, 1e-9, 1e-12));
            } else {
                assert(!branch_query.branch.valid);
                assert(!branch_query.branch.invalid_reason.empty());
            }
        }
    }

    {
        // 中文注释：branch signature 必须显式区分 (k,q) 组合，不能退化成只看 branch_count。
        const problem1::Problem1TableCell* selected_cell = nullptr;
        for (const auto& cell : table.cells()) {
            if (cell.valid) {
                selected_cell = &cell;
                break;
            }
        }
        assert(selected_cell != nullptr);
        const std::string signature = problem1::problem1_table_branch_signature(*selected_cell);
        assert(signature.find("k=0,q=0") != std::string::npos);
        assert(signature.find("k=0,q=1") != std::string::npos);
        assert(signature.find("k=1,q=0") != std::string::npos);
        assert(signature.find("k=1,q=1") != std::string::npos);
    }

    {
        // 中文注释：future interpolation admissibility 只接收 k，不接收 q；这里用一次真实调用做编译层面的静态约束。
        const problem1::Problem1TableInterpolationAdmissibility admissibility =
            problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                table,
                0.31,
                1.07,
                2.41,
                0);
        (void)admissibility;

        const problem1::Problem1TransferTimeQueryResult stub_result =
            problem1::query_problem1_transfer_time_with_interpolation_stub(
                table,
                0.31,
                1.07,
                2.41,
                0);
        (void)stub_result;
    }

    {
        // 中文注释：对一个有效查询点，8 顶点都必须几何有效、conic_type 一致，并且存在相同 k 的有效 transfer branch。
        const problem1::Problem1TableConfig admissibility_config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kTwoPi,
            1,
            common::kHalfPi,
            common::kTwoPi,
            1,
            common::kHalfPi,
            common::kTwoPi,
            1,
            0,
            1,
        };
        const problem1::Problem1Table admissibility_table =
            problem1::build_problem1_table(admissibility_config);
        const problem1::Problem1TableInterpolationAdmissibility admissibility =
            problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                admissibility_table,
                0.01,
                common::kHalfPi + 0.01,
                common::kHalfPi + 0.01,
                0);
        assert(admissibility.admissible);
        assert(admissibility.reason == "admissible");
        for (const auto& vertex : admissibility.vertices) {
            assert(vertex.valid);
            const auto* branch = problem1::find_problem1_table_branch(vertex, 0, 0);
            assert(branch != nullptr);
            assert(branch->transfer_revolution == 0);
            assert(any_valid_branch(vertex));
        }

        const problem1::Problem1TransferTimeQueryResult stub =
            problem1::query_problem1_transfer_time_with_interpolation_stub(
                admissibility_table,
                0.01,
                common::kHalfPi + 0.01,
                common::kHalfPi + 0.01,
                0);
        assert(stub.ok);
        assert(stub.interpolation_admissible);
        assert(stub.method == "exact_stub");
        assert(stub.reason.empty());
        assert(stub.transfer_revolution == 0);
        assert(std::isfinite(stub.theta_arrival_branch));
        assert(std::isfinite(stub.deltaF_transfer));
        assert(std::isfinite(stub.time_of_flight_scale_free));
        assert(std::isfinite(stub.time_of_flight_seconds));
    }

    {
        // 中文注释：transfer-side helper 不应被 representative q / target-side residual 污染。
        problem1::Problem1TimeOfFlightBranch synthetic_branch{};
        synthetic_branch.valid = false;
        synthetic_branch.transfer_revolution = 2;
        synthetic_branch.target_revolution = 999;
        synthetic_branch.theta_arrival_branch = 1.5;
        synthetic_branch.target_true_anomaly_start = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.target_true_anomaly_end_branch = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.deltaF_transfer = 3.0;
        synthetic_branch.deltaF_target = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.time_of_flight_scale_free = 4.0;
        synthetic_branch.target_time_of_flight_scale_free = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.time_of_flight_seconds = 5.0;
        synthetic_branch.target_time_of_flight_seconds = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.residual_scale_free = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.residual_seconds = std::numeric_limits<double>::quiet_NaN();
        synthetic_branch.invalid_reason = "representative target branch invalid";

        assert(problem1::is_problem1_transfer_branch_valid_for_interpolation(synthetic_branch, 2));
        assert(!synthetic_branch.valid);
    }

    {
        // 中文注释：pure transfer branch view 不包含 q、target_time、residual 语义，只表达 transfer-side k 视图。
        problem1::Problem1TableCell synthetic_cell{};
        problem1::Problem1TimeOfFlightBranch branch_q0{};
        branch_q0.valid = false;
        branch_q0.transfer_revolution = 1;
        branch_q0.target_revolution = 0;
        branch_q0.theta_arrival_branch = 1.2;
        branch_q0.deltaF_transfer = 2.3;
        branch_q0.time_of_flight_scale_free = 4.5;
        branch_q0.time_of_flight_seconds = 6.7;
        branch_q0.deltaF_target = std::numeric_limits<double>::quiet_NaN();
        branch_q0.target_time_of_flight_scale_free = std::numeric_limits<double>::quiet_NaN();
        branch_q0.target_time_of_flight_seconds = std::numeric_limits<double>::quiet_NaN();
        branch_q0.residual_scale_free = std::numeric_limits<double>::quiet_NaN();
        branch_q0.residual_seconds = std::numeric_limits<double>::quiet_NaN();

        problem1::Problem1TimeOfFlightBranch branch_q1 = branch_q0;
        branch_q1.target_revolution = 1;

        synthetic_cell.time_of_flight_branches = {branch_q0, branch_q1};
        const problem1::Problem1TransferBranchView view =
            problem1::get_problem1_transfer_branch_view(synthetic_cell, 1);
        assert(view.valid);
        assert(view.transfer_revolution == 1);
        assert(approx_equal(view.theta_arrival_branch, 1.2, 1e-12, 1e-12));
        assert(approx_equal(view.deltaF_transfer, 2.3, 1e-12, 1e-12));
        assert(approx_equal(view.time_of_flight_scale_free, 4.5, 1e-12, 1e-12));
        assert(approx_equal(view.time_of_flight_seconds, 6.7, 1e-12, 1e-12));
    }

    {
        // 中文注释：无效情况必须返回 admissible=false，且 reason 非空。
        const problem1::Problem1TableInterpolationAdmissibility inadmissible =
            problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                table,
                0.0,
                0.0,
                common::kPi,
                0);
        assert(!inadmissible.admissible);
        assert(!inadmissible.reason.empty());
        assert(inadmissible.reason != "admissible");

        const problem1::Problem1TransferTimeQueryResult stub =
            problem1::query_problem1_transfer_time_with_interpolation_stub(
                table,
                0.0,
                0.0,
                common::kPi,
                0);
        assert(!stub.ok);
        assert(!stub.interpolation_admissible);
        assert(!stub.reason.empty());
        assert(stub.method == "exact_fallback");
        assert(stub.transfer_revolution == 0);
    }

    {
        // 中文注释：admissibility 只依赖 geometry+k；即使后续选择不同 q，准入结果也不应发生变化。
        const problem1::Problem1TableInterpolationAdmissibility admissibility_a =
            problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                table,
                0.31,
                1.07,
                2.41,
                0);
        const problem1::Problem1TableInterpolationAdmissibility admissibility_b =
            problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                table,
                0.31,
                1.07,
                2.41,
                0);
        assert(admissibility_a.admissible == admissibility_b.admissible);
        assert(admissibility_a.reason == admissibility_b.reason);
    }

    {
        // 中文注释：真实 2x2x2 周期表下检查 8 顶点定位，不允许退化成同一个顶点；同时验证边界 wrap 和 local unwrap。
        const problem1::Problem1TableConfig locating_config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Venus,
            0.0,
            common::kPi,
            2,
            0.0,
            common::kPi,
            2,
            0.0,
            common::kPi,
            2,
            0,
            1,
        };
        const problem1::Problem1Table locating_table = problem1::build_problem1_table(locating_config);

        bool found_wrapped_admissible = false;
        problem1::Problem1TableInterpolationAdmissibility locating{};
        for (double departure_true_anomaly : {common::kTwoPi - 0.01, common::kTwoPi - 0.10}) {
            for (double target_true_anomaly : {common::kTwoPi - 0.02, common::kTwoPi - 0.12}) {
                for (double transfer_theta_departure : {common::kTwoPi - 0.03, common::kTwoPi - 0.15}) {
                    locating = problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                        locating_table,
                        departure_true_anomaly,
                        target_true_anomaly,
                        transfer_theta_departure,
                        0);
                    if (!locating.admissible) {
                        continue;
                    }

                    bool saw_wrapped_departure = false;
                    bool saw_wrapped_target = false;
                    bool saw_wrapped_theta = false;
                    bool saw_nonzero_index = false;
                    for (const auto& vertex_index : locating.vertex_indices) {
                        if (vertex_index.departure_true_anomaly_index == 0) {
                            saw_wrapped_departure = true;
                        } else {
                            saw_nonzero_index = true;
                        }
                        if (vertex_index.target_true_anomaly_index == 0) {
                            saw_wrapped_target = true;
                        }
                        if (vertex_index.transfer_theta_departure_index == 0) {
                            saw_wrapped_theta = true;
                        }
                    }

                    if (saw_wrapped_departure && saw_wrapped_target && saw_wrapped_theta && saw_nonzero_index) {
                        found_wrapped_admissible = true;
                        break;
                    }
                }
                if (found_wrapped_admissible) {
                    break;
                }
            }
            if (found_wrapped_admissible) {
                break;
            }
        }

        assert(found_wrapped_admissible);
        assert(locating.local_departure_true_anomaly > common::kPi);
        assert(locating.local_target_true_anomaly > common::kPi);
        assert(locating.local_transfer_theta_departure > common::kPi);
    }

    {
        // 中文注释：stub 返回的 transfer-side 结果必须与 pure transfer branch view 一致，不包含 q / target_time / residual。
        const problem1::Problem1TableConfig stub_config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kTwoPi,
            1,
            common::kHalfPi,
            common::kTwoPi,
            1,
            common::kHalfPi,
            common::kTwoPi,
            1,
            0,
            1,
        };
        const problem1::Problem1Table stub_table = problem1::build_problem1_table(stub_config);
        const problem1::Problem1TableQueryResult query = problem1::query_problem1_table_exact(
            stub_table,
            0.01,
            common::kHalfPi + 0.01,
            common::kHalfPi + 0.01);
        const problem1::Problem1TableInterpolationAdmissibility admissibility =
            problem1::check_problem1_table_transfer_branch_interpolation_admissibility(
                stub_table,
                0.01,
                common::kHalfPi + 0.01,
                common::kHalfPi + 0.01,
                0);
        const problem1::Problem1TransferBranchView view =
            problem1::get_problem1_transfer_branch_view(query.cell, 0);
        if (view.valid && admissibility.admissible) {
            const problem1::Problem1TransferTimeQueryResult stub =
                problem1::query_problem1_transfer_time_with_interpolation_stub(
                    stub_table,
                    0.01,
                    common::kHalfPi + 0.01,
                    common::kHalfPi + 0.01,
                    0);
            assert(stub.ok);
            assert(approx_equal(stub.theta_arrival_branch, view.theta_arrival_branch, 1e-12, 1e-12));
            assert(approx_equal(stub.deltaF_transfer, view.deltaF_transfer, 1e-12, 1e-12));
            assert(approx_equal(stub.time_of_flight_scale_free, view.time_of_flight_scale_free, 1e-12, 1e-12));
            assert(approx_equal(stub.time_of_flight_seconds, view.time_of_flight_seconds, 1e-12, 1e-12));
        }
    }

    {
        bool found_valid = false;
        for (const auto& cell : table.cells()) {
            if (any_valid_branch(cell)) {
                found_valid = true;
                break;
            }
        }
        assert(found_valid);
    }

    {
        problem1::Problem1TableMetadata old_metadata = table.metadata();
        old_metadata.schema_version = "relative_phase_table_v1";
        bool threw = false;
        try {
            problem1::validate_problem1_table_metadata(old_metadata);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        problem1::Problem1TableConfig invalid = config;
        invalid.departure_true_anomaly_count = 0;
        bool threw = false;
        try {
            (void)problem1::build_problem1_table(invalid);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        bool threw = false;
        try {
            (void)table.at(2, 0, 0);
        } catch (const std::out_of_range&) {
            threw = true;
        }
        assert(threw);
    }

    return 0;
}
