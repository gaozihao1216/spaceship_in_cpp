/*
 * 文件作用：测试 Problem 1 表格理论约束。
 * 主要工作：检查表格几何、时间分支和有效性标记是否符合预期关系。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_table.hpp"

#include <cassert>
#include <cmath>

namespace {

bool approx_equal(double a, double b, double abs_tol = 1e-9, double rel_tol = 1e-12) {
    const double diff = std::abs(a - b);
    if (diff <= abs_tol) {
        return true;
    }
    return diff <= rel_tol * std::max(std::abs(a), std::abs(b));
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    // Elliptic-orbit counterexample:
    // same delta_lambda, different endpoint anomalies -> different endpoint radii and transfer geometry.
    {
        const double delta_lambda = 0.8;
        const problem1::Problem1TableCell case1 = problem1::evaluate_problem1_table_cell_geometry(
            1.0,
            0.1,
            1.2,
            0.1 + delta_lambda,
            0.4,
            0);
        const problem1::Problem1TableCell case2 = problem1::evaluate_problem1_table_cell_geometry(
            1.4,
            1.3,
            1.1,
            1.3 + delta_lambda,
            0.4,
            0);
        assert(approx_equal(case1.delta_global_angle, case2.delta_global_angle, 1e-12, 1e-12));
        assert(!approx_equal(case1.departure_radius, case2.departure_radius, 1e-6, 1e-12));
        assert(!approx_equal(case1.target_radius, case2.target_radius, 1e-6, 1e-12));
        assert(!approx_equal(case1.transfer_e_raw, case2.transfer_e_raw, 1e-6, 1e-12));
        assert(!approx_equal(case1.transfer_p, case2.transfer_p, 1e-6, 1e-12));
    }

    // Circular-limit degeneration:
    // if both radii are constant, same delta_lambda and theta_A imply the same transfer geometry.
    {
        const double theta_A = 0.6;
        const double delta_lambda = 1.1;
        const problem1::Problem1TableCell case1 = problem1::evaluate_problem1_table_cell_geometry(
            2.0,
            0.2,
            3.0,
            0.2 + delta_lambda,
            theta_A,
            0);
        const problem1::Problem1TableCell case2 = problem1::evaluate_problem1_table_cell_geometry(
            2.0,
            1.7,
            3.0,
            1.7 + delta_lambda,
            theta_A,
            0);
        assert(approx_equal(case1.transfer_e_raw, case2.transfer_e_raw, 1e-12, 1e-12));
        assert(approx_equal(case1.transfer_p, case2.transfer_p, 1e-12, 1e-12));
        assert(approx_equal(case1.transfer_theta_arrival, case2.transfer_theta_arrival, 1e-12, 1e-12));
    }

    // Launch-phase counterexample for q:
    // same nu_A at t and t + T_A, but target planet departure anomaly generally differs,
    // so target-time branches cannot be universal functions of (nu_A, nu_B_arrive, theta_A) alone.
    {
        const double t1 = 0.0;
        const double t2 = planet_params::planet_orbital_period(planet_params::PlanetId::Earth);
        const double nu_A_t1 = planet_params::planet_true_anomaly_at_time(planet_params::PlanetId::Earth, t1);
        const double nu_A_t2 = planet_params::planet_true_anomaly_at_time(planet_params::PlanetId::Earth, t2);
        const double nu_B_depart_t1 = planet_params::planet_true_anomaly_at_time(planet_params::PlanetId::Mars, t1);
        const double nu_B_depart_t2 = planet_params::planet_true_anomaly_at_time(planet_params::PlanetId::Mars, t2);

        assert(approx_equal(nu_A_t1, nu_A_t2, 1e-9, 1e-12));
        assert(!approx_equal(nu_B_depart_t1, nu_B_depart_t2, 1e-6, 1e-12));

        const problem1::Problem1TableConfig config{
            planet_params::PlanetId::Earth,
            planet_params::PlanetId::Mars,
            0.0,
            common::kTwoPi / 8.0,
            8,
            0.0,
            common::kTwoPi / 8.0,
            8,
            0.0,
            common::kTwoPi / 16.0,
            16,
            0,
            0,
        };
        const problem1::Problem1Table table = problem1::build_problem1_table(config);

        bool found_valid_cell = false;
        for (double nu_B_arrive : {0.2, 0.8, 1.4, 2.0, 2.6, 3.2}) {
            for (double theta_A : {0.3, 0.7, 1.1, 1.7, 2.3, 2.9}) {
                const problem1::Problem1TableQueryResult geometry_query =
                    problem1::query_problem1_table_exact(table, nu_A_t1, nu_B_arrive, theta_A);
                if (!geometry_query.cell.valid) {
                    continue;
                }

                const problem1::Problem1TimeOfFlightBranch branch_t1 =
                    problem1::evaluate_problem1_table_branch_with_target_departure_true_anomaly(
                        table,
                        geometry_query.cell,
                        nu_B_depart_t1,
                        0,
                        0);
                const problem1::Problem1TimeOfFlightBranch branch_t2 =
                    problem1::evaluate_problem1_table_branch_with_target_departure_true_anomaly(
                        table,
                        geometry_query.cell,
                        nu_B_depart_t2,
                        0,
                        0);
                if (!branch_t1.valid || !branch_t2.valid) {
                    continue;
                }

                found_valid_cell = true;
                assert(!approx_equal(
                    branch_t1.target_time_of_flight_seconds,
                    branch_t2.target_time_of_flight_seconds,
                    1e-6,
                    1e-12));
                break;
            }
            if (found_valid_cell) {
                break;
            }
        }

        assert(found_valid_cell);
    }

    return 0;
}
