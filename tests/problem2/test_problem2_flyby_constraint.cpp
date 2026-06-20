/*
 * 文件作用：测试 Problem 2 弹弓约束 F 函数、入射缓存与 branch 候选识别。
 */
#include "spaceship_cpp/problem2/problem2_flyby_constraint.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool approx_equal(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

}  // namespace

int main() {
    namespace problem2 = spaceship_cpp::problem2;

    const double phi = 1.0;
    const double e_J = 0.1;
    const double e = 0.3;
    const double theta_global = 0.4;

    {
        const auto F = problem2::evaluate_flyby_constraint_F(e, theta_global, phi, e_J);
        assert(F.valid);
        assert(std::isfinite(F.value));
        assert(F.denominator > 0.0);
    }

    {
        const auto residual = problem2::evaluate_flyby_constraint_residual(
            e,
            theta_global,
            e,
            theta_global,
            phi,
            e_J,
            false);
        assert(residual.valid);
        assert(std::abs(residual.residual) < 1e-12);
        assert(approx_equal(residual.incoming_F, residual.outgoing_F, 1e-12));
    }

    {
        const auto cache = problem2::build_flyby_constraint_incoming_cache(
            e,
            theta_global,
            phi,
            e_J,
            false);
        assert(cache.valid);
        assert(std::isfinite(cache.incoming_F));

        const auto cached_residual = problem2::evaluate_flyby_constraint_residual_from_incoming_cache(
            cache,
            e,
            theta_global,
            false);
        const auto direct_residual = problem2::evaluate_flyby_constraint_residual(
            e,
            theta_global,
            e,
            theta_global,
            phi,
            e_J,
            false);
        assert(cached_residual.valid);
        assert(direct_residual.valid);
        assert(approx_equal(cached_residual.residual, direct_residual.residual, 1e-15));
        assert(approx_equal(cached_residual.incoming_F, cache.incoming_F, 1e-15));
    }

    {
        const double theta_local = phi - theta_global;
        const auto local_residual = problem2::evaluate_flyby_constraint_residual(
            e,
            theta_local,
            e,
            theta_local,
            phi,
            e_J,
            true);
        assert(local_residual.valid);
        assert(std::abs(local_residual.residual) < 1e-12);
    }

    {
        const std::vector<problem2::FlybyThetaPrimeBranchSample> samples{
            {
                .theta_prime = 0.0,
                .branch_id = 0,
                .flyby_constraint_residual = 1.0,
                .flyby_constraint_valid = true,
            },
            {
                .theta_prime = 0.1,
                .branch_id = 0,
                .flyby_constraint_residual = -0.5,
                .flyby_constraint_valid = true,
            },
            {
                .theta_prime = 0.2,
                .branch_id = 0,
                .flyby_constraint_residual = 5e-4,
                .flyby_constraint_valid = true,
            },
        };
        const auto candidates =
            problem2::detect_flyby_theta_prime_candidates_from_branch_samples(samples, 1e-3);
        assert(candidates.size() == 2);

        bool has_sign_change = false;
        bool has_near_zero = false;
        for (const auto& candidate : candidates) {
            if (candidate.type == problem2::FlybyThetaPrimeCandidateType::SignChangeInterval) {
                has_sign_change = true;
                assert(candidate.branch_id == 0);
                assert(approx_equal(candidate.theta_prime_left, 0.0, 1e-12));
                assert(approx_equal(candidate.theta_prime_right, 0.1, 1e-12));
            }
            if (candidate.type == problem2::FlybyThetaPrimeCandidateType::NearZeroNode) {
                has_near_zero = true;
                assert(approx_equal(candidate.theta_prime_node, 0.2, 1e-12));
            }
        }
        assert(has_sign_change);
        assert(has_near_zero);
    }

    std::cout << "problem2_flyby_constraint_ok\n";
    return 0;
}
