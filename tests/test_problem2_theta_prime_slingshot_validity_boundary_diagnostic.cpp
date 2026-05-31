#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_nearest_node_query.hpp"
#include "spaceship_cpp/problem2/problem2_slingshot.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::is_finite;
using spaceship_cpp::common::kTwoPi;
using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kThetaSampleCount = 32;
constexpr double kGeometryDenominatorEpsilon = 1e-12;
constexpr double kBackSubstitutionTolerance = 1e-10;

struct ResidualBranch {
    bool problem1_valid = false;
    bool slingshot_valid = false;
    std::string problem1_invalid_reason;
    std::string slingshot_invalid_reason;
    double theta_prime = 0.0;
    double alpha = 0.0;
    int transfer_revolution = 0;
    int target_revolution = 0;
    double time_of_flight_seconds = 0.0;
    double target_time_seconds = 0.0;
    double problem1_residual_seconds = 0.0;
    double slingshot_residual = 0.0;
    double outgoing_eccentricity = 0.0;
    double outgoing_semi_latus_rectum = 0.0;
};

struct GeometryDiagnostic {
    bool valid_geometry = false;
    std::string invalid_reason;
    double r0 = std::numeric_limits<double>::quiet_NaN();
    double r1 = std::numeric_limits<double>::quiet_NaN();
    double denominator = std::numeric_limits<double>::quiet_NaN();
    double e_prime = std::numeric_limits<double>::quiet_NaN();
    double R_prime = std::numeric_limits<double>::quiet_NaN();
    double B_encounter = std::numeric_limits<double>::quiet_NaN();
    double B_target = std::numeric_limits<double>::quiet_NaN();
    double encounter_reconstructed_radius = std::numeric_limits<double>::quiet_NaN();
    double target_reconstructed_radius = std::numeric_limits<double>::quiet_NaN();
    double encounter_radius_error = std::numeric_limits<double>::quiet_NaN();
    double target_radius_error = std::numeric_limits<double>::quiet_NaN();
};

struct BranchWithGeometry {
    ResidualBranch branch;
    GeometryDiagnostic geometry;
};

struct PairedBranch {
    BranchWithGeometry direct;
    BranchWithGeometry table;
    int rank = 0;
};

struct Sample {
    double theta_prime = 0.0;
    std::vector<BranchWithGeometry> branches;
    bool fallback_used = false;
};

struct Summary {
    int analyzed_interval_count = 0;
    int analyzed_branch_pair_count = 0;
    double direct_table_alpha_max_diff = 0.0;
    double direct_table_time_max_diff = 0.0;
    double direct_table_denominator_max_diff = 0.0;
    double direct_table_R_prime_max_diff = 0.0;
    int near_zero_R_prime_count = 0;
    int near_zero_B_target_count = 0;
    int near_geometry_tolerance_count = 0;
    int validity_flip_count = 0;
    std::map<std::string, int> flip_reason_count;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/tmp/problem1_root_table_2deg_full";
}

bool branch_less(const BranchWithGeometry& lhs, const BranchWithGeometry& rhs) {
    return std::tie(lhs.branch.transfer_revolution, lhs.branch.target_revolution, lhs.branch.time_of_flight_seconds) <
           std::tie(rhs.branch.transfer_revolution, rhs.branch.target_revolution, rhs.branch.time_of_flight_seconds);
}

std::string side_name(int endpoint_index, int interval_left) {
    return endpoint_index == interval_left ? "left" : "right";
}

GeometryDiagnostic compute_geometry_diagnostic(
    double R_J,
    double e_J,
    double R_K,
    double e_K,
    double phi,
    double alpha,
    double theta_prime
) {
    namespace problem2 = spaceship_cpp::problem2;
    GeometryDiagnostic g{};
    g.r0 = problem2::problem2_orbit_radius(R_J, e_J, phi);
    if (!is_finite(g.r0) || !(g.r0 > 0.0)) {
        g.invalid_reason = "invalid_encounter_radius";
        return g;
    }
    g.r1 = problem2::problem2_orbit_radius(R_K, e_K, alpha);
    if (!is_finite(g.r1) || !(g.r1 > 0.0)) {
        g.invalid_reason = "invalid_target_radius";
        return g;
    }

    const double cos_encounter = std::cos(phi - theta_prime);
    const double cos_target = std::cos(alpha - theta_prime);
    g.denominator = g.r0 * cos_encounter - g.r1 * cos_target;
    if (!is_finite(g.denominator) || std::abs(g.denominator) <= kGeometryDenominatorEpsilon) {
        g.invalid_reason = "geometry_denominator_too_small";
        return g;
    }

    g.e_prime = (g.r1 - g.r0) / g.denominator;
    if (!is_finite(g.e_prime)) {
        g.invalid_reason = "non_finite_outgoing_eccentricity";
        return g;
    }

    g.B_encounter = 1.0 + g.e_prime * cos_encounter;
    g.R_prime = g.r0 * g.B_encounter;
    if (!is_finite(g.R_prime) || !(g.R_prime > 0.0)) {
        g.invalid_reason = "non_positive_outgoing_semi_latus_rectum";
        return g;
    }
    if (!is_finite(g.B_encounter) || !(g.B_encounter > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }

    g.B_target = 1.0 + g.e_prime * cos_target;
    if (!is_finite(g.B_target) || !(g.B_target > 0.0)) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }

    g.encounter_reconstructed_radius = g.R_prime / g.B_encounter;
    g.encounter_radius_error = g.encounter_reconstructed_radius - g.r0;
    if (!is_finite(g.encounter_reconstructed_radius) ||
        std::abs(g.encounter_radius_error) > kBackSubstitutionTolerance) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_encounter";
        return g;
    }

    g.target_reconstructed_radius = g.R_prime / g.B_target;
    g.target_radius_error = g.target_reconstructed_radius - g.r1;
    if (!is_finite(g.target_reconstructed_radius) ||
        std::abs(g.target_radius_error) > kBackSubstitutionTolerance) {
        g.invalid_reason = "outgoing_orbit_does_not_pass_target";
        return g;
    }

    g.valid_geometry = true;
    return g;
}

std::vector<BranchWithGeometry> convert_branches(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double phi,
    double incoming_e,
    double incoming_theta,
    double theta_prime,
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem2 = spaceship_cpp::problem2;
    const auto& departure_params = planet_params::get_planet_params(departure_planet);
    const auto& target_params = planet_params::get_planet_params(target_planet);

    std::vector<BranchWithGeometry> out;
    out.reserve(branches.size());
    for (const auto& branch : branches) {
        BranchWithGeometry item{};
        auto& rb = item.branch;
        rb.theta_prime = theta_prime;
        rb.alpha = branch.target_arrival_true_anomaly;
        rb.transfer_revolution = branch.transfer_revolution;
        rb.target_revolution = branch.target_revolution;
        rb.time_of_flight_seconds = branch.time_of_flight_seconds;
        rb.target_time_seconds = branch.target_time_seconds;
        rb.problem1_residual_seconds = branch.residual_seconds;
        item.geometry = compute_geometry_diagnostic(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            phi,
            rb.alpha,
            theta_prime);
        if (!branch.valid) {
            rb.problem1_invalid_reason =
                branch.invalid_reason.empty() ? "problem1_branch_invalid" : branch.invalid_reason;
            out.push_back(item);
            continue;
        }
        rb.problem1_valid = true;
        const auto residual = problem2::evaluate_problem2_slingshot_residual_from_theta_alpha(
            departure_params.orbit.p,
            departure_params.orbit.e,
            target_params.orbit.p,
            target_params.orbit.e,
            phi,
            rb.alpha,
            incoming_e,
            incoming_theta,
            theta_prime);
        if (!residual.valid) {
            rb.slingshot_invalid_reason =
                residual.invalid_reason.empty() ? "slingshot_invalid" : residual.invalid_reason;
            out.push_back(item);
            continue;
        }
        rb.slingshot_valid = true;
        rb.slingshot_residual = residual.slingshot_residual;
        rb.outgoing_eccentricity = residual.outgoing_eccentricity;
        rb.outgoing_semi_latus_rectum = residual.outgoing_semi_latus_rectum;
        out.push_back(item);
    }
    std::sort(out.begin(), out.end(), branch_less);
    return out;
}

Sample build_direct_sample(
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto branches = problem1::solve_problem1_from_departure_anomalies(
        departure_planet, target_planet, phi, beta, theta_A, 1, 1);
    return Sample{theta_prime, convert_branches(departure_planet, target_planet, phi, incoming_e, incoming_theta,
                                                theta_prime, branches), false};
}

Sample build_table_sample(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const spaceship_cpp::problem1::Problem1NearestNodeQueryOptions& options,
    spaceship_cpp::planet_params::PlanetId departure_planet,
    spaceship_cpp::planet_params::PlanetId target_planet,
    double encounter_time,
    double phi,
    double beta,
    double incoming_e,
    double incoming_theta,
    double theta_prime
) {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const double theta_A = normalize_angle_0_2pi(departure_state.theta_global - theta_prime);
    const auto query = problem1::query_problem1_from_2deg_nearest_node(
        loader, departure_planet, target_planet, phi, beta, theta_A, 1, 1, options);
    return Sample{theta_prime, convert_branches(departure_planet, target_planet, phi, incoming_e, incoming_theta,
                                                theta_prime, query.branches), query.used_direct_solve_fallback};
}

std::vector<PairedBranch> pair_direct_table_branches(const Sample& direct, const Sample& table) {
    std::map<std::pair<int, int>, std::vector<BranchWithGeometry>> direct_by_kq;
    std::map<std::pair<int, int>, std::vector<BranchWithGeometry>> table_by_kq;
    for (const auto& branch : direct.branches) {
        if (branch.branch.problem1_valid) {
            direct_by_kq[{branch.branch.transfer_revolution, branch.branch.target_revolution}].push_back(branch);
        }
    }
    for (const auto& branch : table.branches) {
        if (branch.branch.problem1_valid) {
            table_by_kq[{branch.branch.transfer_revolution, branch.branch.target_revolution}].push_back(branch);
        }
    }
    std::vector<PairedBranch> pairs;
    int rank = 0;
    for (auto& [key, direct_group] : direct_by_kq) {
        auto table_it = table_by_kq.find(key);
        if (table_it == table_by_kq.end()) {
            continue;
        }
        auto& table_group = table_it->second;
        std::sort(direct_group.begin(), direct_group.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds;
        });
        std::sort(table_group.begin(), table_group.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.branch.time_of_flight_seconds < rhs.branch.time_of_flight_seconds;
        });
        const int pair_count = std::min<int>(direct_group.size(), table_group.size());
        for (int i = 0; i < pair_count; ++i) {
            pairs.push_back(PairedBranch{
                direct_group[static_cast<std::size_t>(i)],
                table_group[static_cast<std::size_t>(i)],
                rank++});
        }
    }
    return pairs;
}

void print_branch_case(
    int interval_left,
    const char* side,
    const char* source,
    int rank,
    const BranchWithGeometry& item
) {
    const auto& b = item.branch;
    std::cout << "Problem2SlingshotValidityBoundaryBranchCase"
              << " interval=" << interval_left << "-" << (interval_left + 1)
              << " side=" << side
              << " source=" << source
              << " rank=" << rank
              << " k=" << b.transfer_revolution
              << " q=" << b.target_revolution
              << " theta_prime=" << b.theta_prime
              << " alpha=" << b.alpha
              << " time_of_flight_seconds=" << b.time_of_flight_seconds
              << " problem1_residual_seconds=" << b.problem1_residual_seconds
              << " problem1_valid=" << (b.problem1_valid ? 1 : 0)
              << " slingshot_valid=" << (b.slingshot_valid ? 1 : 0)
              << " slingshot_invalid_reason=" << b.slingshot_invalid_reason
              << " slingshot_residual=" << b.slingshot_residual
              << " outgoing_eccentricity=" << b.outgoing_eccentricity
              << " outgoing_semi_latus_rectum=" << b.outgoing_semi_latus_rectum
              << '\n';
}

void print_geometry_case(
    int interval_left,
    const char* side,
    const char* source,
    int rank,
    const BranchWithGeometry& item
) {
    const auto& b = item.branch;
    const auto& g = item.geometry;
    std::cout << "Problem2SlingshotGeometryBoundaryDiagnostic"
              << " interval=" << interval_left << "-" << (interval_left + 1)
              << " side=" << side
              << " source=" << source
              << " rank=" << rank
              << " k=" << b.transfer_revolution
              << " q=" << b.target_revolution
              << " theta_prime=" << b.theta_prime
              << " alpha=" << b.alpha
              << " r0=" << g.r0
              << " r1=" << g.r1
              << " denominator=" << g.denominator
              << " abs_denominator=" << std::abs(g.denominator)
              << " e_prime=" << g.e_prime
              << " R_prime=" << g.R_prime
              << " B_encounter=" << g.B_encounter
              << " B_target=" << g.B_target
              << " encounter_reconstructed_radius=" << g.encounter_reconstructed_radius
              << " target_reconstructed_radius=" << g.target_reconstructed_radius
              << " encounter_radius_error=" << g.encounter_radius_error
              << " target_radius_error=" << g.target_radius_error
              << " abs_encounter_radius_error=" << std::abs(g.encounter_radius_error)
              << " abs_target_radius_error=" << std::abs(g.target_radius_error)
              << " valid_geometry=" << (g.valid_geometry ? 1 : 0)
              << " invalid_reason=" << g.invalid_reason
              << '\n';
}

bool near_zero_R_prime(double R_prime) {
    return is_finite(R_prime) && std::abs(R_prime) <= 1.0e6;
}

bool near_zero_B_target(double B_target) {
    return is_finite(B_target) && std::abs(B_target) <= 1.0e-8;
}

bool near_geometry_tolerance(const GeometryDiagnostic& g) {
    const double encounter_error = std::abs(g.encounter_radius_error);
    const double target_error = std::abs(g.target_radius_error);
    return (is_finite(encounter_error) && encounter_error >= 0.1 * kBackSubstitutionTolerance) ||
           (is_finite(target_error) && target_error >= 0.1 * kBackSubstitutionTolerance);
}

void update_boundary_counts(const BranchWithGeometry& item, Summary* summary) {
    if (near_zero_R_prime(item.geometry.R_prime)) {
        summary->near_zero_R_prime_count += 1;
    }
    if (near_zero_B_target(item.geometry.B_target)) {
        summary->near_zero_B_target_count += 1;
    }
    if (near_geometry_tolerance(item.geometry)) {
        summary->near_geometry_tolerance_count += 1;
    }
}

std::string dominant_reason(const std::map<std::string, int>& reasons) {
    std::string best = "none";
    int best_count = 0;
    for (const auto& [reason, count] : reasons) {
        if (count > best_count) {
            best = reason;
            best_count = count;
        }
    }
    return best;
}

std::string recommendation(const Summary& summary) {
    if (summary.validity_flip_count == 0) {
        return "no_validity_flip_detected_enter_continuation_debug";
    }
    if (summary.direct_table_alpha_max_diff > 1e-6 || summary.direct_table_time_max_diff > 10.0) {
        return "investigate_nearest_node_refinement_accuracy";
    }
    if (summary.near_zero_R_prime_count > 0 ||
        summary.near_zero_B_target_count > 0 ||
        summary.near_geometry_tolerance_count > 0) {
        return "treat_slingshot_geometry_boundary_as_ambiguous_and_skip_for_sign_change";
    }
    return "fix_branch_pairing_metric";
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem2_theta_prime_slingshot_validity_boundary_skipped_missing_table\n";
        return 0;
    }

    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    problem1::Problem1NearestNodeQueryOptions options{};
    options.residual_tolerance_seconds = 1e-2;
    options.max_newton_iterations = 80;
    options.fallback_direct_solve = true;

    const auto departure_planet = planet_params::PlanetId::Earth;
    const auto target_planet = planet_params::PlanetId::Mars;
    const double encounter_time = 0.17 * planet_params::planet_orbital_period(departure_planet);
    const auto departure_state = planet_params::planet_state_at_time(departure_planet, encounter_time);
    const auto target_state = planet_params::planet_state_at_time(target_planet, encounter_time);
    const double phi = departure_state.varphi;
    const double beta = target_state.varphi;
    const double incoming_e = 0.3;
    const double incoming_theta = 0.4;

    Summary summary{};
    summary.analyzed_interval_count = 2;
    for (const int interval_left : std::array<int, 2>{2, 18}) {
        for (const int endpoint_index : std::array<int, 2>{interval_left, interval_left + 1}) {
            const double theta_prime =
                kTwoPi * static_cast<double>(endpoint_index) / static_cast<double>(kThetaSampleCount);
            const Sample direct = build_direct_sample(
                departure_planet, target_planet, encounter_time, phi, beta, incoming_e, incoming_theta, theta_prime);
            const Sample table = build_table_sample(
                loader, options, departure_planet, target_planet, encounter_time, phi, beta, incoming_e,
                incoming_theta, theta_prime);
            const std::string side = side_name(endpoint_index, interval_left);
            const auto pairs = pair_direct_table_branches(direct, table);
            for (const auto& pair : pairs) {
                summary.analyzed_branch_pair_count += 1;
                const auto& d = pair.direct;
                const auto& t = pair.table;
                summary.direct_table_alpha_max_diff = std::max(
                    summary.direct_table_alpha_max_diff,
                    std::abs(normalize_angle_minus_pi_pi(d.branch.alpha - t.branch.alpha)));
                summary.direct_table_time_max_diff = std::max(
                    summary.direct_table_time_max_diff,
                    std::abs(d.branch.time_of_flight_seconds - t.branch.time_of_flight_seconds));
                if (is_finite(d.geometry.denominator) && is_finite(t.geometry.denominator)) {
                    summary.direct_table_denominator_max_diff = std::max(
                        summary.direct_table_denominator_max_diff,
                        std::abs(d.geometry.denominator - t.geometry.denominator));
                }
                if (is_finite(d.geometry.R_prime) && is_finite(t.geometry.R_prime)) {
                    summary.direct_table_R_prime_max_diff = std::max(
                        summary.direct_table_R_prime_max_diff,
                        std::abs(d.geometry.R_prime - t.geometry.R_prime));
                }
                update_boundary_counts(d, &summary);
                update_boundary_counts(t, &summary);
                if (d.branch.slingshot_valid != t.branch.slingshot_valid) {
                    summary.validity_flip_count += 1;
                    const std::string reason = d.branch.slingshot_valid
                        ? t.branch.slingshot_invalid_reason
                        : d.branch.slingshot_invalid_reason;
                    summary.flip_reason_count[reason.empty() ? "unknown" : reason] += 1;
                    std::cout << "Problem2SlingshotValidityBoundaryPairDiff"
                              << " interval=" << interval_left << "-" << (interval_left + 1)
                              << " side=" << side
                              << " rank=" << pair.rank
                              << " k=" << d.branch.transfer_revolution
                              << " q=" << d.branch.target_revolution
                              << " direct_slingshot_valid=" << (d.branch.slingshot_valid ? 1 : 0)
                              << " table_slingshot_valid=" << (t.branch.slingshot_valid ? 1 : 0)
                              << " flip_reason=" << (reason.empty() ? "unknown" : reason)
                              << " alpha_diff=" << std::abs(normalize_angle_minus_pi_pi(d.branch.alpha - t.branch.alpha))
                              << " time_diff_seconds="
                              << std::abs(d.branch.time_of_flight_seconds - t.branch.time_of_flight_seconds)
                              << " denominator_diff=" << std::abs(d.geometry.denominator - t.geometry.denominator)
                              << " R_prime_diff=" << std::abs(d.geometry.R_prime - t.geometry.R_prime)
                              << '\n';
                }

                print_branch_case(interval_left, side.c_str(), "direct", pair.rank, d);
                print_branch_case(interval_left, side.c_str(), "table", pair.rank, t);
                print_geometry_case(interval_left, side.c_str(), "direct", pair.rank, d);
                print_geometry_case(interval_left, side.c_str(), "table", pair.rank, t);
            }
        }
    }

    std::cout << "Problem2SlingshotValidityBoundarySummary\n";
    std::cout << "analyzed_interval_count=" << summary.analyzed_interval_count << '\n';
    std::cout << "analyzed_branch_pair_count=" << summary.analyzed_branch_pair_count << '\n';
    std::cout << "direct_table_alpha_max_diff=" << summary.direct_table_alpha_max_diff << '\n';
    std::cout << "direct_table_time_max_diff=" << summary.direct_table_time_max_diff << '\n';
    std::cout << "direct_table_denominator_max_diff=" << summary.direct_table_denominator_max_diff << '\n';
    std::cout << "direct_table_R_prime_max_diff=" << summary.direct_table_R_prime_max_diff << '\n';
    std::cout << "near_zero_R_prime_count=" << summary.near_zero_R_prime_count << '\n';
    std::cout << "near_zero_B_target_count=" << summary.near_zero_B_target_count << '\n';
    std::cout << "near_geometry_tolerance_count=" << summary.near_geometry_tolerance_count << '\n';
    std::cout << "validity_flip_count=" << summary.validity_flip_count << '\n';
    std::cout << "dominant_flip_reason=" << dominant_reason(summary.flip_reason_count) << '\n';
    std::cout << "recommendation=" << recommendation(summary) << '\n';
    std::cout << "diagnostic_ok=1\n";
    return 0;
}
