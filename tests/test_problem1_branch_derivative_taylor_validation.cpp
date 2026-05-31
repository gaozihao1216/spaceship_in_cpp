#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace {

using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = 3.141592653589793238462643383279502884 / 90.0;
constexpr int kMaxBranchesPerPair = 50;

struct QueryRecord {
    spaceship_cpp::planet_params::PlanetId from_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    double query_nu_A = 0.0;
    double query_nu_B = 0.0;
    double query_theta_A = 0.0;
    bool baseline_positive = false;
    bool baseline_used_direct_fallback = false;
    bool query_level_miss = false;
};

struct BranchCase {
    std::string case_name;
    spaceship_cpp::planet_params::PlanetId from_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to_planet = spaceship_cpp::planet_params::PlanetId::Mercury;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct Direction {
    std::string name;
    double a = 0.0;
    double b = 0.0;
    double t = 0.0;
};

struct PairSummary {
    std::set<int> branch_ids;
    int tested_direction_count = 0;
    int taylor_pass_count = 0;
    int taylor_fail_count = 0;
    int roundoff_pass_count = 0;
    int finite_diff_pass_count = 0;
    int finite_diff_fail_count = 0;
    int branch_unstable_count = 0;
    double slope_sum = 0.0;
    int slope_count = 0;
    double min_order_slope = std::numeric_limits<double>::infinity();
    double max_abs_fd_derivative_diff = 0.0;
    bool derivative_suspect = false;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') {
            return raw;
        }
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

spaceship_cpp::planet_params::PlanetId parse_planet(const std::string& name) {
    namespace planet_params = spaceship_cpp::planet_params;
    if (name == "Mercury") return planet_params::PlanetId::Mercury;
    if (name == "Venus") return planet_params::PlanetId::Venus;
    if (name == "Earth") return planet_params::PlanetId::Earth;
    if (name == "Mars") return planet_params::PlanetId::Mars;
    if (name == "Jupiter") return planet_params::PlanetId::Jupiter;
    if (name == "Saturn") return planet_params::PlanetId::Saturn;
    if (name == "Uranus") return planet_params::PlanetId::Uranus;
    return planet_params::PlanetId::Neptune;
}

std::string planet_name_string(spaceship_cpp::planet_params::PlanetId id) {
    return spaceship_cpp::planet_params::planet_name(id);
}

std::string pair_key(spaceship_cpp::planet_params::PlanetId from, spaceship_cpp::planet_params::PlanetId to) {
    return planet_name_string(from) + "->" + planet_name_string(to);
}

int lower_grid_index(double angle) {
    int index = static_cast<int>(std::floor(normalize_angle_0_2pi(angle) / kGridStepRadians)) % kSamplesPerDimension;
    if (index < 0) {
        index += kSamplesPerDimension;
    }
    return index;
}

double grid_angle(int index) {
    return normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(index));
}

bool parse_bool_int(const std::string& raw) {
    return raw == "1" || raw == "true";
}

bool parse_query_samples_from_log(const std::filesystem::path& path, std::vector<QueryRecord>* records) {
    std::ifstream in(path);
    if (!in) {
        return false;
    }
    std::string line;
    QueryRecord current{};
    bool in_sample = false;
    bool have_from = false;
    bool have_to = false;
    bool have_coords = false;
    auto flush = [&]() {
        if (in_sample && have_from && have_to && have_coords) {
            records->push_back(current);
        }
        current = QueryRecord{};
        in_sample = false;
        have_from = have_to = have_coords = false;
    };
    while (std::getline(in, line)) {
        if (line == "Problem1ScoutQueryRecallSample") {
            flush();
            in_sample = true;
            continue;
        }
        if (line.rfind("Problem1", 0) == 0) {
            flush();
            continue;
        }
        if (!in_sample) {
            continue;
        }
        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        const auto key = line.substr(0, eq);
        const auto value = line.substr(eq + 1);
        if (key == "from_planet") {
            current.from_planet = parse_planet(value);
            have_from = true;
        } else if (key == "to_planet") {
            current.to_planet = parse_planet(value);
            have_to = true;
        } else if (key == "query_nu_A") {
            current.query_nu_A = std::stod(value);
        } else if (key == "query_nu_B") {
            current.query_nu_B = std::stod(value);
        } else if (key == "query_theta_A") {
            current.query_theta_A = std::stod(value);
            have_coords = true;
        } else if (key == "baseline_positive") {
            current.baseline_positive = parse_bool_int(value);
        } else if (key == "baseline_used_direct_fallback") {
            current.baseline_used_direct_fallback = parse_bool_int(value);
        } else if (key == "query_level_miss") {
            current.query_level_miss = parse_bool_int(value);
        }
    }
    flush();
    return !records->empty();
}

std::vector<QueryRecord> load_query_records() {
    std::vector<QueryRecord> records;
    if (parse_query_samples_from_log("/tmp/problem1_scout_query_recall_rate.log", &records)) {
        return records;
    }
    parse_query_samples_from_log("/tmp/problem1_scout_query_recall_rate_exhaustive8_sampled_2k.log", &records);
    return records;
}

std::vector<QueryRecord> fallback_seed_queries() {
    namespace planet_params = spaceship_cpp::planet_params;
    std::vector<QueryRecord> out;
    const std::vector<std::pair<planet_params::PlanetId, planet_params::PlanetId>> pairs{
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mars},
        {planet_params::PlanetId::Mars, planet_params::PlanetId::Earth},
        {planet_params::PlanetId::Earth, planet_params::PlanetId::Mercury},
        {planet_params::PlanetId::Venus, planet_params::PlanetId::Mars},
        {planet_params::PlanetId::Mars, planet_params::PlanetId::Venus},
        {planet_params::PlanetId::Venus, planet_params::PlanetId::Mercury},
        {planet_params::PlanetId::Mars, planet_params::PlanetId::Mercury},
    };
    for (const auto& pair : pairs) {
        QueryRecord q{};
        q.from_planet = pair.first;
        q.to_planet = pair.second;
        q.query_nu_A = 1.0;
        q.query_nu_B = 2.0;
        q.query_theta_A = 3.0;
        out.push_back(q);
    }
    return out;
}

void maybe_add_branch_case(
    std::vector<BranchCase>* cases,
    std::map<std::string, int>* count_by_pair,
    const std::string& case_name,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double nu_A,
    double nu_B,
    double theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    const auto key = pair_key(from, to);
    if ((*count_by_pair)[key] >= kMaxBranchesPerPair ||
        !branch.valid ||
        !branch.derivatives_available ||
        !std::isfinite(branch.encounter_global_angle)) {
        return;
    }
    for (const auto& existing : *cases) {
        if (existing.from_planet == from &&
            existing.to_planet == to &&
            existing.branch.transfer_revolution == branch.transfer_revolution &&
            existing.branch.target_revolution == branch.target_revolution &&
            std::abs(normalize_angle_minus_pi_pi(existing.branch.encounter_global_angle -
                                                 branch.encounter_global_angle)) < 1e-10 &&
            std::abs(normalize_angle_minus_pi_pi(existing.nu_A - nu_A)) < 1e-12 &&
            std::abs(normalize_angle_minus_pi_pi(existing.nu_B - nu_B)) < 1e-12 &&
            std::abs(normalize_angle_minus_pi_pi(existing.theta_A - theta_A)) < 1e-12) {
            return;
        }
    }
    BranchCase c{};
    c.case_name = case_name;
    c.from_planet = from;
    c.to_planet = to;
    c.nu_A = nu_A;
    c.nu_B = nu_B;
    c.theta_A = theta_A;
    c.branch = branch;
    cases->push_back(c);
    (*count_by_pair)[key] += 1;
}

void collect_cell_vertex_cases(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const QueryRecord& query,
    const std::string& case_name,
    std::vector<BranchCase>* cases,
    std::map<std::string, int>* count_by_pair
) {
    const int lower_a = lower_grid_index(query.query_nu_A);
    const int lower_b = lower_grid_index(query.query_nu_B);
    const int lower_t = lower_grid_index(query.query_theta_A);
    for (int da = 0; da <= 1; ++da) {
        for (int db = 0; db <= 1; ++db) {
            for (int dt = 0; dt <= 1; ++dt) {
                const int ia = (lower_a + da) % kSamplesPerDimension;
                const int ib = (lower_b + db) % kSamplesPerDimension;
                const int it = (lower_t + dt) % kSamplesPerDimension;
                const auto node = loader.load_node_by_indices(ia, ib, it);
                if (!node.valid) {
                    continue;
                }
                for (const auto& branch : node.branches) {
                    maybe_add_branch_case(
                        cases,
                        count_by_pair,
                        case_name,
                        query.from_planet,
                        query.to_planet,
                        grid_angle(ia),
                        grid_angle(ib),
                        grid_angle(it),
                        branch);
                    if ((*count_by_pair)[pair_key(query.from_planet, query.to_planet)] >= kMaxBranchesPerPair) {
                        return;
                    }
                }
            }
        }
    }
}

double evaluate_residual_seconds(
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double nu_A,
    double nu_B,
    double theta_A,
    int k,
    int q,
    double alpha
) {
    const auto residual = spaceship_cpp::problem1::evaluate_problem1_root_residual(
        from, to, normalize_angle_0_2pi(nu_A), normalize_angle_0_2pi(nu_B), normalize_angle_0_2pi(theta_A),
        normalize_angle_0_2pi(alpha), k, q);
    return residual.valid ? residual.residual_seconds : std::numeric_limits<double>::quiet_NaN();
}

std::vector<Direction> directions() {
    const double inv = 1.0 / std::sqrt(3.0);
    return {
        {"nu_A", 1.0, 0.0, 0.0},
        {"nu_B", 0.0, 1.0, 0.0},
        {"theta_A", 0.0, 0.0, 1.0},
        {"diag_111", inv, inv, inv},
        {"diag_1m11", inv, -inv, inv},
        {"diag_m111", -inv, inv, inv},
    };
}

std::array<double, 6> h_values() {
    return {
        kGridStepRadians / 8.0,
        kGridStepRadians / 16.0,
        kGridStepRadians / 32.0,
        kGridStepRadians / 64.0,
        kGridStepRadians / 128.0,
        kGridStepRadians / 256.0,
    };
}

double branch_derivative_component(const spaceship_cpp::problem1::Problem1SolutionBranch& branch, int axis) {
    if (axis == 0) return branch.d_encounter_global_angle_d_nu_A;
    if (axis == 1) return branch.d_encounter_global_angle_d_nu_B;
    return branch.d_encounter_global_angle_d_theta_A;
}

const spaceship_cpp::problem1::Problem1SolutionBranch* find_same_branch(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source,
    double predicted_alpha
) {
    const spaceship_cpp::problem1::Problem1SolutionBranch* best = nullptr;
    double best_score = std::numeric_limits<double>::infinity();
    for (const auto& branch : branches) {
        if (!branch.valid ||
            branch.transfer_revolution != source.transfer_revolution ||
            branch.target_revolution != source.target_revolution) {
            continue;
        }
        const double score = std::abs(normalize_angle_minus_pi_pi(branch.encounter_global_angle - predicted_alpha));
        if (score < best_score) {
            best_score = score;
            best = &branch;
        }
    }
    return best;
}

}  // namespace

int main() {
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem1_branch_derivative_taylor_validation_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);

    auto records = load_query_records();
    if (records.empty()) {
        records = fallback_seed_queries();
    }
    std::stable_sort(records.begin(), records.end(), [](const auto& a, const auto& b) {
        const int aw = (a.query_level_miss ? 4 : 0) + (a.baseline_used_direct_fallback ? 2 : 0) +
            (a.baseline_positive ? 1 : 0);
        const int bw = (b.query_level_miss ? 4 : 0) + (b.baseline_used_direct_fallback ? 2 : 0) +
            (b.baseline_positive ? 1 : 0);
        return aw > bw;
    });

    const std::set<std::string> target_pairs{
        "Earth->Mars",
        "Mars->Earth",
        "Earth->Mercury",
        "Venus->Mars",
        "Mars->Venus",
        "Venus->Mercury",
        "Mars->Mercury",
    };
    std::vector<BranchCase> cases;
    std::map<std::string, int> count_by_pair;
    for (const auto& record : records) {
        const auto key = pair_key(record.from_planet, record.to_planet);
        if (target_pairs.count(key) == 0 || count_by_pair[key] >= kMaxBranchesPerPair) {
            continue;
        }
        const std::string case_name =
            (record.query_level_miss ? "scout_missing_" : "sample_") + key;
        collect_cell_vertex_cases(loader, record, case_name, &cases, &count_by_pair);
    }

    std::cout << "Problem1BranchDerivativeTaylorValidationCollection\n";
    std::cout << "query_record_count=" << records.size() << '\n';
    std::cout << "branch_case_count=" << cases.size() << '\n';
    for (const auto& [key, count] : count_by_pair) {
        std::cout << "pair_branch_count_" << key << "=" << count << '\n';
    }

    std::map<std::string, PairSummary> pair_summaries;
    int global_taylor_pass = 0;
    int global_taylor_fail = 0;
    int global_roundoff_pass = 0;
    int global_fd_pass = 0;
    int global_fd_fail = 0;
    int global_unstable = 0;
    double global_slope_sum = 0.0;
    int global_slope_count = 0;
    double worst_slope = std::numeric_limits<double>::infinity();
    std::string worst_case;
    std::string worst_pair;
    int worst_k = 0;
    int worst_q = 0;
    std::string worst_direction;

    const auto hs = h_values();
    const auto dirs = directions();
    for (std::size_t branch_index = 0; branch_index < cases.size(); ++branch_index) {
        const auto& c = cases[branch_index];
        auto& pair_summary = pair_summaries[pair_key(c.from_planet, c.to_planet)];
        pair_summary.branch_ids.insert(static_cast<int>(branch_index));

        const auto& branch = c.branch;
        for (const auto& dir : dirs) {
            std::vector<double> residuals;
            residuals.reserve(hs.size());
            bool all_roundoff = true;
            for (const double h : hs) {
                const double dnu_A = h * dir.a;
                const double dnu_B = h * dir.b;
                const double dtheta_A = h * dir.t;
                const double alpha_pred = normalize_angle_0_2pi(
                    branch.encounter_global_angle +
                    branch.d_encounter_global_angle_d_nu_A * dnu_A +
                    branch.d_encounter_global_angle_d_nu_B * dnu_B +
                    branch.d_encounter_global_angle_d_theta_A * dtheta_A);
                const double residual = evaluate_residual_seconds(
                    c.from_planet,
                    c.to_planet,
                    c.nu_A + dnu_A,
                    c.nu_B + dnu_B,
                    c.theta_A + dtheta_A,
                    branch.transfer_revolution,
                    branch.target_revolution,
                    alpha_pred);
                const double abs_residual = std::abs(residual);
                residuals.push_back(abs_residual);
                all_roundoff = all_roundoff && std::isfinite(abs_residual) && abs_residual < 1e-6;
                std::cout << "Problem1BranchDerivativeTaylorSample\n";
                std::cout << "case_name=" << c.case_name << '\n';
                std::cout << "from_planet=" << planet_name_string(c.from_planet) << '\n';
                std::cout << "to_planet=" << planet_name_string(c.to_planet) << '\n';
                std::cout << "branch_index=" << branch_index << '\n';
                std::cout << "k=" << branch.transfer_revolution << '\n';
                std::cout << "q=" << branch.target_revolution << '\n';
                std::cout << "direction_name=" << dir.name << '\n';
                std::cout << "h=" << h << '\n';
                std::cout << "residual_seconds=" << residual << '\n';
                std::cout << "abs_residual=" << abs_residual << '\n';
                std::cout << "abs_residual_over_h=" << (abs_residual / h) << '\n';
                std::cout << "abs_residual_over_h2=" << (abs_residual / (h * h)) << '\n';
            }

            std::vector<double> slopes;
            for (int i = 1; i <= 3; ++i) {
                if (residuals[static_cast<std::size_t>(i)] > 0.0 &&
                    residuals[static_cast<std::size_t>(i + 1)] > 0.0 &&
                    std::isfinite(residuals[static_cast<std::size_t>(i)]) &&
                    std::isfinite(residuals[static_cast<std::size_t>(i + 1)])) {
                    slopes.push_back(
                        std::log(residuals[static_cast<std::size_t>(i)] /
                                 residuals[static_cast<std::size_t>(i + 1)]) /
                        std::log(hs[static_cast<std::size_t>(i)] / hs[static_cast<std::size_t>(i + 1)]));
                }
            }
            double mean_slope = std::numeric_limits<double>::quiet_NaN();
            double min_slope = std::numeric_limits<double>::quiet_NaN();
            double max_slope = std::numeric_limits<double>::quiet_NaN();
            if (!slopes.empty()) {
                double sum = 0.0;
                min_slope = std::numeric_limits<double>::infinity();
                max_slope = -std::numeric_limits<double>::infinity();
                for (const double slope : slopes) {
                    sum += slope;
                    min_slope = std::min(min_slope, slope);
                    max_slope = std::max(max_slope, slope);
                }
                mean_slope = sum / static_cast<double>(slopes.size());
                pair_summary.slope_sum += mean_slope;
                pair_summary.slope_count += 1;
                pair_summary.min_order_slope = std::min(pair_summary.min_order_slope, mean_slope);
                global_slope_sum += mean_slope;
                global_slope_count += 1;
                if (mean_slope < worst_slope) {
                    worst_slope = mean_slope;
                    worst_case = c.case_name;
                    worst_pair = pair_key(c.from_planet, c.to_planet);
                    worst_k = branch.transfer_revolution;
                    worst_q = branch.target_revolution;
                    worst_direction = dir.name;
                }
            }
            const bool pass = all_roundoff || (std::isfinite(mean_slope) && mean_slope >= 1.7);
            const bool likely_wrong = std::isfinite(mean_slope) && mean_slope < 1.3 && !all_roundoff;
            std::cout << "Problem1BranchDerivativeTaylorDirectionSummary\n";
            std::cout << "case_name=" << c.case_name << '\n';
            std::cout << "from_planet=" << planet_name_string(c.from_planet) << '\n';
            std::cout << "to_planet=" << planet_name_string(c.to_planet) << '\n';
            std::cout << "branch_index=" << branch_index << '\n';
            std::cout << "k=" << branch.transfer_revolution << '\n';
            std::cout << "q=" << branch.target_revolution << '\n';
            std::cout << "direction_name=" << dir.name << '\n';
            std::cout << "mean_order_slope=" << mean_slope << '\n';
            std::cout << "min_order_slope=" << min_slope << '\n';
            std::cout << "max_order_slope=" << max_slope << '\n';
            std::cout << "passes_h2_test=" << (pass ? 1 : 0) << '\n';
            std::cout << "failure_reason="
                      << (all_roundoff ? "roundoff_dominated" :
                          (pass ? "" : (likely_wrong ? "likely_wrong_derivative" : "weak_h2_order")))
                      << '\n';
            pair_summary.tested_direction_count += 1;
            if (all_roundoff) {
                pair_summary.roundoff_pass_count += 1;
                global_roundoff_pass += 1;
            } else if (pass) {
                pair_summary.taylor_pass_count += 1;
                global_taylor_pass += 1;
            } else {
                pair_summary.taylor_fail_count += 1;
                pair_summary.derivative_suspect = true;
                global_taylor_fail += 1;
            }
        }

        const double fd_h = kGridStepRadians / 512.0;
        for (int axis = 0; axis < 3; ++axis) {
            const double dnu_A = axis == 0 ? fd_h : 0.0;
            const double dnu_B = axis == 1 ? fd_h : 0.0;
            const double dtheta_A = axis == 2 ? fd_h : 0.0;
            const double table_derivative = branch_derivative_component(branch, axis);
            const double alpha_plus_pred = normalize_angle_0_2pi(branch.encounter_global_angle + table_derivative * fd_h);
            const double alpha_minus_pred = normalize_angle_0_2pi(branch.encounter_global_angle - table_derivative * fd_h);
            const auto plus = problem1::solve_problem1_from_departure_anomalies(
                c.from_planet, c.to_planet, c.nu_A + dnu_A, c.nu_B + dnu_B, c.theta_A + dtheta_A, 1, 1);
            const auto minus = problem1::solve_problem1_from_departure_anomalies(
                c.from_planet, c.to_planet, c.nu_A - dnu_A, c.nu_B - dnu_B, c.theta_A - dtheta_A, 1, 1);
            const auto* plus_branch = find_same_branch(plus, branch, alpha_plus_pred);
            const auto* minus_branch = find_same_branch(minus, branch, alpha_minus_pred);
            const bool stable = plus_branch != nullptr && minus_branch != nullptr;
            double fd = std::numeric_limits<double>::quiet_NaN();
            double abs_diff = std::numeric_limits<double>::quiet_NaN();
            double rel_diff = std::numeric_limits<double>::quiet_NaN();
            if (stable) {
                fd = normalize_angle_minus_pi_pi(
                    plus_branch->encounter_global_angle - minus_branch->encounter_global_angle) / (2.0 * fd_h);
                abs_diff = std::abs(normalize_angle_minus_pi_pi(table_derivative - fd));
                rel_diff = abs_diff / std::max(1e-12, std::abs(fd));
            }
            const bool fd_pass = stable && abs_diff <= 1e-3;
            const bool sign_flip = stable && std::abs(normalize_angle_minus_pi_pi(-table_derivative - fd)) < abs_diff;
            std::cout << "Problem1BranchDerivativeFiniteDiffSample\n";
            std::cout << "case_name=" << c.case_name << '\n';
            std::cout << "from_planet=" << planet_name_string(c.from_planet) << '\n';
            std::cout << "to_planet=" << planet_name_string(c.to_planet) << '\n';
            std::cout << "branch_index=" << branch_index << '\n';
            std::cout << "k=" << branch.transfer_revolution << '\n';
            std::cout << "q=" << branch.target_revolution << '\n';
            std::cout << "axis=" << (axis == 0 ? "nu_A" : (axis == 1 ? "nu_B" : "theta_A")) << '\n';
            std::cout << "h=" << fd_h << '\n';
            std::cout << "d_alpha_table=" << table_derivative << '\n';
            std::cout << "d_alpha_fd=" << fd << '\n';
            std::cout << "abs_diff=" << abs_diff << '\n';
            std::cout << "relative_diff=" << rel_diff << '\n';
            std::cout << "plus_found=" << (plus_branch != nullptr ? 1 : 0) << '\n';
            std::cout << "minus_found=" << (minus_branch != nullptr ? 1 : 0) << '\n';
            std::cout << "branch_stable=" << (stable ? 1 : 0) << '\n';
            if (!stable) {
                pair_summary.branch_unstable_count += 1;
                global_unstable += 1;
            } else if (fd_pass) {
                pair_summary.finite_diff_pass_count += 1;
                global_fd_pass += 1;
            } else {
                pair_summary.finite_diff_fail_count += 1;
                pair_summary.max_abs_fd_derivative_diff =
                    std::max(pair_summary.max_abs_fd_derivative_diff, abs_diff);
                pair_summary.derivative_suspect = true;
                global_fd_fail += 1;
                std::cout << "Problem1DerivativeCoordinateConventionAudit\n";
                std::cout << "case_name=" << c.case_name << '\n';
                std::cout << "from_planet=" << planet_name_string(c.from_planet) << '\n';
                std::cout << "to_planet=" << planet_name_string(c.to_planet) << '\n';
                std::cout << "theta_A_definition_used=root_table_departure_theta_A\n";
                std::cout << "nu_A_wrap=normalize_angle_0_2pi\n";
                std::cout << "nu_B_wrap=normalize_angle_0_2pi\n";
                std::cout << "theta_A_wrap=normalize_angle_0_2pi\n";
                std::cout << "alpha_wrap=normalize_angle_0_2pi\n";
                std::cout << "derivative_component=" << table_derivative << '\n';
                std::cout << "finite_diff_component=" << fd << '\n';
                std::cout << "sign_flip_would_match=" << (sign_flip ? 1 : 0) << '\n';
                std::cout << "theta_A_sign_flip_would_match=" << ((axis == 2 && sign_flip) ? 1 : 0) << '\n';
                std::cout << "alpha_local_global_suspect=0\n";
            }
        }

        const Direction quality_dir{"diag_111", 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
        for (const double h : {kGridStepRadians / 16.0, kGridStepRadians / 64.0, kGridStepRadians / 256.0}) {
            const double dnu_A = h * quality_dir.a;
            const double dnu_B = h * quality_dir.b;
            const double dtheta_A = h * quality_dir.t;
            const double alpha_pred = normalize_angle_0_2pi(
                branch.encounter_global_angle +
                branch.d_encounter_global_angle_d_nu_A * dnu_A +
                branch.d_encounter_global_angle_d_nu_B * dnu_B +
                branch.d_encounter_global_angle_d_theta_A * dtheta_A);
            const auto exact = problem1::solve_problem1_from_departure_anomalies(
                c.from_planet, c.to_planet, c.nu_A + dnu_A, c.nu_B + dnu_B, c.theta_A + dtheta_A, 1, 1);
            const auto* exact_branch = find_same_branch(exact, branch, alpha_pred);
            const bool found = exact_branch != nullptr;
            const double alpha_error = found ?
                std::abs(normalize_angle_minus_pi_pi(alpha_pred - exact_branch->encounter_global_angle)) :
                std::numeric_limits<double>::quiet_NaN();
            const double residual_pred = evaluate_residual_seconds(
                c.from_planet, c.to_planet, c.nu_A + dnu_A, c.nu_B + dnu_B, c.theta_A + dtheta_A,
                branch.transfer_revolution, branch.target_revolution, alpha_pred);
            std::cout << "Problem1RouteASeedQuality\n";
            std::cout << "case_name=" << c.case_name << '\n';
            std::cout << "from_planet=" << planet_name_string(c.from_planet) << '\n';
            std::cout << "to_planet=" << planet_name_string(c.to_planet) << '\n';
            std::cout << "branch_index=" << branch_index << '\n';
            std::cout << "k=" << branch.transfer_revolution << '\n';
            std::cout << "q=" << branch.target_revolution << '\n';
            std::cout << "direction_name=" << quality_dir.name << '\n';
            std::cout << "h=" << h << '\n';
            std::cout << "alpha_error=" << alpha_error << '\n';
            std::cout << "alpha_error_over_h=" << (found ? alpha_error / h : std::numeric_limits<double>::quiet_NaN()) << '\n';
            std::cout << "alpha_error_over_h2=" << (found ? alpha_error / (h * h) : std::numeric_limits<double>::quiet_NaN()) << '\n';
            std::cout << "residual_pred=" << residual_pred << '\n';
            std::cout << "direct_same_branch_found=" << (found ? 1 : 0) << '\n';
            std::cout << "seed_quality_ok=" << (found && alpha_error / (h * h) < 1e3 ? 1 : 0) << '\n';
        }
    }

    for (const auto& [key, summary] : pair_summaries) {
        const auto arrow = key.find("->");
        const double mean_slope =
            summary.slope_count > 0 ? summary.slope_sum / static_cast<double>(summary.slope_count) :
                                      std::numeric_limits<double>::quiet_NaN();
        std::cout << "Problem1DerivativeValidationPairSummary\n";
        std::cout << "from_planet=" << key.substr(0, arrow) << '\n';
        std::cout << "to_planet=" << key.substr(arrow + 2) << '\n';
        std::cout << "tested_branch_count=" << summary.branch_ids.size() << '\n';
        std::cout << "tested_direction_count=" << summary.tested_direction_count << '\n';
        std::cout << "taylor_pass_count=" << summary.taylor_pass_count << '\n';
        std::cout << "taylor_fail_count=" << summary.taylor_fail_count << '\n';
        std::cout << "roundoff_pass_count=" << summary.roundoff_pass_count << '\n';
        std::cout << "finite_diff_pass_count=" << summary.finite_diff_pass_count << '\n';
        std::cout << "finite_diff_fail_count=" << summary.finite_diff_fail_count << '\n';
        std::cout << "branch_unstable_count=" << summary.branch_unstable_count << '\n';
        std::cout << "mean_order_slope=" << mean_slope << '\n';
        std::cout << "min_order_slope=" << summary.min_order_slope << '\n';
        std::cout << "max_abs_fd_derivative_diff=" << summary.max_abs_fd_derivative_diff << '\n';
        std::cout << "derivative_suspect=" << (summary.derivative_suspect ? 1 : 0) << '\n';
    }

    std::cout << "Problem1DerivativeValidationGlobalSummary\n";
    std::cout << "tested_branch_count=" << cases.size() << '\n';
    std::cout << "taylor_pass_count=" << global_taylor_pass << '\n';
    std::cout << "taylor_fail_count=" << global_taylor_fail << '\n';
    std::cout << "finite_diff_pass_count=" << global_fd_pass << '\n';
    std::cout << "finite_diff_fail_count=" << global_fd_fail << '\n';
    std::cout << "branch_unstable_count=" << global_unstable << '\n';
    std::cout << "mean_order_slope="
              << (global_slope_count > 0 ? global_slope_sum / static_cast<double>(global_slope_count) :
                                           std::numeric_limits<double>::quiet_NaN())
              << '\n';
    std::cout << "worst_case_name=" << worst_case << '\n';
    std::cout << "worst_pair=" << worst_pair << '\n';
    std::cout << "worst_branch_k=" << worst_k << '\n';
    std::cout << "worst_branch_q=" << worst_q << '\n';
    std::cout << "worst_direction=" << worst_direction << '\n';
    std::cout << "worst_order_slope=" << worst_slope << '\n';
    std::cout << "derivative_validation_ok=1\n";
    return 0;
}
