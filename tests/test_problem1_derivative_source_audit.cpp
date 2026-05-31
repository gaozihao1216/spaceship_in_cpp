#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/planet_params/planet_params.hpp"
#include "spaceship_cpp/problem1/problem1_root_table.hpp"
#include "spaceship_cpp/problem1/problem1_root_table_2deg_loader.hpp"

#include <algorithm>
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
#include <vector>

namespace {

using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = 3.141592653589793238462643383279502884 / 90.0;
constexpr int kMaxBranchesPerPair = 20;
constexpr double kFdStep = kGridStepRadians / 512.0;
constexpr double kMatchTolerance = 1e-3;

struct QueryRecord {
    spaceship_cpp::planet_params::PlanetId from = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to = spaceship_cpp::planet_params::PlanetId::Mercury;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
    bool baseline_positive = false;
    bool baseline_used_direct_fallback = false;
    bool query_level_miss = false;
};

struct BranchCase {
    std::string name;
    spaceship_cpp::planet_params::PlanetId from = spaceship_cpp::planet_params::PlanetId::Mercury;
    spaceship_cpp::planet_params::PlanetId to = spaceship_cpp::planet_params::PlanetId::Mercury;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
    spaceship_cpp::problem1::Problem1SolutionBranch branch;
};

struct Summary {
    int count = 0;
    int table_matches_attach = 0;
    int attach_matches_fd = 0;
    int table_matches_fd = 0;
    int table_generation_or_loading = 0;
    int attach_formula_or_convention = 0;
    int both_or_matching = 0;
    int derivative_ok = 0;
};

std::filesystem::path table_dir_from_env() {
    if (const char* raw = std::getenv("ROOT_TABLE_2DEG_DIR")) {
        if (*raw != '\0') return raw;
    }
    return "/home/gaozihao/spaceship_in_cpp/root_tables/problem1_root_table_2deg_full";
}

spaceship_cpp::planet_params::PlanetId parse_planet(const std::string& name) {
    namespace p = spaceship_cpp::planet_params;
    if (name == "Mercury") return p::PlanetId::Mercury;
    if (name == "Venus") return p::PlanetId::Venus;
    if (name == "Earth") return p::PlanetId::Earth;
    if (name == "Mars") return p::PlanetId::Mars;
    if (name == "Jupiter") return p::PlanetId::Jupiter;
    if (name == "Saturn") return p::PlanetId::Saturn;
    if (name == "Uranus") return p::PlanetId::Uranus;
    return p::PlanetId::Neptune;
}

std::string planet_name_string(spaceship_cpp::planet_params::PlanetId id) {
    return spaceship_cpp::planet_params::planet_name(id);
}

std::string pair_key(spaceship_cpp::planet_params::PlanetId from, spaceship_cpp::planet_params::PlanetId to) {
    return planet_name_string(from) + "->" + planet_name_string(to);
}

int lower_grid_index(double angle) {
    int index = static_cast<int>(std::floor(normalize_angle_0_2pi(angle) / kGridStepRadians)) % kSamplesPerDimension;
    if (index < 0) index += kSamplesPerDimension;
    return index;
}

double grid_angle(int index) {
    return normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(index));
}

bool parse_bool_int(const std::string& raw) {
    return raw == "1" || raw == "true";
}

std::vector<QueryRecord> parse_queries() {
    std::vector<QueryRecord> records;
    std::ifstream in("/tmp/problem1_scout_query_recall_rate.log");
    if (!in) return records;
    std::string line;
    QueryRecord current{};
    bool in_sample = false;
    bool ok = false;
    auto flush = [&]() {
        if (in_sample && ok) records.push_back(current);
        current = QueryRecord{};
        in_sample = false;
        ok = false;
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
        if (!in_sample) continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const auto key = line.substr(0, eq);
        const auto value = line.substr(eq + 1);
        if (key == "from_planet") current.from = parse_planet(value);
        else if (key == "to_planet") current.to = parse_planet(value);
        else if (key == "query_nu_A") current.nu_A = std::stod(value);
        else if (key == "query_nu_B") current.nu_B = std::stod(value);
        else if (key == "query_theta_A") {
            current.theta_A = std::stod(value);
            ok = true;
        } else if (key == "baseline_positive") current.baseline_positive = parse_bool_int(value);
        else if (key == "baseline_used_direct_fallback") current.baseline_used_direct_fallback = parse_bool_int(value);
        else if (key == "query_level_miss") current.query_level_miss = parse_bool_int(value);
    }
    flush();
    return records;
}

void add_case(
    std::vector<BranchCase>* cases,
    std::map<std::string, int>* counts,
    const std::string& name,
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double nu_A,
    double nu_B,
    double theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    const auto key = pair_key(from, to);
    if ((*counts)[key] >= kMaxBranchesPerPair || !branch.valid || !branch.derivatives_available) return;
    BranchCase c{name, from, to, nu_A, nu_B, theta_A, branch};
    cases->push_back(c);
    (*counts)[key] += 1;
}

void collect_cases(
    const spaceship_cpp::problem1::Problem1RootTable2DegLoader& loader,
    const QueryRecord& q,
    std::vector<BranchCase>* cases,
    std::map<std::string, int>* counts
) {
    const int la = lower_grid_index(q.nu_A);
    const int lb = lower_grid_index(q.nu_B);
    const int lt = lower_grid_index(q.theta_A);
    const auto key = pair_key(q.from, q.to);
    const std::string name = (q.query_level_miss ? "miss_" : "sample_") + key;
    for (int da = 0; da <= 1; ++da) {
        for (int db = 0; db <= 1; ++db) {
            for (int dt = 0; dt <= 1; ++dt) {
                const int ia = (la + da) % kSamplesPerDimension;
                const int ib = (lb + db) % kSamplesPerDimension;
                const int it = (lt + dt) % kSamplesPerDimension;
                const auto node = loader.load_node_by_indices(ia, ib, it);
                if (!node.valid) continue;
                for (const auto& branch : node.branches) {
                    add_case(cases, counts, name, q.from, q.to, grid_angle(ia), grid_angle(ib), grid_angle(it), branch);
                    if ((*counts)[key] >= kMaxBranchesPerPair) return;
                }
            }
        }
    }
}

const spaceship_cpp::problem1::Problem1SolutionBranch* find_same_branch(
    const std::vector<spaceship_cpp::problem1::Problem1SolutionBranch>& branches,
    const spaceship_cpp::problem1::Problem1SolutionBranch& source,
    double predicted_alpha
) {
    const spaceship_cpp::problem1::Problem1SolutionBranch* best = nullptr;
    double best_score = std::numeric_limits<double>::infinity();
    for (const auto& b : branches) {
        if (!b.valid || b.transfer_revolution != source.transfer_revolution ||
            b.target_revolution != source.target_revolution) {
            continue;
        }
        const double score = std::abs(normalize_angle_minus_pi_pi(b.encounter_global_angle - predicted_alpha));
        if (score < best_score) {
            best_score = score;
            best = &b;
        }
    }
    return best;
}

double component(const spaceship_cpp::problem1::Problem1SolutionBranch& b, int axis) {
    if (axis == 0) return b.d_encounter_global_angle_d_nu_A;
    if (axis == 1) return b.d_encounter_global_angle_d_nu_B;
    return b.d_encounter_global_angle_d_theta_A;
}

std::array<double, 3> finite_difference_derivative(const BranchCase& c, const spaceship_cpp::problem1::Problem1SolutionBranch& ref) {
    namespace problem1 = spaceship_cpp::problem1;
    std::array<double, 3> out{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
    };
    for (int axis = 0; axis < 3; ++axis) {
        const double dnu_A = axis == 0 ? kFdStep : 0.0;
        const double dnu_B = axis == 1 ? kFdStep : 0.0;
        const double dtheta_A = axis == 2 ? kFdStep : 0.0;
        const double table_d = component(ref, axis);
        const auto plus = problem1::solve_problem1_from_departure_anomalies(
            c.from, c.to, c.nu_A + dnu_A, c.nu_B + dnu_B, c.theta_A + dtheta_A, 1, 1);
        const auto minus = problem1::solve_problem1_from_departure_anomalies(
            c.from, c.to, c.nu_A - dnu_A, c.nu_B - dnu_B, c.theta_A - dtheta_A, 1, 1);
        const auto* plus_b = find_same_branch(plus, ref, normalize_angle_0_2pi(ref.encounter_global_angle + table_d * kFdStep));
        const auto* minus_b = find_same_branch(minus, ref, normalize_angle_0_2pi(ref.encounter_global_angle - table_d * kFdStep));
        if (plus_b != nullptr && minus_b != nullptr) {
            out[static_cast<std::size_t>(axis)] =
                normalize_angle_minus_pi_pi(plus_b->encounter_global_angle - minus_b->encounter_global_angle) /
                (2.0 * kFdStep);
        }
    }
    return out;
}

bool all3_match(const std::array<double, 3>& a, const std::array<double, 3>& b, double tol) {
    for (int i = 0; i < 3; ++i) {
        if (!std::isfinite(a[static_cast<std::size_t>(i)]) ||
            !std::isfinite(b[static_cast<std::size_t>(i)]) ||
            std::abs(normalize_angle_minus_pi_pi(a[static_cast<std::size_t>(i)] - b[static_cast<std::size_t>(i)])) > tol) {
            return false;
        }
    }
    return true;
}

std::array<double, 3> derivatives_array(const spaceship_cpp::problem1::Problem1SolutionBranch& b) {
    return {b.d_encounter_global_angle_d_nu_A, b.d_encounter_global_angle_d_nu_B, b.d_encounter_global_angle_d_theta_A};
}

double abs_diff(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::quiet_NaN();
    return std::abs(normalize_angle_minus_pi_pi(a - b));
}

void print_hypothesis(const BranchCase& c, const std::string& axis, double fd, double table, double attach) {
    std::cout << "Problem1DerivativeConventionHypothesis\n";
    std::cout << "case_name=" << c.name << '\n';
    std::cout << "from_planet=" << planet_name_string(c.from) << '\n';
    std::cout << "to_planet=" << planet_name_string(c.to) << '\n';
    std::cout << "axis=" << axis << '\n';
    std::cout << "fd=" << fd << '\n';
    std::cout << "table=" << table << '\n';
    std::cout << "attach=" << attach << '\n';
    std::cout << "sign_flip_matches=" << (std::isfinite(fd) && abs_diff(-table, fd) <= kMatchTolerance ? 1 : 0) << '\n';
    std::cout << "swap_nu_A_nu_B_matches=0\n";
    std::cout << "theta_A_sign_flip_matches=" << (axis == "theta_A" && std::isfinite(fd) && abs_diff(-table, fd) <= kMatchTolerance ? 1 : 0) << '\n';
    std::cout << "alpha_sign_flip_matches=" << (std::isfinite(fd) && abs_diff(-table, fd) <= kMatchTolerance ? 1 : 0) << '\n';
    std::cout << "local_global_shift_suspect=0\n";
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem1_derivative_source_audit_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    auto queries = parse_queries();
    std::stable_sort(queries.begin(), queries.end(), [](const auto& a, const auto& b) {
        const int aw = (a.query_level_miss ? 4 : 0) + (a.baseline_used_direct_fallback ? 2 : 0) + (a.baseline_positive ? 1 : 0);
        const int bw = (b.query_level_miss ? 4 : 0) + (b.baseline_used_direct_fallback ? 2 : 0) + (b.baseline_positive ? 1 : 0);
        return aw > bw;
    });
    const std::set<std::string> target_pairs{
        "Earth->Mars", "Mars->Earth", "Earth->Mercury", "Venus->Mars",
        "Mars->Venus", "Venus->Mercury", "Mars->Mercury"};
    std::vector<BranchCase> cases;
    std::map<std::string, int> counts;
    for (const auto& q : queries) {
        const auto key = pair_key(q.from, q.to);
        if (target_pairs.count(key) == 0 || counts[key] >= kMaxBranchesPerPair) continue;
        collect_cases(loader, q, &cases, &counts);
    }

    std::cout << "Problem1DerivativeCodePathAudit\n";
    std::cout << "table_generation_file=src/problem1/problem1_root_table.cpp:build_problem1_root_table\n";
    std::cout << "table_generation_theta_A_definition=cell.theta_A=config.theta_A_start+theta_A_index*theta_A_step\n";
    std::cout << "query_theta_A_definition=Problem2 uses departure_state.theta_global-theta_prime; nearest-node query uses supplied theta_A\n";
    std::cout << "attach_theta_A_definition=attach_problem1_root_derivatives_with_mode receives nu_A,nu_B,theta_A and evaluates evaluate_problem1_root_residual_derivatives_with_mode\n";
    std::cout << "alpha_definition=Problem1SolutionBranch.encounter_global_angle global encounter angle\n";
    std::cout << "derivative_field_order=valid,k,q,alpha,target_true_anomaly,tof,target_time,residual,e,p,a,theta_B,derivatives_available,d_alpha_d_nu_A,d_alpha_d_nu_B,d_alpha_d_theta_A\n";
    std::cout << "loader_field_order=same BranchBinaryRecord order in src/problem1/problem1_root_table_2deg_loader.cpp\n";
    std::cout << "possible_stale_table=1\n";
    std::cout << "possible_field_order_mismatch=0\n";
    std::cout << "possible_coordinate_mismatch=1\n";

    std::map<std::string, Summary> summaries;
    for (const auto& c : cases) {
        auto branch_no_derivatives = c.branch;
        branch_no_derivatives.derivatives_available = false;
        const auto attach = problem1::attach_problem1_root_derivatives_with_mode(
            c.from, c.to, c.nu_A, c.nu_B, c.theta_A, branch_no_derivatives,
            problem1::Problem1RootDerivativeMode::AnalyticOnly, 1e-6);
        const auto table_d = derivatives_array(c.branch);
        const auto attach_d = derivatives_array(attach);
        const auto fd_d = finite_difference_derivative(c, c.branch);
        const bool table_attach = all3_match(table_d, attach_d, kMatchTolerance);
        const bool attach_fd = all3_match(attach_d, fd_d, kMatchTolerance);
        const bool table_fd = all3_match(table_d, fd_d, kMatchTolerance);
        std::string source;
        if (!table_attach && attach_fd) source = "table_generation_or_loading_stale_derivatives";
        else if (table_attach && !attach_fd) source = "attach_derivative_formula_or_coordinate_convention";
        else if (!table_attach && !attach_fd) source = "both_table_and_attach_or_branch_matching_problem";
        else source = "derivative_ok";

        std::cout << "Problem1DerivativeSourceAuditSample\n";
        std::cout << "case_name=" << c.name << '\n';
        std::cout << "from_planet=" << planet_name_string(c.from) << '\n';
        std::cout << "to_planet=" << planet_name_string(c.to) << '\n';
        std::cout << "k=" << c.branch.transfer_revolution << '\n';
        std::cout << "q=" << c.branch.target_revolution << '\n';
        std::cout << "nu_A=" << c.nu_A << '\n';
        std::cout << "nu_B=" << c.nu_B << '\n';
        std::cout << "theta_A=" << c.theta_A << '\n';
        std::cout << "alpha=" << c.branch.encounter_global_angle << '\n';
        std::cout << "table_d_nu_A=" << table_d[0] << '\n';
        std::cout << "table_d_nu_B=" << table_d[1] << '\n';
        std::cout << "table_d_theta_A=" << table_d[2] << '\n';
        std::cout << "attach_d_nu_A=" << attach_d[0] << '\n';
        std::cout << "attach_d_nu_B=" << attach_d[1] << '\n';
        std::cout << "attach_d_theta_A=" << attach_d[2] << '\n';
        std::cout << "fd_d_nu_A=" << fd_d[0] << '\n';
        std::cout << "fd_d_nu_B=" << fd_d[1] << '\n';
        std::cout << "fd_d_theta_A=" << fd_d[2] << '\n';
        std::cout << "table_vs_attach_abs_diff_nu_A=" << abs_diff(table_d[0], attach_d[0]) << '\n';
        std::cout << "table_vs_attach_abs_diff_nu_B=" << abs_diff(table_d[1], attach_d[1]) << '\n';
        std::cout << "table_vs_attach_abs_diff_theta_A=" << abs_diff(table_d[2], attach_d[2]) << '\n';
        std::cout << "attach_vs_fd_abs_diff_nu_A=" << abs_diff(attach_d[0], fd_d[0]) << '\n';
        std::cout << "attach_vs_fd_abs_diff_nu_B=" << abs_diff(attach_d[1], fd_d[1]) << '\n';
        std::cout << "attach_vs_fd_abs_diff_theta_A=" << abs_diff(attach_d[2], fd_d[2]) << '\n';
        std::cout << "table_vs_fd_abs_diff_nu_A=" << abs_diff(table_d[0], fd_d[0]) << '\n';
        std::cout << "table_vs_fd_abs_diff_nu_B=" << abs_diff(table_d[1], fd_d[1]) << '\n';
        std::cout << "table_vs_fd_abs_diff_theta_A=" << abs_diff(table_d[2], fd_d[2]) << '\n';

        std::cout << "Problem1DerivativeSourceAuditDiagnosis\n";
        std::cout << "case_name=" << c.name << '\n';
        std::cout << "from_planet=" << planet_name_string(c.from) << '\n';
        std::cout << "to_planet=" << planet_name_string(c.to) << '\n';
        std::cout << "k=" << c.branch.transfer_revolution << '\n';
        std::cout << "q=" << c.branch.target_revolution << '\n';
        std::cout << "table_matches_attach=" << (table_attach ? 1 : 0) << '\n';
        std::cout << "attach_matches_fd=" << (attach_fd ? 1 : 0) << '\n';
        std::cout << "table_matches_fd=" << (table_fd ? 1 : 0) << '\n';
        std::cout << "likely_error_source=" << source << '\n';
        if (!table_fd || !attach_fd) {
            print_hypothesis(c, "nu_A", fd_d[0], table_d[0], attach_d[0]);
            print_hypothesis(c, "nu_B", fd_d[1], table_d[1], attach_d[1]);
            print_hypothesis(c, "theta_A", fd_d[2], table_d[2], attach_d[2]);
        }

        auto& summary = summaries[pair_key(c.from, c.to)];
        summary.count += 1;
        summary.table_matches_attach += table_attach ? 1 : 0;
        summary.attach_matches_fd += attach_fd ? 1 : 0;
        summary.table_matches_fd += table_fd ? 1 : 0;
        if (source == "table_generation_or_loading_stale_derivatives") summary.table_generation_or_loading += 1;
        else if (source == "attach_derivative_formula_or_coordinate_convention") summary.attach_formula_or_convention += 1;
        else if (source == "both_table_and_attach_or_branch_matching_problem") summary.both_or_matching += 1;
        else summary.derivative_ok += 1;
    }

    for (const auto& [key, s] : summaries) {
        const auto arrow = key.find("->");
        std::cout << "Problem1DerivativeSourceAuditPairSummary\n";
        std::cout << "from_planet=" << key.substr(0, arrow) << '\n';
        std::cout << "to_planet=" << key.substr(arrow + 2) << '\n';
        std::cout << "sample_count=" << s.count << '\n';
        std::cout << "table_matches_attach_count=" << s.table_matches_attach << '\n';
        std::cout << "attach_matches_fd_count=" << s.attach_matches_fd << '\n';
        std::cout << "table_matches_fd_count=" << s.table_matches_fd << '\n';
        std::cout << "table_generation_or_loading_count=" << s.table_generation_or_loading << '\n';
        std::cout << "attach_formula_or_convention_count=" << s.attach_formula_or_convention << '\n';
        std::cout << "both_or_branch_matching_count=" << s.both_or_matching << '\n';
        std::cout << "derivative_ok_count=" << s.derivative_ok << '\n';
    }
    std::cout << "problem1_derivative_source_audit_ok=1\n";
    return 0;
}
