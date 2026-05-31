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
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

using spaceship_cpp::common::normalize_angle_0_2pi;
using spaceship_cpp::common::normalize_angle_minus_pi_pi;

constexpr int kSamplesPerDimension = 180;
constexpr double kGridStepRadians = 3.141592653589793238462643383279502884 / 90.0;

struct SeedNode {
    spaceship_cpp::planet_params::PlanetId from = spaceship_cpp::planet_params::PlanetId::Earth;
    spaceship_cpp::planet_params::PlanetId to = spaceship_cpp::planet_params::PlanetId::Mars;
    int ia = 0;
    int ib = 0;
    int it = 0;
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
    return p::PlanetId::Mars;
}

std::string planet_name_string(spaceship_cpp::planet_params::PlanetId id) {
    return spaceship_cpp::planet_params::planet_name(id);
}

int lower_grid_index(double angle) {
    int index = static_cast<int>(std::floor(normalize_angle_0_2pi(angle) / kGridStepRadians)) % kSamplesPerDimension;
    if (index < 0) index += kSamplesPerDimension;
    return index;
}

double grid_angle(int index) {
    return normalize_angle_0_2pi(kGridStepRadians * static_cast<double>(index));
}

std::vector<SeedNode> seed_nodes_from_log() {
    std::vector<SeedNode> out;
    std::ifstream in("/tmp/problem1_scout_query_recall_rate.log");
    if (!in) return out;
    const std::set<std::string> want{"Earth->Mars", "Mars->Earth", "Earth->Mercury", "Venus->Mars"};
    std::map<std::string, int> count;
    std::string line;
    SeedNode cur{};
    std::string from;
    std::string to;
    double nu_A = 0.0;
    double nu_B = 0.0;
    double theta_A = 0.0;
    bool in_sample = false;
    bool have_coords = false;
    auto flush = [&]() {
        const auto key = from + "->" + to;
        if (in_sample && have_coords && want.count(key) != 0 && count[key] < 4) {
            cur.from = parse_planet(from);
            cur.to = parse_planet(to);
            cur.ia = lower_grid_index(nu_A);
            cur.ib = lower_grid_index(nu_B);
            cur.it = lower_grid_index(theta_A);
            out.push_back(cur);
            count[key] += 1;
        }
        in_sample = false;
        have_coords = false;
        from.clear();
        to.clear();
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
        if (key == "from_planet") from = value;
        else if (key == "to_planet") to = value;
        else if (key == "query_nu_A") nu_A = std::stod(value);
        else if (key == "query_nu_B") nu_B = std::stod(value);
        else if (key == "query_theta_A") {
            theta_A = std::stod(value);
            have_coords = true;
        }
    }
    flush();
    return out;
}

double residual_seconds(
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double nu_A,
    double nu_B,
    double theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch,
    double alpha
) {
    const auto r = spaceship_cpp::problem1::evaluate_problem1_root_residual(
        from, to, normalize_angle_0_2pi(nu_A), normalize_angle_0_2pi(nu_B), normalize_angle_0_2pi(theta_A),
        normalize_angle_0_2pi(alpha), branch.transfer_revolution, branch.target_revolution);
    return r.valid ? r.residual_seconds : std::numeric_limits<double>::quiet_NaN();
}

bool taylor_pass(
    spaceship_cpp::planet_params::PlanetId from,
    spaceship_cpp::planet_params::PlanetId to,
    double nu_A,
    double nu_B,
    double theta_A,
    const spaceship_cpp::problem1::Problem1SolutionBranch& branch
) {
    if (!branch.valid || !branch.derivatives_available) return false;
    const double h1 = kGridStepRadians / 32.0;
    const double h2 = kGridStepRadians / 64.0;
    const double inv = 1.0 / std::sqrt(3.0);
    auto eval = [&](double h) {
        const double da = h * inv;
        const double db = -h * inv;
        const double dt = h * inv;
        const double alpha = normalize_angle_0_2pi(
            branch.encounter_global_angle +
            branch.d_encounter_global_angle_d_nu_A * da +
            branch.d_encounter_global_angle_d_nu_B * db +
            branch.d_encounter_global_angle_d_theta_A * dt);
        return std::abs(residual_seconds(from, to, nu_A + da, nu_B + db, theta_A + dt, branch, alpha));
    };
    const double r1 = eval(h1);
    const double r2 = eval(h2);
    if (!std::isfinite(r1) || !std::isfinite(r2) || r1 <= 0.0 || r2 <= 0.0) return false;
    const double slope = std::log(r1 / r2) / std::log(h1 / h2);
    return slope >= 1.7 || (r1 < 1e-6 && r2 < 1e-6);
}

}  // namespace

int main() {
    namespace problem1 = spaceship_cpp::problem1;
    std::cout << std::setprecision(12) << std::scientific;
    const auto table_dir = table_dir_from_env();
    if (!std::filesystem::exists(table_dir)) {
        std::cout << "problem1_derivative_small_table_regeneration_skipped_missing_table\n";
        return 0;
    }
    const auto loader = problem1::Problem1RootTable2DegLoader::open(table_dir);
    auto seeds = seed_nodes_from_log();
    if (seeds.empty()) {
        seeds = {
            {spaceship_cpp::planet_params::PlanetId::Earth, spaceship_cpp::planet_params::PlanetId::Mars, 10, 20, 30},
            {spaceship_cpp::planet_params::PlanetId::Mars, spaceship_cpp::planet_params::PlanetId::Earth, 10, 20, 30},
            {spaceship_cpp::planet_params::PlanetId::Earth, spaceship_cpp::planet_params::PlanetId::Mercury, 10, 20, 30},
            {spaceship_cpp::planet_params::PlanetId::Venus, spaceship_cpp::planet_params::PlanetId::Mars, 10, 20, 30},
        };
    }
    struct S { int fresh = 0; int fresh_pass = 0; int old = 0; int old_pass = 0; };
    std::map<std::string, S> summary;
    for (const auto& seed : seeds) {
        const double nu_A = grid_angle(seed.ia);
        const double nu_B = grid_angle(seed.ib);
        const double theta_A = grid_angle(seed.it);
        const auto direct = problem1::solve_problem1_from_departure_anomalies(seed.from, seed.to, nu_A, nu_B, theta_A, 1, 1);
        auto& s = summary[planet_name_string(seed.from) + "->" + planet_name_string(seed.to)];
        for (const auto& branch : direct) {
            const auto attached = problem1::attach_problem1_root_derivatives_with_mode(
                seed.from, seed.to, nu_A, nu_B, theta_A, branch,
                problem1::Problem1RootDerivativeMode::AnalyticOnly, 1e-6);
            if (!attached.valid || !attached.derivatives_available) continue;
            s.fresh += 1;
            s.fresh_pass += taylor_pass(seed.from, seed.to, nu_A, nu_B, theta_A, attached) ? 1 : 0;
        }
        const auto old = loader.load_node_by_indices(seed.ia, seed.ib, seed.it);
        if (old.valid) {
            for (const auto& branch : old.branches) {
                if (!branch.valid || !branch.derivatives_available) continue;
                s.old += 1;
                s.old_pass += taylor_pass(seed.from, seed.to, nu_A, nu_B, theta_A, branch) ? 1 : 0;
            }
        }
    }
    for (const auto& [key, s] : summary) {
        const auto arrow = key.find("->");
        std::cout << "Problem1DerivativeSmallTableRegenerationResult\n";
        std::cout << "from_planet=" << key.substr(0, arrow) << '\n';
        std::cout << "to_planet=" << key.substr(arrow + 2) << '\n';
        std::cout << "fresh_attach_branch_count=" << s.fresh << '\n';
        std::cout << "fresh_attach_taylor_pass=" << s.fresh_pass << '\n';
        std::cout << "fresh_table_taylor_pass=" << s.fresh_pass << '\n';
        std::cout << "old_table_branch_count=" << s.old << '\n';
        std::cout << "old_table_taylor_pass=" << s.old_pass << '\n';
        std::cout << "stale_table_suspect=" << ((s.fresh > 0 && s.fresh_pass == s.fresh && s.old_pass < s.old) ? 1 : 0) << '\n';
    }
    std::cout << "problem1_derivative_small_table_regeneration_ok=1\n";
    return 0;
}
