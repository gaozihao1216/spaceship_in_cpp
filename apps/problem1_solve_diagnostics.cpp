/*
 * 文件作用：运行 Problem 1 直接求解器的命令行诊断程序。
 * 主要工作：读取全局配置、调用残差/求解接口，并输出候选转移轨道的诊断信息。
 */
#include "spaceship_cpp/common/common.hpp"
#include "spaceship_cpp/config/global_config.hpp"
#include "spaceship_cpp/problem1/problem1.hpp"
#include "spaceship_cpp/problem1/problem1_diagnostics.hpp"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

const char* planet_label(spaceship_cpp::planet_params::PlanetId id) {
    return spaceship_cpp::planet_params::planet_name(id);
}

void run_case(
    spaceship_cpp::planet_params::PlanetId departure,
    spaceship_cpp::planet_params::PlanetId target,
    double transfer_perihelion_angle,
    const spaceship_cpp::config::Problem1SolveDefaults& defaults,
    std::vector<spaceship_cpp::problem1::Problem1Candidate>* all_candidates
) {
    const spaceship_cpp::problem1::Problem1SolveInput input =
        spaceship_cpp::config::make_problem1_solve_input(
            departure,
            target,
            0.0,
            transfer_perihelion_angle,
            defaults);

    const auto candidates = spaceship_cpp::problem1::solve_problem1(input);
    const auto summary = spaceship_cpp::problem1::summarize_problem1_candidates(candidates);

    std::cout
        << planet_label(departure) << " -> " << planet_label(target)
        << ", transfer_perihelion_angle=" << transfer_perihelion_angle
        << ", candidate_count=" << summary.candidate_count << '\n';

    for (const auto& candidate : candidates) {
        std::cout
            << "  k=" << candidate.transfer_revolution
            << ", q=" << candidate.target_revolution
            << ", encounter_global_angle=" << candidate.encounter_global_angle
            << ", time_of_flight_days=" << candidate.time_of_flight_seconds / 86400.0
            << ", transfer_e_raw=" << candidate.residual_result.transfer_e_raw
            << ", transfer_e=" << candidate.residual_result.transfer_e
            << ", transfer_p=" << candidate.residual_result.transfer_p
            << ", relative_residual=" << candidate.relative_residual
            << ", root_bracket_width=" << candidate.root_bracket_width
            << ", bisection_iterations=" << candidate.bisection_iterations
            << '\n';
        all_candidates->push_back(candidate);
    }
}

}  // namespace

int main() {
    namespace common = spaceship_cpp::common;
    namespace config = spaceship_cpp::config;
    namespace planet_params = spaceship_cpp::planet_params;
    namespace problem1 = spaceship_cpp::problem1;

    std::filesystem::create_directories("results");

    const config::GlobalConfig& cfg = config::global_config();
    config::Problem1SolveDefaults solve_defaults = cfg.problem1_solve;
    solve_defaults.phi_scan_count = cfg.problem1_diagnostics.solve_phi_scan_count;
    solve_defaults.max_candidate_relative_residual =
        cfg.problem1_diagnostics.max_candidate_relative_residual;

    std::vector<problem1::Problem1Candidate> all_candidates;
    const std::vector<planet_params::PlanetId> targets{
        planet_params::PlanetId::Mars,
        planet_params::PlanetId::Venus,
        planet_params::PlanetId::Mercury,
    };
    const std::vector<double> angles{
        0.0,
        common::kPi / 4.0,
        common::kHalfPi,
        common::kPi,
        1.5 * common::kPi,
    };

    for (planet_params::PlanetId target : targets) {
        for (double angle : angles) {
            run_case(
                planet_params::PlanetId::Earth,
                target,
                angle,
                solve_defaults,
                &all_candidates);
        }
    }

    problem1::write_problem1_candidates_csv(all_candidates, "results/problem1_solve_diagnostics.csv");
    return 0;
}
