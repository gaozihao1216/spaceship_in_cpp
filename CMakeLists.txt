cmake_minimum_required(VERSION 3.20)

project(spaceship_in_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(spaceship_cpp STATIC
    src/config/global_config.cpp
    src/core_types/core_types.cpp
    src/planet_params/planet_params.cpp
    src/common/common.cpp
    src/common/orbit_math.cpp
    src/common/time_utils.cpp
    src/problem1/problem1.cpp
    src/problem1/problem1_diagnostics.cpp
    src/problem1/problem1_root_table.cpp
    src/problem1/problem1_root_table_2deg_loader.cpp
    src/problem1/problem1_root_table_nearest_node_query.cpp
    src/problem1/problem1_table.cpp
    src/problem2/problem2.cpp
    src/problem2/problem2_gravity_assist_expansion.cpp
    src/problem2/problem2_gravity_assist_solver.cpp
    src/problem2/problem2_gravity_assist_transition.cpp
    src/problem2/problem2_slingshot.cpp
    src/trajectory/flyby_physics.cpp
    src/trajectory/orbit_velocity.cpp
    src/bfs/bfs.cpp
    src/bfs/fixed_launch_theta_search.cpp
    src/bfs/fixed_launch_time_theta_scan.cpp
    src/bfs/fixed_sequence_launch_time_scan.cpp
    src/bfs/fixed_sequence_theta_refinement.cpp
    src/bfs/fixed_sequence_theta_scan.cpp
    src/bfs/fixed_sequence_theta_search.cpp
    src/bfs/gravity_assist_step_expansion.cpp
    src/bfs/initial_launch_expansion.cpp
    src/bfs/problem2_angle_frame_adapter.cpp
    src/bfs/theta_candidate_selector.cpp
    src/bfs/theta_launch_feasibility_scout.cpp
)

target_include_directories(spaceship_cpp PUBLIC include)

add_executable(problem1_solve_diagnostics
    apps/problem1_solve_diagnostics.cpp
)
target_link_libraries(problem1_solve_diagnostics PRIVATE spaceship_cpp)

add_executable(problem1_table_diagnostics
    apps/problem1_table_diagnostics.cpp
)
target_link_libraries(problem1_table_diagnostics PRIVATE spaceship_cpp)

add_executable(build_problem1_root_table_2deg
    tools/build_problem1_root_table_2deg.cpp
)
target_link_libraries(build_problem1_root_table_2deg PRIVATE spaceship_cpp)
set_target_properties(build_problem1_root_table_2deg PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools
)

add_executable(validate_problem1_root_table_2deg_chunk
    tools/validate_problem1_root_table_2deg_chunk.cpp
)
target_link_libraries(validate_problem1_root_table_2deg_chunk PRIVATE spaceship_cpp)
set_target_properties(validate_problem1_root_table_2deg_chunk PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools
)

add_executable(scan_problem1_root_table_2deg_derivative_failures
    tools/scan_problem1_root_table_2deg_derivative_failures.cpp
)
target_link_libraries(scan_problem1_root_table_2deg_derivative_failures PRIVATE spaceship_cpp)
set_target_properties(scan_problem1_root_table_2deg_derivative_failures PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools
)

add_executable(patch_problem1_root_table_2deg_derivative_failures
    tools/patch_problem1_root_table_2deg_derivative_failures.cpp
)
target_link_libraries(patch_problem1_root_table_2deg_derivative_failures PRIVATE spaceship_cpp)
set_target_properties(patch_problem1_root_table_2deg_derivative_failures PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools
)

enable_testing()
add_subdirectory(tests)
