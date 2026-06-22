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
    src/problem2/flyby/problem2_slingshot.cpp
    src/problem2/flyby/problem2_flyby_constraint.cpp
    src/problem2/flyby/problem2_flyby_G_search.cpp
    src/problem2/flyby/problem2_flyby_solve.cpp
    src/problem2/flyby/problem2_theta_prime_scan.cpp
    src/problem2/flyby/problem2_theta_prime_route_a.cpp
    src/trajectory/flyby_physics.cpp
    src/trajectory/orbit_velocity.cpp
    src/bfs/bfs.cpp
    src/bfs/fixed_sequence_bfs.cpp
    src/bfs/leg0_theta_feasibility.cpp
    src/bfs/trajectory_search_config.cpp
    src/bfs/free_path_bfs.cpp
    src/bfs/step3_top_k_sequences.cpp
    src/bfs/step4_fixed_sequence_fine_search.cpp
    src/bfs/adapters/problem2_angle_frame_adapter.cpp
)

target_include_directories(spaceship_cpp PUBLIC include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(spaceship_cpp PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(trajectory_search
    apps/trajectory_search.cpp
)
target_link_libraries(trajectory_search PRIVATE spaceship_cpp)

enable_testing()
add_subdirectory(tests)
