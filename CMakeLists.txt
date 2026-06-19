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
    src/problem1/diagnostics/problem1_diagnostics.cpp
    src/problem1/tables/problem1_table.cpp
    src/problem2/problem2.cpp
    src/problem2/flyby/problem2_slingshot.cpp
    src/trajectory/flyby_physics.cpp
    src/trajectory/orbit_velocity.cpp
    src/bfs/bfs.cpp
    src/bfs/adapters/problem2_angle_frame_adapter.cpp
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

enable_testing()
add_subdirectory(tests)
