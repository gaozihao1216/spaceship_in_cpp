/*
 * 文件作用：声明 BFS 与 Problem 2 之间的角度坐标系适配接口。
 * 主要工作：在全局近日点角和飞掠行星局部角之间做双向转换。
 */
#pragma once

#include "spaceship_cpp/planet_params/planet_params.hpp"

namespace spaceship_cpp::bfs {

// 将全局近日点角转成 Problem 2 以飞掠行星近日点为零点的局部角。
double global_periapsis_angle_to_problem2_local(
    planet_params::PlanetId flyby_planet,
    double global_periapsis_angle
);

// 将 Problem 2 局部近日点角转回全局坐标，供轨迹搜索统一保存。
double problem2_local_periapsis_angle_to_global(
    planet_params::PlanetId flyby_planet,
    double local_periapsis_angle
);

}  // namespace spaceship_cpp::bfs
