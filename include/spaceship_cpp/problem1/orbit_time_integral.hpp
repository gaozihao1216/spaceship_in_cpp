/*
 * 文件作用：声明 Problem 1 使用的轨道飞行时间积分辅助函数。
 * 主要工作：根据轨道形态和异常角差计算无量纲飞行时间。
 */
#pragma once

#include "spaceship_cpp/common/orbit_math.hpp"

namespace spaceship_cpp::problem1 {

using spaceship_cpp::common::ConicType;
using spaceship_cpp::common::OrbitFResult;
using spaceship_cpp::common::classify_conic;
using spaceship_cpp::common::is_orbit_F_defined;
using spaceship_cpp::common::orbit_F_integrand;
using spaceship_cpp::common::try_orbit_F;
using spaceship_cpp::common::orbit_F;
using spaceship_cpp::common::try_orbit_F_xi_derivative;
using spaceship_cpp::common::orbit_F_xi_derivative;

}  // namespace spaceship_cpp::problem1
