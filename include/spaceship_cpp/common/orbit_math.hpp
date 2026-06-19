/*
 * 文件作用：声明轨道数学基础工具。
 * 主要工作：提供椭圆轨道半径、异常角转换和轨道时间积分相关函数。
 */
#pragma once

#include <string>

namespace spaceship_cpp::common {

enum class ConicType {
    Elliptic,
    Parabolic,
    Hyperbolic,
};

struct OrbitFResult {
    bool ok;
    double value;
    std::string message;
};

ConicType classify_conic(double e, double eccentricity_tolerance = 1e-12);

bool is_orbit_F_defined(double e, double xi, double tolerance = 1e-12);

double orbit_F_integrand(double e, double xi);

OrbitFResult try_orbit_F(double e, double xi, double tolerance = 1e-12);

double orbit_F(double e, double xi, double tolerance = 1e-12);

OrbitFResult try_orbit_F_xi_derivative(double e, double xi, double tolerance = 1e-12);

double orbit_F_xi_derivative(double e, double xi, double tolerance = 1e-12);

OrbitFResult try_orbit_F_theta_derivative(double e, double theta, double tolerance = 1e-12);

double orbit_F_theta_derivative(double e, double theta, double tolerance = 1e-12);

OrbitFResult try_orbit_F_e_derivative(double e, double theta, double tolerance = 1e-12);

double orbit_F_e_derivative(double e, double theta, double tolerance = 1e-12);

}  // namespace spaceship_cpp::common
