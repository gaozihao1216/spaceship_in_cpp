/*
 * 文件作用：声明轨道数学基础工具。
 * 主要工作：提供椭圆轨道半径、异常角转换和轨道时间积分相关函数。
 */
#pragma once

#include <string>

namespace spaceship_cpp::common {

enum class ConicType {
    // e < 1，对应闭合椭圆轨道。
    Elliptic,
    // e = 1 的边界轨道。
    Parabolic,
    // e > 1，对应双曲线轨道。
    Hyperbolic,
};

struct OrbitFResult {
    // ok 表示积分函数在当前 e/xi 下有定义。
    bool ok;
    // value 是积分值或导数值；ok=false 时不可使用。
    double value;
    // 失败原因，便于诊断是奇点、非法偏心率还是分母退化。
    std::string message;
};

// 根据偏心率分类轨道类型，解决不同时间积分公式的分派问题。
ConicType classify_conic(double e, double eccentricity_tolerance = 1e-12);

// 判断 F(e, xi) 是否避开 1 + e cos(xi) = 0 的奇点。
bool is_orbit_F_defined(double e, double xi, double tolerance = 1e-12);

// 轨道时间积分的被积函数。
double orbit_F_integrand(double e, double xi);

// 安全计算轨道时间积分，失败时返回原因而不是抛异常。
OrbitFResult try_orbit_F(double e, double xi, double tolerance = 1e-12);

// 直接计算轨道时间积分，适合调用方已经保证输入有效的路径。
double orbit_F(double e, double xi, double tolerance = 1e-12);

// 计算 F 对 xi/theta/e 的导数，用于敏感性分析和表格诊断。
OrbitFResult try_orbit_F_xi_derivative(double e, double xi, double tolerance = 1e-12);

double orbit_F_xi_derivative(double e, double xi, double tolerance = 1e-12);

OrbitFResult try_orbit_F_theta_derivative(double e, double theta, double tolerance = 1e-12);

double orbit_F_theta_derivative(double e, double theta, double tolerance = 1e-12);

OrbitFResult try_orbit_F_e_derivative(double e, double theta, double tolerance = 1e-12);

double orbit_F_e_derivative(double e, double theta, double tolerance = 1e-12);

}  // namespace spaceship_cpp::common
