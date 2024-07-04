#include <iostream>
#include <cmath>
#include "frankaFK.h"

using namespace Eigen;

// 正解函数
Eigen::Matrix4d franka_FK(const Eigen::Matrix<double, 7, 1> &q)
{
    Eigen::Matrix4d T = Matrix4d::Identity(); // 初始化为单位矩阵

    // 定义DH参数
    double d[] = {0.333, 0, 0.316, 0, 0.384, 0, 0.107};                                  // d参数
    double a[] = {0, 0, 0, 0.0825, -0.0825, 0, 0.088};                                    // a参数
    double alpha[] = {0, -M_PI / 2, M_PI / 2, M_PI / 2, -M_PI / 2, M_PI / 2, M_PI / 2}; // alpha参数

    for (int i = 0; i < 7; i++)
    {
        double theta_i = q[i];

        Eigen::Matrix4d A_i;
        A_i << cos(theta_i), -sin(theta_i), 0, a[i],
            sin(theta_i) * cos(alpha[i]), cos(theta_i) * cos(alpha[i]), -sin(alpha[i]), -sin(alpha[i]) * d[i],
            sin(theta_i) * sin(alpha[i]), cos(theta_i) * sin(alpha[i]), cos(alpha[i]), cos(alpha[i]) * d[i],
            0, 0, 0, 1;
        T = T * A_i;
        // std::cout << "--------" << std::endl;
        // std::cout << "q" << q[i] << std::endl;
        // std::cout << A_i << std::endl;
        // std::cout << T << std::endl;
        // std::cout << "--------" << std::endl;
    }
    // std::cout << "+++++" << q.transpose() << std::endl;
    return T;
}