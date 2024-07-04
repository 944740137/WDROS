#ifndef FRANKA_FK
#define FRANKA_FK

#include <Eigen/Dense> 
Eigen::Matrix4d franka_FK(const Eigen::Matrix<double, 7, 1> &q);


#endif 