// Analytical Franka inverse kinematics using q7 as redundant parameter
// - Yanhao He, February 2020

#ifndef FRANKA_IK_HE_H
#define FRANKA_IK_HE_H

#include <array>
#include <cmath>
#include "Eigen/Dense"

// inverse kinematics w.r.t. End Effector Frame (using Franka Hand data)
std::array< std::array<double, 7>, 4 > franka_IK_EE ( std::array<double, 16> O_T_EE_array,
                                                      double q7,
                                                      std::array<double, 7> q_actual_array );

// "Case-Consistent" inverse kinematics w.r.t. End Effector Frame (using Franka Hand data)
std::array<double, 7> franka_IK_EE_CC(std::array<double, 16> O_T_EE_array,
                                      double q7,
                                      std::array<double, 7> q_actual_array);
void NewFunction();
#endif // FRANKA_IK_HE_HPP