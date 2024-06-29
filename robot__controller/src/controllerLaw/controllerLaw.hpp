#include "communication/communication.h"
#include "robot/robot.hpp"
#include "planner/jogPlanner.h"
#include <vector>
#include <queue>
struct ControllerParam
{
    char paramName[2] = {0};
    double value = 0.0;
};

template <int _Dofs>
struct ControllerParamBase
{
public:
    ControllerParam jointParam1[_Dofs];
    ControllerParam jointParam2[_Dofs];
    ControllerParam jointParam3[_Dofs];
    ControllerParam cartesianParam1[6];
    ControllerParam cartesianParam2[6];
    ControllerParam cartesianParam3[6];
};

template <int _Dofs>
void dynamicSetParameterTempl(TaskSpace taskSpace, const ControllerParamBase<_Dofs> &config, unsigned int time,
                              Eigen::Matrix<double, 6, 6> &cartesianK1, Eigen::Matrix<double, 6, 6> &cartesianK2,
                              Eigen::Matrix<double, 6, 6> &cartesianK1_d, Eigen::Matrix<double, 6, 6> &cartesianK2_d,
                              Eigen::Matrix<double, _Dofs, _Dofs> &jointK1, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2,
                              Eigen::Matrix<double, _Dofs, _Dofs> &jointK1_d, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2_d)
{
    if (taskSpace == TaskSpace::jointSpace)
    {
        if (time == 0)
        {
            for (int i = 0; i < _Dofs; i++)
            {
                jointK1(i, i) = config.jointParam1[i].value;
                jointK2(i, i) = config.jointParam2[i].value;
            }
        }
        for (int i = 0; i < _Dofs; i++)
        {
            jointK1_d(i, i) = config.jointParam1[i].value;
            jointK2_d(i, i) = config.jointParam2[i].value;
        }
    }
    else
    {
        if (time == 0)
        {
            for (int i = 0; i < 6; i++)
            {
                cartesianK1(i, i) = config.cartesianParam1[i].value;
                cartesianK2(i, i) = config.cartesianParam2[i].value;
            }
        }
        for (int i = 0; i < 6; i++)
        {
            cartesianK1_d(i, i) = config.cartesianParam1[i].value;
            cartesianK2_d(i, i) = config.cartesianParam2[i].value;
        }
    }
}
template <int _Dofs>
void controllerParamRenewTempl(double filterParams,
                               Eigen::Matrix<double, 6, 6> &cartesianK1, Eigen::Matrix<double, 6, 6> &cartesianK2,
                               Eigen::Matrix<double, 6, 6> &cartesianK1_d, Eigen::Matrix<double, 6, 6> &cartesianK2_d,
                               Eigen::Matrix<double, _Dofs, _Dofs> &jointK1, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2,
                               Eigen::Matrix<double, _Dofs, _Dofs> &jointK1_d, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2_d)
{
    jointK1 = filterParams * jointK1_d + (1.0 - filterParams) * jointK1;
    jointK2 = filterParams * jointK2_d + (1.0 - filterParams) * jointK2;
    cartesianK1 = filterParams * cartesianK1_d + (1.0 - filterParams) * cartesianK1;
    cartesianK2 = filterParams * cartesianK2_d + (1.0 - filterParams) * cartesianK2;
}
//
// 控制律基类
template <int _Dofs>
class ControllerLaw
{
public:
    TaskSpace taskSpace;
    std::string controllerLawName;

    // 当前时刻误差
    Eigen::Matrix<double, _Dofs, 1> jointError;
    Eigen::Matrix<double, _Dofs, 1> djointError;
    Eigen::Matrix<double, 6, 1> cartesianError;
    Eigen::Matrix<double, 6, 1> dcartesianError;
    Eigen::Matrix<double, 6, 1> cartesianOldError;

    // 当前时刻期望 关节空间
    Eigen::Matrix<double, _Dofs, 1> q_d;
    Eigen::Matrix<double, _Dofs, 1> dq_d;
    Eigen::Matrix<double, _Dofs, 1> ddq_d;

    // 当前时刻期望 笛卡尔空间
    Eigen::Vector3d position_d;
    Eigen::Quaterniond orientation_d;
    Eigen::Vector3d dposition_d;
    Eigen::Quaterniond dorientation_d;
    Eigen::Matrix<double, 6, 1> ddX_d; // position+orientation

    // controllerLaw
    Eigen::Matrix<double, _Dofs, 1> tau_d;
    Eigen::Matrix<double, _Dofs, 1> qc;

public:
    ControllerLaw(const ControllerLaw &) = delete;
    void operator=(const ControllerLaw &) = delete;

    ControllerLaw() = delete;
    ControllerLaw(TaskSpace createTaskSpace, std::string createControllerLawName);
    virtual ~ControllerLaw();

    void calError(my_robot::Robot<_Dofs> *robot, double cycleTime);
    void calWaitDesireNext(Eigen::Matrix<double, _Dofs, 1> &q_hold, Eigen::Vector3d &position_hold, Eigen::Quaterniond &orientation_hold);
    void calRunStopDesireNext(std::vector<std::queue<double>> &q_dQueue, std::vector<std::queue<double>> &dq_dQueue,
                              std::vector<std::queue<double>> &ddq_dQueue, std::vector<std::queue<double>> &x_dQueue,
                              std::vector<std::queue<double>> &dx_dQueue, std::vector<std::queue<double>> &ddx_dQueue);
    void calJogMove(my_robot::Robot<_Dofs> *robot, bool &jogMoveFlag, double jogSpeed,
                                          double cycleTime, int jogDir, int jogNum);
    void calJogStop(my_robot::Robot<_Dofs> *robot, bool &jogStopFlag, double cycleTime, int jogNum);

    virtual void setTaskSpace(TaskSpace newTaskSpace);
    virtual void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d) = 0;
    virtual void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time) = 0;
    virtual void controllerParamRenew(double filterParams) = 0;
};
template <int _Dofs>
ControllerLaw<_Dofs>::~ControllerLaw()
{
}
template <int _Dofs>
ControllerLaw<_Dofs>::ControllerLaw(TaskSpace createTaskSpace, std::string createControllerLawName) : jointError(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      djointError(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      cartesianError(Eigen::Matrix<double, 6, 1>::Zero()),
                                                                                                      dcartesianError(Eigen::Matrix<double, 6, 1>::Zero()),
                                                                                                      cartesianOldError(Eigen::Matrix<double, 6, 1>::Zero()),
                                                                                                      q_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      dq_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      ddq_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      position_d(Eigen::Matrix<double, 3, 1>::Zero()),
                                                                                                      orientation_d(Eigen::Quaterniond::Identity()),
                                                                                                      dposition_d(Eigen::Matrix<double, 3, 1>::Zero()),
                                                                                                      dorientation_d(Eigen::Quaterniond::Identity()),
                                                                                                      ddX_d(Eigen::Matrix<double, 6, 1>::Zero()),
                                                                                                      tau_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      qc(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                                                      taskSpace(createTaskSpace),
                                                                                                      controllerLawName(createControllerLawName)

{
    // 全部初始化为0
}

template <int _Dofs>
void ControllerLaw<_Dofs>::setTaskSpace(TaskSpace newTaskSpace)
{
    this->taskSpace = newTaskSpace;
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calError(my_robot::Robot<_Dofs> *robot, double cycleTime)
{
    switch (this->taskSpace)
    {
    case TaskSpace::jointSpace:
    {
        // 关节误差与误差导数
        this->jointError = this->q_d - robot->getq();
        this->djointError = this->dq_d - robot->getdq();
        break;
    }
    case TaskSpace::cartesianSpace:
    {
        // 笛卡尔位姿误差
        Eigen::Quaterniond orientation = robot->getOrientation();
        this->cartesianError.head(3) = this->position_d - robot->getPosition();
        if (this->orientation_d.coeffs().dot(orientation.coeffs()) < 0.0)
        {
            orientation.coeffs() << -orientation.coeffs();
        }
        Eigen::Quaterniond error_quaternion(orientation.inverse() * this->orientation_d);
        this->cartesianError.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
        this->cartesianError.tail(3) << robot->getT().rotation() * this->cartesianError.tail(3);

        this->dcartesianError = (this->cartesianError - this->cartesianOldError) / cycleTime;
        this->cartesianOldError = this->cartesianError;
        break;
    }
    }
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calWaitDesireNext(Eigen::Matrix<double, _Dofs, 1> &q_hold, Eigen::Vector3d &position_hold, Eigen::Quaterniond &orientation_hold)
{
    if (this->taskSpace == jointSpace)
    {
        this->q_d = q_hold;
        this->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
        this->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
    }
    else
    {
        this->position_d = position_hold;
        this->orientation_d = orientation_hold;
    }
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calRunStopDesireNext(std::vector<std::queue<double>> &q_dQueue, std::vector<std::queue<double>> &dq_dQueue,
                                                std::vector<std::queue<double>> &ddq_dQueue, std::vector<std::queue<double>> &x_dQueue,
                                                std::vector<std::queue<double>> &dx_dQueue, std::vector<std::queue<double>> &ddx_dQueue)
{
    if (this->taskSpace == jointSpace)
    {
        for (int i = 0; i < _Dofs; i++)
        {
            this->q_d[i] = q_dQueue[i].front();
            this->dq_d[i] = dq_dQueue[i].front();
            this->ddq_d[i] = ddq_dQueue[i].front();
            q_dQueue[i].pop();
            dq_dQueue[i].pop();
            ddq_dQueue[i].pop();
        }
    }
    else
    {
        for (int i = 0; i < 2; i++)
        {
            this->position_d[i] = x_dQueue[i].front();
            this->dposition_d[i] = dx_dQueue[i].front();
            x_dQueue[i].pop();
            dx_dQueue[i].pop();
        }
        for (int i = 3; i < 6; i++)
        {
            // this->orientation_d = x_dRunQueue[i].front();
            // this->dorientation_d = dx_dRunQueue[i].front();
            // x_dRunQueue[i].pop();
            // dx_dRunQueue[i].pop();
        }
        for (int i = 0; i < 6; i++)
        {
            this->ddX_d[i] = ddx_dQueue[i].front();
            ddx_dQueue[i].pop();
        }
    }
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calJogMove(my_robot::Robot<_Dofs> *robot, bool &jogMoveFlag, double jogSpeed,
                                      double cycleTime, int jogDir, int jogNum)
{
    if (this->taskSpace == jointSpace)
    {
        Eigen::Matrix<double, _Dofs, 1> dddqLimit = robot->getdddqLimit();
        Eigen::Matrix<double, _Dofs, 1> ddqLimit = robot->getddqLimit();
        Eigen::Matrix<double, _Dofs, 1> dqLimit = jogSpeed * robot->getdqLimit();
        if (jogMoveFlag)
        {
            calJogMovePlan(jogMoveFlag, cycleTime, jogDir, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1], dqLimit[jogNum - 1],
                           this->q_d[jogNum - 1], this->dq_d[jogNum - 1], this->ddq_d[jogNum - 1]);
            jogMoveFlag = false;
        }
        else
        {
            calJogMovePlan(jogMoveFlag, cycleTime, jogDir, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1], dqLimit[jogNum - 1],
                           this->q_d[jogNum - 1], this->dq_d[jogNum - 1], this->ddq_d[jogNum - 1]);
        }
    }
    else
    {
    }
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calJogStop(my_robot::Robot<_Dofs> *robot, bool &jogStopFlag, double cycleTime, int jogNum)
{
    if (this->taskSpace == jointSpace)
    {
        Eigen::Matrix<double, _Dofs, 1> dddqLimit = robot->getdddqLimit();
        Eigen::Matrix<double, _Dofs, 1> ddqLimit = robot->getddqLimit();
        if (jogStopFlag)
        {
            calJogStopPlan(jogStopFlag, cycleTime, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1],
                           this->q_d[jogNum - 1], this->dq_d[jogNum - 1], this->ddq_d[jogNum - 1]);
            jogStopFlag = false;
        }
        else
        {
            calJogStopPlan(jogStopFlag, cycleTime, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1],
                           this->q_d[jogNum - 1], this->dq_d[jogNum - 1], this->ddq_d[jogNum - 1]);
        }
    }
    else
    {
    }
}
//
//
// 计算力矩控制器
template <int _Dofs>
class ComputedTorqueMethod : public ControllerLaw<_Dofs>
{
public:
    // 关节空间
    Eigen::Matrix<double, _Dofs, _Dofs> jointKv;
    Eigen::Matrix<double, _Dofs, _Dofs> jointKp;

    Eigen::Matrix<double, _Dofs, _Dofs> jointKv_d;
    Eigen::Matrix<double, _Dofs, _Dofs> jointKp_d;

    // 笛卡尔空间
    Eigen::Matrix<double, 6, 6> cartesianKp;
    Eigen::Matrix<double, 6, 6> cartesianKv;

    Eigen::Matrix<double, 6, 6> cartesianKp_d;
    Eigen::Matrix<double, 6, 6> cartesianKv_d;

public:
    ComputedTorqueMethod(const ComputedTorqueMethod &) = delete;
    void operator=(const ComputedTorqueMethod &) = delete;
    ComputedTorqueMethod() = delete;

    ~ComputedTorqueMethod();
    explicit ComputedTorqueMethod(TaskSpace taskSpace);

    void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

    void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time);
    void controllerParamRenew(double filterParams);
};
template <int _Dofs>
ComputedTorqueMethod<_Dofs>::~ComputedTorqueMethod()
{
}
template <int _Dofs>
ComputedTorqueMethod<_Dofs>::ComputedTorqueMethod(TaskSpace taskSpace) : ControllerLaw<_Dofs>(taskSpace, "ComputedTorqueMethod"),
                                                                         jointKv(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                                         jointKp(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                                         jointKv_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                                         jointKp_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                                         cartesianKp(Eigen::Matrix<double, 6, 6>::Zero()),
                                                                         cartesianKv(Eigen::Matrix<double, 6, 6>::Zero()),
                                                                         cartesianKp_d(Eigen::Matrix<double, 6, 6>::Zero()),
                                                                         cartesianKv_d(Eigen::Matrix<double, 6, 6>::Zero())
{
    std::cout << "[robotController] 设置控制律: " << this->controllerLawName << std::endl;
}
template <int _Dofs>
void ComputedTorqueMethod<_Dofs>::setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d_in)
{
    if (this->taskSpace == jointSpace)
    {
        this->qc = this->ddq_d + jointKp * this->jointError + jointKv * this->djointError;
        this->tau_d << robot->getM() * (this->qc) + robot->getC() * robot->getdq() /* + G */;
        tau_d_in = this->tau_d;
    }
    else
    {
        this->tau_d << robot->getM() * (robot->getJ_inv() * (this->ddX_d - this->cartesianKp_d * this->cartesianError - this->cartesianKv_d * this->dcartesianError - robot->getExternJ() * robot->getdq())) + robot->getC() * robot->getdq() /* + G */;
        tau_d_in = this->tau_d;
    }
}
template <int _Dofs>
void ComputedTorqueMethod<_Dofs>::dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time)
{
    dynamicSetParameterTempl<_Dofs>(this->taskSpace, config, time, this->cartesianKp, this->cartesianKv, this->cartesianKp_d,
                                    this->cartesianKv_d, this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}
template <int _Dofs>
void ComputedTorqueMethod<_Dofs>::controllerParamRenew(double filterParams)
{
    controllerParamRenewTempl<_Dofs>(filterParams, this->cartesianKp, this->cartesianKv, this->cartesianKp_d, this->cartesianKv_d,
                                     this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}

//
//
// 反步控制器
template <int _Dofs>
class Backstepping : public ControllerLaw<_Dofs>
{
    // 关节空间
    Eigen::Matrix<double, _Dofs, _Dofs> jointK1;
    Eigen::Matrix<double, _Dofs, _Dofs> jointK2;

    Eigen::Matrix<double, _Dofs, _Dofs> jointK1_d;
    Eigen::Matrix<double, _Dofs, _Dofs> jointK2_d;

    // 笛卡尔空间
    Eigen::Matrix<double, 6, 6> cartesianK1;
    Eigen::Matrix<double, 6, 6> cartesianK2;

    Eigen::Matrix<double, 6, 6> cartesianK1_d;
    Eigen::Matrix<double, 6, 6> cartesianK2_d;

    // 误差中间变量
    Eigen::Matrix<double, _Dofs, 1> e1;
    Eigen::Matrix<double, _Dofs, 1> e2;
    Eigen::Matrix<double, _Dofs, 1> r;
    Eigen::Matrix<double, _Dofs, 1> dr;

public:
    Backstepping(const Backstepping &) = delete;
    void operator=(const Backstepping &) = delete;
    Backstepping() = delete;

    ~Backstepping();
    explicit Backstepping(TaskSpace taskSpace);

    void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

    void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time);
    void controllerParamRenew(double filterParams);
};
template <int _Dofs>
Backstepping<_Dofs>::~Backstepping()
{
}
template <int _Dofs>
Backstepping<_Dofs>::Backstepping(TaskSpace taskSpace) : ControllerLaw<_Dofs>(taskSpace, "Backstepping"),
                                                         jointK1(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                         jointK2(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                         jointK1_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                         jointK2_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                         cartesianK1(Eigen::Matrix<double, 6, 6>::Zero()),
                                                         cartesianK2(Eigen::Matrix<double, 6, 6>::Zero()),
                                                         cartesianK1_d(Eigen::Matrix<double, 6, 6>::Zero()),
                                                         cartesianK2_d(Eigen::Matrix<double, 6, 6>::Zero()),
                                                         e1(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                         e2(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                         r(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                         dr(Eigen::Matrix<double, _Dofs, 1>::Zero())
{
    std::cout << "[robotController] 设置控制律: " << this->controllerLawName << std::endl;
}
template <int _Dofs>
void Backstepping<_Dofs>::setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d_in)
{
    // note: 父类是抽象模板类，子类使用其成员需显式指定：this->或base::
    if (this->taskSpace == jointSpace)
    {
        this->e1 = this->jointError;
        this->e2 = this->djointError + this->jointK1 * e1;
        this->r = this->dq_d + this->jointK1 * e1;
        this->dr = this->ddq_d + this->jointK1 * this->djointError;
        this->tau_d << robot->getM() * (this->dr) + robot->getC() * (this->dr) /* + G */ + this->jointK2 * this->e2 + this->e1;
        tau_d_in = this->tau_d;
    }
    else
    {
        this->tau_d << robot->getM() * (robot->getJ_inv() * (this->ddX_d - this->cartesianK1 * this->cartesianError - this->cartesianK2 * this->dcartesianError - robot->getExternJ() * robot->getdq())) + robot->getC() * robot->getdq() /* + G */;
        tau_d_in = this->tau_d;
    }
}

template <int _Dofs>
void Backstepping<_Dofs>::dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time)
{
    dynamicSetParameterTempl<_Dofs>(this->taskSpace, config, time, this->cartesianK1, this->cartesianK2, this->cartesianK1_d,
                                    this->cartesianK2_d, this->jointK1, this->jointK2, this->jointK1_d, this->jointK2_d);
}
template <int _Dofs>
void Backstepping<_Dofs>::controllerParamRenew(double filterParams)
{
    controllerParamRenewTempl<_Dofs>(filterParams, this->cartesianK1, this->cartesianK2, this->cartesianK1_d, this->cartesianK2_d,
                                     this->jointK1, this->jointK2, this->jointK1_d, this->jointK2_d);
}

//
//
// PD+重力补偿
template <int _Dofs>
class PD : public ControllerLaw<_Dofs>
{
public:
    // 关节空间
    Eigen::Matrix<double, _Dofs, _Dofs> jointKv;
    Eigen::Matrix<double, _Dofs, _Dofs> jointKp;

    Eigen::Matrix<double, _Dofs, _Dofs> jointKv_d;
    Eigen::Matrix<double, _Dofs, _Dofs> jointKp_d;

    // 笛卡尔空间
    Eigen::Matrix<double, 6, 6> cartesianKp;
    Eigen::Matrix<double, 6, 6> cartesianKv;

    Eigen::Matrix<double, 6, 6> cartesianKp_d;
    Eigen::Matrix<double, 6, 6> cartesianKv_d;

public:
    PD(const PD &) = delete;
    void operator=(const PD &) = delete;
    PD() = delete;

    ~PD();
    explicit PD(TaskSpace taskSpace);

    void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

    void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time);
    void controllerParamRenew(double filterParams);
};
template <int _Dofs>
PD<_Dofs>::~PD()
{
}
template <int _Dofs>
PD<_Dofs>::PD(TaskSpace taskSpace) : ControllerLaw<_Dofs>(taskSpace, "PD"),
                                     jointKv(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                     jointKp(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                     jointKv_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                     jointKp_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                     cartesianKp(Eigen::Matrix<double, 6, 6>::Zero()),
                                     cartesianKv(Eigen::Matrix<double, 6, 6>::Zero()),
                                     cartesianKp_d(Eigen::Matrix<double, 6, 6>::Zero()),
                                     cartesianKv_d(Eigen::Matrix<double, 6, 6>::Zero())
{
    std::cout << "[robotController] 设置控制律: " << this->controllerLawName << std::endl;
}
template <int _Dofs>
void PD<_Dofs>::setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d_in)
{
    if (this->taskSpace == jointSpace)
    {
        this->tau_d << this->ddq_d + jointKp * this->jointError + jointKv * this->djointError /* + G */;
        tau_d_in = this->tau_d;
    }
    else
    {
        this->tau_d << robot->getJ_inv() * (this->ddX_d - this->cartesianKp_d * this->cartesianError - this->cartesianKv_d * this->dcartesianError - robot->getExternJ() * robot->getdq()) + robot->getC() * robot->getdq() /* + G */;
        tau_d_in = this->tau_d;
    }
}
template <int _Dofs>
void PD<_Dofs>::dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time)
{
    dynamicSetParameterTempl<_Dofs>(this->taskSpace, config, time, this->cartesianKp, this->cartesianKv, this->cartesianKp_d,
                                    this->cartesianKv_d, this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}
template <int _Dofs>
void PD<_Dofs>::controllerParamRenew(double filterParams)
{
    controllerParamRenewTempl<_Dofs>(filterParams, this->cartesianKp, this->cartesianKv, this->cartesianKp_d, this->cartesianKv_d,
                                     this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}

// 工厂
template <int _Dofs>
bool newControllerLaw(std::unique_ptr<ControllerLaw<_Dofs>> &controllerLaw, ControllerLawType controllerLawType, TaskSpace taskSpace)
{
    switch (controllerLawType)
    {
    case ControllerLawType::ComputedTorqueMethod_:
        controllerLaw = std::make_unique<ComputedTorqueMethod<_Dofs>>(taskSpace);
        break;
    case ControllerLawType::Backstepping_:
        controllerLaw = std::make_unique<Backstepping<_Dofs>>(taskSpace);
        break;
    case ControllerLawType::PD_:
        controllerLaw = std::make_unique<PD<_Dofs>>(taskSpace);
        break;
    default:
        controllerLaw = std::make_unique<ComputedTorqueMethod<_Dofs>>(TaskSpace::jointSpace);
        return false;
        break;
    }
    return true;
}