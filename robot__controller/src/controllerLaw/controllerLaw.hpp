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
void dynamicSetParameterTempl(const ControllerParamBase<_Dofs> &config, unsigned int time,
                              Eigen::Matrix<double, _Dofs, _Dofs> &jointK1, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2,
                              Eigen::Matrix<double, _Dofs, _Dofs> &jointK1_d, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2_d)
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
template <int _Dofs>
void controllerParamRenewTempl(double filterParams,
                               Eigen::Matrix<double, _Dofs, _Dofs> &jointK1, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2,
                               Eigen::Matrix<double, _Dofs, _Dofs> &jointK1_d, Eigen::Matrix<double, _Dofs, _Dofs> &jointK2_d)
{
    jointK1 = filterParams * jointK1_d + (1.0 - filterParams) * jointK1;
    jointK2 = filterParams * jointK2_d + (1.0 - filterParams) * jointK2;
}
//
// 控制律基类
template <int _Dofs>
class ControllerLaw
{
public:
    std::string controllerLawName;

    // 当前时刻误差
    Eigen::Matrix<double, _Dofs, 1> jointError;
    Eigen::Matrix<double, _Dofs, 1> djointError;

    // 当前时刻期望 关节空间
    Eigen::Matrix<double, _Dofs, 1> q_d;
    Eigen::Matrix<double, _Dofs, 1> dq_d;
    Eigen::Matrix<double, _Dofs, 1> ddq_d;

    // controllerLaw
    Eigen::Matrix<double, _Dofs, 1> tau_d;
    Eigen::Matrix<double, _Dofs, 1> qc;

public:
    ControllerLaw(const ControllerLaw &) = delete;
    void operator=(const ControllerLaw &) = delete;

    ControllerLaw(std::string createControllerLawName);
    virtual ~ControllerLaw();

    void calError(my_robot::Robot<_Dofs> *robot);
    void calWaitDesireNext(Eigen::Matrix<double, _Dofs, 1> &q_hold);
    void calRunStopDesireNext(std::vector<std::queue<double>> &q_dQueue, std::vector<std::queue<double>> &dq_dQueue,
                              std::vector<std::queue<double>> &ddq_dQueue);
    void setNullSpaceDesire(Eigen::Matrix<double, _Dofs, 1> &q_ns);

    virtual void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d) = 0;
    virtual void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time) = 0;
    virtual void controllerParamRenew(double filterParams) = 0;
};
template <int _Dofs>
ControllerLaw<_Dofs>::~ControllerLaw()
{
}
template <int _Dofs>
ControllerLaw<_Dofs>::ControllerLaw(std::string createControllerLawName) : jointError(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           djointError(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           q_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           dq_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           ddq_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           tau_d(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           qc(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                                           controllerLawName(createControllerLawName)
{
    // 全部初始化为0
}

template <int _Dofs>
void ControllerLaw<_Dofs>::calError(my_robot::Robot<_Dofs> *robot)
{
    // 关节误差与误差导数
    this->jointError = this->q_d - robot->getq();
    this->djointError = this->dq_d - robot->getdq();
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calWaitDesireNext(Eigen::Matrix<double, _Dofs, 1> &q_hold)
{
    this->q_d = q_hold;
    this->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
    this->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
}
template <int _Dofs>
void ControllerLaw<_Dofs>::calRunStopDesireNext(std::vector<std::queue<double>> &q_dQueue, std::vector<std::queue<double>> &dq_dQueue,
                                                std::vector<std::queue<double>> &ddq_dQueue)
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



template <int _Dofs>
void ControllerLaw<_Dofs>::setNullSpaceDesire(Eigen::Matrix<double, _Dofs, 1> &q_ns)
{
    this->q_nullSapce = q_ns;
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

public:
    ComputedTorqueMethod(const ComputedTorqueMethod &) = delete;
    void operator=(const ComputedTorqueMethod &) = delete;

    ~ComputedTorqueMethod();
    explicit ComputedTorqueMethod();

    void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

    void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time);
    void controllerParamRenew(double filterParams);
};
template <int _Dofs>
ComputedTorqueMethod<_Dofs>::~ComputedTorqueMethod()
{
}
template <int _Dofs>
ComputedTorqueMethod<_Dofs>::ComputedTorqueMethod() : ControllerLaw<_Dofs>("ComputedTorqueMethod"),
                                                      jointKv(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                      jointKp(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                      jointKv_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                                      jointKp_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero())
{
    std::cout << "[robotController] 设置控制律: " << this->controllerLawName << std::endl;
}
template <int _Dofs>
void ComputedTorqueMethod<_Dofs>::setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d_in)
{
    this->qc = this->ddq_d + jointKp * this->jointError + jointKv * this->djointError;
    this->tau_d << robot->getM() * (this->qc) + robot->getC() * robot->getdq() /* + G */;
    tau_d_in = this->tau_d;
}
template <int _Dofs>
void ComputedTorqueMethod<_Dofs>::dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time)
{
    dynamicSetParameterTempl<_Dofs>(config, time, this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}
template <int _Dofs>
void ComputedTorqueMethod<_Dofs>::controllerParamRenew(double filterParams)
{
    controllerParamRenewTempl<_Dofs>(filterParams, this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
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

    // 误差中间变量
    Eigen::Matrix<double, _Dofs, 1> e1;
    Eigen::Matrix<double, _Dofs, 1> e2;
    Eigen::Matrix<double, _Dofs, 1> r;
    Eigen::Matrix<double, _Dofs, 1> dr;

public:
    Backstepping(const Backstepping &) = delete;
    void operator=(const Backstepping &) = delete;

    ~Backstepping();
    explicit Backstepping();

    void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

    void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time);
    void controllerParamRenew(double filterParams);
};
template <int _Dofs>
Backstepping<_Dofs>::~Backstepping()
{
}
template <int _Dofs>
Backstepping<_Dofs>::Backstepping() : ControllerLaw<_Dofs>("Backstepping"),
                                      jointK1(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                      jointK2(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                      jointK1_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                                      jointK2_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
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
    this->e1 = this->jointError;
    this->e2 = this->djointError + this->jointK1 * e1;
    this->r = this->dq_d + this->jointK1 * e1;
    this->dr = this->ddq_d + this->jointK1 * this->djointError;
    this->tau_d << robot->getM() * (this->dr) + robot->getC() * (this->dr) /* + G */ + this->jointK2 * this->e2 + this->e1;
    tau_d_in = this->tau_d;
}

template <int _Dofs>
void Backstepping<_Dofs>::dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time)
{
    dynamicSetParameterTempl<_Dofs>(config, time, this->jointK1, this->jointK2, this->jointK1_d, this->jointK2_d);
}
template <int _Dofs>
void Backstepping<_Dofs>::controllerParamRenew(double filterParams)
{
    controllerParamRenewTempl<_Dofs>(filterParams, this->jointK1, this->jointK2, this->jointK1_d, this->jointK2_d);
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

public:
    PD(const PD &) = delete;
    void operator=(const PD &) = delete;

    ~PD();
    explicit PD();

    void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

    void dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time);
    void controllerParamRenew(double filterParams);
};
template <int _Dofs>
PD<_Dofs>::~PD()
{
}
template <int _Dofs>
PD<_Dofs>::PD() : ControllerLaw<_Dofs>("PD"),
                  jointKv(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                  jointKp(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                  jointKv_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                  jointKp_d(Eigen::Matrix<double, _Dofs, _Dofs>::Zero())
{
    std::cout << "[robotController] 设置控制律: " << this->controllerLawName << std::endl;
}
template <int _Dofs>
void PD<_Dofs>::setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d_in)
{
    this->tau_d << this->ddq_d + jointKp * this->jointError + jointKv * this->djointError /* + G */;
    tau_d_in = this->tau_d;
}
template <int _Dofs>
void PD<_Dofs>::dynamicSetParameter(const ControllerParamBase<_Dofs> &config, unsigned int time)
{
    dynamicSetParameterTempl<_Dofs>(config, time, this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}
template <int _Dofs>
void PD<_Dofs>::controllerParamRenew(double filterParams)
{
    controllerParamRenewTempl<_Dofs>(filterParams, this->jointKp, this->jointKv, this->jointKp_d, this->jointKv_d);
}

// 工厂
template <int _Dofs>
bool newControllerLaw(std::unique_ptr<ControllerLaw<_Dofs>> &controllerLaw, ControllerLawType controllerLawType)
{
    switch (controllerLawType)
    {
    case ControllerLawType::ComputedTorqueMethod_:
        controllerLaw = std::make_unique<ComputedTorqueMethod<_Dofs>>();
        break;
    case ControllerLawType::Backstepping_:
        controllerLaw = std::make_unique<Backstepping<_Dofs>>();
        break;
    case ControllerLawType::PD_:
        controllerLaw = std::make_unique<PD<_Dofs>>();
        break;
    default:
        controllerLaw = std::make_unique<ComputedTorqueMethod<_Dofs>>();
        return false;
        break;
    }
    return true;
}