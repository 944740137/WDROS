// #pragma once
#include "controllerLaw/controllerLaw.hpp"
#include "planner/runPlanner.hpp"
#include "planner/jogPlanner.h"
#include "algorithm/frankaIK.h"
#include "algorithm/frankaFK.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "math.h"

static int debug = 0;
// constexpr double param = 0.0000001;
namespace robot_controller
{
    // _Dofs自由度机器人控制器
    template <int _Dofs, typename pubDataType>
    class Controller
    {

    public:
        // debug，绘图
        int recordPeriod = 1;           // 数据记录周期
        unsigned int time = 0;          // 当前时刻
        std::ofstream myfile;           // 记录文件io对象
        double filterParams = 0.005;    // 调参滤波参数
        const double cycleTime = 0.001; // 运行周期0.001s

        // 点动参数
        bool jogSign = false;
        int jogNum = 0;      // 点动轴
        int jogDir = 0;      // 点动方向
        double jogSpeed = 0; // 0-1
        double jogSpeed_d = 0;
        bool jogMoveFlag = false;
        bool jogStopFlag = false;

        // 运行参数
        bool newPlan = false;
        Eigen::Matrix<double, _Dofs, 1> q_calQueue;
        Eigen::Matrix<double, 6, 1> x_calQueue;
        Eigen::Matrix<double, 3, 1> rotationAxis;
        double rotationAngle = 0;
        double runSpeed = 0;
        double runSpeed_d = 0;
        Eigen::Matrix<double, 3, 3> R1;

        // 坐标系
        TaskSpace runTaskSpace = TaskSpace::jointSpace;

        // 急停参数
        bool newStop = false;

        // 当前运行状态
        RunStatus nowControllerStatus = RunStatus::wait_; // 当前状态

        // 运行队列
        std::vector<std::queue<double>> q_dRunQueue{_Dofs};
        std::vector<std::queue<double>> dq_dRunQueue{_Dofs};
        std::vector<std::queue<double>> ddq_dRunQueue{_Dofs};
        std::vector<std::queue<double>> x_dRunQueue{4};
        std::vector<std::queue<double>> dx_dRunQueue{4};
        std::vector<std::queue<double>> ddx_dRunQueue{4};

        // 急停队列
        std::vector<std::queue<double>> q_dStopQueue{_Dofs};
        std::vector<std::queue<double>> dq_dStopQueue{_Dofs};
        std::vector<std::queue<double>> ddq_dStopQueue{_Dofs};
        std::vector<std::queue<double>> x_dStopQueue{4};
        std::vector<std::queue<double>> dx_dStopQueue{4};
        std::vector<std::queue<double>> ddx_dStopQueue{4};

        // 进程通讯
        bool connectStatus = false;
        Communication communicationModel;
        struct RobotData *robotDataBuff = nullptr;
        struct ControllerCommand *controllerCommandBUff = nullptr;
        struct ControllerState *controllerStateBUff = nullptr;

        //  控制律
        ControllerLawType controllerLawType = ControllerLawType::PD_;
        ControllerLawType controllerLawType_d = ControllerLawType::PD_;
        std::unique_ptr<ControllerLaw<_Dofs>> controllerLaw;
        ControllerParamBase<_Dofs> controllerParam;

        // 规划器
        std::unique_ptr<Planner<_Dofs>> jointPlanner;
        std::unique_ptr<Planner<4>> taskPlanner;
        PlannerType plannerType = PlannerType::Quintic_;
        PlannerType plannerType_d = PlannerType::Quintic_;

    public:
        Controller(const Controller &) = delete;
        void operator=(const Controller &) = delete;

        Controller();
        virtual ~Controller();

        // api
        void setRecord(int recordPeriod);
        const std::string &getControllerLawName();
        void dynamicSetParameter();
        void changeControllerLaw(ControllerLawType type, my_robot::Robot<_Dofs> *robot);
        void changePlanner(PlannerType type);
        void initStateToMaster();
        void calJointRunQueue(my_robot::Robot<_Dofs> *robot);
        void calJointStopQueue(my_robot::Robot<_Dofs> *robot);
        void calCartesianRunQueue(my_robot::Robot<_Dofs> *robot);
        void calCartesianStopQueue(my_robot::Robot<_Dofs> *robot);
        void calJointJogMove(my_robot::Robot<_Dofs> *robot);
        void calJointJogStop(my_robot::Robot<_Dofs> *robot);
        void calCartesianJogMove(my_robot::Robot<_Dofs> *robot);
        void calCartesianJogStop(my_robot::Robot<_Dofs> *robot);
        bool RtoAxisAngle(Eigen::Matrix3d &R, Eigen::Matrix<double, 3, 1> &axis, double &angle);
        bool AxisAngletoR(Eigen::Matrix3d &R, Eigen::Matrix<double, 3, 1> &axis, double &angle);
        void startMotion();
        void stopMotion();
        void clearMotionQueue(std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);
        bool checkMotionQueue(std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);
        void HandleiKFailed();
        // init
        void init(int recordPeriod, my_robot::Robot<_Dofs> *robot);

        // run
        void updateTime();
        void controllerParamRenew();
        void communication(my_robot::Robot<_Dofs> *robot);
        void updateStatus(my_robot::Robot<_Dofs> *robot);
        void calDesireNext(my_robot::Robot<_Dofs> *robot); // 计算下一个周期的期望
        void calError(my_robot::Robot<_Dofs> *robot);
        void setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

        //
        void recordData(my_robot::Robot<_Dofs> *robot);
        void pubData(pubDataType &param_debug, my_robot::Robot<_Dofs> *robot);
    };

    //****************************************************************default****************************************************************//
    // default
    template <int _Dofs, typename pubDataType>
    Controller<_Dofs, pubDataType>::~Controller()
    {
    }
    template <int _Dofs, typename pubDataType>
    Controller<_Dofs, pubDataType>::Controller() : q_calQueue(Eigen::Matrix<double, _Dofs, 1>::Zero())
    {
    }

    //******************************************************************API******************************************************************//
    // API
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::setRecord(int recordPeriod)
    {
        this->recordPeriod = recordPeriod;
    }
    template <int _Dofs, typename pubDataType>
    const std::string &Controller<_Dofs, pubDataType>::getControllerLawName()
    {
        if (controllerLaw != nullptr)
            return this->controllerLaw->controllerLawName;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::changeControllerLaw(ControllerLawType type, my_robot::Robot<_Dofs> *robot)
    {
        if (!newControllerLaw(controllerLaw, type))
            printf("ControllerLaw create Error\n");
        dynamicSetParameter();
        this->controllerLaw->initDesire(robot->getq(), robot->getPosition(), robot->getOrientation());
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::changePlanner(PlannerType type)
    {
        if (!newPlanner(this->jointPlanner, type))
            printf("jointPlanner create Error\n");
        if (!newPlanner(this->taskPlanner, type))
            printf("taskPlanner create Error\n");
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::dynamicSetParameter()
    {
        if (controllerLaw.get() != nullptr)
            controllerLaw->dynamicSetParameter(this->controllerParam, this->time);
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::initStateToMaster()
    {
        this->controllerStateBUff->robotDof = _Dofs;
        strcpy(this->controllerStateBUff->name, "panda\0\0\0");
        this->controllerStateBUff->controllerStatus = this->nowControllerStatus;
        this->controllerCommandBUff->plannerTaskSpace = this->runTaskSpace;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::startMotion()
    {
        this->nowControllerStatus = RunStatus::run_; // 开始运动
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::stopMotion()
    {
        this->nowControllerStatus = RunStatus::stop_; // 开始停止
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::clearMotionQueue(std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
    {
        for (int i = 0; i < q_d.size(); i++)
        {
            std::queue<double>().swap(q_d[i]);
            std::queue<double>().swap(dq_d[i]);
            std::queue<double>().swap(ddq_d[i]);
        }
    }
    template <int _Dofs, typename pubDataType>
    bool Controller<_Dofs, pubDataType>::checkMotionQueue(std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
    {
        int size = 0;
        for (int i = 0; i < q_d.size(); i++)
        {
            size = size + q_d[i].size();
        }
        if ((size / q_d.size()) != q_d[0].size())
        {
            this->clearMotionQueue(q_d, dq_d, ddq_d);
            std::cout << "[robotController] 队列长度不一致 清空处理" << std::endl;
            return false;
        }
        return true;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calJointRunQueue(my_robot::Robot<_Dofs> *robot)
    {
        Eigen::Matrix<double, _Dofs, 1> velLimit = this->runSpeed * robot->getdqLimit(); // 这里的velLimit是速度
        Eigen::Matrix<double, _Dofs, 1> accLimit = robot->getddqLimit();
        Eigen::Matrix<double, _Dofs, 1> jerkLimit = robot->getdddqLimit();

        // note: 显示指定模板类型
        if (!this->jointPlanner->calPlanQueue(true, this->cycleTime, velLimit, accLimit, jerkLimit, robot->getq(), this->q_calQueue, robot->getdq(),
                                              this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue))
        {
            this->clearMotionQueue(this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue);
            return;
        }
        if (!this->checkMotionQueue(this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue))
            this->clearMotionQueue(this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue);
        else
            this->startMotion();
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calJointStopQueue(my_robot::Robot<_Dofs> *robot)
    {
        Eigen::Matrix<double, _Dofs, 1> maxAcc = robot->getddqLimit();
        Eigen::Matrix<double, _Dofs, 1> maxJerk = robot->getdddqLimit();
        Eigen::Matrix<double, _Dofs, 1> q = Eigen::Matrix<double, _Dofs, 1>::Zero();
        Eigen::Matrix<double, _Dofs, 1> dq = Eigen::Matrix<double, _Dofs, 1>::Zero();
        Eigen::Matrix<double, _Dofs, 1> ddq = Eigen::Matrix<double, _Dofs, 1>::Zero();
        for (int i = 0; i < _Dofs; i++)
        {
            q[i] = this->controllerLaw->q_d[i];
            dq[i] = this->controllerLaw->dq_d[i];
            ddq[i] = this->controllerLaw->ddq_d[i];
        }
        if (!calStopPlan<_Dofs>(true, this->cycleTime, maxJerk, maxAcc, q, dq, ddq, this->q_dStopQueue, this->dq_dStopQueue, this->ddq_dStopQueue))
        {
            this->clearMotionQueue(this->q_dStopQueue, this->dq_dStopQueue, this->ddq_dStopQueue);
            return;
        }

        if (!this->checkMotionQueue(this->q_dStopQueue, this->dq_dStopQueue, this->ddq_dStopQueue))
            this->clearMotionQueue(this->q_dStopQueue, this->dq_dStopQueue, this->ddq_dStopQueue);
        else
        {
            this->clearMotionQueue(this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue);
            this->stopMotion();
        }
    }
    template <int _Dofs, typename pubDataType>
    bool Controller<_Dofs, pubDataType>::RtoAxisAngle(Eigen::Matrix3d &R, Eigen::Matrix<double, 3, 1> &axis, double &angle)
    {
        double r = (R(0, 0) + R(1, 1) + R(2, 2) - 1) / 2;
        if (std::fabs(r - 1) < param) // 0
            return false;

        if (std::fabs(r + 1) < param) // pi
        {
            angle = M_PI;
            auto compute_axis = [&](int i, int j, int k)
            {
                double val = (R(i, i) + 1) / 2;
                if (val < param)
                    val = 0;
                axis[i] = std::sqrt(val);
                axis[j] = R(i, j) / (2 * axis[i]);
                axis[k] = R(i, k) / (2 * axis[i]);
            };
            if (std::fabs(R(0, 0)) >= std::max(std::fabs(R(1, 1)), std::fabs(R(2, 2))) && std::fabs(R(0, 0) + 1) > param)
                compute_axis(0, 1, 2);
            else if (std::fabs(R(1, 1)) >= std::max(std::fabs(R(0, 0)), std::fabs(R(2, 2))) && std::fabs(R(1, 1) + 1) > param)
                compute_axis(1, 0, 2);
            else if (std::fabs(R(2, 2)) >= std::max(std::fabs(R(0, 0)), std::fabs(R(1, 1))) && std::fabs(R(2, 2) + 1) > param)
                compute_axis(2, 0, 1);
        }
        else
        {
            angle = std::acos(r);
            double d = 1 / (2 * sin(angle));
            axis[0] = d * (R(2, 1) - R(1, 2));
            axis[1] = d * (R(0, 2) - R(2, 0));
            axis[2] = d * (R(1, 0) - R(0, 1));
        }
        return true;
    }
    template <int _Dofs, typename pubDataType>
    bool Controller<_Dofs, pubDataType>::AxisAngletoR(Eigen::Matrix3d &R, Eigen::Matrix<double, 3, 1> &axis, double &angle)
    {
        double x = axis[0];
        double y = axis[1];
        double z = axis[2];
        R(0, 0) = x * x * (1 - std::cos(angle)) + std::cos(angle);
        R(0, 1) = x * y * (1 - std::cos(angle)) - z * std::sin(angle);
        R(0, 2) = x * z * (1 - std::cos(angle)) + y * std::sin(angle);

        R(1, 0) = x * y * (1 - std::cos(angle)) + z * std::sin(angle);
        R(1, 1) = y * y * (1 - std::cos(angle)) + std::cos(angle);
        R(1, 2) = y * z * (1 - std::cos(angle)) - x * std::sin(angle);

        R(2, 0) = x * z * (1 - std::cos(angle)) - y * std::sin(angle);
        R(2, 1) = y * z * (1 - std::cos(angle)) + x * std::sin(angle);
        R(2, 2) = z * z * (1 - std::cos(angle)) + std::cos(angle);
        return true;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calCartesianRunQueue(my_robot::Robot<_Dofs> *robot)
    {
        double dddposLimit = robot->dddposLimit;
        double ddposLimit = robot->ddposLimit;
        double dposLimit = this->runSpeed * robot->dposLimit;
        double dddoriLimit = robot->dddoriLimit;
        double ddoriLimit = robot->ddoriLimit;
        double doriLimit = this->runSpeed * robot->doriLimit;

        Eigen::Matrix<double, 4, 1> jerkLimit(dddposLimit, dddposLimit, dddposLimit, dddoriLimit);
        Eigen::Matrix<double, 4, 1> accLimit(ddposLimit, ddposLimit, ddposLimit, ddoriLimit);
        Eigen::Matrix<double, 4, 1> velLimit(dposLimit, dposLimit, dposLimit, doriLimit);
        Eigen::Matrix<double, 4, 1> pose = Eigen::Matrix<double, 4, 1>::Zero();
        Eigen::Matrix<double, 4, 1> dpose = Eigen::Matrix<double, 4, 1>::Zero();
        Eigen::Matrix<double, 4, 1> pose_calQueue;
        pose.block(0, 0, 3, 1) = robot->getPosition();
        pose_calQueue.block(0, 0, 3, 1) = x_calQueue.block(0, 0, 3, 1);
        pose[3] = 0;
        this->R1 = robot->getOrientation().toRotationMatrix();
        Eigen::AngleAxisd Rx(Eigen::AngleAxisd(this->x_calQueue(3), Eigen::Vector3d::UnitX()));
        Eigen::AngleAxisd Ry(Eigen::AngleAxisd(this->x_calQueue(4), Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd Rz(Eigen::AngleAxisd(this->x_calQueue(5), Eigen::Vector3d::UnitZ()));
        Eigen::Matrix3d R2;
        R2 = Rx * Ry * Rz;
        Eigen::Matrix3d dR = R1.transpose() * R2;
        if (this->RtoAxisAngle(dR, this->rotationAxis, this->rotationAngle))
            pose_calQueue[3] = this->rotationAngle;
        else
            pose_calQueue[3] = 0;
        if (!this->taskPlanner->calPlanQueue(true, this->cycleTime, velLimit, accLimit, jerkLimit, pose, pose_calQueue, dpose,
                                             this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue))
        {
            this->clearMotionQueue(this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue);
            return;
        }
        if (!this->checkMotionQueue(this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue))
            this->clearMotionQueue(this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue);
        else
            this->startMotion();
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calCartesianStopQueue(my_robot::Robot<_Dofs> *robot)
    {
        double dddposLimit = robot->dddposLimit;
        double ddposLimit = robot->ddposLimit;
        double dddoriLimit = robot->dddoriLimit;
        double ddoriLimit = robot->ddoriLimit;
        Eigen::Matrix<double, 4, 1> maxJerk(dddposLimit, dddposLimit, dddposLimit, dddoriLimit);
        Eigen::Matrix<double, 4, 1> maxAcc(ddposLimit, ddposLimit, ddposLimit, ddoriLimit);
        Eigen::Matrix<double, 4, 1> x = Eigen::Matrix<double, 4, 1>::Zero();
        Eigen::Matrix<double, 4, 1> dx = Eigen::Matrix<double, 4, 1>::Zero();
        Eigen::Matrix<double, 4, 1> ddx = Eigen::Matrix<double, 4, 1>::Zero();
        for (int i = 0; i < 3; i++)
        {
            x[i] = this->x_dRunQueue[i].front();
            dx[i] = this->dx_dRunQueue[i].front();
            ddx[i] = this->ddx_dRunQueue[i].front();
        }
        if (!calStopPlan<4>(true, this->cycleTime, maxJerk, maxAcc, x, dx, ddx, this->x_dStopQueue, this->dx_dStopQueue, this->ddx_dStopQueue))
        {
            this->clearMotionQueue(this->x_dStopQueue, this->dx_dStopQueue, this->ddx_dStopQueue);
            return;
        }

        if (!this->checkMotionQueue(this->x_dStopQueue, this->dx_dStopQueue, this->ddx_dStopQueue))
            this->clearMotionQueue(this->x_dStopQueue, this->dx_dStopQueue, this->ddx_dStopQueue);
        else
        {
            this->clearMotionQueue(this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue);
            this->stopMotion();
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calJointJogMove(my_robot::Robot<_Dofs> *robot)
    {
        Eigen::Matrix<double, _Dofs, 1> dddqLimit = robot->getdddqLimit();
        Eigen::Matrix<double, _Dofs, 1> ddqLimit = robot->getddqLimit();
        Eigen::Matrix<double, _Dofs, 1> dqLimit = this->jogSpeed * robot->getdqLimit();
        if (this->jogMoveFlag)
        {
            calJogMovePlan(jogMoveFlag, cycleTime, this->jogDir, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1], dqLimit[jogNum - 1],
                           this->controllerLaw->q_d[jogNum - 1], this->controllerLaw->dq_d[jogNum - 1], this->controllerLaw->ddq_d[jogNum - 1]);
            this->jogMoveFlag = false;
        }
        else
        {
            calJogMovePlan(jogMoveFlag, cycleTime, this->jogDir, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1], dqLimit[jogNum - 1],
                           this->controllerLaw->q_d[jogNum - 1], this->controllerLaw->dq_d[jogNum - 1], this->controllerLaw->ddq_d[jogNum - 1]);
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calJointJogStop(my_robot::Robot<_Dofs> *robot)
    {
        Eigen::Matrix<double, _Dofs, 1> dddqLimit = robot->getdddqLimit();
        Eigen::Matrix<double, _Dofs, 1> ddqLimit = robot->getddqLimit();
        if (this->jogStopFlag)
        {
            calJogStopPlan(jogStopFlag, cycleTime, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1],
                           this->controllerLaw->q_d[jogNum - 1], this->controllerLaw->dq_d[jogNum - 1], this->controllerLaw->ddq_d[jogNum - 1]);
            this->jogStopFlag = false;
        }
        else
        {
            calJogStopPlan(jogStopFlag, cycleTime, dddqLimit[jogNum - 1], ddqLimit[jogNum - 1],
                           this->controllerLaw->q_d[jogNum - 1], this->controllerLaw->dq_d[jogNum - 1], this->controllerLaw->ddq_d[jogNum - 1]);
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calCartesianJogMove(my_robot::Robot<_Dofs> *robot)
    {
        if (this->jogNum <= 3)
        {
            double dddposLimit = robot->dddposLimit;
            double ddposLimit = robot->ddposLimit;
            double dposLimit = this->jogSpeed * robot->dposLimit;
            if (this->jogMoveFlag)
            {
                calJogMovePlan(this->jogMoveFlag, this->cycleTime, this->jogDir, dddposLimit, ddposLimit, dposLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
                this->jogMoveFlag = false;
            }
            else
            {
                calJogMovePlan(this->jogMoveFlag, this->cycleTime, this->jogDir, dddposLimit, ddposLimit, dposLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
            }
            Eigen::Matrix<double, 4, 4> TO2EE = robot->getT().matrix();
            TO2EE(jogNum - 1, 3) = this->controllerLaw->x_d[jogNum - 1]; // x y z 12 13 14
            std::array<double, 7> IKresult = franka_IK_EE_CC(TO2EE, robot->getq()[6], robot->getqArray());
            if (IKresult[0] == NAN || IKresult[1] == NAN || IKresult[2] == NAN ||
                IKresult[3] == NAN || IKresult[4] == NAN || IKresult[5] == NAN)
            {
                this->HandleiKFailed();
                return;
            }
            this->controllerLaw->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->q_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(IKresult.data());
        }
        else
        {
            double dddoriLimit = robot->dddoriLimit;
            double ddoriLimit = robot->ddoriLimit;
            double doriLimit = this->jogSpeed * robot->doriLimit;
            if (this->jogMoveFlag)
            {
                calJogMovePlan(this->jogMoveFlag, this->cycleTime, this->jogDir, dddoriLimit, ddoriLimit, doriLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
                this->jogMoveFlag = false;
            }
            else
            {
                calJogMovePlan(this->jogMoveFlag, this->cycleTime, this->jogDir, dddoriLimit, ddoriLimit, doriLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
            }
            Eigen::Matrix<double, 4, 4> TO2EE = robot->getT().matrix();
            Eigen::AngleAxisd Rx(Eigen::AngleAxisd(this->controllerLaw->x_d(3), Eigen::Vector3d::UnitX()));
            Eigen::AngleAxisd Ry(Eigen::AngleAxisd(this->controllerLaw->x_d(4), Eigen::Vector3d::UnitY()));
            Eigen::AngleAxisd Rz(Eigen::AngleAxisd(this->controllerLaw->x_d(5), Eigen::Vector3d::UnitZ()));
            Eigen::Matrix3d R;
            R = Rx * Ry * Rz;
            TO2EE.block(0, 0, 3, 3) = R;
            std::array<double, 7> IKresult = franka_IK_EE_CC(TO2EE, robot->getq()[6], robot->getqArray());
            this->controllerLaw->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->q_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(IKresult.data());
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calCartesianJogStop(my_robot::Robot<_Dofs> *robot)
    {
        if (this->jogNum <= 3)
        {
            double dddposLimit = robot->dddposLimit;
            double ddposLimit = robot->ddposLimit;
            if (this->jogStopFlag)
            {
                calJogStopPlan(jogStopFlag, cycleTime, dddposLimit, ddposLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
                this->jogStopFlag = false;
            }
            else
            {
                calJogStopPlan(jogStopFlag, cycleTime, dddposLimit, ddposLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
            }
            Eigen::Matrix<double, 4, 4> TO2EE = robot->getT().matrix();
            TO2EE(jogNum - 1, 3) = this->controllerLaw->x_d[jogNum - 1]; // x y z 12 13 14
            std::array<double, 7> IKresult = franka_IK_EE_CC(TO2EE, robot->getq()[6], robot->getqArray());
            if (IKresult[0] == NAN || IKresult[1] == NAN || IKresult[2] == NAN ||
                IKresult[3] == NAN || IKresult[4] == NAN || IKresult[5] == NAN)
            {
                this->HandleiKFailed();
                return;
            }
            this->controllerLaw->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->q_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(IKresult.data());
        }
        else
        {
            double dddoriLimit = robot->dddoriLimit;
            double ddoriLimit = robot->ddoriLimit;
            if (this->jogStopFlag)
            {
                calJogStopPlan(jogStopFlag, cycleTime, dddoriLimit, ddoriLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
                this->jogStopFlag = false;
            }
            else
            {
                calJogStopPlan(jogStopFlag, cycleTime, dddoriLimit, ddoriLimit, this->controllerLaw->x_d[jogNum - 1],
                               this->controllerLaw->dx_d[jogNum - 1], this->controllerLaw->ddx_d[jogNum - 1]);
            }
            Eigen::Matrix<double, 4, 4> TO2EE = robot->getT().matrix();
            Eigen::AngleAxisd Rx(Eigen::AngleAxisd(this->controllerLaw->x_d(3), Eigen::Vector3d::UnitX()));
            Eigen::AngleAxisd Ry(Eigen::AngleAxisd(this->controllerLaw->x_d(4), Eigen::Vector3d::UnitY()));
            Eigen::AngleAxisd Rz(Eigen::AngleAxisd(this->controllerLaw->x_d(5), Eigen::Vector3d::UnitZ()));
            Eigen::Matrix3d R;
            R = Rx * Ry * Rz;
            TO2EE.block(0, 0, 3, 3) = R;
            std::array<double, 7> IKresult = franka_IK_EE_CC(TO2EE, robot->getq()[6], robot->getqArray());
            this->controllerLaw->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
            this->controllerLaw->q_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(IKresult.data());
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::HandleiKFailed()
    {
        std::cout << "[robotController] 逆解失败，停止运动" << std::endl;
        clearMotionQueue(this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue);
        clearMotionQueue(this->x_dStopQueue, this->dx_dStopQueue, this->ddx_dStopQueue);
        this->nowControllerStatus = RunStatus::wait_;
    }
    //******************************************************************init******************************************************************//
    // init
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::init(int recordPeriod, my_robot::Robot<_Dofs> *robot)
    {
        // 设置数据记录周期
        setRecord(this->recordPeriod);

        // tmp
        for (int i = 0; i < _Dofs; i++)
        {
            controllerParam.jointParam1[i].value = 800;
            controllerParam.jointParam2[i].value = 40;
        }

        // 建立通信 建立数据映射
        if (this->communicationModel.createConnect((key_t)SM_ID, (key_t)MS_ID, this->robotDataBuff,
                                                   this->controllerCommandBUff, this->controllerStateBUff))
            printf("通信模型建立成功\n");

        // init state
        this->initStateToMaster();
        // 丢弃数据
        this->controllerCommandBUff->stopSign = false;
        this->controllerCommandBUff->runSign = false;
        this->controllerCommandBUff->jogSign = false;
        this->controllerCommandBUff->newLimit = false;
        changeControllerLaw(this->controllerLawType, robot);
        changePlanner(this->plannerType);
    }

    //******************************************************************run******************************************************************//
    // run
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::updateTime()
    {
        this->time++;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::controllerParamRenew()
    {
        if (this->controllerLaw != nullptr)
            this->controllerLaw->controllerParamRenew(this->filterParams);
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::communication(my_robot::Robot<_Dofs> *robot)
    {
        bool isConnect = false;
        this->connectStatus = this->communicationModel.comSendMessage(isConnect);

        if (!this->connectStatus)
            return;

        // write
        for (int i = 0; i < _Dofs; i++)
        {
            this->robotDataBuff->q_d[i] = this->controllerLaw->q_d[i] * 180.0 / M_PI;
            this->robotDataBuff->q[i] = robot->getq()[i] * 180.0 / M_PI;
            this->robotDataBuff->dq[i] = robot->getdq()[i] * 180.0 / M_PI;
        }
        // note: matrix->数组
        std::copy(robot->getTorque().data(), robot->getTorque().data() + _Dofs, this->robotDataBuff->tau);
        for (int i = 0; i < 3; i++)
        {
            this->robotDataBuff->position[i] = robot->getPosition()[i];
            this->robotDataBuff->orientation[i] = robot->getOrientation().toRotationMatrix().eulerAngles(0, 1, 2)[i] * 180.0 / M_PI;
        }
        this->controllerStateBUff->controllerStatus = this->nowControllerStatus;

        // read
        this->jogSpeed_d = (double)this->controllerCommandBUff->jogSpeed_d / 100.0;
        this->runSpeed_d = (double)this->controllerCommandBUff->runSpeed_d / 100.0;
        this->controllerLawType_d = this->controllerCommandBUff->controllerLawType_d;
        this->plannerType_d = this->controllerCommandBUff->plannerType_d;
        if (this->controllerCommandBUff->newLimit) // 新的限位设置
        {                                          // note: 数组->matrix
            for (int i = 0; i < _Dofs; i++)
            {
                robot->qMax[i] = this->controllerCommandBUff->qMax[i] * M_PI / 180.0;
                robot->qMin[i] = this->controllerCommandBUff->qMin[i] * M_PI / 180.0;
                robot->dqLimit[i] = this->controllerCommandBUff->dqLimit[i] * M_PI / 180.0;
                robot->ddqLimit[i] = this->controllerCommandBUff->ddqLimit[i] * M_PI / 180.0;
                robot->dddqLimit[i] = this->controllerCommandBUff->dddqLimit[i] * M_PI / 180.0;
            }
            this->controllerCommandBUff->newLimit = false;
        }
        if (this->controllerCommandBUff->runSign) // 新的规划任务
        {
            for (int i = 0; i < _Dofs; i++)
            {
                this->q_calQueue[i] = this->controllerCommandBUff->q_final[i] * M_PI / 180.0;
            }
            for (int i = 0; i < 3; i++)
            {
                this->x_calQueue[i] = this->controllerCommandBUff->x_final[i];
            }
            for (int i = 3; i < 6; i++)
            {
                this->x_calQueue[i] = this->controllerCommandBUff->x_final[i] * M_PI / 180.0;
            }
            this->controllerCommandBUff->runSign = false;
            this->newPlan = true;
        }
        if (this->controllerCommandBUff->stopSign) // 新的急停任务
        {
            this->controllerCommandBUff->stopSign = false;
            this->newStop = true;
        }
        // 直接赋值
        this->runTaskSpace = this->controllerCommandBUff->plannerTaskSpace;
        this->jogSign = this->controllerCommandBUff->jogSign;
        this->jogNum = this->controllerCommandBUff->jogNum;
        this->jogDir = this->controllerCommandBUff->jogDir;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::updateStatus(my_robot::Robot<_Dofs> *robot)
    {
        if (!this->connectStatus)
            return;

        if (this->controllerLawType != this->controllerLawType_d) // 切换控制器
        {
            changeControllerLaw(this->controllerLawType_d, robot);
            this->controllerLawType = this->controllerLawType_d;
        }
        if (this->plannerType != this->plannerType_d) // 切换规划器
        {
            changePlanner(this->plannerType_d);
            this->plannerType = this->plannerType_d;
        }
        if (this->jogSpeed != this->jogSpeed_d) // 更改点动速度
        {
            this->jogSpeed_d = std::max(0.01, std::min(1.0, this->jogSpeed_d)); // 1%和100%
            // recal
            this->jogSpeed = this->jogSpeed_d;
        }
        if (this->runSpeed != this->runSpeed_d) // 更改运行速度
        {
            this->runSpeed_d = std::max(0.01, std::min(1.0, this->runSpeed_d));
            // recal
            this->runSpeed = this->runSpeed_d;
        }
        // run stop
        if (this->newPlan && this->nowControllerStatus == RunStatus::wait_) // 新的规划
        {
            if (this->runTaskSpace == TaskSpace::jointSpace)
                this->calJointRunQueue(robot);
            else
                this->calCartesianRunQueue(robot);
        }
        if (this->newStop && this->nowControllerStatus == RunStatus::run_) // 新的急停规划
        {
            if (this->runTaskSpace == TaskSpace::jointSpace)
                this->calJointStopQueue(robot);
            else
                this->calCartesianStopQueue(robot);
        }
        this->newPlan = false;
        this->newStop = false;
        // jog jogStop
        if (this->jogSign && (this->nowControllerStatus == RunStatus::wait_)) // 点动
        {
            this->jogMoveFlag = true;
            this->nowControllerStatus = RunStatus::jog_;
        }
        if (!this->jogSign && (this->nowControllerStatus == RunStatus::jog_)) // 点动停止
        {
            this->jogStopFlag = true;
            this->nowControllerStatus = RunStatus::jogStop_;
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calDesireNext(my_robot::Robot<_Dofs> *robot)
    {
        auto updateX_d = [this]()
        {
            Eigen::Matrix4d T = franka_FK(this->controllerLaw->q_d);
            this->controllerLaw->x_d.head(3) = T.block<3, 1>(0, 3);
            Eigen::Matrix3d R = T.block<3, 3>(0, 0);
            Eigen::Vector3d eulerAngles = R.eulerAngles(0, 1, 2);
            this->controllerLaw->x_d.tail(3) = eulerAngles;
        };
        switch (this->nowControllerStatus)
        {
        case RunStatus::wait_:
            this->controllerLaw->calWaitDesireNext();
            break;
        case RunStatus::run_:
            if ((q_dRunQueue[0].empty() && this->runTaskSpace == TaskSpace::jointSpace) ||
                (x_dRunQueue[0].empty() && this->runTaskSpace == TaskSpace::cartesianSpace))
            {
                if (this->runTaskSpace == TaskSpace::jointSpace)
                    updateX_d();
                this->nowControllerStatus = RunStatus::wait_;
                break;
            }
            if (this->runTaskSpace == TaskSpace::cartesianSpace)
            {
                Eigen::Matrix<double, 4, 4> TO2EE;
                Eigen::Matrix3d dR;
                this->AxisAngletoR(dR, this->rotationAxis, this->x_dRunQueue[3].front());
                this->x_dRunQueue[3].pop();
                TO2EE.block(0, 0, 3, 3) = this->R1 * dR;
                for (size_t i = 0; i < 4; i++)
                {
                    this->dx_dRunQueue[i].pop();
                    this->ddx_dRunQueue[i].pop();
                    if (i == 3)
                        break;
                    TO2EE(i, 3) = this->x_dRunQueue[i].front(); // x y z 12 13 14
                    this->x_dRunQueue[i].pop();
                }
                std::array<double, 7> IKresult = franka_IK_EE_CC(TO2EE, robot->getq()[6], robot->getqArray());
                if (IKresult[0] == NAN || IKresult[1] == NAN || IKresult[2] == NAN ||
                    IKresult[3] == NAN || IKresult[4] == NAN || IKresult[5] == NAN)
                {
                    this->HandleiKFailed();
                    return;
                }
                this->controllerLaw->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
                this->controllerLaw->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
                this->controllerLaw->q_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(IKresult.data());
            }
            else
            {
                this->controllerLaw->calRunStopDesireNext(this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue);
            }
            break;
        case RunStatus::stop_:
            if ((q_dStopQueue[0].empty() && this->runTaskSpace == TaskSpace::jointSpace) ||
                (x_dStopQueue[0].empty() && this->runTaskSpace == TaskSpace::cartesianSpace))
            {
                if (this->runTaskSpace == TaskSpace::jointSpace)
                    updateX_d();
                this->nowControllerStatus = RunStatus::wait_;
                break;
            }
            if (this->runTaskSpace == TaskSpace::cartesianSpace)
            {
                Eigen::Matrix<double, 4, 4> TO2EE;
                Eigen::Matrix3d dR;
                this->AxisAngletoR(dR, this->rotationAxis, this->x_dStopQueue[3].front());
                this->x_dStopQueue[3].pop();
                TO2EE.block(0, 0, 3, 3) = this->R1 * dR;
                for (size_t i = 0; i < 4; i++)
                {
                    this->dx_dStopQueue[i].pop();
                    this->ddx_dStopQueue[i].pop();
                    if (i == 3)
                        break;
                    TO2EE(i, 3) = this->x_dStopQueue[i].front(); // x y z 12 13 14
                    this->x_dStopQueue[i].pop();
                }
                std::array<double, 7> IKresult = franka_IK_EE_CC(TO2EE, robot->getq()[6], robot->getqArray());
                if (IKresult[0] == NAN || IKresult[1] == NAN || IKresult[2] == NAN ||
                    IKresult[3] == NAN || IKresult[4] == NAN || IKresult[5] == NAN)
                {
                    this->HandleiKFailed();
                    return;
                }
                this->controllerLaw->ddq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
                this->controllerLaw->dq_d = Eigen::Matrix<double, _Dofs, 1>::Zero();
                this->controllerLaw->q_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(IKresult.data());
            }
            else
            {
                this->controllerLaw->calRunStopDesireNext(this->q_dStopQueue, this->dq_dStopQueue, this->ddq_dStopQueue);
            }
            break;
        case RunStatus::jog_:
            if (this->runTaskSpace == TaskSpace::jointSpace)
                this->calJointJogMove(robot);
            else
                this->calCartesianJogMove(robot);
            break;
        case RunStatus::jogStop_:
            if ((this->controllerLaw->dq_d[jogNum - 1] == 0 && this->runTaskSpace == TaskSpace::jointSpace) ||
                (this->controllerLaw->dx_d[jogNum - 1] == 0 && this->runTaskSpace == TaskSpace::cartesianSpace))
            {
                if (this->runTaskSpace == TaskSpace::jointSpace)
                    updateX_d();
                this->nowControllerStatus = RunStatus::wait_;
                break;
            }
            if (this->runTaskSpace == TaskSpace::jointSpace)
                this->calJointJogStop(robot);
            else
                this->calCartesianJogStop(robot);
            break;
        default:
            printf("calDesireNext error\n");
            break;
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calError(my_robot::Robot<_Dofs> *robot)
    {
        if (controllerLaw.get() != nullptr)
            this->controllerLaw->calError(robot);
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::setU(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d)
    {
        if (controllerLaw.get() != nullptr)
            controllerLaw->setU(robot, tau_d);
    }

    //******************************************************************ROS******************************************************************//
    // ROS
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::recordData(my_robot::Robot<_Dofs> *robot)
    {
        if (controllerLaw.get() == nullptr)
            return;
        if (recordPeriod == 0)
            return;
        if (time == 1)
        {
            this->myfile.open("/home/wd/log/franka/main/pandaController.txt"); // todo
            this->myfile << "pandaController"
                         << "\n";
            this->myfile << this->controllerLaw->controllerLawName << "\n";
            this->myfile << "程序编译日期:" << __DATE__ << "\n";
            this->myfile << "程序编译时刻:" << __TIME__ << std::endl;
        }
        if (time % recordPeriod != 0)
            return;
        // static const char *n = "\n";
        // std::array<double, 7> tmp;
        // tmp = franka_IK_EE_CC(robot->getT02EEArray(), robot->getq()[6], robot->getqArray());
        // this->myfile << "time: " << this->time << "_" << n;
        // this->myfile << "q: " << robot->getq().transpose() << "\n";
        // this->myfile << "tmp: " << tmp[0] << "  " << tmp[1] << "  " << tmp[2] << "  " << tmp[3]
        //              << " " << tmp[4] << "  " << tmp[5] << "  " << tmp[6] << "\n";
        // this->myfile << "q0: " << robot->getq0().transpose() << "\n";

        // this->myfile << "dq: " << robot->getdq().transpose() << "\n";
        // this->myfile << "q_d: " << this->controllerLaw->q_d.transpose() << "\n";
        // this->myfile << "dq_d: " << this->controllerLaw->dq_d.transpose() << "\n";
        // this->myfile << "ddq_d: " << this->controllerLaw->ddq_d.transpose() << "\n";
        // this->myfile << "Position0: " << robot->getPosition0().transpose() << "\n";
        // this->myfile << "Orientation0: " << robot->getOrientation0().toRotationMatrix().eulerAngles(0, 1, 2).transpose() << "\n";
        // this->myfile << "this->controllerLaw->x_d: " << this->controllerLaw->x_d.transpose() << "\n";
        // this->myfile << "Position: " << robot->getPosition().transpose() << "\n";
        // this->myfile << "Orientation: " << robot->getOrientation().toRotationMatrix().eulerAngles(0, 1, 2).transpose() << "\n";
        // this->myfile << "orientation_hold: " << this->orientation_hold.toRotationMatrix().eulerAngles(0, 1, 2).transpose() << "\n";
        // this->myfile << "Position: " << robot->getdPosition().transpose() << "\n";
        // this->myfile << "Orientation: " << robot->getdOrientation().toRotationMatrix().eulerAngles(0, 1, 2).transpose() << "\n";
        // this->myfile << "T:" << n;
        // this->myfile << robot->getT().matrix() << "\n";
        // this->myfile << "M:" << n;
        // this->myfile << robot->getM() << "\n";
        // this->myfile << "ExternM:" << n;
        // this->myfile << robot->getExternM() << "\n";
        // this->myfile << "C: " << n;
        // this->myfile << robot->getC() * robot->getdq() << "\n";
        // this->myfile << "Externc: " << n;
        // this->myfile << robot->getExternc() << "\n";
        // this->myfile << "G: " << n;
        // this->myfile << robot->getG() << n;
        // this->myfile << "ExternG: " << n;
        // this->myfile << robot->getExternG() << n;
        // this->myfile << "J: " << n;
        // this->myfile << robot->getJ() << "\n";
        // this->myfile << "ExternJ: " << n;
        // this->myfile << robot->getExternJ() << "\n";
        // this->myfile << "getdJ: " << n;
        // this->myfile << robot->getdJ() << "\n";
        // this->myfile << "getTorque: " << robot->getTorque().transpose() << "\n";
        // this->myfile << "tau_d: " << this->controllerLaw->tau_d.transpose() << "\n";
        // this->myfile << "-------------------" << std::endl;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::pubData(pubDataType &param_debug, my_robot::Robot<_Dofs> *robot)
    {
        static bool tmp = false;
        if (controllerLaw.get() == nullptr)
            return;
        for (int i = 0; i < _Dofs; i++)
        {
            param_debug.jointError[i] = this->controllerLaw->jointError[i] * 180.0 / M_PI;
            param_debug.q[i] = robot->getq()[i] * 180.0 / M_PI;
            param_debug.dq[i] = robot->getdq()[i] * 180.0 / M_PI;

            param_debug.q_d[i] = this->controllerLaw->q_d[i] * 180.0 / M_PI;
            param_debug.dq_d[i] = this->controllerLaw->dq_d[i] * 180.0 / M_PI;
            param_debug.ddq_d[i] = this->controllerLaw->ddq_d[i] * 180.0 / M_PI;

            param_debug.tau_d[i] = this->controllerLaw->tau_d[i];
        }
        for (int i = 0; i < 3; i++)
        {
            param_debug.position[i] = robot->getPosition()[i];
            param_debug.dposition[i] = robot->getdPosition()[i];
            param_debug.position_d[i] = this->controllerLaw->x_d[i];

            param_debug.orientation[i] = robot->getOrientation().toRotationMatrix().eulerAngles(0, 1, 2)[i];
            param_debug.dorientation[i] = robot->getdOrientation().toRotationMatrix().eulerAngles(0, 1, 2)[i];
            param_debug.orientation_d[i] = this->controllerLaw->x_d[i + 3];
        }

        // for (int i = 0; i < 6; i++)
        // {
        //     param_debug.cartesianError[i] = this->controllerLaw->cartesianError[i];
        // }
    }
};
