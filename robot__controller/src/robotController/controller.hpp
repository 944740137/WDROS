// #pragma once
#include "controllerLaw/controllerLaw.hpp"
#include "planner/planner.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "math.h"

namespace robot_controller
{
    // _Dofs自由度机器人控制器
    template <int _Dofs, typename pubDataType>
    class Controller
    {

    public:
        // debug，绘图
        int recordPeriod = 1;                   // 数据记录周期
        unsigned int time = 0;                  // 当前时刻
        std::ofstream myfile;                   // 记录文件io对象
        double filterParams = 0.005;            // 调参滤波参数
        const double cycleTime = 0.001;         // 运行周期0.001s
        Eigen::Matrix<double, _Dofs, 1> q_hold; // 维持位置
        Eigen::Vector3d position_hold;
        Eigen::Quaterniond orientation_hold;

        // 点动参数
        bool jogSign = false;
        int jogNum = 0;      // 点动轴
        int jogDir = 0;      // 点动轴
        double jogSpeed = 0; // 0-1
        double jogSpeed_d = 0;

        bool jogMoveFlag = false;
        bool jogStopFlag = false;

        // 运行参数
        bool newPlan = false;
        Eigen::Matrix<double, _Dofs, 1> q_calQueue;
        double runSpeed = 0;
        double runSpeed_d = 0;

        // 坐标系
        TaskSpace runTaskSpace = TaskSpace::jointSpace;
        TaskSpace runTaskSpace_d = TaskSpace::jointSpace;

        // 急停参数
        bool newStop = false;

        // 当前运行状态
        RunStatus nowControllerStatus = RunStatus::wait_; // 当前状态

        // 运行队列
        std::vector<std::queue<double>> q_dRunQueue{_Dofs};
        std::vector<std::queue<double>> dq_dRunQueue{_Dofs};
        std::vector<std::queue<double>> ddq_dRunQueue{_Dofs};
        std::vector<std::queue<double>> x_dRunQueue{6};
        std::vector<std::queue<double>> dx_dRunQueue{6};
        std::vector<std::queue<double>> ddx_dRunQueue{6};

        // 急停队列
        std::vector<std::queue<double>> q_dStopQueue{_Dofs};
        std::vector<std::queue<double>> dq_dStopQueue{_Dofs};
        std::vector<std::queue<double>> ddq_dStopQueue{_Dofs};
        std::vector<std::queue<double>> x_dStopQueue{6};
        std::vector<std::queue<double>> dx_dStopQueue{6};
        std::vector<std::queue<double>> ddx_dStopQueue{6};

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
        std::unique_ptr<Planner<6>> taskPlanner;
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
        void changeControllerLaw(ControllerLawType type);
        void changePlanner(PlannerType type);
        void changeTaskSpace();
        void initStateToMaster();
        void calRunQueue(my_robot::Robot<_Dofs> *robot);
        void calStopQueue(my_robot::Robot<_Dofs> *robot);
        void startMotion();
        void stopMotion();
        void clearMotionQueue(std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);
        bool checkMotionQueue(std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);

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
    Controller<_Dofs, pubDataType>::Controller() : q_calQueue(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                                                   q_hold(Eigen::Matrix<double, _Dofs, 1>::Zero())
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
    void Controller<_Dofs, pubDataType>::changeControllerLaw(ControllerLawType type)
    {
        if (!newControllerLaw(controllerLaw, type, runTaskSpace))
            printf("ControllerLaw create Error\n");
        dynamicSetParameter();
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::changePlanner(PlannerType type)
    {
        if (!newPlanner(this->jointPlanner, type))
            printf("Planner create Error\n");
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
        for (int i = 0; i < _Dofs; i++)
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
        for (int i = 0; i < _Dofs; i++)
        {
            size = size + q_d[i].size();
        }
        if ((size / _Dofs) != q_d[0].size())
        {
            this->clearMotionQueue(q_d, dq_d, ddq_d);
            std::cout << "[robotController] 队列长度不一致 清空处理" << std::endl;
            return false;
        }
        return true;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calRunQueue(my_robot::Robot<_Dofs> *robot)
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
    void Controller<_Dofs, pubDataType>::calStopQueue(my_robot::Robot<_Dofs> *robot)
    {
        Eigen::Matrix<double, _Dofs, 1> maxAcc = robot->getddqLimit();
        Eigen::Matrix<double, _Dofs, 1> maxJerk = robot->getdddqLimit();
        Eigen::Matrix<double, _Dofs, 1> q = Eigen::Matrix<double, _Dofs, 1>::Zero();
        Eigen::Matrix<double, _Dofs, 1> dq = Eigen::Matrix<double, _Dofs, 1>::Zero();
        Eigen::Matrix<double, _Dofs, 1> ddq = Eigen::Matrix<double, _Dofs, 1>::Zero();
        for (int i = 0; i < _Dofs; i++)
        {
            q[i] = this->q_dRunQueue[i].front();
            dq[i] = this->dq_dRunQueue[i].front();
            ddq[i] = this->ddq_dRunQueue[i].front();
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
    void Controller<_Dofs, pubDataType>::changeTaskSpace()
    {
        if (this->controllerLaw.get() != nullptr)
            this->controllerLaw->setTaskSpace(this->runTaskSpace);
    }

    //******************************************************************init******************************************************************//
    // init
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::init(int recordPeriod, my_robot::Robot<_Dofs> *robot)
    {
        // 设置数据记录周期
        setRecord(this->recordPeriod);
        this->q_hold = robot->getq();

        // tmp
        for (int i = 0; i < _Dofs; i++)
        {
            controllerParam.jointParam1[i].value = 60;
            controllerParam.jointParam2[i].value = 5;
        }
        for (int i = 0; i < 6; i++)
        {
            controllerParam.cartesianParam1[i].value = 40;
            controllerParam.cartesianParam2[i].value = 5;
        }

        // 建立通信 建立数据映射
        if (this->communicationModel.createConnect((key_t)SM_ID, (key_t)MS_ID, this->robotDataBuff,
                                                   this->controllerCommandBUff, this->controllerStateBUff))
        {
            printf("通信模型建立成功\n");
        }

        // init state
        this->initStateToMaster();
        // 丢弃数据
        this->controllerCommandBUff->stopSign = false;
        this->controllerCommandBUff->runSign = false;
        this->controllerCommandBUff->newLimit = false;
        changeControllerLaw(this->controllerLawType);
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
            this->robotDataBuff->orientation[i] = robot->getOrientation().toRotationMatrix().eulerAngles(2, 1, 0)[i] * 180.0 / M_PI;
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
            this->runTaskSpace = this->controllerCommandBUff->plannerTaskSpace;
            for (int i = 0; i < _Dofs; i++)
            {
                this->q_calQueue[i] = this->controllerCommandBUff->q_final[i] * M_PI / 180.0;
            }
            this->controllerCommandBUff->runSign = false;
            this->newPlan = true;
        }
        if (this->controllerCommandBUff->stopSign) // 新的急停任务
        {
            this->controllerCommandBUff->stopSign = false;
            this->newStop = true;
        }
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
            changeControllerLaw(this->controllerLawType_d);
            this->controllerLawType = this->controllerLawType_d;
        }
        if (this->plannerType != this->plannerType_d) // 切换规划器
        {
            changePlanner(this->plannerType_d);
            this->plannerType = this->plannerType_d;
        }
        if (this->runTaskSpace != this->runTaskSpace_d) // 切换坐标系：关节/笛卡尔
        {
            if (this->controllerLaw.get() != nullptr)
            {
                this->controllerLaw->setTaskSpace(this->runTaskSpace_d);
                this->runTaskSpace = this->runTaskSpace_d;
            }
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
            this->calRunQueue(robot);
        }
        if (this->newStop && this->nowControllerStatus == RunStatus::run_) // 新的急停规划
        {
            this->calStopQueue(robot);
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
        switch (this->nowControllerStatus)
        {
        case RunStatus::wait_:
            this->controllerLaw->calWaitDesireNext(this->q_hold, this->position_hold, this->orientation_hold);
            break;
        case RunStatus::run_:
            if ((q_dRunQueue[0].empty() && runTaskSpace == TaskSpace::jointSpace) ||
                (x_dRunQueue[0].empty() && runTaskSpace == TaskSpace::cartesianSpace))
            {
                this->q_hold = robot->getq();
                this->position_hold = robot->getPosition();
                this->orientation_hold = robot->getOrientation();
                this->nowControllerStatus = RunStatus::wait_;
                break;
            }
            this->controllerLaw->calRunStopDesireNext(this->q_dRunQueue, this->dq_dRunQueue, this->ddq_dRunQueue, this->x_dRunQueue, this->dx_dRunQueue, this->ddx_dRunQueue);
            break;
        case RunStatus::stop_:
            if ((q_dStopQueue[0].empty() && runTaskSpace == TaskSpace::jointSpace) ||
                (x_dStopQueue[0].empty() && runTaskSpace == TaskSpace::cartesianSpace))
            {
                this->q_hold = robot->getq();
                this->position_hold = robot->getPosition();
                this->orientation_hold = robot->getOrientation();
                this->nowControllerStatus = RunStatus::wait_;
                break;
            }
            this->controllerLaw->calRunStopDesireNext(this->q_dStopQueue, this->dq_dStopQueue, this->ddq_dStopQueue, this->x_dStopQueue, this->dx_dStopQueue, this->ddx_dStopQueue);
            break;
        case RunStatus::jog_:
            this->controllerLaw->calJogMove(robot, this->jogMoveFlag, this->jogSpeed, this->cycleTime, this->jogDir, this->jogNum);
            break;
        case RunStatus::jogStop_:
            if ((this->controllerLaw->dq_d[jogNum - 1] == 0 && runTaskSpace == TaskSpace::jointSpace) /* ||
                (this->controllerLaw->dposition_d[jogNum - 1] == 0 && runTaskSpace == TaskSpace::cartesianSpace) */
            )
            {
                this->q_hold = robot->getq();
                this->position_hold = robot->getPosition();
                this->orientation_hold = robot->getOrientation();
                this->nowControllerStatus = RunStatus::wait_;
                break;
            }
            this->controllerLaw->calJogStop(robot, this->jogStopFlag, this->cycleTime, this->jogNum);
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
            this->controllerLaw->calError(robot, this->cycleTime);
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

        static const char *n = "\n";
        // this->myfile << "time: " << this->time << "_" << n;
        // this->myfile << "q0: " << robot->getq0().transpose() << "\n";
        // this->myfile << "q: " << robot->getq().transpose() << "\n";
        // this->myfile << "dq: " << robot->getdq().transpose() << "\n";
        // this->myfile << "q_d: " << this->controllerLaw->q_d.transpose() << "\n";
        // this->myfile << "dq_d: " << this->controllerLaw->dq_d.transpose() << "\n";
        // this->myfile << "ddq_d: " << this->controllerLaw->ddq_d.transpose() << "\n";
        // this->myfile << "Position0: " << robot->getPosition0().transpose() << "\n";
        // this->myfile << "Orientation0: " << robot->getOrientation0().toRotationMatrix().eulerAngles(2, 1, 0).transpose() << "\n";
        // this->myfile << "Position: " << robot->getPosition().transpose() << "\n";
        // this->myfile << "Orientation: " << robot->getOrientation().toRotationMatrix().eulerAngles(2, 1, 0).transpose() << "\n";
        // this->myfile << "Position: " << robot->getdPosition().transpose() << "\n";
        // this->myfile << "Orientation: " << robot->getdOrientation().toRotationMatrix().eulerAngles(2, 1, 0).transpose() << "\n";
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
        for (int i = 0; i < 6; i++)
        {
            param_debug.cartesianError[i] = this->controllerLaw->cartesianError[i];
        }
        for (int i = 0; i < 3; i++)
        {
            param_debug.position[i] = robot->getPosition()[i];
            param_debug.position_d[i] = this->controllerLaw->position_d[i];
            param_debug.orientation[i] = robot->getOrientation().toRotationMatrix().eulerAngles(2, 1, 0)[i] * 180.0 / M_PI;
            param_debug.orientation_d[i] = this->controllerLaw->orientation_d.toRotationMatrix().eulerAngles(2, 1, 0)[i] * 180.0 / M_PI;
        }
    }
};
