// #pragma once
#include "planner/planner.hpp"
#include "controllerLaw/controllerLaw.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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

        // 点动参数
        unsigned int jogSign = 0;
        double jogSpeed = 0;
        double jogSpeed_d = 0;

        // 规划参数
        unsigned int nowPlanTaskNum = 0; // 当前规划任务号
        bool newPlan = false;
        TaskSpace plannerTaskSpace = TaskSpace::jointSpace;
        Eigen::Matrix<double, _Dofs, 1> q_calQueue;
        double runSpeed = 0;
        double runSpeed_d = 0;

        // 急停参数
        unsigned int nowStopTaskNum = 0; // 当前急停任务号
        bool newStop = false;

        // 当前运行状态
        RunStatus nowControllerStatus = RunStatus::wait_; // 当前状态

        // 运行队列
        std::vector<std::queue<double>> q_dQueue{_Dofs};
        std::vector<std::queue<double>> dq_dQueue{_Dofs};
        std::vector<std::queue<double>> ddq_dQueue{_Dofs};

        // 急停队列
        std::vector<std::queue<double>> q_stopQueue{_Dofs};
        std::vector<std::queue<double>> dq_stopQueue{_Dofs};
        std::vector<std::queue<double>> ddq_stopQueue{_Dofs};

        // 通讯
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
        void initParamFromMaster(my_robot::Robot<_Dofs> *robot);
        void initParamToMaster();
        void calRunQueue(my_robot::Robot<_Dofs> *robot);
        void calStopQueue(my_robot::Robot<_Dofs> *robot);
        void startMotion();
        void clearRunQueue();

        // init
        void init(int recordPeriod, my_robot::Robot<_Dofs> *robot);

        // run
        void updateTime();
        void controllerParamRenew();
        void communication(my_robot::Robot<_Dofs> *robot);
        void updateStatus(my_robot::Robot<_Dofs> *robot);
        void calDesireNext(my_robot::Robot<_Dofs> *robot); // 计算下一个周期的期望
        void calError(my_robot::Robot<_Dofs> *robot);
        void setControllerLaw(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d);

        //
        void recordData(my_robot::Robot<_Dofs> *robot);
        void pubData(pubDataType &param_debug, my_robot::Robot<_Dofs> *robot);
    };

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
        if (!newControllerLaw(controllerLaw, type, plannerTaskSpace))
            printf("ControllerLaw create Error\n");
        dynamicSetParameter();
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::dynamicSetParameter()
    {
        if (controllerLaw.get() != nullptr)
            controllerLaw->dynamicSetParameter(this->controllerParam, this->time);
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::initParamToMaster()
    {
        this->controllerStateBUff->robotDof = _Dofs;
        this->controllerStateBUff->controllerStatus = this->nowControllerStatus;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::initParamFromMaster(my_robot::Robot<_Dofs> *robot)
    {
        for (int i = 0; i < _Dofs; i++)
        {
            robot->qMax[i] = this->controllerCommandBUff->qMax[i];
            robot->qMin[i] = this->controllerCommandBUff->qMin[i];
            robot->dqLimit[i] = this->controllerCommandBUff->dqLimit[i];
            robot->ddqLimit[i] = this->controllerCommandBUff->ddqLimit[i];
            robot->dddqLimit[i] = this->controllerCommandBUff->dddqLimit[i];
            printf("%f %f %f %f %f\n", robot->qMax[i], robot->qMin[i], robot->dqLimit[i], robot->ddqLimit[i], robot->dddqLimit[i]);
        }
        this->nowStopTaskNum = this->controllerCommandBUff->stopTaskNum;
        this->nowPlanTaskNum = this->controllerCommandBUff->planTaskNum;
        this->jogSpeed = this->controllerCommandBUff->jogSpeed_d;
        this->jogSpeed_d = this->controllerCommandBUff->jogSpeed_d;
        this->runSpeed = this->controllerCommandBUff->runSpeed_d;
        this->runSpeed_d = this->controllerCommandBUff->runSpeed_d;
        this->controllerLawType = this->controllerCommandBUff->controllerLawType_d;
        this->controllerLawType_d = this->controllerCommandBUff->controllerLawType_d;
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::startMotion()
    {
        this->nowControllerStatus = RunStatus::run_; // 开始运动
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calRunQueue(my_robot::Robot<_Dofs> *robot)
    {
        double velLimit[_Dofs] = {0.0};
        double accLimit[_Dofs] = {0.0};
        for (int i = 0; i < _Dofs; i++)
        {
            velLimit[i] = this->runSpeed * robot->getdqLimit()[i];
            accLimit[i] = robot->getddqLimit()[i];
        }
        calQuinticPlan(true, this->cycleTime, velLimit, accLimit, robot->getq(), this->q_calQueue,
                       this->q_dQueue, this->dq_dQueue, this->ddq_dQueue);
        this->startMotion();
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calStopQueue(my_robot::Robot<_Dofs> *robot)
    {
        double maxJerk[_Dofs] = {0.0};
        double maxAcc[_Dofs] = {0.0};
        double q[_Dofs] = {0.0};
        double dq[_Dofs] = {0.0};
        double ddq[_Dofs] = {0.0};
        for (int i = 0; i < _Dofs; i++)
        {
            std::queue<double> empty[3];
            maxJerk[i] = robot->getdddqLimit()[i];
            maxAcc[i] = robot->getddqLimit()[i];
            q[i] = this->q_dQueue[i].front();
            dq[i] = this->dq_dQueue[i].front();
            ddq[i] = this->ddq_dQueue[i].front();
            swap(empty[0], q_dQueue[i]);
            swap(empty[1], dq_dQueue[i]);
            swap(empty[2], ddq_dQueue[i]);
        }
        calStopPlan<_Dofs>(true, this->cycleTime, maxJerk, maxAcc, q, dq, ddq, q_dQueue, dq_dQueue, ddq_dQueue);
    }
    
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
            controllerParam.jointParam1[i].value = 50;
            controllerParam.jointParam2[i].value = 5;
        }
        for (int i = 0; i < 6; i++)
        {
            controllerParam.cartesianParam1[i].value = 50;
            controllerParam.cartesianParam2[i].value = 5;
        }

        // 建立通信 建立数据映射
        if (this->communicationModel.createConnect((key_t)SM_ID, (key_t)MS_ID, this->robotDataBuff,
                                                   this->controllerCommandBUff, this->controllerStateBUff))
        {
            printf("通信模型建立成功\n");
        }

        // read
        this->initParamToMaster();
        this->initParamFromMaster(robot);
        changeControllerLaw(this->controllerLawType);
    }

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
        if (!this->connectStatus && this->communicationModel.comSendMessage())
        {
            this->initParamFromMaster(robot);
        }
        this->connectStatus = this->communicationModel.comSendMessage();

        if (!this->connectStatus)
            return;
        std::copy(this->controllerLaw->q_d.data(), this->controllerLaw->q_d.data() + _Dofs, this->robotDataBuff->q_d);
        std::copy(robot->getq().data(), robot->getq().data() + _Dofs, this->robotDataBuff->q);
        std::copy(robot->getdq().data(), robot->getdq().data() + _Dofs, this->robotDataBuff->dq);
        std::copy(robot->getTorque().data(), robot->getTorque().data() + _Dofs, this->robotDataBuff->tau);

        for (int i = 0; i < 3; i++)
        {
            this->robotDataBuff->position[i] = robot->getdOrientation().toRotationMatrix().eulerAngles(2, 1, 0)[i];
            this->robotDataBuff->orientation[i] = robot->getdPosition()[i];
        }
        this->controllerStateBUff->controllerStatus = this->nowControllerStatus;

        // read
        this->jogSpeed_d = this->controllerCommandBUff->jogSpeed_d;
        this->runSpeed_d = this->controllerCommandBUff->runSpeed_d;
        this->controllerLawType_d = this->controllerCommandBUff->controllerLawType_d;

        this->jogSign = this->controllerCommandBUff->jogSign;

        if (this->nowPlanTaskNum != this->controllerCommandBUff->planTaskNum) // 新的规划任务
        {
            this->plannerTaskSpace = this->controllerCommandBUff->plannerTaskSpace;
            this->q_calQueue = Eigen::Map<Eigen::MatrixXd>(this->controllerCommandBUff->q_final, _Dofs, 1);
            this->newPlan = true;
            this->nowPlanTaskNum = this->controllerCommandBUff->planTaskNum;
        }
        if (this->nowStopTaskNum != this->controllerCommandBUff->stopTaskNum) // 新的急停任务
        {
            this->newStop = true;
            this->nowStopTaskNum = this->controllerCommandBUff->stopTaskNum;
        }
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
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calDesireNext(my_robot::Robot<_Dofs> *robot)
    {
        for (int i = 0; i < _Dofs; i++)
        {
            switch (this->nowControllerStatus)
            {
            case RunStatus::wait_:
                this->controllerLaw->q_d[i] = this->q_hold[i];
                this->controllerLaw->dq_d[i] = 0;
                this->controllerLaw->ddq_d[i] = 0;
                break;
            case RunStatus::run_:
                this->controllerLaw->q_d[i] = q_dQueue[i].front();
                this->controllerLaw->dq_d[i] = dq_dQueue[i].front();
                this->controllerLaw->ddq_d[i] = ddq_dQueue[i].front();
                q_dQueue[i].pop();
                dq_dQueue[i].pop();
                ddq_dQueue[i].pop();
                if (q_dQueue[i].empty())
                {
                    this->q_hold = this->controllerLaw->q_d;
                    this->nowControllerStatus = RunStatus::wait_;
                }
                break;
            default:
                printf("calDesireNext error\n");
                break;
            }
        }
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::calError(my_robot::Robot<_Dofs> *robot)
    {

        if (controllerLaw.get() == nullptr)
            return;
        // 关节误差与误差导数
        this->controllerLaw->jointError = this->controllerLaw->q_d - robot->getq();
        this->controllerLaw->djointError = this->controllerLaw->dq_d - robot->getdq();

        // 笛卡尔位姿误差
        Eigen::Quaterniond orientation = robot->getOrientation();
        this->controllerLaw->cartesianError.head(3) = this->controllerLaw->position_d - robot->getPosition();
        if (this->controllerLaw->orientation_d.coeffs().dot(orientation.coeffs()) < 0.0)
        {
            orientation.coeffs() << -orientation.coeffs();
        }
        Eigen::Quaterniond error_quaternion(orientation.inverse() * this->controllerLaw->orientation_d);
        this->controllerLaw->cartesianError.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
        this->controllerLaw->cartesianError.tail(3) << robot->getT().rotation() * this->controllerLaw->cartesianError.tail(3);

        // 笛卡尔位姿误差导数
        Eigen::Quaterniond dorientation = robot->getdOrientation();
        this->controllerLaw->dcartesianError.head(3) = this->controllerLaw->dposition_d - robot->getdPosition();
        if (this->controllerLaw->dorientation_d.coeffs().dot(dorientation.coeffs()) < 0.0)
        {
            dorientation.coeffs() << -dorientation.coeffs();
        }
        Eigen::Quaterniond derror_quaternion(dorientation.inverse() * this->controllerLaw->dorientation_d);
        this->controllerLaw->cartesianError.tail(3) << derror_quaternion.x(), derror_quaternion.y(), derror_quaternion.z();
        this->controllerLaw->cartesianError.tail(3) << robot->getT().rotation() * this->controllerLaw->cartesianError.tail(3);
    }
    template <int _Dofs, typename pubDataType>
    void Controller<_Dofs, pubDataType>::setControllerLaw(my_robot::Robot<_Dofs> *robot, Eigen::Matrix<double, _Dofs, 1> &tau_d)
    {
        if (controllerLaw.get() != nullptr)
            controllerLaw->setControllerLaw(robot, tau_d);
    }

    //
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

        this->myfile << "M:" << n;
        this->myfile << robot->getM() << "\n";
        this->myfile << "ExternM:" << n;
        this->myfile << robot->getExternM() << "\n";
        this->myfile << "C: " << n;
        this->myfile << robot->getC() * robot->getdq() << "\n";
        this->myfile << "Externc: " << n;
        this->myfile << robot->getExternc() << "\n";
        this->myfile << "G: " << n;
        this->myfile << robot->getG() << n;
        this->myfile << "ExternG: " << n;
        this->myfile << robot->getExternG() << n;
        this->myfile << "J: " << n;
        this->myfile << robot->getJ() << "\n";
        this->myfile << "ExternJ: " << n;
        this->myfile << robot->getExternJ() << "\n";

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
            param_debug.jointError[i] = this->controllerLaw->jointError[i];
            param_debug.q[i] = robot->getq()[i];
            param_debug.dq[i] = robot->getdq()[i];

            param_debug.q_d[i] = this->controllerLaw->q_d[i];
            param_debug.dq_d[i] = this->controllerLaw->dq_d[i];
            param_debug.ddq_d[i] = this->controllerLaw->ddq_d[i];

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
            param_debug.orientation[i] = robot->getOrientation().toRotationMatrix().eulerAngles(2, 1, 0)[i];
            param_debug.orientation_d[i] = this->controllerLaw->orientation_d.toRotationMatrix().eulerAngles(2, 1, 0)[i];
        }
    }
};
