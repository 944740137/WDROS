#pragma once

#include <algorithm/pandaDynLibManager.h>
#include <algorithm/calPseudoInverseMatrix.hpp>
// extern pinLibInteractive *pinInteractive;

namespace my_robot
{
    template <int _Dofs = 7>
    class Robot
    {

    private:
        // limit
        Eigen::Matrix<double, 3, 1> dposMax; /* = {1.7, 1.7, 1.7} */
        Eigen::Matrix<double, 3, 1> doriMax; /*  = {2.5, 2.5, 2.5} */

        // sensor
        Eigen::Matrix<double, _Dofs, 1> q0;
        Eigen::Vector3d position0;
        Eigen::Quaterniond orientation0;

        Eigen::Matrix<double, _Dofs, 1> theta;
        Eigen::Matrix<double, _Dofs, 1> q;
        Eigen::Matrix<double, _Dofs, 1> dq;

        Eigen::Vector3d position;
        Eigen::Quaterniond orientation; // todo 换成欧拉角
        Eigen::Vector3d dposition;
        Eigen::Quaterniond dorientation;

        Eigen::Affine3d TO2E;
        std::array<double, 16> T02EEArray;
        std::array<double, 7> qArray;
        Eigen::Matrix<double, _Dofs, 1> tau;

        // pino
        Eigen::Matrix<double, _Dofs, _Dofs> C;
        Eigen::Matrix<double, _Dofs, _Dofs> M;
        Eigen::Matrix<double, _Dofs, 1> G;
        // Eigen::Matrix<double, _Dofs * 10, _Dofs> Y;
        Eigen::Matrix<double, 6, _Dofs> J;
        Eigen::Matrix<double, 6, _Dofs> dJ;
        Eigen::Matrix<double, _Dofs, 6> J_inv;

        // panda
        Eigen::Matrix<double, _Dofs, 1> externc; // 科氏项,非矩阵
        Eigen::Matrix<double, _Dofs, 7> externM;
        Eigen::Matrix<double, _Dofs, 1> externG;
        Eigen::Matrix<double, 6, _Dofs> externJ;
        Eigen::Matrix<double, 6, _Dofs> externdJ;
        Eigen::Matrix<double, _Dofs, 6> externJ_inv;

    public:
        Eigen::Matrix<double, _Dofs, 1> qMax;      /* = {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973}; */
        Eigen::Matrix<double, _Dofs, 1> qMin;      /* = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973}; */
        Eigen::Matrix<double, _Dofs, 1> dqLimit;   /* = {2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100}; */
        Eigen::Matrix<double, _Dofs, 1> ddqLimit;  /* = {15, 7.5, 10, 12.5, 15, 20, 20}; */
        Eigen::Matrix<double, _Dofs, 1> dddqLimit; /* = {1500, 750, 1000, 1250, 1500, 2000, 2000}; */
        double dposLimit = 0.4;
        double ddposLimit = 4;
        double dddposLimit = 20;
        double doriLimit = 0.4;
        double ddoriLimit = 4;
        double dddoriLimit = 12;

        Robot(const Robot &) = delete;
        void operator=(const Robot &) = delete;

        Robot();
        virtual ~Robot();

        // gei limit
        const Eigen::Matrix<double, 3, 1> &getdposMax() const;
        const Eigen::Matrix<double, 3, 1> &getdoriMax() const;
        const Eigen::Matrix<double, _Dofs, 1> &getqMax() const;
        const Eigen::Matrix<double, _Dofs, 1> &getqMin() const;
        const Eigen::Matrix<double, _Dofs, 1> &getdqLimit() const;
        const Eigen::Matrix<double, _Dofs, 1> &getddqLimit() const;
        const Eigen::Matrix<double, _Dofs, 1> &getdddqLimit() const;

        // get kinematics
        const Eigen::Matrix<double, _Dofs, 1> &getq0();
        const Eigen::Matrix<double, _Dofs, 1> &getq();
        const std::array<double, 7> &getqArray();
        const Eigen::Matrix<double, _Dofs, 1> &getdq();

        const Eigen::Vector3d &getPosition0();
        const Eigen::Quaterniond &getOrientation0();
        const Eigen::Vector3d &getPosition();
        const Eigen::Quaterniond &getOrientation();
        const Eigen::Vector3d &getdPosition();
        const Eigen::Quaterniond &getdOrientation();
        const Eigen::Affine3d &getT();
        const std::array<double, 16> &getT02EEArray();
        const Eigen::Matrix<double, _Dofs, 1> &getTorque();

        // set kinematics
        void setq0(const Eigen::Matrix<double, _Dofs, 1> &q);
        void setPosAndOri0(const Eigen::Vector3d &position, const Eigen::Quaterniond &orientation);
        void updateJointData(const Eigen::Matrix<double, _Dofs, 1> &q, const std::array<double, 7> &qArray, const Eigen::Matrix<double, _Dofs, 1> &theta, const Eigen::Matrix<double, _Dofs, 1> &dq, const Eigen::Matrix<double, _Dofs, 1> &tau);
        void updateEndeffectorData(const Eigen::Vector3d &position, const Eigen::Quaterniond &orientation, const Eigen::Affine3d &TO2E, const std::array<double, 16> &T02EEArray);

        // other dyn api
        void setExternM(const Eigen::Matrix<double, _Dofs, _Dofs> &externM);
        void setExternc(const Eigen::Matrix<double, _Dofs, 1> &externc);
        void setExternG(const Eigen::Matrix<double, _Dofs, 1> &externG);
        void setExternJ(const Eigen::Matrix<double, 6, _Dofs> &externJ);
        const Eigen::Matrix<double, _Dofs, _Dofs> &getExternM();
        const Eigen::Matrix<double, _Dofs, 1> &getExternc(); // 科氏项
        const Eigen::Matrix<double, _Dofs, 1> &getExternG();
        const Eigen::Matrix<double, 6, _Dofs> &getExternJ();
        const Eigen::Matrix<double, _Dofs, 6> &getExternJ_inv();

        // pinocchino dyn api
        void calculation(Eigen::Matrix<double, _Dofs, 1> ddq_d);
        const Eigen::Matrix<double, 6, _Dofs> &getJ();
        const Eigen::Matrix<double, 6, _Dofs> &getdJ();
        const Eigen::Matrix<double, _Dofs, 6> &getJ_inv();
        const Eigen::Matrix<double, _Dofs, _Dofs> &getM();
        const Eigen::Matrix<double, _Dofs, _Dofs> &getC();
        const Eigen::Matrix<double, _Dofs, 1> &getG();
    };

}

namespace my_robot
{
    template <int _Dofs>
    Robot<_Dofs>::~Robot()
    {
    }
    template <int _Dofs>
    Robot<_Dofs>::Robot() : q0(Eigen::Matrix<double, _Dofs, 1>::Zero()), position0(Eigen::Matrix<double, 3, 1>::Zero()),
                            orientation0(Eigen::Quaterniond::Identity()), theta(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                            q(Eigen::Matrix<double, _Dofs, 1>::Zero()), dq(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                            position(Eigen::Matrix<double, 3, 1>::Zero()), orientation(Eigen::Quaterniond::Identity()),
                            dposition(Eigen::Matrix<double, 3, 1>::Zero()), dorientation(Eigen::Quaterniond::Identity()),
                            TO2E(Eigen::Affine3d::Identity()), tau(Eigen::Matrix<double, _Dofs, 1>::Zero()),
                            C(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()), M(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                            G(Eigen::Matrix<double, _Dofs, 1>::Zero()), J(Eigen::Matrix<double, 6, _Dofs>::Zero()),
                            dJ(Eigen::Matrix<double, 6, _Dofs>::Zero()), J_inv(Eigen::Matrix<double, _Dofs, 6>::Zero()),
                            externc(Eigen::Matrix<double, _Dofs, 1>::Zero()), externM(Eigen::Matrix<double, _Dofs, _Dofs>::Zero()),
                            externG(Eigen::Matrix<double, _Dofs, 1>::Zero()), externJ(Eigen::Matrix<double, 6, _Dofs>::Zero()),
                            externdJ(Eigen::Matrix<double, 6, _Dofs>::Zero()), externJ_inv(Eigen::Matrix<double, _Dofs, 6>::Zero())
    {
        // 全部数据初始化为0
    }

    template <int _Dofs>
    const Eigen::Matrix<double, 3, 1> &Robot<_Dofs>::getdposMax() const
    {
        return this->dposMax;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, 3, 1> &Robot<_Dofs>::getdoriMax() const
    {
        return this->doriMax;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getqMax() const
    {
        return this->qMax;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getqMin() const
    {
        return this->qMin;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getdqLimit() const
    {
        return this->dqLimit;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getddqLimit() const
    {
        return this->ddqLimit;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getdddqLimit() const
    {
        return this->dddqLimit;
    }
    // 初值设置和获取
    template <int _Dofs>
    void Robot<_Dofs>::setq0(const Eigen::Matrix<double, _Dofs, 1> &q)
    {
        this->q0 = q;
    }
    template <int _Dofs>
    void Robot<_Dofs>::setPosAndOri0(const Eigen::Vector3d &position0, const Eigen::Quaterniond &orientation0)
    {
        this->position0 = position0;
        this->orientation0 = orientation0;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getq0()
    {
        return this->q0;
    }
    template <int _Dofs>
    const Eigen::Vector3d &Robot<_Dofs>::getPosition0()
    {
        return this->position0;
    }
    template <int _Dofs>
    const Eigen::Quaterniond &Robot<_Dofs>::getOrientation0()
    {
        return this->orientation0;
    }

    // 获取传感器数据
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getq()
    {
        return this->q;
    }
    template <int _Dofs>
    const std::array<double, 7> &Robot<_Dofs>::getqArray()
    {
        return this->qArray;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getdq()
    {
        return this->dq;
    }
    template <int _Dofs>
    const Eigen::Vector3d &Robot<_Dofs>::getPosition()
    {
        return this->position;
    }
    template <int _Dofs>
    const Eigen::Quaterniond &Robot<_Dofs>::getOrientation()
    {
        return this->orientation;
    }
    template <int _Dofs>
    const Eigen::Vector3d &Robot<_Dofs>::getdPosition()
    {
        return this->dposition;
    }
    template <int _Dofs>
    const Eigen::Quaterniond &Robot<_Dofs>::getdOrientation()
    {
        return this->dorientation;
    }
    template <int _Dofs>
    const Eigen::Affine3d &Robot<_Dofs>::getT()
    {
        return this->TO2E;
    }
    template <int _Dofs>
    const std::array<double, 16> &Robot<_Dofs>::getT02EEArray()
    {
        return this->T02EEArray;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getTorque()
    {
        return this->tau;
    }

    // 传感器数据更新
    template <int _Dofs>
    void Robot<_Dofs>::updateJointData(const Eigen::Matrix<double, _Dofs, 1> &q, const std::array<double, 7> &qArray, const Eigen::Matrix<double, _Dofs, 1> &theta, const Eigen::Matrix<double, _Dofs, 1> &dq, const Eigen::Matrix<double, _Dofs, 1> &tau)
    {
        this->q = q;
        this->qArray = qArray;
        this->theta = theta;
        this->dq = dq;
        this->tau = tau;
    }
    template <int _Dofs>
    void Robot<_Dofs>::updateEndeffectorData(const Eigen::Vector3d &position, const Eigen::Quaterniond &orientation, const Eigen::Affine3d &TO2E, const std::array<double, 16> &T02EEArray)
    {
        this->position = position;
        this->orientation = orientation;

        this->dposition = (this->externJ * dq).head(3);
        Eigen::Matrix<double, 3, 1> tmp = (this->externJ * dq).tail(3); // 欧拉角
        this->dorientation = Eigen::AngleAxisd(tmp[0], Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(tmp[1], Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(tmp[2], Eigen::Vector3d::UnitZ());

        this->TO2E = TO2E;
        this->T02EEArray = T02EEArray;
    }

    // pino dyn
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, _Dofs> &Robot<_Dofs>::getM()
    {
        return this->M;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, _Dofs> &Robot<_Dofs>::getC()
    {
        return this->C;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getG()
    {
        return this->G;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, 6, _Dofs> &Robot<_Dofs>::getJ()
    {
        return this->J;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, 6, _Dofs> &Robot<_Dofs>::getdJ()
    {
        return this->dJ;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 6> &Robot<_Dofs>::getJ_inv()
    {
        return this->J_inv;
    }
    template <int _Dofs>
    void Robot<_Dofs>::calculation(Eigen::Matrix<double, _Dofs, 1> ddq_d)
    {
        pPandaDynLibManager->upDataModel(this->q);
        pPandaDynLibManager->computeKinData(this->J, this->dJ, this->q, this->dq);
        pPandaDynLibManager->computeDynData(this->M, this->C, this->G, this->q, this->dq);
        weightedPseudoInverse(this->J, this->J_inv, this->M);
    }

    // other
    template <int _Dofs>
    void Robot<_Dofs>::setExternM(const Eigen::Matrix<double, _Dofs, _Dofs> &externM)
    {
        this->externM = externM;
    }
    template <int _Dofs>
    void Robot<_Dofs>::setExternc(const Eigen::Matrix<double, _Dofs, 1> &externc)
    {
        this->externc = externc;
    }
    template <int _Dofs>
    void Robot<_Dofs>::setExternG(const Eigen::Matrix<double, _Dofs, 1> &externG)
    {
        this->externG = externG;
    }
    template <int _Dofs>
    void Robot<_Dofs>::setExternJ(const Eigen::Matrix<double, 6, _Dofs> &externJ)
    {
        this->externJ = externJ;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, _Dofs> &Robot<_Dofs>::getExternM()
    {
        return this->externM;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getExternc()
    {
        return this->externc;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 1> &Robot<_Dofs>::getExternG()
    {
        return this->externG;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, 6, _Dofs> &Robot<_Dofs>::getExternJ()
    {
        return this->externJ;
    }
    template <int _Dofs>
    const Eigen::Matrix<double, _Dofs, 6> &Robot<_Dofs>::getExternJ_inv()
    {
        weightedPseudoInverse(this->externJ, this->externJ_inv, this->M);
        return this->externJ_inv;
    }
}
