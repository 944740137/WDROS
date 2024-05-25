#include <cmath>
#include <queue>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex.h>
// 关节规划器
template <int _Dofs>
class Planner
{
public:
public:
    Planner(const Planner &) = delete;
    void operator=(const Planner &) = delete;

    Planner();
    virtual ~Planner();

    virtual void calPlanQueue(bool isCoordinated, double deltaT,
                              const Eigen::Matrix<double, _Dofs, 1> &maxVel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &jerk,
                              const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                              std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d) = 0;
};
template <int _Dofs>
Planner<_Dofs>::~Planner()
{
}
template <int _Dofs>
Planner<_Dofs>::Planner()
{
    // 全部初始化为0
}

// 五次多项式规划：计算时间
template <int _Dofs>
void calQuinticPlanTime(bool isCoordinated, double *T,
                        const Eigen::Matrix<double, _Dofs, 1> &maxVel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc,
                        const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf)
{
    double Tf = 0.0;
    for (int i = 0; i < _Dofs; i++)
    {
        double error = std::fabs(qf[i] - q0[i]);
        double t1 = (15.0 / 8.0 * error) / (maxVel[i]); // note 15.0/8.0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double t2 = sqrt(5.7735 * error / maxAcc[i]);   // 10/3^0.5
        T[i] = std::max(t1, t2);
        Tf = std::max(Tf, T[i]);
    }
    if (isCoordinated)
    {
        for (int i = 0; i < _Dofs; i++)
        {
            T[i] = Tf;
        }
    }
}
// 五次多项式规划：计算队列
template <int _Dofs>
void calQuinticPlan(bool isCoordinated, double deltaT,
                    const Eigen::Matrix<double, _Dofs, 1> &maxVel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc,
                    const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf,
                    std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d,
                    std::vector<std::queue<double>> &ddq_d)
{
    // isCoordinated为false未测试
    double T[_Dofs] = {0.0};
    double dq0[_Dofs] = {0.0};
    double ddq0[_Dofs] = {0.0};
    double dqf[_Dofs] = {0.0};
    double ddqf[_Dofs] = {0.0};
    int pointNum[_Dofs] = {0};

    calQuinticPlanTime<_Dofs>(isCoordinated, T, maxVel, maxAcc, q0, qf);

    double t = 0;
    for (int i = 0; i < _Dofs; i++)
    {
        pointNum[i] = static_cast<int>(T[i] / deltaT) + 1;
        double a0 = q0[i];
        double a1 = 0;
        double a2 = 0;
        double a3 = (20 * qf[i] - 20 * q0[i] - (8 * dqf[i] + 12 * dq0[i]) * T[i] - (3 * ddq0[i] - ddqf[i]) * pow(T[i], 2)) / (2 * pow(T[i], 3));
        double a4 = (30 * q0[i] - 30 * qf[i] + (14 * dqf[i] + 16 * dq0[i]) * T[i] + (3 * ddq0[i] - 2 * ddqf[i]) * pow(T[i], 2)) / (2 * pow(T[i], 4));
        double a5 = (12 * qf[i] - 12 * q0[i] - (6 * dqf[i] + 6 * dq0[i]) * T[i] - (ddq0[i] - ddqf[i]) * pow(T[i], 2)) / (2 * pow(T[i], 5));

        for (int j = 1; j <= pointNum[i]; j++)
        {
            t = j * deltaT;
            if (j == pointNum[i])
                q_d[i].push(qf[i]);
            else
                q_d[i].push(a0 + a1 * t + a2 * pow(t, 2) + a3 * pow(t, 3) + a4 * pow(t, 4) + a5 * pow(t, 5));
            dq_d[i].push(a1 + 2 * a2 * t + 3 * a3 * pow(t, 2) + 4 * a4 * pow(t, 3) + 5 * a5 * pow(t, 4));
            ddq_d[i].push(2 * a2 + 6 * a3 * t + 12 * a4 * pow(t, 2) + 20 * a5 * pow(t, 3));
        }
    }
}

//***********************************************************************************************************************************//

// 急停规划：最短时间规划
template <int _Dofs>
bool calStopPlanParam(double deltaT, double jerk, double dq0, double ddq0,
                      double &jerk1, double &jerk2, double &maxAcc, double &planTime)
{
    double S = -dq0;
    double S1, S2, S3;
    if (S < 0)
        maxAcc = -maxAcc;

    if (ddq0 <= maxAcc)
        jerk1 = jerk;
    else
        jerk1 = -jerk;
    if (0 <= maxAcc)
        jerk2 = -jerk;
    else
        jerk2 = jerk;

    if ((maxAcc * ddq0 > 0) && std::fabs(S) <= (std::fabs(ddq0 * ddq0) / (2 * jerk1)))
    {
        return false;
    }

    if (std::fabs((maxAcc * maxAcc - ddq0 * ddq0) / (2 * jerk1) + (0 - maxAcc * maxAcc) / (2 * jerk2)) < (std::fabs(S)))
    {
        S1 = (maxAcc * maxAcc - ddq0 * ddq0) / (2 * jerk1);
        S3 = (0 - maxAcc * maxAcc) / (2 * (jerk2));
        S2 = S - (S1 + S3);
    }
    else
    {
        maxAcc = sqrt((2 * jerk1 * S + ddq0 * ddq0) / 2);
        if (S < 0)
            maxAcc = -maxAcc;
        S1 = (maxAcc * maxAcc - ddq0 * ddq0) / (2 * jerk1);
        S3 = (0 - maxAcc * maxAcc) / (2 * (jerk2));
        S2 = 0;
    }

    double T1 = std::fabs((maxAcc - ddq0) / (jerk1));
    double T2 = std::fabs(S2 / maxAcc);
    double T3 = std::fabs((0 - maxAcc) / (jerk2));

    planTime = T1 + T2 + T3;
    return true;
}
// 急停规划：固定时间规划
template <int _Dofs>
void calStopPlanQueue(double deltaT, double q0, double dq0, double ddq0, double maxPlanTime,
                      double jerk1, double jerk2, double maxAcc,
                      std::queue<double> &q_d, std::queue<double> &dq_d, std::queue<double> &ddq_d)
{
    double S = -dq0;
    if (ddq0 * maxAcc > 0 && std::fabs(maxAcc) > std::fabs(ddq0)) // 确定jerk1的符号
    {
        double T = std::fabs(ddq0 / jerk1) + std::fabs((S - ((-ddq0 * ddq0) / (2 * jerk2))) / ddq0);
        if (maxPlanTime <= T)
        {
            jerk1 = jerk1;
        }
        else
        {
            jerk1 = -jerk1;
        }
    }
    double A = jerk1 - jerk2;
    double B = (2 * jerk1 * jerk2 * maxPlanTime + 2 * jerk2 * ddq0);
    double C = -2 * jerk1 * jerk2 * S - jerk2 * ddq0 * ddq0;
    double Acc;
    if (A == 0)
    {
        Acc = -C / B;
    }
    else
    {
        Acc = (-B - std::sqrt(B * B - 4 * A * C)) / (2 * A);
    }

    double S1 = (Acc * Acc - ddq0 * ddq0) / (2 * jerk1);
    double S3 = (0 - Acc * Acc) / (2 * jerk2);
    double S2 = S - (S1 + S3);

    double T1 = std::fabs((Acc - ddq0) / jerk1);
    double T2 = std::fabs(S2 / Acc);
    double T3 = std::fabs((0 - Acc) / jerk2);
    double q = q0;
    int totalNum = std::ceil((T1 + T2 + T3) / deltaT);
    for (int num = 0; num <= totalNum; ++num)
    {
        double t = num * deltaT;
        double detladq;
        double detladdq;
        if (t >= 0 && t < T1)
        {
            double nowt = t;
            detladdq = jerk1 * nowt;
            detladq = ddq0 * nowt + 0.5 * jerk1 * nowt * nowt;
            q = q + (dq0 + detladq) * deltaT;
            q_d.push(q);
            dq_d.push(dq0 + detladq);
            ddq_d.push(ddq0 + detladdq);
        }
        else if (t >= T1 && t < (T1 + T2))
        {
            double nowt = t - T1;
            detladq = S1 + nowt * Acc;
            q = q + (dq0 + detladq) * deltaT;
            q_d.push(q);
            dq_d.push(dq0 + detladq);
            ddq_d.push(Acc);
        }
        else if (t >= (T1 + T2) && t < (T1 + T2 + T3))
        {
            double nowt = t - T1 - T2;
            detladdq = jerk2 * nowt;
            detladq = S1 + S2 + Acc * nowt + 0.5 * jerk2 * nowt * nowt;
            q = q + (dq0 + detladq) * deltaT;
            q_d.push(q);
            dq_d.push(dq0 + detladq);
            ddq_d.push(Acc + detladdq);
        }
        else if (num == totalNum)
        {
            q_d.push(q);
            dq_d.push(0.0);
            ddq_d.push(0.0);
        }
    }
}
// 急停规划：计算队列
template <int _Dofs>
bool calStopPlan(bool isCoordinated, double deltaT, Eigen::Matrix<double, _Dofs, 1> &dddq, Eigen::Matrix<double, _Dofs, 1> &ddq,
                 Eigen::Matrix<double, _Dofs, 1> &q0, Eigen::Matrix<double, _Dofs, 1> &dq0, Eigen::Matrix<double, _Dofs, 1> &ddq0,
                 std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
{
    double jerk1[_Dofs] = {0.0};
    double jerk2[_Dofs] = {0.0};
    double maxAcc[_Dofs] = {0.0};
    double planTime = 0;
    double maxPlanTime = 0;
    for (int i = 0; i < _Dofs; i++)
    {
        maxAcc[i] = ddq[i];
        if (!calStopPlanParam<_Dofs>(deltaT, dddq[i], dq0[i], ddq0[i], jerk1[i], jerk2[i], maxAcc[i], planTime))
            return false;
        if (i == 0 || maxPlanTime < planTime)
            maxPlanTime = planTime;
    }
    for (int i = 0; i < _Dofs; i++)
    {
        calStopPlanQueue<_Dofs>(deltaT, q0[i], dq0[i], ddq0[i], maxPlanTime, jerk1[i], jerk2[i], maxAcc[i], q_d[i], dq_d[i], ddq_d[i]);
    }
    return true;
}

//***********************************************************************************************************************************//

// TVP规划：最短时间规划
template <int _Dofs>
bool calTVPPlanParam(double deltaT, double q0, double qf, double dq0, double maxAcc,
                     double &vel, double &acc1, double &acc2, double &planTime)
{
    double S = qf - q0;
    double S1, S2, S3;
    if (S < 0)
        vel = -vel;

    if (dq0 <= vel)
        acc1 = maxAcc;
    else
        acc1 = -maxAcc;
    if (0 <= vel)
        acc2 = -maxAcc;
    else
        acc2 = maxAcc;

    if ((vel * dq0 > 0) && std::fabs(S) <= (std::fabs(dq0 * dq0) / (2 * acc1)))
    {
        return false;
    }

    if (std::fabs((vel * vel - dq0 * dq0) / (2 * acc1) + (0 - vel * vel) / (2 * acc2)) < (std::fabs(S)))
    {
        S1 = (vel * vel - dq0 * dq0) / (2 * acc1);
        S3 = (0 - vel * vel) / (2 * (acc2));
        S2 = S - (S1 + S3);
    }
    else
    {
        vel = sqrt((2 * acc1 * S + dq0 * dq0) / 2);
        if (S < 0)
            vel = -vel;
        S1 = (vel * vel - dq0 * dq0) / (2 * acc1);
        S3 = (0 - vel * vel) / (2 * (acc2));
        S2 = 0;
    }

    double T1 = std::fabs((vel - dq0) / (acc1));
    double T2 = std::fabs(S2 / vel);
    double T3 = std::fabs((0 - vel) / (acc2));
    planTime = T1 + T2 + T3;
    return true;
}
// TVP规划：固定时间规划
template <int _Dofs>
void calTVPPlanQueue(double deltaT, double q0, double qf, double dq0, double maxPlanTime,
                     double acc1, double acc2, double vel,
                     std::queue<double> &q_d, std::queue<double> &dq_d, std::queue<double> &ddq_d)
{
    double S = qf - q0;
    if (dq0 * vel > 0 && std::fabs(vel) > std::fabs(dq0)) // 确定acc1的符号
    {
        double T = std::fabs(dq0 / acc1) + std::fabs((S - ((-dq0 * dq0) / (2 * acc2))) / dq0);
        if (maxPlanTime <= T)
            acc1 = acc1;
        else
            acc1 = -acc1;
    }
    double A = acc1 - acc2;
    double B = (2 * acc1 * acc2 * maxPlanTime + 2 * acc2 * dq0);
    double C = -2 * acc1 * acc2 * S - acc2 * dq0 * dq0;
    if (A == 0)
    {
        vel = -C / B;
    }
    else
    {
        vel = (-B - std::sqrt(B * B - 4 * A * C)) / (2 * A);
    }

    double S1 = (vel * vel - dq0 * dq0) / (2 * acc1);
    double S3 = (0 - vel * vel) / (2 * acc2);
    double S2 = S - (S1 + S3);

    double T1 = std::fabs((vel - dq0) / acc1);
    double T2 = std::fabs(S2 / vel);
    double T3 = std::fabs((0 - vel) / acc2);

    int totalNum = std::ceil((T1 + T2 + T3) / deltaT);
    for (int num = 0; num <= totalNum; ++num)
    {
        double t = num * deltaT;
        double detlaq;
        double detladq;
        if (t >= 0 && t < T1)
        {
            double nowt = t;
            detladq = acc1 * nowt;
            detlaq = dq0 * nowt + 0.5 * acc1 * nowt * nowt;
            q_d.push(q0 + detlaq);
            dq_d.push(dq0 + detladq);
            ddq_d.push(acc1);
        }
        else if (t >= T1 && t < (T1 + T2))
        {
            double nowt = t - T1;
            detlaq = S1 + nowt * vel;
            q_d.push(q0 + detlaq);
            dq_d.push(vel);
            ddq_d.push(0);
        }
        else if (t >= (T1 + T2) && t < (T1 + T2 + T3))
        {
            double nowt = t - T1 - T2;
            detladq = acc2 * nowt;
            detlaq = S1 + S2 + vel * nowt + 0.5 * acc2 * nowt * nowt;
            q_d.push(q0 + detlaq);
            dq_d.push(vel + detladq);
            ddq_d.push(acc2);
        }
        else if (num == totalNum)
        {
            q_d.push(qf);
            dq_d.push(0.0);
            ddq_d.push(0.0);
        }
    }
}
// TVP规划：计算队列
template <int _Dofs>
bool calTVPPlan(bool isCoordinated, double deltaT,
                const Eigen::Matrix<double, _Dofs, 1> &vel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc,
                const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
{
    double acc1[_Dofs] = {0.0};
    double acc2[_Dofs] = {0.0};
    double planTime = 0;
    double maxPlanTime = 0;

    Eigen::Matrix<double, _Dofs, 1> atuVel = vel;
    for (int i = 0; i < _Dofs; i++)
    {
        if (!calTVPPlanParam<_Dofs>(deltaT, q0[i], qf[i], dq0[i], maxAcc[i], atuVel[i], acc1[i], acc2[i], planTime))
            return false;
        if (i == 0 || maxPlanTime < planTime)
            maxPlanTime = planTime;
    }
    for (int i = 0; i < _Dofs; i++)
    {
        calTVPPlanQueue<_Dofs>(deltaT, q0[i], qf[i], dq0[i], maxPlanTime, acc1[i], acc2[i], atuVel[i], q_d[i], dq_d[i], ddq_d[i]);
    }
    return true;
}

// SS规划：最短时间规划
template <int _Dofs>
bool calSSPlanParam(double deltaT, double q0, double qf, double vel, double maxAcc, double maxJerk,
                    double &acc1, double &acc2, double &j1, double &j2, double &j3, double &j4, double &planTime)
{
    double S = qf - q0;
    double T1, T2, T3, T4, T5, T6, T7;
    double S1, S2, S3, S4, S5, S6, S7;
    if (S < 0)
    {
        vel = -vel;
        acc1 = -maxAcc; // 加速段最大加速度
        acc2 = maxAcc;  // 减速段最大加速度
        j1 = -maxJerk;
        j2 = maxJerk;
        j3 = maxJerk;
        j4 = -maxJerk;
    }
    else
    {
        acc1 = maxAcc;
        acc2 = -maxAcc;
        j1 = maxJerk;
        j2 = -maxJerk;
        j3 = -maxJerk;
        j4 = maxJerk;
    }

    // 根据Vm判断有无匀加减速段
    if (std::fabs((std::pow(acc1, 2) - 0) / (2 * j1) + (0 - std::pow(acc1, 2)) / (2 * j2)) < std::fabs(vel - 0))
    {
        T1 = std::fabs(acc1 / j1);
        T3 = std::fabs(-acc1 / j2);
        T2 = (vel - (std::pow(acc1, 2) - 0) / (2 * j1) - (0 - std::pow(acc1, 2)) / (2 * j2)) / acc1;

        // 判断有无匀速段
        S1 = (1.0 / 6.0) * j1 * std::pow(T1, 3);
        S2 = (1.0 / 2.0) * j1 * std::pow(T1, 2) * T2 + (1.0 / 2.0) * acc1 * std::pow(T2, 2);
        S3 = ((1.0 / 2.0) * j1 * std::pow(T1, 2) + acc1 * T2) * T3 + (1.0 / 2.0) * acc1 * std::pow(T3, 2) + (1.0 / 6.0) * j2 * std::pow(T1, 3);

        if (2 * std::fabs(S1 + S2 + S3) < std::fabs(S))
        {
            // 有匀加减速段,有匀速段 1
            T4 = (S - 2 * (S1 + S2 + S3)) / vel;
            T5 = T3;
            T6 = T2;
            T7 = T1;
        }
        else
        {
            S1 = (1.0 / 6.0) * j1 * std::pow(T1, 3);
            S3 = (1.0 / 2.0) * j1 * std::pow(T1, 2) * T3 + (1.0 / 2.0) * acc1 * std::pow(T3, 2) + (1.0 / 6.0) * j2 * std::pow(T1, 3);
            if (2 * std::fabs(S1 + S3) < std::fabs(S))
            {
                // 有匀加减速段,无匀速段 2
                double Se = (S - (S1 + S3)) / 2.0;
                double A = (1.0 / 2.0) * acc1;
                double B = (1.0 / 2.0) * j1 * std::pow(T1, 2) + T3 * acc1;
                double C = (1.0 / 2.0) * j1 * std::pow(T1, 2) * T3 - Se;
                T2 = std::fabs(std::max((-B - std::sqrt(std::pow(B, 2) - 4 * A * C)) / (2 * A), (-B + std::sqrt(std::pow(B, 2) - 4 * A * C)) / (2 * A)));
                T4 = 0;
                T5 = T3;
                T6 = T2;
                T7 = T1;
            }
            else
            {
                // 无匀加减速段,无匀速段 3
                T1 = std::fabs(std::pow(S / (2 * ((1.0 / 6.0) * j1 + (1.0 / 6.0) * j2 + j1)), 1.0 / 3.0));
                acc1 = j1 * T1;
                acc2 = -acc1;
                T2 = 0;
                T3 = T1;
                T4 = 0;
                T5 = T1;
                T6 = 0;
                T7 = T1;
            }
        }
    }
    else
    {
        acc1 = (vel >= 0 ? 1 : -1) * std::sqrt(vel * j1);
        T1 = std::fabs(acc1 / j1);
        T3 = T1;
        S1 = (1.0 / 6.0) * j1 * std::pow(T1, 3);
        S3 = (1.0 / 2.0) * j1 * std::pow(T1, 2) * T3 + (1.0 / 2.0) * acc1 * std::pow(T3, 2) + (1.0 / 6.0) * j2 * std::pow(T1, 3);

        if (2 * std::fabs(S1 + S3) < std::fabs(S))
        {
            // 无匀加减速段,,有匀速段 4
            acc2 = -acc1;
            T4 = std::fabs((S - 2 * (S1 + S3)) / vel);
            T2 = 0;
            T5 = T3;
            T6 = 0;
            T7 = T1;
        }
        else
        {
            // 无匀加减速段,无匀速段 5
            T1 = std::fabs(std::pow(S / (2 * ((1.0 / 6.0) * j1 + (1.0 / 6.0) * j2 + j1)), 1.0 / 3.0));
            acc1 = j1 * T1;
            acc2 = -acc1;
            T2 = 0;
            T3 = T1;
            T4 = 0;
            T5 = T1;
            T6 = 0;
            T7 = T1;
        }
    }
    planTime = T1 + T2 + T3 + T4 + T5 + T6 + T7;
    return true;
}
// SS规划：固定时间规划
template <int _Dofs>
void calSSPlanQueue(double deltaT, double q0, double qf, double vel, double maxAcc, double maxJerk,
                    double &acc1, double &acc2, double &j1, double &j2, double &j3, double &j4, double &planTime,
                    std::queue<double> &q_d, std::queue<double> &dq_d, std::queue<double> &ddq_d)
{
    double S = qf - q0;
    double T1, T2, T3, T4, T5, T6, T7;
    double S1, S2, S3, S4, S5, S6, S7;

    T1 = std::fabs(acc1 / j1);
    double v = 2 * (1 / 2.0) * j1 * T1 * T1;
    double S13 = j1 * std::pow(T1, 3);

    auto findCubicRoot = [&]() mutable
    {
        double a = -2 * j1;
        double b = j1 * planTime;
        double c = 0;
        double d = -S;
        double p = (3 * a * c - b * b) / (3 * a * a);
        double q = (2 * std::pow(b, 3) - 9 * a * b * c + 27 * a * a * d) / (27 * std::pow(a, 3));

        std::complex<double> theta[3], t[3], x[3];
        for (int i = 0; i < 3; ++i)
        {
            theta[i] = (1.0 / 3.0) * std::acos((3 * q / (2 * p)) * std::sqrt(-3.0 / p)) - (2 * i / 3.0) * M_PI;
            t[i] = 2 * std::sqrt(-p / 3.0) * std::cos(theta[i]);
            x[i] = t[i] - b / (3.0 * a);
        }
        T1 = 0;
        for (int i = 0; i < 3; ++i)
        {
            if (std::fabs(x[i].imag()) <= 0.000000001 && x[i].real() > 0)
            {
                if ((T1 != 0 && x[i].real() < T1) || T1 == 0)
                {
                    T1 = x[i].real();
                }
            }
        }
        acc1 = j1 * T1;
        acc2 = -acc1;
        T2 = 0;
        T3 = T1;
        T4 = planTime - 4 * T1;
        T5 = T1;
        T6 = 0;
        T7 = T1;
    };

    if ((std::fabs(S) - std::fabs(2 * S13)) > 0)
    {
        if ((planTime + deltaT) > 4 * T1)
        {
            T4 = std::fabs((S - 2 * S13) / v);
            if (4 * T1 + T4 > (planTime + deltaT))
            {
                // std::cout << "有匀加减速段,有匀速段 1 " << std::endl;
                S13 = S13 - (0.5) * j1 * T1 * T1 * T1;
                double Se = S - 2 * S13;
                double Te = planTime - 4 * T1;
                double A = -acc1;
                double B = j1 * std::pow(T1, 2) + acc1 * Te;
                double C = j1 * std::pow(T1, 3) + acc1 * T1 * Te - Se;
                T2 = std::min((-B - std::sqrt(std::pow(B, 2) - 4 * A * C)) / (2 * A), (-B + std::sqrt(std::pow(B, 2) - 4 * A * C)) / (2 * A));
                T3 = T1;
                T4 = planTime - 4 * T1 - 2 * T2;
                T5 = T1;
                T6 = T2;
                T7 = T1;
            }
            else
            {
                // std::cout << "无加减速段,有匀速段 2 " << std::endl;
                findCubicRoot();
            }
        }
        else
        {
            std::cout << "SS 求解错误 1" << std::endl;
        }
    }
    else
    {
        if ((4 * T1) < (planTime + deltaT))
        {
            // std::cout << "无加减速段,有匀速段 3 " << std::endl;
            findCubicRoot();
        }
        else
        {
            std::cout << "SS 求解错误 2" << std::endl;
        }
    }

    int totalNum = std::ceil((T1 + T2 + T3 + T4 + T5 + T6 + T7) / deltaT);
    double q = q0;
    double dq = 0;
    double ddq = 0;
    for (int num = 0; num <= totalNum; ++num)
    {
        double t = num * deltaT;
        double detlaq;
        double detladq;
        if ((t >= 0) && (t < T1))
        {
            double nowt = t;
            ddq = j1 * nowt;
        }
        else if ((t >= T1) && (t < (T1 + T2)))
        {
            ddq = acc1;
        }
        else if ((t >= (T1 + T2)) && (t < (T1 + T2 + T3)))
        {
            double nowt = t - T1 - T2;
            ddq = acc1 + j2 * nowt;
        }
        else if ((t >= (T1 + T2 + T3)) && (t < (T1 + T2 + T3 + T4)))
        {
            ddq = 0;
        }
        else if ((t >= (T1 + T2 + T3 + T4)) && (t < (T1 + T2 + T3 + T4 + T5)))
        {
            double nowt = t - T1 - T2 - T3 - T4;
            ddq = j3 * nowt;
        }
        else if ((t >= (T1 + T2 + T3 + T4 + T5)) && (t < (T1 + T2 + T3 + T4 + T5 + T6)))
        {
            ddq = acc2;
        }
        else if ((t >= (T1 + T2 + T3 + T4 + T5 + T6)) && (t < (T1 + T2 + T3 + T4 + T5 + T6 + T7)))
        {
            double nowt = t - T1 - T2 - T3 - T4 - T5 - T6;
            ddq = acc2 + j4 * nowt;
        }
        if (num == totalNum)
        {
            q_d.push(qf);
            dq_d.push(0);
            ddq_d.push(0);
        }
        else
        {
            dq = dq + ddq * deltaT;
            q = q + dq * deltaT;
            q_d.push(q);
            dq_d.push(dq);
            ddq_d.push(ddq);
        }
    }
}
// SS规划：计算队列
template <int _Dofs>
bool calSSPlan(bool isCoordinated, double deltaT,
               const Eigen::Matrix<double, _Dofs, 1> &vel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &maxJerk,
               const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
               std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
{
    double acc1[_Dofs] = {0.0};
    double acc2[_Dofs] = {0.0};
    double j1[_Dofs] = {0.0};
    double j2[_Dofs] = {0.0};
    double j3[_Dofs] = {0.0};
    double j4[_Dofs] = {0.0};
    double planTime = 0;
    double maxPlanTime = 0;

    Eigen::Matrix<double, _Dofs, 1> atuVel = vel;
    for (int i = 0; i < _Dofs; i++)
    {
        if (!calSSPlanParam<_Dofs>(deltaT, q0[i], qf[i], atuVel[i], maxAcc[i], maxJerk[i], acc1[i], acc2[i], j1[i], j2[i], j3[i], j4[i], planTime))
            return false;
        if (i == 0 || maxPlanTime < planTime)
            maxPlanTime = planTime;
    }
    std::cout << "maxPlanTime  " << maxPlanTime << std::endl;

    for (int i = 0; i < _Dofs; i++)
    {
        calSSPlanQueue<_Dofs>(deltaT, q0[i], qf[i], atuVel[i], maxAcc[i], maxJerk[i], acc1[i], acc2[i], j1[i], j2[i], j3[i], j4[i], maxPlanTime, q_d[i], dq_d[i], ddq_d[i]);
    }
    std::cout << "--------------------------------------------" << std::endl;

    return true;
}

// 五次多项式规划器
template <int _Dofs>
class QuinticPlanner : public Planner<_Dofs>
{
public:
public:
    QuinticPlanner(const QuinticPlanner &) = delete;
    void operator=(const QuinticPlanner &) = delete;

    QuinticPlanner();
    ~QuinticPlanner();

    void calPlanQueue(bool isCoordinated, double deltaT,
                      const Eigen::Matrix<double, _Dofs, 1> &maxVel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &jerk,
                      const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                      std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);
};
template <int _Dofs>
QuinticPlanner<_Dofs>::~QuinticPlanner()
{
}
template <int _Dofs>
QuinticPlanner<_Dofs>::QuinticPlanner()
{
    std::cout << "[robotController] 设置规划器: QuinticPlanner" << std::endl;
}
template <int _Dofs>
void QuinticPlanner<_Dofs>::calPlanQueue(bool isCoordinated, double deltaT, const Eigen::Matrix<double, _Dofs, 1> &maxVel,
                                         const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &jerk,
                                         const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                                         std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
{
    calQuinticPlan(isCoordinated, deltaT, maxVel, maxAcc, q0, qf, q_d, dq_d, ddq_d);
}

// TVP规划器
template <int _Dofs>
class TVPPlanner : public Planner<_Dofs>
{
public:
public:
    TVPPlanner(const TVPPlanner &) = delete;
    void operator=(const TVPPlanner &) = delete;

    TVPPlanner();
    ~TVPPlanner();

    void calPlanQueue(bool isCoordinated, double deltaT,
                      const Eigen::Matrix<double, _Dofs, 1> &maxVel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &jerk,
                      const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                      std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);
};
template <int _Dofs>
TVPPlanner<_Dofs>::~TVPPlanner()
{
}
template <int _Dofs>
TVPPlanner<_Dofs>::TVPPlanner()
{
    std::cout << "[robotController] 设置规划器: TVPPlanner" << std::endl;
}
template <int _Dofs>
void TVPPlanner<_Dofs>::calPlanQueue(bool isCoordinated, double deltaT, const Eigen::Matrix<double, _Dofs, 1> &maxVel,
                                     const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &jerk,
                                     const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                                     std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
{
    calTVPPlan(isCoordinated, deltaT, maxVel, maxAcc, q0, qf, dq0, q_d, dq_d, ddq_d);
}

// SS规划器
template <int _Dofs>
class SSPlanner : public Planner<_Dofs>
{
public:
public:
    SSPlanner(const SSPlanner &) = delete;
    void operator=(const SSPlanner &) = delete;

    SSPlanner();
    ~SSPlanner();

    void calPlanQueue(bool isCoordinated, double deltaT,
                      const Eigen::Matrix<double, _Dofs, 1> &maxVel, const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &maxJerk,
                      const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                      std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d);
};
template <int _Dofs>
SSPlanner<_Dofs>::~SSPlanner()
{
}
template <int _Dofs>
SSPlanner<_Dofs>::SSPlanner()
{
    std::cout << "[robotController] 设置规划器: SSPlanner" << std::endl;
}
template <int _Dofs>
void SSPlanner<_Dofs>::calPlanQueue(bool isCoordinated, double deltaT, const Eigen::Matrix<double, _Dofs, 1> &maxVel,
                                    const Eigen::Matrix<double, _Dofs, 1> &maxAcc, const Eigen::Matrix<double, _Dofs, 1> &maxJerk,
                                    const Eigen::Matrix<double, _Dofs, 1> &q0, const Eigen::Matrix<double, _Dofs, 1> &qf, const Eigen::Matrix<double, _Dofs, 1> &dq0,
                                    std::vector<std::queue<double>> &q_d, std::vector<std::queue<double>> &dq_d, std::vector<std::queue<double>> &ddq_d)
{
    calSSPlan(isCoordinated, deltaT, maxVel, maxAcc, maxJerk, q0, qf, dq0, q_d, dq_d, ddq_d);
}

// 工厂
template <int _Dofs>
bool newPlanner(std::unique_ptr<Planner<_Dofs>> &Planner, PlannerType plannerType)
{
    switch (plannerType)
    {
    case PlannerType::Quintic_:
        Planner = std::make_unique<QuinticPlanner<_Dofs>>();
        break;
    case PlannerType::TVP_:
        Planner = std::make_unique<TVPPlanner<_Dofs>>();
        break;
    case PlannerType::SS_:
        Planner = std::make_unique<SSPlanner<_Dofs>>();
        break;
    default:
        Planner = std::make_unique<TVPPlanner<_Dofs>>();
        return false;
        break;
    }
    return true;
}