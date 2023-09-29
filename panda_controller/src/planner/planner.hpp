#include <cmath>
#include <queue>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
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
