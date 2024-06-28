#include "jogPlanner.h"
#include <cmath>
#include <iostream>

void calStopTAP(double vel, double maxAcc, double maxJerk, double &t1, double &t2, double &t3)
{
    t1 = maxAcc / maxJerk; // j/a
    t2 = (std::fabs(vel) - t1 * t1 * maxJerk) / maxAcc;
    if (t2 < 0)
    {
        double vj = std::fabs(vel) / maxJerk;
        if (vj < 0)
            vj = 0;
        t2 = 0;
        t1 = std::sqrt(vj);
    }
    t3 = t1;
}
bool calJogMovePlan(bool flag, double deltaT, int dir,
                    double maxJerk, double maxAcc, double maxVel, // vel是带百分比系数的
                    double &q_d, double &dq_d, double &ddq_d)
{
    static double jogTime = 0;
    static double t1 = 0;
    static double t2 = 0;
    static double t3 = 0;
    if (flag)
    {
        jogTime = 0;
        calStopTAP(maxVel, maxAcc, maxJerk, t1, t2, t3);
    }
    jogTime = jogTime + deltaT;

    if (jogTime <= t1)
    {
        ddq_d = ddq_d + dir * maxJerk * deltaT;
    }
    else if (jogTime > t1 && jogTime <= t1 + t2)
    {
        // 不操作
    }
    else if (jogTime > t1 + t2 && jogTime <= t1 + t2 + t3 - deltaT)
    {
        ddq_d = ddq_d - dir * maxJerk * deltaT;
    }
    else
    {
        ddq_d = 0;
    }
    dq_d = dq_d + ddq_d * deltaT;
    q_d = q_d + dq_d * deltaT;
    return true;
}

bool calJogStopPlan(bool flag, double deltaT, double dq0, double ddq0,
                    double maxJerk, double maxAcc, // vel是带百分比系数的
                    double &q_d, double &dq_d, double &ddq_d)
{
    static double jogTime = 0;
    static bool AccFlag = false;
    static double t0 = 0;
    static double t1 = 0;
    static double t2 = 0;
    static double t3 = 0;
    static int sign = 0;
    if (flag)
    {
        jogTime = 0;
        sign = (dq_d > 0) ? 1 : -1;
        if (ddq0 == 0)
        {
            t0 = 0;
            calStopTAP(dq0, maxAcc, maxJerk, t1, t2, t3);
            AccFlag = false;
        }
        else
        {
            t0 = std::fabs(ddq0 / maxJerk);
            AccFlag = true;
        }
    }
    jogTime = jogTime + deltaT;
    if (AccFlag && jogTime > t0)
    {
        AccFlag = false;
        calStopTAP(dq_d, maxAcc, maxJerk, t1, t2, t3);
    }

    if (jogTime <= t0)
    {
        ddq_d = ddq_d - sign * maxJerk * deltaT;
    }
    else if (jogTime > t0 && jogTime <= t0 + t1)
    {
        ddq_d = ddq_d - sign * maxJerk * deltaT;
    }
    else if (jogTime > t0 + t1 && jogTime <= t0 + t1 + t2)
    {
        // 不操作
    }
    else if (jogTime > t0 + t1 + t2 && jogTime <= t0 + t1 + t2 + t3 - deltaT)
    {
        ddq_d = ddq_d + sign * maxJerk * deltaT;
    }
    else
    {
        ddq_d = 0;
        dq_d = 0;
    }
    dq_d = dq_d + ddq_d * deltaT;
    q_d = q_d + dq_d * deltaT;
    return true;
}