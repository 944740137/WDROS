bool calJogMovePlan(bool flag, double deltaT, int dir,
                    double maxJerk, double maxAcc, double maxVel, // vel是带百分比系数的
                    double &q_d, double &dq_d, double &ddq_d);
bool calJogStopPlan(bool flag, double deltaT, double maxJerk, double maxAcc, // vel是带百分比系数的
                    double &q_d, double &dq_d, double &ddq_d);