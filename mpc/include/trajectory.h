//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_TRAJECTORY_H
#define BILEVEL_GAIT_GEN_TRAJECTORY_H

#include <Eigen/Core>

#include "inputs.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;

    class Trajectory {
    public:
        Trajectory(int len, int state_size, int num_joints,
                   const std::vector<std::vector<double>>& switching_times, double node_dt);

//        Trajectory(const Trajectory& traj);

        std::vector<vector_t> GetStates() const;

        const Inputs& GetInputs() const;

        const std::vector<std::array<Spline, 3>>& GetPositions() const;

        Eigen::Matrix<double, 3, 4> GetPositionsPolyVarsLin(int end_effector, double time) const;

        /**
         * Resets all states_ and inputs to 0
         */
        void Reset();

        void SetState(int idx, const vector_t& state);
        void SetInput(const Inputs& input);
    protected:
    private:
        std::vector<vector_t> states_;
        Inputs inputs_;
        std::vector<std::array<Spline, 3>> end_effector_pos_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_TRAJECTORY_H
