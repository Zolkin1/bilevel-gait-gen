//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_TRAJECTORY_H
#define BILEVEL_GAIT_GEN_TRAJECTORY_H

#include <Eigen/Core>
#include "mpc.h"

namespace mpc {
    class Trajectory {
    public:
        Trajectory(int len, int state_size, int input_size);

        Trajectory(const Trajectory& traj);

        std::vector<vector_t> GetStates() const;

        std::vector<vector_t> GetInputs() const;

        /**
         * Resets all states_ and inputs to 0
         */
        void Reset();

        void SetState(int idx, const vector_t& state);
        void SetInput(int idx, const vector_t& input);
    protected:
    private:
        std::vector<vector_t> states_;
        std::vector<vector_t> inputs_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_TRAJECTORY_H
