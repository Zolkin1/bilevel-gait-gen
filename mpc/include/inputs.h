//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_INPUTS_H
#define BILEVEL_GAIT_GEN_INPUTS_H

#include <Eigen/Core>

#include "spline.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;

    class Inputs {
    public:
        Inputs(const std::vector<std::vector<double>>& switching_times, int num_joints, int num_nodes, double node_dt);

        Eigen::Vector3d GetForce(int end_effector, double time) const;

        vector_t GetVels(double time) const;
    protected:
    private:
        std::vector<std::array<Spline, 3>> forces_;    // We need a spline for each coordinate of each end effector
        std::vector<std::array<Spline, 3>> positions_;
        std::vector<vector_t> joint_vels_;      // These are discrete, ZOH velocities.
        double node_dt_;

        static int constexpr NUM_POLY = 2;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_INPUTS_H
