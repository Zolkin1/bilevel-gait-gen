//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "inputs.h"

namespace mpc {
    Inputs::Inputs(const std::vector<std::vector<double>>& switching_times, int num_joints, int num_nodes, double node_dt) :
            node_dt_(node_dt) {
        for (auto & switching_time : switching_times) {
            // TODO: Clean up
            std::array<Spline, 3> end_effector_force =
                    {Spline(2, switching_time, true), Spline(2, switching_time, true), Spline(2, switching_time, true)};
            std::array<Spline, 3> end_effector_pos =
                    {Spline(2, switching_time, false), Spline(2, switching_time, false), Spline(2, switching_time, false)};
            forces_.emplace_back(end_effector_force);
            positions_.emplace_back(end_effector_pos);
        }

        for (int i = 0; i < num_nodes; i++) {
            joint_vels_.emplace_back(vector_t::Zero(num_joints));
        }
    }

    Eigen::Vector3d Inputs::GetForce(int end_effector, double time) const {
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for (int i = 0; i < 3; i++) {
            force(i) = forces_.at(end_effector).at(i).ValueAt(time);
        }
        return force;
    }

    vector_t Inputs::GetVels(double time) const {
        int idx = floor(time/node_dt_); // TODO: Check
        return joint_vels_.at(idx);
    }
}