//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "include/trajectory.h"

namespace mpc {

    Trajectory::Trajectory(int len, int state_size, int num_joints,
                           const std::vector<std::vector<double>>& switching_times, double node_dt) :
            inputs_(switching_times, num_joints, len, node_dt) {
        for (int i = 0; i < len; i++) {
            states_.push_back(vector_t::Zero(state_size));
        }

        for (const auto& switching_time : switching_times) {
            std::array<Spline, 3> end_effector_pos =
                    {Spline(2, switching_time, false), Spline(2, switching_time, false),
                     Spline(2, switching_time, false)};
            end_effector_pos_.emplace_back(end_effector_pos);
        }
    }

    std::vector<vector_t> Trajectory::GetStates() const {
        return states_;
    }

    const Inputs& Trajectory::GetInputs() const {
        return inputs_;
    }

    const std::vector<std::array<Spline, 3>>& Trajectory::GetPositions() const {
        return end_effector_pos_;
    }

    Eigen::Matrix<double, 3, 4> Trajectory::GetPositionsPolyVarsLin(int end_effector, double time) const {
        Eigen::Matrix<double, 3, 4> vars = Eigen::Matrix<double, 3, 4>::Zero();
        for (int coord = 0; coord < 3; coord++) {
            vars.row(coord) = end_effector_pos_.at(end_effector).at(coord).GetPolyVarsLin(time);
        }

        return vars;
    }

    void Trajectory::Reset() {
//        assert(states_.size() == inputs_.size() + 1);
//        for (int i = 0; i < states_.size(); i++) {
//            states_.at(i) = vector_t::Zero(states_.at(i).size());
//            if (i != states_.size() - 1) {
//                inputs_.at(i) = vector_t::Zero(inputs_.at(i).size());
//            }
//        }
    }

    void Trajectory::SetState(int idx, const vector_t& state) {
        assert(state.size() == states_.at(0).size());
        states_.at(idx) = state;
    }

    void Trajectory::SetInput(const Inputs& input) {
        inputs_ = input;
    }

} // mpc