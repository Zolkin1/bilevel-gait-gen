//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "include/trajectory.h"

namespace mpc {

    Trajectory::Trajectory(int len, int state_size, int input_size) {
        for (int i = 0; i < len; i++) {
            states_.push_back(vector_t::Zero(state_size));
            if (i != len - 1) {
                inputs_.push_back(vector_t::Zero(input_size));
            }
        }
    }

    // copy constructor
    Trajectory::Trajectory(const mpc::Trajectory &traj) {
        std::vector<vector_t> states = traj.GetStates();
        std::vector<vector_t> inputs = traj.GetInputs();
        for (int i = 0; i < states.size(); i++) {
            states_.push_back(states.at(i));
            if (i != states.size() - 1) {
                inputs_.push_back(inputs.at(i));
            }
        }
    }

    std::vector<vector_t> Trajectory::GetStates() const {
        return states_;
    }

    std::vector<vector_t> Trajectory::GetInputs() const {
        return inputs_;
    }

    void Trajectory::Reset() {
        assert(states_.size() == inputs_.size() + 1);
        for (int i = 0; i < states_.size(); i++) {
            states_.at(i) = vector_t::Zero(states_.at(i).size());
            if (i != states_.size() - 1) {
                inputs_.at(i) = vector_t::Zero(inputs_.at(i).size());
            }
        }
    }

    void Trajectory::SetState(int idx, const vector_t& state) {
        assert(state.size() == states_.at(0).size());
        states_.at(idx) = state;
    }

    void Trajectory::SetInput(int idx, const vector_t& input) {
        assert(input.size() == inputs_.at(0).size());
        inputs_.at(idx) = input;
    }

} // mpc