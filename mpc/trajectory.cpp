//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <fstream>

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
                    {Spline(3, switching_time, false), Spline(3, switching_time, false),
                     Spline(3, switching_time, false)};
            end_effector_pos_.emplace_back(end_effector_pos);
        }

        UpdatePosSplineVarsCount();
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

//    void Trajectory::UpdatePosition(int end_effector, int coord,
//                                    const std::vector<std::array<double, Spline::POLY_ORDER>>& vars) {
//        end_effector_pos_.at(end_effector).at(coord).SetAllSplineVars(vars);
//    }
//
//    void Trajectory::UpdateForce(int end_effector, int coord,
//                                 const std::vector<std::array<double, Spline::POLY_ORDER>>& vars) {
//        inputs_.UpdateForce(end_effector, coord, vars);
//    }

    void Trajectory::SetState(int idx, const vector_t& state) {
        assert(state.size() == states_.at(idx).size());
        states_.at(idx) = state;
    }

    void Trajectory::SetInput(const Inputs& input) {
        inputs_ = input;
    }

    void Trajectory::SetInputVels(int idx, const vector_t& joint_vels) {
        inputs_.SetJointVelsNoTime(joint_vels, idx);
    }

    int Trajectory::GetTotalPosSplineVars() const {
        return pos_spline_vars_;
    }

    void Trajectory::UpdateForceSpline(int end_effector, int coord, const vector_t& vars) {
        int idx = 0;
        for (int i = 0; i < inputs_.GetForces().at(end_effector).at(coord).GetNumPolyTimes(); i++) {
            int num_vars = inputs_.GetForces().at(end_effector).at(coord).GetNumPolyVars(i);
            inputs_.UpdateForcePoly(end_effector, coord, i, vars.segment(idx, num_vars));
            if (num_vars == 1 && i+1 < end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes() &&
            inputs_.GetForces().at(end_effector).at(coord).GetNumPolyVars(i+1) == 1) {
                // Set then skip the additional constant term
                // TODO: In theory we can skip this
                inputs_.UpdateForcePoly(end_effector, coord, i+1, vars.segment(idx, num_vars));
                i++;
            }
            idx += num_vars;
        }
    }

    void Trajectory::UpdatePositionSpline(int end_effector, int coord, const vector_t& vars) {
        int idx = 0;
        for (int i = 0; i < end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes(); i++) {
            int num_vars = end_effector_pos_.at(end_effector).at(coord).GetNumPolyVars(i);

            std::vector<double> temp(num_vars);
            for (int j = 0; j < num_vars; j++) {
                temp.at(j) = vars.segment(idx, num_vars)(j);
            }
            end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i, temp);

            if (num_vars == 1 && i+1 < end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes() &&
            end_effector_pos_.at(end_effector).at(coord).GetNumPolyVars(i+1) == 1) {
                // Set then skip the additional constant term
                end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i+1, temp);
                i++;
            }
            idx += num_vars;
        }
    }

    std::pair<int, int> Trajectory::GetPositionSplineIndex(int end_effector, double time, int coord) const {
        int num_spline_vars_before = 0;
        for (int ee = 0; ee < end_effector; ee++) {     // TODO: start at 1 or 0?
            for (int j = 0; j < POS_VARS; j++) {
                num_spline_vars_before += end_effector_pos_.at(ee).at(j).GetTotalPolyVars();
            }
        }

        int idx_into_ee_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            idx_into_ee_coord_spline_vars += end_effector_pos_.at(end_effector).at(j).GetTotalPolyVars();
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = end_effector_pos_.at(end_effector).at(coord).GetVarsIndexEnd(time);

        return std::make_pair(num_spline_vars_before + idx_into_ee_coord_spline_vars + vars_idx, vars_affecting);
    }

    void Trajectory::SetPositionsForAllTime(int ee, const std::array<double, POS_VARS>& ee_pos) {
        for (int coord = 0; coord < POS_VARS; coord++) {
            end_effector_pos_.at(ee).at(coord).SetAllPositions(ee_pos.at(coord));
        }
    }

    void Trajectory::PrintTrajectoryToFile(const std::string& file_name) const {
        std::ofstream file;
        file.open(file_name);

        matrix_t states = matrix_t::Zero(states_.size(), states_.at(0).size());
        for (int i = 0; i < states_.size(); i++) {
            states.row(i) = states_.at(i);
        }

        file << "states: " << std::endl;
        file << states << std::endl;

//        file << "node: [states], [joint vels]" << std::endl;
//        for (int node = 0; node < states_.size(); node++) {
//            if (node == 0) {
//                file << states_.at(node).transpose() << std::endl;
//            } else {
//                file << states_.at(node).transpose() << ", " << inputs_.GetVel(node-1).transpose() << std::endl;
//            }
//        }

        file << "force spline: " << std::endl;
        for (const auto& force : inputs_.GetForces()) {
            for (int coord = 0; coord < 3; coord++) {
                for (int i = 0; i < force.at(coord).GetPolyVars().size(); i++) {
                    for (int j = 0; j < force.at(coord).GetPolyVars().at(i).size(); j++) {
                        file << force.at(coord).GetPolyVars().at(i).at(j);
                        if (j != force.at(coord).GetPolyVars().at(i).size()-1) {
                            file << " ";
                        }
                    }
                    file << ", ";
                }
                file << std::endl;
            }
        }

        file << "position spline: " << std::endl;
        for (const auto& pos : end_effector_pos_) {
            for (int coord = 0; coord < 3; coord++) {
                for (int i = 0; i < pos.at(coord).GetPolyVars().size(); i++) {
                    for (int j = 0; j < pos.at(coord).GetPolyVars().at(i).size(); j++) {
                        file << pos.at(coord).GetPolyVars().at(i).at(j);
                        if (j != pos.at(coord).GetPolyVars().at(i).size()-1) {
                            file << " ";
                        }
                    }
                    file << ", ";
                }
                file << std::endl;
            }
        }

        matrix_t joint_vels = matrix_t::Zero(inputs_.GetAllVels().size(), inputs_.GetAllVels().at(0).size());
        for (int i = 0; i < inputs_.GetAllVels().size(); i++) {
            joint_vels.row(i) = inputs_.GetVel(i).transpose();
        }

        file << "timings: " << std::endl;
        for (int ee = 0; ee < end_effector_pos_.size(); ee++) {
            file << "end effector #" << ee << ": " << std::endl;
            const auto& times = end_effector_pos_.at(ee).at(0).GetPolyTimes();
            const auto& forces = inputs_.GetForces().at(ee).at(0).GetPolyTimes();
            for (double time : times) {
                file << time << ", ";
            }
            file << std::endl;

            // Sanity check that they are the same
            for (int i = 0; i < times.size(); i++) {
                assert(times.at(i) == forces.at(i));
            }
        }

        file << "input joint vels: " << std::endl;
        file << joint_vels << std::endl;

        file.close();
    }

    int Trajectory::GetTotalPosConstants() const {
        int num_constants = 0;
        for (const auto& pos: end_effector_pos_) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                num_constants += pos.at(coord).GetNumConstant();
            }
        }

        return num_constants;
    }

    double Trajectory::GetTotalTime() const {
        assert(end_effector_pos_.at(0).at(0).GetTotalTime() == inputs_.GetForces().at(0).at(0).GetTotalTime());
        return end_effector_pos_.at(0).at(0).GetTotalTime();
    }

    void Trajectory::AddPolys(double final_time) {
        if (GetTotalTime() < final_time) {
            for (auto &end_effector_po: end_effector_pos_) {
                for (int coord = 0; coord < POS_VARS; coord++) {
                    end_effector_po.at(coord).AddPoly(0.5);
                }
            }
            inputs_.AddPolys(0.5);
        }

        UpdatePosSplineVarsCount();
    }

    void Trajectory::RemoveUnusedPolys(double init_time) {
        for (auto &end_effector_po: end_effector_pos_) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                end_effector_po.at(coord).RemoveUnused(init_time);
            }
        }

        inputs_.RemoveUnusedPolys(init_time);

        UpdatePosSplineVarsCount();
    }

    void Trajectory::SetInitTime(double time) {
        inputs_.SetInitTime(time);
    }

    void Trajectory::UpdatePosSplineVarsCount() {
        pos_spline_vars_ = 0;
        for (const auto& end_effector_po : end_effector_pos_) {
            for (int coord = 0; coord < 3; coord++) {
                pos_spline_vars_ += end_effector_po.at(coord).GetTotalPolyVars();
            }
        }
    }

    void Trajectory::SetEndEffectorSplines(int ee, const Spline& spline) {
        for (int coord = 0; coord < POS_VARS; coord++) {
            end_effector_pos_.at(ee).at(coord) = spline;
            inputs_.SetForceSpline(ee, coord, spline);
        }
    }

} // mpc