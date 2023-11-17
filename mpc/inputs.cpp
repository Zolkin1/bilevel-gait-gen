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
//            positions_.emplace_back(end_effector_pos);        NOTE: Removed positions from the inputs
        }

        for (int i = 0; i < num_nodes; i++) {
            joint_vels_.emplace_back(vector_t::Zero(num_joints));
        }

        force_spline_vars_ = 0;
        for (const auto& force : forces_) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                force_spline_vars_ += force.at(coord).GetTotalPolyVars();
            }
        }
    }

    Inputs::Inputs(const Inputs& input) {
        const std::vector<std::array<Spline, 3>>& force = input.GetForces();
        const std::vector<std::array<Spline, 3>>& pos = input.GetPositions();
        force_spline_vars_ = 0;
        for (int i = 0; i < force.size(); i++) {
            forces_.push_back(force.at(i));
            for (int coord = 0; coord < POS_VARS; coord++) {
                force_spline_vars_ += forces_.at(i).at(coord).GetTotalPolyVars();
            }
        }

        joint_vels_ = input.GetAllVels();

        node_dt_ = input.GetNodeDt();
    }

    Inputs Inputs::operator=(const mpc::Inputs &input) {
        const std::vector<std::array<Spline, 3>>& force = input.GetForces();
        const std::vector<std::array<Spline, 3>>& pos = input.GetPositions();

        forces_.erase(forces_.begin(), forces_.end());
        force_spline_vars_ = 0;
        for (int i = 0; i < force.size(); i++) {
            forces_.push_back(force.at(i));
            for (int coord = 0; coord < POS_VARS; coord++) {
                force_spline_vars_ += forces_.at(i).at(coord).GetTotalPolyVars();
            }
        }

        node_dt_ = input.GetNodeDt();

        return *this;
//        Inputs inp(input);
//        return inp;
    }

    Eigen::Vector3d Inputs::GetForce(int end_effector, double time) const {
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for (int i = 0; i < 3; i++) {
            force(i) = forces_.at(end_effector).at(i).ValueAt(time);
        }
        return force;
    }

    Eigen::Vector3d Inputs::GetPosition(int end_effector, double time) const {
        Eigen::Vector3d position = Eigen::Vector3d::Zero();
        for (int i = 0; i < 3; i++) {
            position(i) = positions_.at(end_effector).at(i).ValueAt(time);
        }

        return position;
    }

    vector_t Inputs::GetVels(double time) const {
        assert(time >= 0);
        assert(time <= joint_vels_.size()*node_dt_);

        int idx = floor(time/node_dt_); // TODO: Check
        if (idx == joint_vels_.size()) {
            idx -= 1;
        }
        return joint_vels_.at(idx);
    }

    const std::vector<std::array<Spline, 3>>& Inputs::GetForces() const {
        return forces_;
    }

    const std::vector<std::array<Spline, 3>>& Inputs::GetPositions() const {
        return positions_;
    }

    double Inputs::GetNodeDt() const {
        return node_dt_;
    }

    // TODO: Check that its valid
    void Inputs::SetEndEffectorForce(int end_effector, const std::array<Spline, 3>& force) {
        forces_.at(end_effector) = force;
    }

    void Inputs::SetEndEffectorPosition(int end_effector, const std::array<Spline, 3> &position) {
        positions_.at(end_effector) = position;
    }

    void Inputs::SetJointVels(const mpc::vector_t &vels, double time) {
        int idx = floor(time/node_dt_); // TODO: Check
        joint_vels_.at(idx) = vels;
    }

    void Inputs::SetJointVelsNoTime(const mpc::vector_t& vels, int index) {
        joint_vels_.at(index) = vels;
    }

    int Inputs::GetNumInputs() const {
        return 3*forces_.size()*forces_.at(0).at(0).GetTotalPolyVars() + joint_vels_.at(0).size();
    }

    int Inputs::GetNumForces() const {
        return 3*forces_.size();
    }

    int Inputs::GetNumPositions() const {
        return 3*positions_.size();
    }

    vector_t Inputs::GetInputVector(double time) const {
        const int POLY_ORDER = 4;
        vector_t input = vector_t::Zero(GetNumInputs());

        int poly_vars_per_spline = forces_.at(0).at(0).GetTotalPolyVars();

        vector_t force_vec = vector_t::Zero(GetNumForces()*poly_vars_per_spline);
        int j = 0;
        for (const auto & force : forces_) {
            vector_t fx = force.at(0).GetAllPolyVars();
            vector_t fy = force.at(1).GetAllPolyVars();
            vector_t fz = force.at(2).GetAllPolyVars();

            force_vec.segment(j*poly_vars_per_spline, poly_vars_per_spline) = fx;
            force_vec.segment((j+1)*poly_vars_per_spline, poly_vars_per_spline) = fy;
            force_vec.segment((j+2)*poly_vars_per_spline, poly_vars_per_spline) = fz;

            j += 3;
        }
        input << force_vec, GetVels(time);

        return input;
    }

    std::vector<vector_t> Inputs::GetAllVels() const {
        return joint_vels_;
    }

    int Inputs::GetForcePolyIdx(double time) const {
        return forces_.at(0).at(0).GetPolyIdx(time);
    }

//    void Inputs::UpdateForce(int end_effector, int coord,
//                             const std::vector<std::array<double, Spline::POLY_ORDER>>& vars) {
//        forces_.at(end_effector).at(coord).SetAllSplineVars(vars);
//    }

    // Always returns the index of the end of the end point of the polynomial
    std::pair<int, int> Inputs::GetForceSplineIndex(int end_effector, double time, int coord) const {
        int num_spline_vars_before = 0;
        for (int ee = 1; ee < end_effector; ee++) {
            for (int j = 0; j < POS_VARS; j++) {
                num_spline_vars_before += forces_.at(ee).at(j).GetTotalPolyVars();
            }
        }

        int idx_into_ee_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            idx_into_ee_coord_spline_vars += forces_.at(end_effector).at(j).GetTotalPolyVars();
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = forces_.at(end_effector).at(coord).GetVarsIndexEnd(time);

        return std::make_pair(num_spline_vars_before + idx_into_ee_coord_spline_vars +
                              vars_idx, vars_affecting);
    }

    int Inputs::GetTotalForceSplineVars() const {
        return force_spline_vars_;
    }

    void Inputs::UpdateForcePoly(int end_effector, int coord, int idx, const vector_t& vars) {
        std::vector<double> temp(vars.size());

        // TODO: how to convert eigen to std::vector
        for (int i = 0; i < vars.size(); i++) {
            temp.at(i) = vars(i);
        }

        forces_.at(end_effector).at(coord).UpdatePolyVar(idx, temp);
    }
}