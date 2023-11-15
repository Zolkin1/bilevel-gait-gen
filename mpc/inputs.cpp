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
    }

    Inputs::Inputs(const Inputs& input) {
        const std::vector<std::array<Spline, 3>>& force = input.GetForces();
        const std::vector<std::array<Spline, 3>>& pos = input.GetPositions();
        for (int i = 0; i < force.size(); i++) {
            forces_.push_back(force.at(i));
//            positions_.push_back(pos.at(i));
        }

        joint_vels_ = input.GetAllVels();

        node_dt_ = input.GetNodeDt();
    }

    Inputs Inputs::operator=(const mpc::Inputs &input) {
        const std::vector<std::array<Spline, 3>>& force = input.GetForces();
        const std::vector<std::array<Spline, 3>>& pos = input.GetPositions();

        forces_.erase(forces_.begin(), forces_.end());
//        positions_.erase(positions_.begin(), positions_.end());

        for (int i = 0; i < force.size(); i++) {
            forces_.push_back(force.at(i));
//            positions_.push_back(pos.at(i));
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

    Eigen::Matrix<double, 3, 4> Inputs::GetForcePolyVarsLin(int end_effector, double time) const {
        Eigen::Matrix<double, 3, 4> vars = Eigen::Matrix<double, 3, 4>::Zero();
        for (int coord = 0; coord < 3; coord++) {
            vars.row(coord) = forces_.at(end_effector).at(coord).GetPolyVarsLin(time);
        }

        return vars;
    }

    int Inputs::GetForcePolyIdx(double time) const {
        return forces_.at(0).at(0).GetPolyIdx(time);
    }

}