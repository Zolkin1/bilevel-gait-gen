//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <fstream>

#include "trajectory.h"
#include "models/model.h"

namespace mpc {
    Trajectory::Trajectory(int len, int state_size, bool using_joints,
                           const std::vector<std::vector<double>>& switching_times, double node_dt,
                           double swing_height, double foot_offset) :
            swing_height_(swing_height),
            foot_offset_(foot_offset),
            fk_traj_(5),
            using_joints_(using_joints),
            node_dt_(node_dt) { // TODO: Make not hard coded
        for (int i = 0; i < len; i++) {
            states_.push_back(vector_t::Zero(state_size));
//            full_velocity_.push_back(vector_t::Zero(state_size));
        }

        for (int i = 0; i < switching_times.size(); i++) {
            bool in_contact = false;
            if (i == 1 || i == 2) {
                in_contact = true;
            }
            assert(switching_times.at(i).at(0) == 0);

            ee_splines_.emplace_back(switching_times.at(i).size(), switching_times.at(i),
                                     in_contact, 3); // TODO: Make this not hard coded
        }

        for (int i = 0; i < full_config_.size(); i++) {
            full_config_.at(i) = vector_t::Zero(state_size);
            full_velocity_.at(i) = vector_t::Zero(state_size-1);
        }

        UpdateSplineVarsCount();
        SetSwingPosZ();

        for (auto& ee_traj : fk_traj_) {
            ee_traj.resize(len);
        }
    }

    Trajectory& Trajectory::operator=(const Trajectory& traj) {
        if (this == &traj) {
            return *this;
        }

        this->pos_spline_vars_ = traj.pos_spline_vars_;
        this->force_spline_vars_ = traj.force_spline_vars_;
        this->swing_height_ = traj.swing_height_;
        this->foot_offset_ = traj.foot_offset_;
        this->states_ = traj.states_;
        this->fk_traj_ = traj.fk_traj_;
        this->node_dt_ = traj.node_dt_;
        this->full_config_ = traj.full_config_;
        this->full_velocity_ = traj.full_velocity_;
        this->init_time_ = traj.init_time_;
        this->ee_splines_ = traj.ee_splines_;

        return *this;
    }

    std::vector<vector_t> Trajectory::GetStates() const {
        return states_;
    }

    void Trajectory::SetState(int idx, const vector_t& state) {
        assert(state.size() == states_.at(idx).size());
        states_.at(idx) = state;
    }

    int Trajectory::GetTotalPosSplineVars() const {
        return pos_spline_vars_;
    }

    void Trajectory::UpdateForceSpline(int end_effector, int coord, const vector_t& vars) {
        int idx = 0;

        const std::vector<int> mut_nodes = ee_splines_.at(end_effector).GetMutableNodes(EndEffectorSplines::Force, coord);
        assert(vars.size() == ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Force, coord));
        for(auto& node : mut_nodes) {
            vector_2t temp;
            for (int j = 0; j < 2; j++) {
                temp(j) = vars(idx + j);
            }

            ee_splines_.at(end_effector).SetVars(EndEffectorSplines::Force, coord, node, temp);
            idx += 2;
        }
    }

    void Trajectory::UpdatePositionSpline(int end_effector, int coord, const vector_t& vars) {
        int idx = 0;

        const std::vector<int> mut_nodes = ee_splines_.at(end_effector).GetMutableNodes(EndEffectorSplines::Position, coord);
        assert(mut_nodes.size() == vars.size());
        for(auto& node : mut_nodes) {
            vector_2t temp;
            temp << vars(idx), 0;

            ee_splines_.at(end_effector).SetVars(EndEffectorSplines::Position, coord, node, temp);
            idx++;
        }
    }

    std::pair<int, int> Trajectory::GetPositionSplineIndex(int end_effector, double time, int coord) const {
        if (coord == 2) {
            throw std::runtime_error("The chosen spline is not mutable and thus does not provide a index.");
        }

        int num_spline_vars_before = 0;
        for (int ee = 0; ee < end_effector; ee++) {
            num_spline_vars_before += 2*ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Position, coord);
        }

        int idx_into_ee_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            idx_into_ee_coord_spline_vars += ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Position, coord);
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = ee_splines_.at(end_effector).GetVarsIdx(
                EndEffectorSplines::Position, coord, time);

        return std::make_pair(num_spline_vars_before + idx_into_ee_coord_spline_vars + vars_idx, vars_affecting);
    }

//    void Trajectory::SetPositionsForAllTime(int ee, const vector_3t& ee_pos) {
//        for (int coord = 0; coord < 2; coord++) {
//            const std::vector<int> mut_nodes = ee_splines_.at(ee).GetMutableNodes(EndEffectorSplines::Position, coord);
//            for (auto& node : mut_nodes) {
//                vector_2t vars = {ee_pos(coord), 0};
//                ee_splines_.at(ee).SetVars(EndEffectorSplines::Position, coord, node, vars);
//            }
//        }
//    }

    // TODO: Fix this
    void Trajectory::PrintTrajectoryToFile(const std::string& file_name) const {
        std::ofstream file;
        file.open(file_name);

        matrix_t states = matrix_t::Zero(states_.size(), states_.at(0).size());
        for (int i = 0; i < states_.size(); i++) {
            states.row(i) = states_.at(i);
        }

        file << "states: " << std::endl;
        file << states << std::endl;

        file << "force spline: " << std::endl;
//        for (const auto& force : forces_) {
//            for (int coord = 0; coord < 3; coord++) {
//                for (int i = 0; i < force.at(coord).GetPolyVars().size(); i++) {
//                    for (int j = 0; j < force.at(coord).GetPolyVars().at(i).size(); j++) {
//                        file << force.at(coord).GetPolyVars().at(i).at(j);
//                        if (j != force.at(coord).GetPolyVars().at(i).size()-1) {
//                            file << " ";
//                        }
//                    }
//                    file << ", ";
//                }
//                file << std::endl;
//            }
//        }

        file << "position spline: " << std::endl;
//        for (const auto& pos : end_effector_pos_) {
//            for (int coord = 0; coord < 3; coord++) {
//                for (int i = 0; i < pos.at(coord).GetPolyVars().size(); i++) {
//                    for (int j = 0; j < pos.at(coord).GetPolyVars().at(i).size(); j++) {
//                        file << pos.at(coord).GetPolyVars().at(i).at(j);
//                        if (j != pos.at(coord).GetPolyVars().at(i).size()-1) {
//                            file << " ";
//                        }
//                    }
//                    file << ", ";
//                }
//                file << std::endl;
//            }
//        }

        if (using_joints_) {
//            matrix_t joint_vels = matrix_t::Zero(inputs_.GetAllVels().size(), inputs_.GetAllVels().at(0).size());
//            for (int i = 0; i < inputs_.GetAllVels().size(); i++) {
//                joint_vels.row(i) = inputs_.GetVel(i).transpose();
//            }
        }

        file << "timings: " << std::endl;
//        for (int ee = 0; ee < end_effector_pos_.size(); ee++) {
//            file << "end effector #" << ee << ": " << std::endl;
//            const auto& pos_times = end_effector_pos_.at(ee).at(0).GetPolyTimes();
//            const auto& force_times = forces_.at(ee).at(0).GetPolyTimes();
//            file << "force: ";
//            for (double force_time : force_times) {
//                file << force_time << ", ";
//            }
//            file << std::endl;
//
//            file << "positions: ";
//            for (double pos_time : pos_times) {
//                file << pos_time << ", ";
//            }
//            file << std::endl;
//        }

        if (using_joints_) {
//            file << "input joint vels: " << std::endl;
//            file << joint_vels << std::endl;
        }

        file << "spline vec: \n" << SplinesAsVec() << std::endl;

        file.close();
    }

    void Trajectory::AddPolys(double final_time) {
        for (auto & ee_spline : ee_splines_) {
            while (ee_spline.GetEndTime() < final_time) {
                ee_spline.AddPoly(0.2);    // TODO: Make not hard coded, 0.2
            }
        }
        SetSwingPosZ();
        UpdateSplineVarsCount();
    }

    void Trajectory::RemoveUnusedPolys(double init_time) {
        for (auto & ee_spline : ee_splines_) {
            ee_spline.RemovePoly(init_time);
        }

        UpdateSplineVarsCount();
    }

    void Trajectory::SetInitTime(double time) {
        init_time_ = time;
    }

    void Trajectory::UpdateSplineVarsCount() {
        pos_spline_vars_ = 0;
        force_spline_vars_ = 0;
        for (const auto& ee_spline : ee_splines_) {
            // Note the coordinates don't matter here because we don't sure the z position
            pos_spline_vars_ += 2*ee_spline.GetTotalPolyVars(EndEffectorSplines::Position, 0);
            force_spline_vars_ += 3*ee_spline.GetTotalPolyVars(EndEffectorSplines::Force, 0);

            assert(ee_spline.GetTotalPolyVars(EndEffectorSplines::Position, 0) == ee_spline.GetTotalPolyVars(EndEffectorSplines::Position, 1));
            assert(ee_spline.GetTotalPolyVars(EndEffectorSplines::Force, 0) == ee_spline.GetTotalPolyVars(EndEffectorSplines::Force, 1));
            assert(ee_spline.GetTotalPolyVars(EndEffectorSplines::Force, 1) == ee_spline.GetTotalPolyVars(EndEffectorSplines::Force, 2));
        }

    }

//    vector_t Trajectory::PositionAsQPVector() const {
//        vector_t pos_vec = vector_t::Zero(pos_spline_vars_);
//        int idx = 0;
//        for (int ee = 0; ee < end_effector_pos_.size(); ee++) {
//            for (int coord = 0; coord < POS_VARS; coord++) {
//                if (mut_flags_.at(ee).at(coord)) {
//                    pos_vec.segment(idx, end_effector_pos_.at(ee).at(coord).GetTotalPolyVars()) =
//                            end_effector_pos_.at(ee).at(coord).GetAllPolyVars();
//                    idx += end_effector_pos_.at(ee).at(coord).GetTotalPolyVars();
//                }
//            }
//        }
//
//        return pos_vec;
//    }

    vector_t Trajectory::GetState(int node) const {
        return states_.at(node);
    }

    std::vector<bool> Trajectory::GetContacts(double time) const {
        std::vector<bool> in_contact(ee_splines_.size());
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = ee_splines_.at(ee).GetVarsIdx(
                    EndEffectorSplines::Position, 0, time);
            if (vars_affecting == 1) {
                in_contact.at(ee) = true;
            } else {
                in_contact.at(ee) = false;
            }
        }

        return in_contact;
    }

    void Trajectory::SetSwingPosZ() {
        const int coord = 2;
        for (auto& ee_spline : ee_splines_) {
            const std::vector<int> mut_nodes = ee_spline.GetMutableNodes(EndEffectorSplines::Position, coord);
            for (auto& node: mut_nodes) {
                if (ee_spline.GetNodeType(EndEffectorSplines::Position, coord, node) == NodeType::FullDeriv) {
                    ee_spline.SetVars(EndEffectorSplines::Position, coord,
                                      node, {swing_height_, 0});
                } else {
                    ee_spline.SetVars(EndEffectorSplines::Position, coord,
                                      node, {foot_offset_, 0});
                }
            }
        }
    }

    void Trajectory::UpdateFullVelocity(int node, const vector_t& vel) {
        full_velocity_.at(node) = vel;
    }

    void Trajectory::UpdateFullConfig(int node, const mpc::vector_t& q) {
        full_config_.at(node) = q;
    }

    vector_t Trajectory::GetFullVelocity(int node) const {
        return full_velocity_.at(node);
    }

    vector_t Trajectory::GetFullConfig(int node) const {
        return full_config_.at(node);
    }

    int Trajectory::GetNumContactNodes(int ee) const {
        return ee_splines_.at(ee).GetNumContacts();
    }

    std::vector<time_v> Trajectory::GetContactTimes() const {
        std::vector<time_v> contact_times(ee_splines_.size());
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            contact_times.at(ee) = ee_splines_.at(ee).GetContactTimes();
        }

        return contact_times;
    }

    vector_t Trajectory::GetSplineLin(const mpc::Trajectory::SplineTypes& spline_type, int ee, int coord, double time) const {
        // TODO: Change this function to be two or change the spline type function
        switch (spline_type) {
            case SplineTypes::Force:
                return ee_splines_.at(ee).GetPolyVarsLin(EndEffectorSplines::Force, coord, time);
            case SplineTypes::Position:
                if (coord == 2) {
                    throw std::runtime_error("You cannot request spline linearizations for position z axis.");
                }
                return ee_splines_.at(ee).GetPolyVarsLin(EndEffectorSplines::Position, coord, time);
            default:
                throw std::runtime_error("Spline type not implemented.");
        }
    }

    std::pair<int, int> Trajectory::GetForceSplineIndex(int end_effector, double time, int coord) const {
        int num_spline_vars_before = 0;
        for (int ee = 0; ee < end_effector; ee++) {
            num_spline_vars_before += 3*ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Force, coord);
        }

        int idx_into_force_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            idx_into_force_coord_spline_vars += ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Force, coord);
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = ee_splines_.at(end_effector).GetVarsIdx(EndEffectorSplines::Force, coord, time);

        return std::make_pair(num_spline_vars_before + idx_into_force_coord_spline_vars + vars_idx, vars_affecting);
    }

    int Trajectory::GetTotalPolyVars(const mpc::Trajectory::SplineTypes& spline_type, int end_effector, int coord) {
        // TODO: Change this switch statement/function to be better
        switch (spline_type) {
            case SplineTypes::Force:
                return ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Force, coord);
            case SplineTypes::Position:
                if (coord == 2) {
                    throw std::runtime_error("You cannot request spline poly vars count for position z axis.");
                }
                return ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Position, coord);
            default:
                throw std::runtime_error("Spline type not implemented.");
        }
    }

    Eigen::Vector3d Trajectory::GetForce(int end_effector, double time) const {
        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        for (int coord = 0; coord < POS_VARS; coord++) {
            force(coord) = ee_splines_.at(end_effector).ValueAt(EndEffectorSplines::Force, coord, time);
        }

        return force;
    }

    Eigen::Vector3d Trajectory::GetEndEffectorLocation(int end_effector, double time) const {
        Eigen::Vector3d ee_location = Eigen::Vector3d::Zero();
        for (int coord = 0; coord < POS_VARS; coord++) {
            ee_location(coord) = ee_splines_.at(end_effector).ValueAt(EndEffectorSplines::Position, coord, time);
        }

        return ee_location;
    }

    double Trajectory::GetTime(int node) const {
        return init_time_ + node_dt_*node;
    }

    int Trajectory::GetTotalForceSplineVars() const {
        return force_spline_vars_;
    }

    bool Trajectory::IsForceMutable(int ee, double time) const {
        return ee_splines_.at(ee).IsForceMutable(time);
    }

    int Trajectory::GetTotalVariables() const {
        return force_spline_vars_ + pos_spline_vars_ + states_.size()*(states_.at(0).size()-1);
    }

    vector_t Trajectory::SplinesAsVec() const {
        // TODO: DMA
        vector_t qp_vec = vector_t::Zero(GetTotalForceSplineVars() + GetTotalPosSplineVars());

        int force_idx = 0;
        int pos_idx = GetTotalForceSplineVars();
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                qp_vec.segment(force_idx, ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Force, coord)) =
                        ee_splines_.at(ee).GetSplineAsQPVec(EndEffectorSplines::Force, coord);
                force_idx += ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Force, coord);

                if (coord < 2) {
                    qp_vec.segment(pos_idx, ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Position, coord)) =
                            ee_splines_.at(ee).GetSplineAsQPVec(EndEffectorSplines::Position, coord);
                    pos_idx += ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Position, coord);
                }
            }
        }

        assert(force_idx == GetTotalForceSplineVars());
        assert(pos_idx == GetTotalForceSplineVars() + GetTotalPosSplineVars());
        return qp_vec;
    }

    int Trajectory::GetNode(double time) const {
        return std::ceil((time - init_time_)/node_dt_);
    }

    controller::Contact Trajectory::GetDesiredContacts(double time) const {
        controller::Contact contact(ee_splines_.size());
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            if (GetForce(ee, time)(2) <= 0.0) {
                contact.in_contact_.at(ee) = false;
            } else {
                contact.in_contact_.at(ee) = true;
            }
        }
        // TODO: Get the contact frames

        return contact;
    }


    vector_3t Trajectory::GetForcePartialWrtContactTime(int end_effector, double time, int contact_idx) const {
        vector_3t force_partials = vector_3t::Zero();

        for (int coord = 0; coord < POS_VARS; coord++) {
            force_partials(coord) = ee_splines_.at(end_effector).
                    ComputePartialWrtTime(EndEffectorSplines::Force, coord, time, contact_idx);
        }

        return force_partials;
    }

    vector_3t Trajectory::GetPositionPartialWrtContactTime(int end_effector, double time, int contact_idx) const {
        vector_3t position_partials = vector_3t::Zero();

        for (int coord = 0; coord < 2; coord++) {
            position_partials(coord) = ee_splines_.at(end_effector).
                    ComputePartialWrtTime(EndEffectorSplines::Position, coord, time, contact_idx);
        }

        return position_partials;
    }

    vector_t Trajectory::GetForceCoefPartialsWrtContactTime(int end_effector, int coord, double time,
                                                            int contact_idx) const {
        return ee_splines_.at(end_effector).ComputeCoefPartialWrtTime(EndEffectorSplines::Force,
                                                                      coord, time, contact_idx);
    }

    vector_t Trajectory::GetPositionCoefPartialsWrtContactTime(int end_effector, int coord, double time,
                                                               int contact_idx) const {
        return ee_splines_.at(end_effector).ComputeCoefPartialWrtTime(EndEffectorSplines::Position,
                                                                      coord, time, contact_idx);
    }

    void Trajectory::UpdateContactTimes(const std::vector<time_v>& contact_times) {
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            ee_splines_.at(ee).SetContactTimes(contact_times.at(ee));
        }
    }

} // mpc