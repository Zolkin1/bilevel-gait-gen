//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <fstream>

#include "trajectory.h"
#include "model.h"

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

//    const std::vector<std::array<Spline, 3>>& Trajectory::GetPositions() const {
//        return end_effector_pos_;
//    }

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

//    void Trajectory::SetInputVels(int idx, const vector_t& joint_vels) {
//        if (!using_joints_) {
//            throw std::runtime_error("Cannot set input velocities");
//        }
//    }

    int Trajectory::GetTotalPosSplineVars() const {
        return pos_spline_vars_;
    }

    void Trajectory::UpdateForceSpline(int end_effector, int coord, const vector_t& vars) {
        int idx = 0;

        const std::vector<int> mut_nodes = ee_splines_.at(end_effector).GetMutableNodes(EndEffectorSplines::Force);
        for(auto& node : mut_nodes) {
            vector_2t temp;
            for (int j = 0; j < 2; j++) {
                temp(j) = vars(idx + j);
            }

            ee_splines_.at(end_effector).SetVars(EndEffectorSplines::Force, coord, node, temp);
            idx += 2;
        }

//        for (int i = 0; i < forces_.at(end_effector).at(coord).GetNumPolyTimes(); i++) {
//            if (IsForceMutable(end_effector, coord, i)) {
//
//                // TODO: DMA
//                std::vector<double> temp(2);
//                for (int j = 0; j < 2; j++) {
//                    temp.at(j) = vars(idx + j);
//                }
//
//                forces_.at(end_effector).at(coord).UpdatePolyVar(i, temp);
//                idx += 2;
//            }
//        }
    }

    void Trajectory::UpdatePositionSpline(int end_effector, int coord, const vector_t& vars) {
        int idx = 0;

        const std::vector<int> mut_nodes = ee_splines_.at(end_effector).GetMutableNodes(EndEffectorSplines::Position);
        for(auto& node : mut_nodes) {
            vector_2t temp;
            temp << vars(idx), 0;
//            for (int j = 0; j < 2; j++) {
//                temp(j) = vars(idx + j);
//            }

            ee_splines_.at(end_effector).SetVars(EndEffectorSplines::Position, coord, node, temp);
            idx++;
        }
//
//        for (int i = 0; i < end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes(); i++) {
//            // TODO: The below only works with constant splines (see commented code for otherwise)
//            std::vector<double> temp(1);
//            temp.at(0) = vars.segment(idx, 1)(0);
//
//            if (i == 0) {
//                end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i, temp);
//                if (end_effector_pos_.at(end_effector).at(coord).IsStartPairConstant()) {
//                    end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i + 1, temp);
//                    i++;
//                }
//            } else if (i == end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes()-2) {
//                end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i, temp);
//                if (end_effector_pos_.at(end_effector).at(coord).IsEndPairConstant()) {
//                    end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i + 1, temp);
//                    i++;
//                }
//            } else {
//                end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i, temp);
//                if (i != end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes()-1) {
//                    end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i + 1, temp);
//                    i++;
//                }
//            }
//
//            idx++;


//            int num_vars = end_effector_pos_.at(end_effector).at(coord).GetNumPolyVars(i);
//
//            std::vector<double> temp(num_vars);
//            for (int j = 0; j < num_vars; j++) {
//                temp.at(j) = vars.segment(idx, num_vars)(j);
//            }
//            end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i, temp);

//            if (num_vars == 1 && i+1 < end_effector_pos_.at(end_effector).at(coord).GetNumPolyTimes() &&
//            end_effector_pos_.at(end_effector).at(coord).GetNumPolyVars(i+1) == 1) {
//                // Set then skip the additional constant term
//                end_effector_pos_.at(end_effector).at(coord).UpdatePolyVar(i+1, temp);
//                i++;
//            }
//            idx += num_vars;
//
//        }
    }

    std::pair<int, int> Trajectory::GetPositionSplineIndex(int end_effector, double time, int coord) const {
        if (coord == 2) {
            throw std::runtime_error("The chosen spline is not mutable and thus does not provide a index.");
        }

        int num_spline_vars_before = 0;
        for (int ee = 0; ee < end_effector; ee++) {
            num_spline_vars_before += 2*ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Position);
        }

        int idx_into_ee_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            idx_into_ee_coord_spline_vars += ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Position);
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = ee_splines_.at(end_effector).GetVarsIdx(
                EndEffectorSplines::Position, coord, time);

        return std::make_pair(num_spline_vars_before + idx_into_ee_coord_spline_vars + vars_idx, vars_affecting);
    }

    void Trajectory::SetPositionsForAllTime(int ee, const vector_3t& ee_pos) {
        for (int coord = 0; coord < 2; coord++) {
            const std::vector<int> mut_nodes = ee_splines_.at(ee).GetMutableNodes(EndEffectorSplines::Position);
            for (auto& node : mut_nodes) {
                vector_2t vars = {ee_pos(coord), 0};
                ee_splines_.at(ee).SetVars(EndEffectorSplines::Position, coord, node, vars);
            }
        }
    }

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

//    int Trajectory::GetTotalPosConstantsZ() const {
//        int num_constants = 0;
//        for (const auto& pos: end_effector_pos_) {
//            num_constants += pos.at(2).GetNumConstant();
//        }
//
//        return num_constants;
//    }

//    double Trajectory::GetTotalTime() const {
//        return std::min(end_effector_pos_.at(0).at(0).GetEndTime(), forces_.at(0).at(0).GetEndTime());  // In theory this is redundant
//    }

    void Trajectory::AddPolys(double final_time) {
        for (auto & ee_spline : ee_splines_) {
            while (ee_spline.GetEndTime() < final_time) {
                ee_spline.AddPoly(0.2);    // TODO: Make not hard coded
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
            pos_spline_vars_ += 2*ee_spline.GetTotalPolyVars(EndEffectorSplines::Position);
            force_spline_vars_ += 3*ee_spline.GetTotalPolyVars(EndEffectorSplines::Force);
        }
    }

//    void Trajectory::SetEndEffectorSplines(int ee, const Spline& force_spline, const Spline& pos_spline) {
//        for (int coord = 0; coord < POS_VARS; coord++) {
//            if (coord == 2) {
//                std::vector<double> times(pos_spline.GetNumPolyTimes()-1);
//                for (int i = 0; i < times.size(); i++) {
//                    times.at(i) = pos_spline.GetPolyTimes().at(i+1);
//                }
//                Spline position1(2, times, !pos_spline.IsStartPairConstant(), Spline::Normal);
//                end_effector_pos_.at(ee).at(coord) = position1;
//            } else {
//                end_effector_pos_.at(ee).at(coord) = pos_spline;
//            }
//            forces_.at(ee).at(coord) = force_spline;
//        }
//
//        SetSwingPosZ();
//        UpdateSplineVarsCount();
//    }

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

//    Eigen::Vector3d Trajectory::GetPosition(int ee, double time) const {
//        Eigen::Vector3d position;
//        for (int coord = 0; coord < POS_VARS; coord++) {
//            position(coord) = end_effector_pos_.at(ee).at(coord).ValueAt(time);
//        }
//
//        return position;
//    }

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

//    bool Trajectory::PosIsConstant(int ee, int coord, double time) const {
//        return end_effector_pos_.at(ee).at(coord).IsConstant(time);
//    }

//    int Trajectory::GetTotalPosNonConstantZ() const {
//        int non_const = 0;
//        const int coord = 2;
//        for (const auto& end_effector_pos : end_effector_pos_) {
//            const auto& poly_vars = end_effector_pos.at(coord).GetPolyVars();
//            for (int i = 0; i < poly_vars.size(); i++) {
//                if (poly_vars.size() != 1) {
//                    non_const++;
//                }
//            }
//        }
//
//        return non_const;
//    }

//    bool Trajectory::IsSplineMutable(int ee, int coord) const {
//        return mut_flags_.at(ee).at(coord);
//    }

// TODO: Fix this
    void Trajectory::SetSwingPosZ() {
        const int coord = 2;
        for (auto& ee_spline : ee_splines_) {
            const std::vector<int> mut_nodes = ee_spline.GetMutableNodes(EndEffectorSplines::Position);
            for (auto& node: mut_nodes) {
                if (ee_spline.GetNodeType(EndEffectorSplines::Force, node) == NodeType::FullDeriv) {
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

//    vector_t Trajectory::GetAcc(int node, double dt) {
//        return -(full_velocities_.at(node+1) - full_velocities_.at(node))/dt;
//    }

//    std::vector<std::vector<Eigen::Vector3d>> Trajectory::CreateVizData(const Model* model) {
//        for (int ee = 0; ee < 5; ee++) {
//            for (int node = 0; node < states_.size(); node++) {
//                if (ee == 4) {
//                    fk_traj_.at(ee).at(node) = model->GetCOMPosition(GetState(node));
//                } else {
//                    // TODO: How to deal with this
//                    fk_traj_.at(ee).at(node) = GetEndEffectorLocation(ee, GetTime(node));
////                            model.GetEndEffectorLocationCOMFrame(GetState(node),
////                                                                 model.GetEndEffectorFrame(ee))
////                            + model->GetCOMPosition(GetState(node));
//                }
//            }
//        }
//
//        return fk_traj_;
//    }

    int Trajectory::GetNumContactNodes(int ee) const {
        return ee_splines_.at(ee).GetMutableNodes(EndEffectorSplines::Position).size();
    }

    std::vector<std::vector<double>> Trajectory::GetContactTimes() const {
        std::vector<std::vector<double>> contact_times(ee_splines_.size());
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            std::vector<int> mut_nodes = ee_splines_.at(ee).GetMutableNodes(EndEffectorSplines::Position);
            for (int mut_node : mut_nodes) {
                contact_times.at(ee).push_back(ee_splines_.at(ee).GetTimes().at(mut_node));
            }
        }

        return contact_times;
    }

    vector_t Trajectory::GetSplineLin(const mpc::Trajectory::SplineTypes& spline_type, int ee, int coord, double time) const {
        // TODO: Change this function to be two or change the spline type function
        switch (spline_type) {
            case SplineTypes::Force:
                return ee_splines_.at(ee).GetPolyVarsLin(EndEffectorSplines::Force, coord, time);
            case SplineTypes::Position:
                return ee_splines_.at(ee).GetPolyVarsLin(EndEffectorSplines::Position, coord, time);
            default:
                throw std::runtime_error("Spline type not implemented.");
        }
    }

    std::pair<int, int> Trajectory::GetForceSplineIndex(int end_effector, double time, int coord) const {
        int num_spline_vars_before = 0;
        for (int ee = 0; ee < end_effector; ee++) {
            num_spline_vars_before += 3*ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Force);
        }

        int idx_into_force_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            idx_into_force_coord_spline_vars += ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Force);
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = ee_splines_.at(end_effector).GetVarsIdx(EndEffectorSplines::Force, coord, time);

        return std::make_pair(num_spline_vars_before + idx_into_force_coord_spline_vars + vars_idx, vars_affecting);
    }

    int Trajectory::GetTotalPolyVars(const mpc::Trajectory::SplineTypes& spline_type, int end_effector, int coord) {
        // TODO: Change this switch statement/function to be better
        switch (spline_type) {
            case SplineTypes::Force:
                return ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Force);
            case SplineTypes::Position:
                return ee_splines_.at(end_effector).GetTotalPolyVars(EndEffectorSplines::Position);
            default:
                throw std::runtime_error("Spline type not impelemented.");
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

//    bool Trajectory::IsForceMutable(int ee, int coord, int idx) const {
//        return forces_.at(ee).at(coord).IsMutable(idx);
//    }

    int Trajectory::GetTotalVariables() const {
        // TODO: Account for joint potentially
        // TODO: Do I want states in here?
        return force_spline_vars_ + pos_spline_vars_ + states_.size()*(states_.at(0).size()-1);
    }

    vector_t Trajectory::SplinesAsVec() const {
        // TODO: DMA
        vector_t qp_vec = vector_t::Zero(GetTotalForceSplineVars() + GetTotalPosSplineVars());

        // Force spline
        int force_idx = 0;
        int pos_idx = GetTotalForceSplineVars();
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                qp_vec.segment(force_idx, ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Force)) =
                        ee_splines_.at(ee).GetSplineAsQPVec(EndEffectorSplines::Force, coord);
                force_idx += ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Force);

                if (coord < 2) {
                    qp_vec.segment(pos_idx, ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Position)) =
                            ee_splines_.at(ee).GetSplineAsQPVec(EndEffectorSplines::Position, coord);
                    pos_idx += ee_splines_.at(ee).GetTotalPolyVars(EndEffectorSplines::Position);
                }
            }
        }

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

    void Trajectory::UpdateContactTimes(const std::vector<std::vector<double>>& contact_times) {
        for (int ee = 0; ee < ee_splines_.size(); ee++) {
            ee_splines_.at(ee).SetContactTimes(contact_times.at(ee));
        }
    }

} // mpc