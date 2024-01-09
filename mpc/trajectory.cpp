//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <fstream>

#include "trajectory.h"
#include "centroidal_model.h"

namespace mpc {

    Trajectory::Trajectory(int len, int state_size, int num_joints,
                           const std::vector<std::vector<double>>& switching_times, double node_dt,
                           double swing_height, double foot_offset) :
            inputs_(switching_times, num_joints, len, node_dt),
            swing_height_(swing_height),
            foot_offset_(foot_offset),
            fk_traj_(5)  { // TODO: Make not hard coded
        for (int i = 0; i < len; i++) {
            states_.push_back(vector_t::Zero(state_size));
            full_velocities_.push_back(vector_t::Zero(state_size));
        }

        for (const auto& switching_time : switching_times) {
            std::array<Spline, 3> end_effector_pos =
                    {Spline(2, switching_time, false, Spline::Normal),
                     Spline(2, switching_time, false, Spline::Normal),
                     Spline(2, switching_time, false, Spline::Normal)};
            end_effector_pos_.emplace_back(end_effector_pos);
            const std::array<bool, 3> mut_arr = {true, true, false}; // z pos is always a constant trajectory
            mut_flags_.emplace_back(mut_arr);
        }

        UpdatePosSplineVarsCount();
        SetSwingPosZ();

        for (auto& ee_traj : fk_traj_) {
            ee_traj.resize(len);
        }
    }

    // TODO: Implement an equality operator
    Trajectory& Trajectory::operator=(const Trajectory& traj) {
        if (this == &traj) {
            return *this;
        }

        this->pos_spline_vars_ = traj.pos_spline_vars_;
        this->swing_height_ = traj.swing_height_;
        this->foot_offset_ = traj.foot_offset_;
        this->states_ = traj.states_;
        this->inputs_ = traj.inputs_;
        this->end_effector_pos_ = traj.end_effector_pos_;
        this->mut_flags_ = traj.mut_flags_;
        this->full_velocities_ = traj.full_velocities_;
        this->fk_traj_ = traj.fk_traj_;

        return *this;
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
            if (inputs_.IsForceMutable(end_effector, coord, i)) {
                inputs_.UpdateForcePoly(end_effector, coord, i, vars.segment(idx, 2));
                idx += 2;
            }
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
        if (!mut_flags_.at(end_effector).at(coord)) {
            throw std::runtime_error("The chosen spline is not mutable and thus does not provide a index.");
        }

        int num_spline_vars_before = 0;
        for (int ee = 0; ee < end_effector; ee++) {
            for (int j = 0; j < POS_VARS; j++) {
                if (mut_flags_.at(ee).at(j)) {
                    num_spline_vars_before += end_effector_pos_.at(ee).at(j).GetTotalPolyVars();
                }
            }
        }

        int idx_into_ee_coord_spline_vars = 0;
        for (int j = 0; j < coord; j++) {
            if (mut_flags_.at(end_effector).at(j)) {
                idx_into_ee_coord_spline_vars += end_effector_pos_.at(end_effector).at(j).GetTotalPolyVars();
            }
        }

        int vars_idx, vars_affecting;
        std::tie(vars_idx, vars_affecting) = end_effector_pos_.at(end_effector).at(coord).GetVarsIndexEnd(time);

        return std::make_pair(num_spline_vars_before + idx_into_ee_coord_spline_vars + vars_idx, vars_affecting);
    }

    void Trajectory::SetPositionsForAllTime(int ee, const std::array<double, POS_VARS>& ee_pos) {
        for (int coord = 0; coord < POS_VARS; coord++) {
            if (mut_flags_.at(ee).at(coord)) {
                end_effector_pos_.at(ee).at(coord).SetAllPositions(ee_pos.at(coord));
            }
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
            const auto& pos_times = end_effector_pos_.at(ee).at(0).GetPolyTimes();
            const auto& force_times = inputs_.GetForces().at(ee).at(0).GetPolyTimes();
            file << "force: ";
            for (double force_time : force_times) {
                file << force_time << ", ";
            }
            file << std::endl;

            file << "positions: ";
            for (double pos_time : pos_times) {
                file << pos_time << ", ";
            }
            file << std::endl;
        }

        file << "input joint vels: " << std::endl;
        file << joint_vels << std::endl;

        file.close();
    }

    int Trajectory::GetTotalPosConstantsZ() const {
        int num_constants = 0;
        for (const auto& pos: end_effector_pos_) {
            num_constants += pos.at(2).GetNumConstant();
        }

        return num_constants;
    }

    // TODO: think about this more
    double Trajectory::GetTotalTime() const {
//        assert(end_effector_pos_.at(0).at(0).GetEndTime() == inputs_.GetForces().at(0).at(0).GetEndTime());
        return std::min(end_effector_pos_.at(0).at(0).GetEndTime(), inputs_.GetForces().at(0).at(0).GetEndTime());
    }

    void Trajectory::AddPolys(double final_time) {
        while (GetTotalTime() < final_time) {
            for (auto &end_effector_pos: end_effector_pos_) {
                for (int coord = 0; coord < POS_VARS; coord++) {
                    end_effector_pos.at(coord).AddPoly(0.2);
                }
            }
            inputs_.AddPolys(0.2);
        }
        SetSwingPosZ();
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
        for (int ee = 0; ee < end_effector_pos_.size(); ee++) {
            for (int coord = 0; coord < 3; coord++) {
                if (mut_flags_.at(ee).at(coord)) {
                    pos_spline_vars_ += end_effector_pos_.at(ee).at(coord).GetTotalPolyVars();
                }
            }
        }
    }

    void Trajectory::SetEndEffectorSplines(int ee, const Spline& force_spline, const Spline& pos_spline) {
        for (int coord = 0; coord < POS_VARS; coord++) {
            end_effector_pos_.at(ee).at(coord) = pos_spline;
            inputs_.SetForceSpline(ee, coord, force_spline);
        }
        // TODO: Consider doing this differently
        SetSwingPosZ();
    }

    vector_t Trajectory::ConvertToQPVector() const {
        vector_t qp_vec = vector_t::Zero(pos_spline_vars_ + inputs_.GetTotalForceSplineVars()
                + inputs_.GetAllVels().at(0).size()*(inputs_.GetAllVels().size()-1) // TODO: Why?
                + (states_.at(0).size()-1)*states_.size());

        for (int i = 0; i < states_.size(); i++) {
            qp_vec.segment(i*(states_.at(0).size()-1), states_.at(i).size()-1) =
                    CentroidalModel::ConvertManifoldStateToAlgebraState(states_.at(i), states_.at(0));
        }

        qp_vec.segment( (states_.at(0).size()-1)*states_.size(), inputs_.GetTotalForceSplineVars()
                            + inputs_.GetAllVels().at(0).size()*inputs_.GetAllVels().size()) = inputs_.AsQPVector();
        qp_vec.tail(pos_spline_vars_) = PositionAsQPVector();

        return qp_vec;
    }

    vector_t Trajectory::PositionAsQPVector() const {
        vector_t pos_vec = vector_t::Zero(pos_spline_vars_);
        int idx = 0;
        for (int ee = 0; ee < end_effector_pos_.size(); ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                if (mut_flags_.at(ee).at(coord)) {
                    pos_vec.segment(idx, end_effector_pos_.at(ee).at(coord).GetTotalPolyVars()) =
                            end_effector_pos_.at(ee).at(coord).GetAllPolyVars();
                    idx += end_effector_pos_.at(ee).at(coord).GetTotalPolyVars();
                }
            }
        }

        return pos_vec;
    }

    vector_t Trajectory::GetState(int node) const {
        return states_.at(node);
    }

    Eigen::Vector3d Trajectory::GetPosition(int ee, double time) const {
        Eigen::Vector3d position;
        for (int coord = 0; coord < POS_VARS; coord++) {
            position(coord) = end_effector_pos_.at(ee).at(coord).ValueAt(time);
        }

        return position;
    }

    std::vector<bool> Trajectory::GetContacts(double time) const {
        std::vector<bool> in_contact(end_effector_pos_.size());
        for (int ee = 0; ee < end_effector_pos_.size(); ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = end_effector_pos_.at(ee).at(coord).GetVarsIndexEnd(time);
                if (vars_affecting == 1) {
                    in_contact.at(ee) = true;
                } else {
                    in_contact.at(ee) = false;
                }
            }
        }

        return in_contact;
    }

    bool Trajectory::PosIsConstant(int ee, int coord, double time) const {
        return end_effector_pos_.at(ee).at(coord).IsConstant(time);
    }

    int Trajectory::GetTotalPosNonConstantZ() const {
        int non_const = 0;
        const int coord = 2;
        for (const auto& end_effector_pos : end_effector_pos_) {
            const auto& poly_vars = end_effector_pos.at(coord).GetPolyVars();
            for (int i = 0; i < poly_vars.size(); i++) {
                if (poly_vars.size() != 1) {
                    non_const++;
                }
            }
        }

        return non_const;
    }

    bool Trajectory::IsSplineMutable(int ee, int coord) const {
        return mut_flags_.at(ee).at(coord);
    }

    void Trajectory::SetSwingPosZ() {
        const int coord = 2;
        for (auto& end_effector_pos : end_effector_pos_) {
            const auto& poly_vars = end_effector_pos.at(coord).GetPolyVars();
            for (int i = 0; i < poly_vars.size(); i++) {
                if (poly_vars.at(i).size() == 2) {
                    end_effector_pos.at(coord).SetPolyVars(i, {swing_height_, 0});
                } else {
                    end_effector_pos.at(coord).SetPolyVars(i, {foot_offset_});
                }
            }
        }
    }

    void Trajectory::UpdateFullVelocity(int node, const vector_t& vel) {
        full_velocities_.at(node) = vel;
    }

    vector_t Trajectory::GetFullVelocity(int node) {
        return full_velocities_.at(node);
    }

    vector_t Trajectory::GetAcc(int node, double dt) {
        return -0.5*(full_velocities_.at(node+1) - full_velocities_.at(node))/dt;
    }

    std::vector<std::vector<Eigen::Vector3d>> Trajectory::CreateVizData(const CentroidalModel& model) {
        for (int ee = 0; ee < 5; ee++) {
            for (int node = 0; node < states_.size(); node++) {
                if (ee == 4) {
                    fk_traj_.at(ee).at(node) = model.GetCOMPosition(GetState(node));
                } else {
                    fk_traj_.at(ee).at(node) =
                            model.GetEndEffectorLocationCOMFrame(GetState(node),
                                                                 model.GetEndEffectorFrame(ee))
                            + model.GetCOMPosition(GetState(node));
                }
            }
        }


        return fk_traj_;
    }

} // mpc