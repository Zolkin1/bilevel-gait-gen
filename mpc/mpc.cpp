//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <array>
#include <Eigen/SparseCore>

#include "mpc.h"
#include "qp/osqp_interface.h"

namespace mpc {

    using triplet_t = Eigen::Triplet<double>;

    MPCInfo::MPCInfo() {}
    MPCInfo::MPCInfo(const MPCInfo& info) {
        num_nodes = info.num_nodes;
        num_qp_iterations = info.num_qp_iterations;
        num_contacts = info.num_contacts;
        friction_coef = info.friction_coef;
        vel_bounds = info.vel_bounds;
        joint_bounds_lb = info.joint_bounds_lb;
        joint_bounds_ub = info.joint_bounds_ub;
        ee_frames = info.ee_frames;
        discretization_steps = info.discretization_steps;
        num_switches = info.num_switches;
        integrator_dt = info.integrator_dt;
        force_bound = info.force_bound;
        swing_height = info.swing_height;
        foot_offset = info.foot_offset;
        nom_state = info.nom_state;
        ee_box_size = info.ee_box_size;
        real_time_iters = info.real_time_iters;
        verbose = info.verbose;
    }

    MPC::MPC(const MPCInfo& info, const std::string& robot_urdf) :
        info_(info),
        model_(robot_urdf, info.ee_frames, info.discretization_steps, info.integrator_dt, info.nom_state),
        num_states_(model_.GetNumTangentStates()),
        num_ee_(model_.GetNumEndEffectors()),
        prev_traj_(info.num_nodes+1, model_.GetNumManifoldStates(), model_.UsesJoints(),
               CreateDefaultSwitchingTimes(info.num_switches, num_ee_,
                    info.integrator_dt*(info.num_nodes)),
                    info.integrator_dt, info.swing_height, info.foot_offset),
        data_(false, 25000, 2000, model_.GetApplicableConstraints()),
        integrator_(info.integrator_dt),
        using_clarabel_(true) {

        assert(info_.ee_frames.size() == num_ee_);

        init_time_ = 0;

        UpdateNumInputs();

        SetFrictionPyramid();

        Phi_ = matrix_t::Zero(num_states_, num_states_);
        Phi_w_ = vector_t::Zero(num_states_);

//        prev_qp_sol = vector_t::Zero(data_.num_decision_vars);

//        line_search_res_ = vector_t::Zero(data_.num_decision_vars);
        mu_ = 5000;
        run_num_ = 0;
        constraint_idx_ = 0;

        in_real_time_ = false;

        data_.using_clarabel_ = using_clarabel_;
    }

    Trajectory MPC::CreateInitialRun(const mpc::vector_t& state, const std::vector<vector_3t>& ee_start_locations) {
        bool converged = false;
        int num_iter = 0;
        in_real_time_ = false;

        qp_solver->ConfigureForInitialRun();

        while (!converged && num_iter < 10) {
            Solve(state, 0, ee_start_locations);
            num_iter++;
        }
        return prev_traj_;
    }

    Trajectory MPC::GetRealTimeUpdate(const vector_t& state, double init_time,
                                      const std::vector<vector_3t>& ee_start_locations,
                                      bool high_quality) {
//        if (high_quality) {
//            qp_solver->SetSolveTolerances(1e-4, 1e-4);
//        } else {
//            qp_solver->SetSolveTolerances(1e-4, 4e-4);
//        }

        if (in_real_time_) {
            return Solve(state, init_time, ee_start_locations);
        } else {
            qp_solver->ConfigureForRealTime(info_.real_time_iters);
            in_real_time_ = true;
            return Solve(state, init_time, ee_start_locations);
        }
    }

    void MPC::SetWarmStartTrajectory(const mpc::Trajectory &trajectory) {
        prev_traj_ = trajectory;

        UpdateNumInputs();

        // TODO: Might want to put this back
//        prev_qp_sol = ConvertTrajToQPVec(prev_traj_);

        init_time_ = trajectory.GetTime(0);
    }

    void MPC::SetQuadraticCostTerm(const matrix_t& Q) {
        if (Q.rows() != num_states_ || Q.cols() != num_states_) {
            throw std::runtime_error("Supplied quadratic cost term is the wrong size.");
        }

        Q_ = Q;
    }

    void MPC::SetLinearCostTerm(const vector_t &w) {
        if (w.rows() != num_states_) {
            throw std::runtime_error("Supplied linear cost term is the wrong size.");
        }

        w_ = w;
    }

    void MPC::SetQuadraticFinalCost(const mpc::matrix_t& Phi) {
        if (Phi.rows() != num_states_ || Phi.cols() != num_states_) {
            throw std::runtime_error("Wrong sized quadratic final cost.");
        }

        Phi_ = Phi;
    }

    void MPC::SetLinearFinalCost(const vector_t& w) {
        if (w.size() != num_states_) {
            throw std::runtime_error("Wrong sized linear final cost term.");
        }

        Phi_w_ = w;
    }

    void MPC::SetFrictionPyramid() {
        // Assumes flat ground
        Eigen::Vector3d h = {1, 0, 0};
        Eigen::Vector3d l = {0, 1, 0};
        Eigen::Vector3d n = {0, 0, 1};

        friction_pyramid_ << (h - n*info_.friction_coef).transpose(),
                            -(h + n*info_.friction_coef).transpose(),
                            (l - n*info_.friction_coef).transpose(),
                            -(l + n*info_.friction_coef).transpose();
    }

    // TODO: Examine what this looks like on areas when the force is NOT mutable
    void MPC::AddFrictionConeConstraints() {
        int force_idx = GetForceSplineStartIdx();
        int row_idx = 0;

        const std::vector<time_v> contact_times = prev_traj_.GetContactTimes();
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int time_idx = 0; time_idx < contact_times.at(ee).size() - 1; time_idx++) {
                if (contact_times.at(ee).at(time_idx).GetType() == TouchDown) {
                    for (int i = 0; i < FB_PER_FORCE; i++) { // for each cone constraint
                        for (int coord = 0; coord < POS_VARS; coord++) {
                            double lower_time = contact_times.at(ee).at(time_idx).GetTime();
                            double upper_time = contact_times.at(ee).at(time_idx + 1).GetTime();
                            double time =
                                    (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) *
                                    (upper_time - lower_time) +
                                    lower_time;
                            assert(prev_traj_.IsForceMutable(ee, time));


                            int vars_index, vars_affecting;
                            std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee, time, coord);
                            const vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Force, ee,
                                                                              coord,
                                                                              time);

                            for (int fric_con = 0; fric_con < 4; fric_con++) {
                                // All the friction constraints are effected by all 3 coordinates of the end effectors
                                data_.constraint_mat_.SetMatrix(friction_pyramid_(fric_con, coord) * vars_lin.transpose(),
                                                                constraint_idx_ + row_idx + fric_con,
                                                                force_idx + vars_index);
                                data_.friction_cone_ub_(row_idx + fric_con) = 0;
                                data_.friction_cone_lb_(row_idx + fric_con) = -1e30;//qp_solver->GetInfinity(1)(0);
                            }
                        }
                        row_idx += 4;
                    }
                }
            }
        }


        assert(row_idx == data_.num_cone_constraints_);
        constraint_idx_ += row_idx;

        // -------------------- Friction pyramid -------------------- //
//        int force_offset = GetForceSplineStartIdx();
//        int idx = 0;
//        for (int node = 0; node < info_.num_nodes+1; node++) {
//            double time = GetTime(node);
//            for (int ee = 0; ee < num_ee_; ee++) {
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    if (prev_traj_.IsForceMutable(ee, time)) {
//                        vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Force, ee, coord, time);
//
//                        int vars_index, vars_affecting;
//                        std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee,time,coord);
//
//                        for (int fric_con = 0; fric_con < 4; fric_con++) {
//                            // All the friction constraints are effected by all 3 coordinates of the end effectors
//                            data_.constraint_mat_.SetMatrix(friction_pyramid_(fric_con, coord) * vars_lin.transpose(),
//                                                            constraint_idx_ + idx + fric_con,
//                                                            force_offset + vars_index);
//                            data_.friction_cone_ub_(idx + fric_con) = 0;
//                            data_.friction_cone_lb_(idx + fric_con) = -1e30;//qp_solver->GetInfinity(1)(0);
//                        }
//                    }
//                }
//                idx += 4;
//            }
//        }
//
//        constraint_idx_ += idx;
    }

    void MPC::AddFrictionConeConstraintPartials(utils::SparseMatrixBuilder& builder, int contact_idx,
                                                int start_idx, int ee) {

        int force_idx = GetForceSplineStartIdx();

        const std::vector<time_v> contact_times = prev_traj_.GetContactTimes();
        int row_idx = 0;

        for (int i = 0; i < ee; i++) {
            for (int contact_time = 0; contact_time < contact_times.at(i).size()-1; contact_time++) {
                if (contact_times.at(i).at(contact_time).GetType() == TouchDown) {
                    row_idx += 4*FB_PER_FORCE;
                }
            }
        }

        for (int contact_time = 0; contact_time < contact_idx; contact_time++) {
            if (contact_times.at(ee).at(contact_time).GetType() == TouchDown) {
                row_idx += 4*FB_PER_FORCE;
            }
        }

        if (contact_times.at(ee).at(contact_idx).GetType() == LiftOff && contact_idx > 0) {
            row_idx -= 4*FB_PER_FORCE;
        }


        if (contact_times.at(ee).at(contact_idx).GetType() == TouchDown
                && contact_idx < contact_times.at(ee).size()-1) {
            for (int i = 0; i < FB_PER_FORCE; i++) { // for each cone constraint
                for (int coord = 0; coord < POS_VARS; coord++) {
                    double lower_time = contact_times.at(ee).at(contact_idx).GetTime();
                    double upper_time = contact_times.at(ee).at(contact_idx + 1).GetTime();
                    double time =
                            (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) *
                            (upper_time - lower_time) +
                            lower_time;
                    assert(prev_traj_.IsForceMutable(ee, time));

                    double dtimedth = -(static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) + 1.0;

                    int vars_index, vars_affecting;
                    std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee, time, coord);
                    const vector_t vars_partial = prev_traj_.GetForceCoefPartialsWrtContactTime(
                            ee,coord,time, contact_idx, dtimedth);



                    for (int fric_con = 0; fric_con < 4; fric_con++) {
                        // All the friction constraints are effected by all 3 coordinates of the end effectors
                        builder.SetMatrix(friction_pyramid_(fric_con, coord) * vars_partial.transpose(),
                                                        start_idx + row_idx + fric_con,
                                                        force_idx + vars_index);
                    }
                }
                row_idx += 4;
            }
        } else if (contact_idx > 0 && contact_times.at(ee).at(contact_idx).GetType() == LiftOff) {
            for (int i = 0; i < FB_PER_FORCE; i++) { // for each cone constraint
                for (int coord = 0; coord < POS_VARS; coord++) {
                    double lower_time = contact_times.at(ee).at(contact_idx - 1).GetTime();
                    double upper_time = contact_times.at(ee).at(contact_idx).GetTime();
                    double time =
                            (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) *
                            (upper_time - lower_time) +
                            lower_time;
                    assert(prev_traj_.IsForceMutable(ee, time));

                    double dtimedth = (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE));

                    int vars_index, vars_affecting;
                    std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee, time, coord);
                    const vector_t vars_partial = prev_traj_.GetForceCoefPartialsWrtContactTime(
                            ee, coord, time, contact_idx, dtimedth);


                    for (int fric_con = 0; fric_con < 4; fric_con++) {
                        // All the friction constraints are effected by all 3 coordinates of the end effectors
                        builder.SetMatrix(friction_pyramid_(fric_con, coord) * vars_partial.transpose(),
                                                        start_idx + row_idx + fric_con,
                                                        force_idx + vars_index);
                    }
                }
                row_idx += 4;
            }
        }



//        int force_offset = GetForceSplineStartIdx();
//        int idx = 0;
//        for (int node = 0; node < info_.num_nodes+1; node++) {
//            double time = GetTime(node);
//            for (int coord = 0; coord < POS_VARS; coord++) {
//                if (prev_traj_.IsForceMutable(ee, time)) {
//                    vector_t vars_coef = prev_traj_.GetForceCoefPartialsWrtContactTime( ee, coord, time, contact_idx);
//
//                    int vars_index, vars_affecting;
//                    std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee,time,coord);
//
//                    for (int fric_con = 0; fric_con < 4; fric_con++) {
//                        // All the friction constraints are effected by all 3 coordinates of the end effectors
//                        builder.SetMatrix(friction_pyramid_(fric_con, coord) * vars_coef.transpose(),
//                                                        start_idx + idx + fric_con,
//                                                        force_offset + vars_index);
//                    }
//                }
//            }
//            idx += 4;
//        }
    }

    void MPC::AddForceBoxConstraints() {
        int force_idx = GetForceSplineStartIdx();
        int row_idx = 0;

        const int coord = 2;


        int extra_runs = 1;
        if (using_clarabel_) {
            extra_runs = 2;
        }

        constexpr double FB_LB = 0.0;

        const std::vector<time_v> contact_times = prev_traj_.GetContactTimes();
        for (int j = 0; j < extra_runs; j++) {
            for (int ee = 0; ee < num_ee_; ee++) {
                for (int time_idx = 0; time_idx < contact_times.at(ee).size() - 1; time_idx++) {
                    if (contact_times.at(ee).at(time_idx).GetType() == TouchDown) {
                        for (int i = 0; i < FB_PER_FORCE; i++) {
                            double lower_time = contact_times.at(ee).at(time_idx).GetTime();
                            double upper_time = contact_times.at(ee).at(time_idx + 1).GetTime();
                            double time =
                                    (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) * (upper_time - lower_time) +
                                    lower_time;
                            assert(prev_traj_.IsForceMutable(ee, time));


                            int vars_index, vars_affecting;
                            std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee, time, coord);
                            const vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Force, ee, coord,
                                                                              time);

                            if (j == 0) {
                                data_.constraint_mat_.SetMatrix(vars_lin.transpose(),
                                                                constraint_idx_ + row_idx,
                                                                force_idx + vars_index);
                            } else {
                                data_.constraint_mat_.SetMatrix(-vars_lin.transpose(),
                                                                constraint_idx_ + row_idx,
                                                                force_idx + vars_index);
                            }

                            if (j == 0) {
                                data_.force_box_lb_(row_idx) = FB_LB;
                                data_.force_box_ub_(row_idx) = info_.force_bound;
                            }

                            row_idx++;

                        }
                    }
                }
            }
        }

        assert(row_idx == data_.num_force_box_constraints_);
        constraint_idx_ += row_idx;
    }

    void MPC::AddForceBoxConstraintPartials(utils::SparseMatrixBuilder& builder, int contact_idx,
                                            int start_idx, int ee) {
        int force_idx = GetForceSplineStartIdx();

        const std::vector<time_v> contact_times = prev_traj_.GetContactTimes();
        int row_idx = 0;

        for (int i = 0; i < ee; i++) {
            for (int contact_time = 0; contact_time < contact_times.at(i).size()-1; contact_time++) {
                if (contact_times.at(i).at(contact_time).GetType() == TouchDown) {
                    row_idx += FB_PER_FORCE;
                }
            }
        }

        for (int contact_time = 0; contact_time < contact_idx; contact_time++) {
            if (contact_times.at(ee).at(contact_time).GetType() == TouchDown) {
                row_idx += FB_PER_FORCE;
            }
        }

        if (contact_times.at(ee).at(contact_idx).GetType() == LiftOff && contact_idx > 0) {
            row_idx -= FB_PER_FORCE;
        }

        const int coord = 2;


        int extra_runs = 1;
        if (using_clarabel_) {
            extra_runs = 2;
        }

        for (int j = 0; j < extra_runs; j++) {
            vector_t vars_partials;
            int vars_index, vars_affecting;

            if (contact_times.at(ee).at(contact_idx).GetType() == TouchDown && contact_idx < contact_times.at(ee).size()-1) { // TODO: Account for the change in lift off (look backwards in time)
                double lower_time = contact_times.at(ee).at(contact_idx).GetTime();
                double upper_time = contact_times.at(ee).at(contact_idx + 1).GetTime();

                for (int i = 0; i < FB_PER_FORCE; i++) {
                    double time =
                            (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) * (upper_time - lower_time) +
                            lower_time;

                    assert(prev_traj_.IsForceMutable(ee, time));

                    double dtimedth = -(static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) + 1.0;

                    std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee, time, coord);
                    vars_partials = prev_traj_.GetForceCoefPartialsWrtContactTime(ee, coord,
                                                                      time, contact_idx, dtimedth);

                    if (j == 0) {
                        builder.SetMatrix(vars_partials.transpose(),
                                          start_idx + row_idx,
                                          force_idx + vars_index);
                    } else {
                        builder.SetMatrix(-vars_partials.transpose(),
                                          start_idx + row_idx,
                                          force_idx + vars_index);
                    }

                    row_idx++;
                }
            } else if (contact_times.at(ee).at(contact_idx).GetType() == LiftOff && contact_idx > 0) { // Liftoff - so effect the derivatives behind in time
                double lower_time = contact_times.at(ee).at(contact_idx-1).GetTime();
                double upper_time = contact_times.at(ee).at(contact_idx).GetTime();

                for (int i = 0; i < FB_PER_FORCE; i++) {
                    double time =
                            (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE)) * (upper_time - lower_time) +
                            lower_time;

                    assert(prev_traj_.IsForceMutable(ee, time));

                    double dtimedth = (static_cast<double>(i) / static_cast<double>(FB_PER_FORCE));


                    std::tie(vars_index, vars_affecting) = prev_traj_.GetForceSplineIndex(ee, time, coord);
                    vars_partials = prev_traj_.GetForceCoefPartialsWrtContactTime(ee, coord,
                                                                                                 time, contact_idx, dtimedth);

                    if (j == 0) {
                        builder.SetMatrix(vars_partials.transpose(),
                                          start_idx + row_idx,
                                          force_idx + vars_index);
                    } else {
                        builder.SetMatrix(-vars_partials.transpose(),
                                          start_idx + row_idx,
                                          force_idx + vars_index);
                    }
                    row_idx++;

                }
            }


//            int skip = 0;
//            for (int end_effector = 0; end_effector < num_ee_; end_effector++) {
//                if (end_effector != ee) {
//                    for (int time_idx = 0; time_idx < contact_times.at(end_effector).size()-1; time_idx++) {
//                        if (contact_times.at(end_effector).at(time_idx).GetType() == TouchDown) {
//                            skip += FB_PER_FORCE;
//                        }
//                    }
//                }
//            }

            row_idx += data_.num_force_box_constraints_/2 - FB_PER_FORCE;

        }
//        assert(row_idx == data_.num_force_box_constraints_); // this is tough to do because I need to know the number of force segments for that end effector
//        start_idx += row_idx;
    }

    void MPC::AddQuadraticTrackingCost(const vector_t& state_des, const matrix_t& Q) {
        if (Q.rows() != num_states_ || Q.cols() != num_states_) {
            throw std::runtime_error("Supplied quadratic cost term is the wrong size.");
        }

        Q_ = Q;
        w_ = -1*Q*state_des;
    }

    void MPC::AddHessianApproxCost() {
        for (int node = 0; node < info_.num_nodes; node++) {
            data_.cost_mat_.SetMatrix(Q_, node*num_states_, node*num_states_);
        }

        if (Q_forces_.size() > 0) {
            data_.cost_mat_.SetMatrix(Q_forces_, GetForceSplineStartIdx(), GetForceSplineStartIdx());
        }

//        data_.cost_mat_.SetDiagonalMatrix(10, GetPosSplineStartIdx(), GetPosSplineStartIdx(), prev_traj_.GetTotalPosSplineVars());
    }

    void MPC::AddGradientCost() {
        for (int node = 0; node < info_.num_nodes; node++) {
            data_.cost_linear.segment(node * num_states_, num_states_) = w_;
        }
    }

    void MPC::AddFinalCost() {
        data_.cost_mat_.SetMatrix(Phi_, info_.num_nodes*num_states_, info_.num_nodes*num_states_);

        data_.cost_linear.segment(info_.num_nodes*num_states_, num_states_) = Phi_w_;
    }

    std::vector<std::vector<double>> MPC::CreateDefaultSwitchingTimes(int num_switches, int num_ee, double horizon) {
        std::vector<std::vector<double>> switching_times;
        for (int i = 0; i < num_ee; i++) {
            std::vector<double> times;
//            times.push_back(0);
//            for (int j = 0; j < num_switches; j++) {
//                times.push_back((j+1)*horizon/num_switches);
//            }
//            for (int ee = 0; ee < num_ee_; ee++) {

//            times.push_back(0);
//            times.push_back(0.15);
//            times.push_back(0.3);
//            times.push_back(0.45);
//            times.push_back(0.6);
//            times.push_back(0.75);
//            times.push_back(0.9);
//            times.push_back(1.05);

            times.push_back(0);
            times.push_back(0.3);
            times.push_back(0.6);
            times.push_back(0.9);
            times.push_back(1.2);

//            times.push_back(0);
//            times.push_back(0.2);
//            times.push_back(0.4);
//            times.push_back(0.65);
//            times.push_back(0.8);
//            times.push_back(1);
//            }

            switching_times.push_back(times);
        }

        return switching_times;
    }

    void MPC::UpdateQPSizes() {
        UpdateNumInputs();
        data_.num_decision_vars = (info_.num_nodes+1)*num_states_ + num_inputs_;


        // TODO: May want to put in a check for this, but it should always be active
        if (using_clarabel_) {
            data_.num_force_box_constraints_ = GetNumForceBoxConstraints();
            data_.num_cone_constraints_ = GetNumFricConeConstraints();
            data_.num_td_pos_constraints_ = GetNumTDConstraints();
        } else {
            data_.num_force_box_constraints_ = GetNodeIntersectMutableForces(); // TODO: Change
            // TODO: Add other constraint
        }
    }

    void MPC::SetDefaultGaitTrajectory(Gaits gait, int num_polys, const std::vector<vector_3t>& ee_pos) {
        if (num_ee_ != 4) {
            throw std::runtime_error("Default gaits have only been implemented for 4 end effectors!");
        }

        std::vector<std::vector<double>> switching_times;

        switch (gait) {
            case Trot: {
                std::vector<std::vector<double>> times(num_ee_);

                // TODO: Make this a yaml configuration
                for (int ee = 0; ee < num_ee_; ee++) {
                    times.at(ee).push_back(0);
                    times.at(ee).push_back(0.2);
//                    times.at(ee).push_back(0.4);
//                    times.at(ee).push_back(0.6);
//                    times.at(ee).push_back(0.8);
//                    times.at(ee).push_back(1);
                }

//                prev_traj_.UpdateContactTimes(times);

//                Spline position1(1, times, true, Spline::Constants);
//                Spline position2(1, times, false, Spline::Constants);
//
//                Spline force1(num_polys, times, false, Spline::Force);
//                Spline force2(num_polys, times, true, Spline::Force);
//                for (int j = 0; j < force1.GetNumPolyTimes(); j++) {
//                    if (force1.IsMutable(j)) {
//                        std::vector<double> vars(2);
//                        vars.at(0) = 100;
//                        vars.at(1) = 0;
//                        force1.SetPolyVars(j, vars);
//                    }
//                }
//                for (int j = 0; j < force2.GetNumPolyTimes(); j++) {
//                    if (force2.IsMutable(j)) {
//                        std::vector<double> vars(2);
//                        vars.at(0) = 100;
//                        vars.at(1) = 0;
//                        force2.SetPolyVars(j, vars);
//                    }
//                }

//
//                for (int ee = 0; ee < num_ee_; ee++) {
//                    if (ee == 0 || ee == 3) {
//                        prev_traj_.SetEndEffectorSplines(ee, force1, position1);
//                    } else {
//                        prev_traj_.SetEndEffectorSplines(ee, force2, position2);
//                    }
//                    for (int coord = 0; coord < POS_VARS; coord++) {
//                        prev_traj_.SetPositionsForAllTime(ee, ee_pos.at(ee));
//                    }
//                }

                prev_traj_.PrintTrajectoryToFile("trot_test.txt");
                break;
            }
            case Amble: {
                throw std::runtime_error("Amble not implemented yet!");
                break;
            }
            case Static_Walk: {
                throw std::runtime_error("Static Walk not implemented yet!");
                break;
            }
            default:
                throw std::runtime_error("Unsupported gait.");
        }

    }

    void MPC::SetStateTrajectoryWarmStart(const std::vector<vector_t>& states) {
        assert(states.size() == info_.num_nodes + 1 && states.at(0).size() == prev_traj_.GetStates().at(0).size());

        for (int node = 0; node < info_.num_nodes + 1; node++) {
            prev_traj_.SetState(node, states.at(node));
        }
    }

    vector_t MPC::GetTargetConfig(double time) const {
        int node = floor((time - init_time_)/info_.integrator_dt);
        return prev_traj_.GetStates().at(node).tail(num_states_ - 5);
    }

    vector_t MPC::GetForceTarget(double time) const {
        controller::Contact contacts = GetDesiredContacts(time);
        vector_t forces(contacts.GetNumContacts()*3);
        int idx = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            if (contacts.in_contact_.at(ee)) {
                forces.segment(idx, POS_VARS) = prev_traj_.GetForce(ee, time);
                idx += 3;
            }
        }
        return forces;
    }

    const Model* MPC::GetModel() const {
        return &model_;
    }

    double MPC::LineSearch(const vector_t& direction, const vector_t& init_state) {
        double alpha = 1;

       // Note: for now just using the equality constraints on the merit function
        double merit = GetMeritValue(prev_qp_sol, mu_, init_state);
        double merit_step = GetMeritValue(alpha*direction + prev_qp_sol, mu_, init_state);
        double merit_directional = GetMeritGradient(prev_qp_sol, direction, mu_, init_state);

        int i = 0;
        while ((merit - merit_step) < -0.00001 * alpha * merit_directional && i < 5) {
            alpha *= 0.5;
            vector_t temp = (alpha * direction) + prev_qp_sol;
            merit_step = GetMeritValue(temp, mu_, init_state);
            i++;
        }

        return alpha;
    }

    double MPC::GetMeritValue(const vector_t& x, double mu, const vector_t& init_state) {
        const Trajectory temp_traj = ConvertQPSolToTrajectory(x, init_state);

        return mu*GetEqualityConstraintValues(temp_traj, init_state).lpNorm<1>() + GetCostValue(x);
    }

    double MPC::GetMeritValue(const Trajectory& traj, double mu, const mpc::vector_t& init_state) {
        return mu*GetEqualityConstraintValues(traj, init_state).lpNorm<1>() + GetCostValue(ConvertTrajToQPVec(traj));
    }

    double MPC::GetCostValue(const vector_t& x) const {
        return 0.5*x.dot(data_.sparse_cost_*(x)) + data_.cost_linear.dot(x);
    }

    // TODO: Optimize this and make it so it is called less
    vector_t MPC::GetEqualityConstraintValues(const Trajectory& traj, const vector_t& init_state) {
        // TODO: DMA
        vector_t eq_constraints = vector_t::Zero(data_.num_dynamics_constraints - num_states_); // data_.num_fk_constraints_

        for (int node = 0; node < info_.num_nodes; node++) {
            eq_constraints.segment(node*num_states_, num_states_) =
                    model_.ConvertManifoldStateToTangentState(traj.GetState(node+1), init_state)
                    - integrator_.CalcIntegral(
                            model_.ConvertManifoldStateToTangentState(traj.GetState(node), init_state),
                                               traj, GetTime(node), 1, model_, init_state);
        }
        return eq_constraints;
    }

    // TODO: Consider replacing with the trajectory GetTime function
    double MPC::GetTime(int node) const {
        return node * info_.integrator_dt + init_time_;
    }

    double MPC::GetMeritGradient(const vector_t& x, const vector_t& p, double mu, const vector_t& init_state) {
        const Trajectory temp_traj = ConvertQPSolToTrajectory(x, init_state);

        return (data_.sparse_cost_*x + data_.cost_linear).dot(p)
        - mu * GetEqualityConstraintValues(temp_traj, init_state).lpNorm<1>();
    }


    // TODO: Make this weight only on force values, so we don't hit the derivatives
    void MPC::AddForceCost(double weight) {
        force_weight_ = weight;
        const int num_forces = prev_traj_.GetTotalForceSplineVars();
        Q_forces_.resize(num_forces, num_forces);
        Q_forces_.setZero();        // TODO: Note. Without this my performance was totally shot and occasional errors
        for (int i = 0; i < num_forces; i++) {
            double weight1 = weight;
            if (i < prev_traj_.GetTotalForceSplineVars()/num_ee_) {
                weight1 = weight1 * 1;
            }
            Q_forces_(i,i) = weight1;
        }
    }

    void MPC::RecordStats(double alpha, const vector_t& direction, const SolveQuality& solve_type,
                          const vector_t& ref_state, double solve_time, double cost) {
        equality_constraint_violations_.push_back(GetEqualityConstraintValues(prev_traj_, ref_state).lpNorm<1>());
        step_norm_.push_back(direction.norm());
        alpha_.push_back(alpha);
        cost_result_.push_back(GetCostValue(prev_qp_sol));
        merit_result_.push_back(GetMeritValue(prev_traj_, mu_, ref_state));
        merit_directional_deriv_.push_back(GetMeritGradient(prev_qp_sol - alpha*direction, direction, mu_, ref_state));
        solve_type_.push_back(solve_type);
        ref_state_.push_back(ref_state);
        solve_time_.push_back(solve_time);
        cost_.push_back(cost);
    }

    void MPC::PrintStats() {
        using std::setw;
        using std::setfill;

        const int col_width = 15;
        const int table_width = 10*col_width;

        std::cout << setfill('-') << setw(table_width) << "" << std::endl;
        std::cout << std::left << setfill(' ') << setw(table_width/2 - 7) << "" << "MPC Statistics" << std::endl;
        std::cout << setfill('-') << setw(table_width) << "" << std::endl;

        std::cout << setfill(' ');
        std::cout << setw(col_width) << "Solve #"
                << setw(col_width) << "Time (ms)"
                << setw(col_width) << "Constraints"
                << setw(col_width) << "Step Norm"
                << setw(col_width) << "Alpha"
                << setw(col_width) << "Cost"
                << setw(col_width) << "Merit"
                << setw(col_width) << "Merit dd"
                << setw(col_width) << "Solve Type"
                << setw(col_width) << "QP Cost" << std::endl;
        std::cout << std::setfill('-') << setw(table_width) << "" << std::endl;
        std::cout << setfill(' ');
        for (int i = 0; i < alpha_.size(); i++) {

            for (int j = 0; j < contact_sched_change_.size(); j++) {
                if (i == contact_sched_change_.at(j)) {
                    std::cout << "Contact schedule changed " << setw(table_width - 25) << std::endl;
                }
            }

            std::string solve_type;
            switch (solve_type_.at(i)) {
                case Solved:
                    solve_type = "Solved";
                    break;
                case SolvedInacc:
                    solve_type = "Solved Inacc";
                    break;
                case PrimalInfeasible:
                    solve_type = "P - Infeasible";
                    break;
                case DualInfeasible:
                    solve_type = "D - Infeasible";
                    break;
                case PrimalInfeasibleInacc:
                    solve_type = "P - Infeasible Inacc";
                    break;
                case DualInfeasibleInacc:
                    solve_type = "D - Infeasible Inacc";
                    break;
                case MaxIter:
                    solve_type = "Max Iter";
                    break;
                case Unsolved:
                    solve_type = "Unsolved";
                    break;
                default:
                    solve_type = "Other";
            }

            std::cout << setw(col_width) << i
                      << setw(col_width) << solve_time_.at(i)
                      << setw(col_width) << equality_constraint_violations_.at(i)
                      << setw(col_width) << step_norm_.at(i)
                      << setw(col_width) << alpha_.at(i)
                      << setw(col_width) << cost_result_.at(i)
                      << setw(col_width) << merit_result_.at(i)
                      << setw(col_width) << merit_directional_deriv_.at(i)
                      << setw(col_width) << solve_type
                      << setw(col_width) << cost_.at(i) << std::endl;
        }
        std::cout << std::endl;
    }

    controller::Contact MPC::GetDesiredContacts(double time) const {
        std::vector<bool> in_contact = prev_traj_.GetContacts(time);
        std::vector<int> contact_frames = model_.GetContactFrames();
        controller::Contact contact;
        contact.in_contact_ = in_contact;
        contact.contact_frames_ = contact_frames;
        return contact;
    }

    int MPC::GetNodeIntersectMutableForces() const {
        int count = 0;
        for (int node = 0; node < info_.num_nodes+1; node++) {
            for (int ee = 0; ee < num_ee_; ee++) {
                const int coord = 2;
                if (prev_traj_.IsForceMutable(ee, GetTime(node))) {
                    count++;
                }
            }
        }

        return count;
    }

    Trajectory MPC::GetTrajectory() const {
        return prev_traj_;
    }

    int MPC::GetNode(double time) const {
        return std::ceil((time - init_time_)/info_.integrator_dt);
    }

    int MPC::GetNumDecisionVars() const {
        return data_.num_decision_vars;
    }

    int MPC::GetNumConstraints() const {
        return data_.GetTotalNumConstraints();
    }

    vector_t MPC::Getdx() {
        if (solve_type_.at(solve_type_.size()-1) == Solved) {
            return qp_solver->Getdx();
        } else {
            return vector_t::Zero(data_.num_decision_vars); // TODO: Change this return
        }
    }

    bool MPC::ComputeDerivativeTerms() {
        if (solve_type_.at(solve_type_.size()-1) == Solved) {
            vector_t dx = qp_solver->Computedx(data_.sparse_cost_, data_.cost_linear, prev_qp_sol);
            vector_t dy_lu = vector_t::Zero(data_.GetTotalNumConstraints());
            qp_solver->SetupDerivativeCalcs(dx, dy_lu, dy_lu, data_);
            return true;
        } else {
            return false;
        }
    }

    bool MPC::GetQPPartials(QPPartialsDense& partials) const {
        if (solve_type_.at(solve_type_.size()-1) == Solved) {
            partials.SetZero();

            qp_solver->CalcDerivativeWrtMats(partials.dP, partials.dA, partials.dG);
            qp_solver->CalcDerivativeWrtVecs(partials.dq, partials.db, partials.dh);

            return true;
        } else {
            return false;
        }
    }


    vector_t MPC::GetQPSolution() const {
        return prev_qp_sol;
    }

    int MPC::UpdateNumInputs() {
        num_inputs_ = prev_traj_.GetTotalPosSplineVars() + prev_traj_.GetTotalForceSplineVars();
        if (model_.UsesJoints()) {
            num_inputs_ += info_.num_nodes * model_.GetNumJoints();
        }

        return num_inputs_;
    }

    void MPC::UpdateContactTimes(std::vector<time_v>& contact_times) {
        prev_traj_.UpdateContactTimes(contact_times);
        contact_sched_change_.push_back(run_num_);
    }

    void MPC::AddDiagonalCost() {
        // Adds a positive small cost to every value so that Q is PD
        // TODO: 1e-2 seemed to work well!
        // TODO: The gradient descent appears to be better with this more PD!
        data_.cost_mat_.SetDiagonalMatrix(1e-3,0,0,data_.num_decision_vars); // 1e-3
    }

    const QPData& MPC::GetQPData() const {
        return data_;
    }

    int MPC::GetNumForceBoxConstraints() const {
        const std::vector<time_v> contact_times = prev_traj_.GetContactTimes();
        int num_constraints = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int i = 0; i < contact_times.at(ee).size() - 1; i++) {
                if (contact_times.at(ee).at(i).GetType() == TouchDown) {
                    num_constraints += 2*FB_PER_FORCE;
                }
            }
        }

        return num_constraints;
    }

    int MPC::GetNumFricConeConstraints() const {
        const std::vector<time_v> contact_times = prev_traj_.GetContactTimes();
        int num_constraints = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int i = 0; i < contact_times.at(ee).size() - 1; i++) {
                if (contact_times.at(ee).at(i).GetType() == TouchDown) {
                    num_constraints += 4*FB_PER_FORCE;
                }
            }
        }

        return num_constraints;
    }

    void MPC::SetQPData(const mpc::QPData& data) {
        data_ = data;
    }

    MPC& MPC::operator=(const mpc::MPC& mpc) {
        if (this != &mpc) {
            model_ = mpc.model_;
            data_ = mpc.data_;
            info_ = mpc.info_;
            num_states_ = mpc.num_states_;
            num_ee_ = mpc.num_ee_;
            num_inputs_ = mpc.num_inputs_;
            friction_pyramid_ = mpc.friction_pyramid_;
            qp_solver = std::make_unique<ClarabelInterface>(*mpc.qp_solver);
            prev_traj_ = mpc.prev_traj_;
            Q_ = mpc.Q_;
            w_ = mpc.w_;
            Phi_ = mpc.Phi_;
            Phi_w_ = mpc.Phi_w_;
            Q_forces_ = mpc.Q_forces_;
            prev_qp_sol = mpc.prev_qp_sol;
            init_time_ = mpc.init_time_;
            force_weight_ = mpc.force_weight_;
            run_num_ = mpc.run_num_;
            line_search_res_ = mpc.line_search_res_;
            prev_dual_sol_ = mpc.prev_dual_sol_;
            mu_ = mpc.mu_;

            // TODO: Do the recording stuff

            in_real_time_ = mpc.in_real_time_;
            A_ = mpc.A_;
            B_ = mpc.B_;
            C_ = mpc.C_;
            C2_ = mpc.C2_;
            using_clarabel_ = mpc.using_clarabel_;
            integrator_ = mpc.integrator_;
            constraint_idx_ = mpc.constraint_idx_;
            num_ineq_fk_ = mpc.num_ineq_fk_;
        }

        return *this;
    }

    MPC::MPC(const mpc::MPC& other) : model_(other.model_), prev_traj_(other.prev_traj_),
        integrator_(other.integrator_) {
        *this = other;

        UpdateNumInputs();

    }

    void MPC::SetVerbosityLevel(mpc::MPCVerbosityLevel verbosity) {
        info_.verbose = verbosity;

        if (qp_solver != nullptr) {
            if (info_.verbose == Optimization || info_.verbose == All) {
                qp_solver->SetVerbosity(true);
            } else {
                qp_solver->SetVerbosity(false);
            }
        }
    }

    void MPC::AdjustForCurrentContacts(double time, const controller::Contact& contact) {
        const controller::Contact& traj_contact = prev_traj_.GetDesiredContacts(time);
        for (int ee = 0; ee < num_ee_; ee++) {
            if (contact.in_contact_.at(ee) && !traj_contact.in_contact_.at(ee) &&
            std::abs(prev_traj_.GetNextContactTime(ee, time) - time) < 7e-2) {
                prev_traj_.SetEEInContact(ee, time);
            }
        }
    }

    int MPC::GetNumTDConstraints() const {
        int num_constraints = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            if (prev_traj_.GetNextContactTime(ee, init_time_) - init_time_ < prev_traj_.GetCurrentSwingTime(ee)/2) {
                num_constraints += 2;
            }
        }

        return num_constraints;
    }
} // mpc