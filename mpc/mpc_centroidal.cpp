//
// Created by zolkin on 1/19/24.
//


#include <array>
#include <Eigen/SparseCore>

#include "mpc.h"
#include "qp/osqp_interface.h"
#include "mpc_centroidal.h"

namespace mpc {

//    using triplet_t = Eigen::Triplet<double>;
//
//    MPCInfo::MPCInfo() {}
//
//    MPCInfo::MPCInfo(const MPCInfo& info) {
//        num_nodes = info.num_nodes;
//        num_qp_iterations = info.num_qp_iterations;
//        num_contacts = info.num_contacts;
//        friction_coef = info.friction_coef;
//        vel_bounds = info.vel_bounds;
//        joint_bounds_lb = info.joint_bounds_lb;
//        joint_bounds_ub = info.joint_bounds_ub;
//        ee_frames = info.ee_frames;
//        discretization_steps = info.discretization_steps;
//        num_switches = info.num_switches;
//        integrator_dt = info.integrator_dt;
//        force_bound = info.force_bound;
//        swing_height = info.swing_height;
//        foot_offset = info.foot_offset;
//    }
//
//    MPC::MPC(const MPCInfo& info, const std::string& robot_urdf) :
//            info_(info),
//            model_(robot_urdf, info.ee_frames, info.discretization_steps, info.integrator_dt),
//            num_ee_(model_.GetNumEndEffectors()),
//            prev_traj_(info.num_nodes + 1, num_states_ + 1, model_.UsesJoints(),
//                       CreateDefaultSwitchingTimes(info.num_switches, num_ee_,
//                                                   info.integrator_dt * (info.num_nodes)),
//                       info.integrator_dt, info.swing_height, info.foot_offset),
//            constraint_projection_(false),
//            data_(false, 25000, 2000) {
//
//        assert(info_.ee_frames.size() == num_ee_);
//
////        if (info_.vel_bounds.size() != num_joints_ || info_.joint_bounds_lb.size() != num_joints_
////        || info_.joint_bounds_ub.size() != num_joints_) {
////            throw std::runtime_error("Velocity or joint bounds do not match the number of joints on the robot.");
////        }
//
////        num_ineq_fk_ = 2;
//
//        // --------------------------------- //
//        // k = # of nodes
//        // QP Decision vector:
//        // [s_{n1}, ..., s_{nk}, u_{splines}, u_{j1}, ..., u_{jk}, pos_{splines}]
//        // ^ in words: the state at each node, then the splines (forces), then the joint velocities at each node,
//        // then the end effector position spline variables.
//        // --------------------------------- //
//
//        // state vector: [lcom, kcom, qb, q0, ..., qnj]
//        // input vector: [f1, ..., fnee, v0, ..., vnj]. Note each f is a set of 3 splines.
//
//        init_time_ = 0;
//
//        SetInitQPSizes();
//
//        qp_solver = std::make_unique<OSQPInterface>(data_, false);
//
//        SetFrictionPyramid();
//
//        Phi_ = matrix_t::Zero(num_states_, num_states_);
//        Phi_w_ = vector_t::Zero(num_states_);
//
//        prev_qp_sol = vector_t::Zero(data_.num_decision_vars);
//
//        line_search_res_ = vector_t::Zero(data_.num_decision_vars);
//        mu_ = 5000;
//        run_num_ = 0;
//        constraint_idx_ = 0;
//
//        in_real_time_ = false;
//
//        // Set everything to 0
//        data_.InitQPMats();
//    }
//
//    Trajectory MPC::CreateInitialRun(const mpc::vector_t& state) {
//        bool converged = false;
//        int num_iter = 0;
//        in_real_time_ = false;
//
//        qp_solver->ConfigureForInitialRun();
//
//        while (!converged && num_iter < 10) {
//            Solve(state, 0);
//            num_iter++;
//        }
//        return prev_traj_;
//    }
//
//    Trajectory MPC::GetRealTimeUpdate(double run_time_iters, const vector_t& state, double init_time) {
//        if (in_real_time_) {
//            return Solve(state, init_time);
//        } else {
//            qp_solver->ConfigureForRealTime(run_time_iters);
//            in_real_time_ = true;
//            return Solve(state, init_time);
//        }
//    }
//
//    Trajectory MPC::Solve(const vector_t& state, double init_time) {
//        utils::Timer solve_timer("mpc solve");
//        solve_timer.StartTimer();
//
//        assert(state.size() == num_states_ + 1);
//
//        init_time_ = init_time;
//
//        utils::Timer poly_update_timer("spline update");
//        poly_update_timer.StartTimer();
//        prev_traj_.SetInitTime(init_time);
//        prev_traj_.AddPolys(info_.integrator_dt * info_.num_nodes + init_time);
//        prev_traj_.RemoveUnusedPolys(init_time);
//        poly_update_timer.StopTimer();
//
//        utils::Timer data_update_timer("data update");
//        data_update_timer.StartTimer();
//        UpdateQPSizes();
//        data_.InitQPMats();
//
//        prev_traj_.SetState(0, state);
//        prev_qp_sol = prev_traj_.ConvertToQPVector();
//        data_update_timer.StopTimer();
//
//        assert(prev_qp_sol.size() == data_.num_decision_vars);
//
//        utils::Timer constraint_costs_timer("constraints and costs");
//        constraint_costs_timer.StartTimer();
//        // ----------------------- Costs ------------------------- //
//        AddHessianApproxCost();
//        AddGradientCost();
//        AddFinalCost();
//
//        // -------------------- Constraints ---------------------- //
//        constraint_idx_ = 0;
//        utils::Timer dynamics_timer("dynamics constraints");
//        dynamics_timer.StartTimer();
//        AddDynamicsConstraints(state);  // TODO: Speed up
//        dynamics_timer.StopTimer();
//
////        dynamics_timer.PrintElapsedTime();
//        if (constraint_projection_) {
//            AddFKConstraintsProjection(state);
//        } else {
//            AddFKConstraints(state);
//        }
//
//        AddFrictionConeConstraints();
//        AddBoxConstraints();
//        AddForceBoxConstraints();
//        constraint_costs_timer.StopTimer();
//
//        data_.ConstructSparseMats();
//        data_.ConstructVectors();
//        if (constraint_projection_) {
//            data_.ApplyProjection();
//        }
//
//        // ----------------------- Solve ------------------------- //
//        qp_solver->SetupQP(data_, prev_qp_sol);
//        utils::Timer qp_solve_timer("QP solve");
//        qp_solve_timer.StartTimer();
//        const vector_t sol = qp_solver->Solve(data_);
//        qp_solve_timer.StopTimer();
//
//        if (qp_solver->GetSolveQuality() != "Solved Inaccurate" && qp_solver->GetSolveQuality() != "Solved"
//            && qp_solver->GetSolveQuality() != "Max Iter Reached") {
//            std::cerr << "Warning: " << qp_solver->GetSolveQuality() << std::endl;
////            throw std::runtime_error("Bad solve.");
//        }
//
//        std::cout << "Solve type: " << qp_solver->GetSolveQuality() << std::endl;
//
//        const vector_t p = sol - prev_qp_sol;
//        double max = 0;
//        int max_idx = 0;
//        for (int i = 0; i < p.size(); i++) {
//            if (std::abs(sol(i) - prev_qp_sol(i)) > max) {
//                max = std::abs(sol(i) - prev_qp_sol(i));
//                max_idx = i;
//            }
//        }
//
////        std::cout << "max difference is: " << max << " and occurs at index: " << max_idx << std::endl;
////        std::cout << "dynamics variables: " << num_states_*(info_.num_nodes+1) << std::endl;
////        std::cout << "force variables: " << GetVelocityIndex(0) << std::endl;
//
//        utils::Timer line_search_timer("line search");
//        double alpha = 0;
//        if (sol.size() == prev_qp_sol.size()) {
//            line_search_timer.StartTimer();
//            alpha = LineSearch(p, state);
//            line_search_timer.StopTimer();
//        }
//
////        if (run_num_ > 0){
////            std::cout << (alpha*(qp_solver->GetDualSolution() - prev_dual_sol_) + prev_dual_sol_).lpNorm<Eigen::Infinity>() << std::endl;
////        }
//
//        prev_dual_sol_ = qp_solver->GetDualSolution();
//
//        prev_qp_sol = ((alpha * p) + prev_qp_sol).eval();
//
//        prev_traj_ = ConvertQPSolToTrajectory(prev_qp_sol, state);
//        for (int node = 0; node < info_.num_nodes; node++) {
//            prev_traj_.UpdateFullVelocity(node, GetTargetVelocity(GetTime(node)));
//        }
//
//        solve_timer.StopTimer();
//
//        utils::Timer stats_timer("recording stats");
//        stats_timer.StartTimer();
//        RecordStats(alpha, p, qp_solver->GetSolveQuality(), state,
//                    solve_timer.GetElapsedTimeMilliseconds());
//        stats_timer.StopTimer();
//
//        run_num_++;
//
//        constraint_costs_timer.PrintElapsedTime();
//        line_search_timer.PrintElapsedTime();
//        poly_update_timer.PrintElapsedTime();
//        data_update_timer.PrintElapsedTime();
//        qp_solve_timer.PrintElapsedTime();
////        stats_timer.PrintElapsedTime();
//        solve_timer.PrintElapsedTime();
//        std::cout << "-----------" << std::endl;
//
////        prev_traj_.PrintTrajectoryToFile("mpc_demo_traj.txt");
//
//        return prev_traj_;
//    }
//
//    void MPC::SetWarmStartTrajectory(const mpc::Trajectory& trajectory) {
//        prev_traj_ = trajectory;
//
//        prev_qp_sol = prev_traj_.ConvertToQPVector();
//    }
//
//    void MPC::SetQuadraticCostTerm(const matrix_t& Q) {
//        if (Q.rows() != num_states_ || Q.cols() != num_states_) {
//            throw std::runtime_error("Supplied quadratic cost term is the wrong size.");
//        }
//
//        Q_ = Q;
//    }
//
//    void MPC::SetLinearCostTerm(const vector_t& w) {
//        if (w.rows() != num_states_) {
//            throw std::runtime_error("Supplied linear cost term is the wrong size.");
//        }
//
//        w_ = w;
//    }
//
//    void MPC::SetQuadraticFinalCost(const mpc::matrix_t& Phi) {
//        if (Phi.rows() != num_states_ || Phi.cols() != num_states_) {
//            throw std::runtime_error("Wrong sized quadratic final cost.");
//        }
//
//        Phi_ = Phi;
//    }
//
//    void MPC::SetLinearFinalCost(const vector_t& w) {
//        if (w.size() != num_states_) {
//            throw std::runtime_error("Wrong sized linear final cost term.");
//        }
//
//        Phi_w_ = w;
//    }
//
//    void MPC::SetFrictionPyramid() {
//        // Assumes flat ground
//        Eigen::Vector3d h = {1, 0, 0};
//        Eigen::Vector3d l = {0, 1, 0};
//        Eigen::Vector3d n = {0, 0, 1};
//
//        friction_pyramid_ << (h - n * info_.friction_coef).transpose(),
//                -(h + n * info_.friction_coef).transpose(),
//                (l - n * info_.friction_coef).transpose(),
//                -(l + n * info_.friction_coef).transpose();
//    }
//
//    void MPC::AddDynamicsConstraints(const vector_t& state) {
//
//        // initial condition
//        data_.constraint_mat_.SetDiagonalMatrix(-1, 0, 0, num_states_);
//        data_.dynamics_constants.head(num_states_) =
//                -CentroidalModel::ConvertManifoldStateToAlgebraState(state, state);
//
//        const int force_spline_vars = prev_traj_.GetInputs().GetTotalForceSplineVars();
//
//        A_.resize(num_states_, num_states_);
//        B_.resize(num_states_, force_spline_vars + num_joints_);
//        C_.resize(num_states_);
//
//        utils::Timer get_dyn_timer("linear discrete dynamics");
//
//        for (int node = 0; node < info_.num_nodes; node++) {
//            double time = GetTime(node);
//            // A can go in normally, B needs to be split
//            get_dyn_timer.StartTimer();
//            model_.GetLinearDiscreteDynamics(prev_traj_.GetStates().at(node),
//                                             state, prev_traj_.GetInputs(), time, A_, B_, C_);
//            get_dyn_timer.StopTimer();
////            get_dyn_timer.PrintElapsedTime();
//            data_.constraint_mat_.SetMatrix(A_, (node + 1) * num_states_, node * num_states_);
//
//            data_.constraint_mat_.SetDiagonalMatrix(-1, (node + 1) * num_states_,
//                                                    (node + 1) * num_states_, num_states_);
//
//            const matrix_t B_spline = B_.leftCols(force_spline_vars);
//            const matrix_t B_vels = B_.rightCols(num_joints_);
//
//            data_.constraint_mat_.SetMatrix(B_vels, (node + 1) * num_states_, GetVelocityIndex(node));
//            data_.constraint_mat_.SetMatrix(B_spline, (node + 1) * num_states_, GetForceSplineStartIdx());
//            data_.dynamics_constants.segment((node + 1) * num_states_, num_states_) = -C_;
//        }
//        constraint_idx_ = (info_.num_nodes + 1) * num_states_;
//    }
//
//    void MPC::AddFKConstraints(const vector_t& state) {
//        matrix_t G;
//        vector_t g;
//
//        const int spline_offset = GetPosSplineStartIdx();
//        int idx_eq = 0;
//        int idx_ineq = 0;
//        const double margin = .01;
//
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            double time = GetTime(node);
//            vector_t state1 = prev_traj_.GetStates().at(node);
//            if (node == 0) {
//                state1 = state;
//            }
//
//            for (int ee = 0; ee < num_ee_; ee++) {
//                model_.GetFKLinearization(state1, state, prev_traj_.GetInputs(), ee, G, g);
//                if (node >= num_ineq_fk_) {
//                    data_.constraint_mat_.SetMatrix(G,
//                                                    constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_,
//                                                    node * num_states_ + CentroidalModel::MOMENTUM_OFFSET);
//
//                    data_.fk_constants_.segment(idx_eq, POS_VARS) = -g;
//                } else {
//                    // We only need the margin in the z direction, in x and y we can treat it as an equality
//                    data_.constraint_mat_.SetMatrix(G.bottomRows<1>(),
//                                                    constraint_idx_ + idx_ineq,
//                                                    node * num_states_ + CentroidalModel::MOMENTUM_OFFSET);
//
//                    data_.fk_lb_(idx_ineq) = -g(2) - ((num_ineq_fk_ - node) * margin);
//                    data_.fk_ub_(idx_ineq) = -g(2) + ((num_ineq_fk_ - node) * margin);
//
//
//                    data_.constraint_mat_.SetMatrix(G.topRows<2>(),
//                                                    constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_,
//                                                    node * num_states_ + CentroidalModel::MOMENTUM_OFFSET);
//
//                    data_.fk_constants_.segment(idx_eq, 2) = -g.head<2>();
//                }
//
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    if (prev_traj_.IsSplineMutable(ee, coord)) {
//                        int vars_index, vars_affecting;
//                        std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
//
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            data_.constraint_mat_.SetMatrix(
//                                    -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose(),
//                                    constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_,
//                                    spline_offset + vars_index - vars_affecting);
//                            idx_eq++;
//                        } else {
//                            data_.constraint_mat_.SetMatrix(
//                                    -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose(),
//                                    constraint_idx_ + idx_ineq,
//                                    spline_offset + vars_index - vars_affecting);
//                            idx_ineq++;
//                        }
//                    } else {
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            data_.fk_constants_(idx_eq) +=
//                                    prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
//                            idx_eq++;
//                        } else {
//                            data_.fk_lb_(idx_ineq) +=
//                                    prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
//                            data_.fk_ub_(idx_ineq) +=
//                                    prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
//                            idx_ineq++;
//                        }
//                    }
//                }
//            }
//        }
//        constraint_idx_ = constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_;
//    }
//
//    void MPC::AddFKConstraintsProjection(const vector_t& state) {
//        matrix_t G;
//        vector_t g;
//
//        const int spline_offset = GetPosSplineStartIdx();
//        int idx_eq = 0;
//        int idx_ineq = 0;
//        const double margin = .01;
//
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            double time = GetTime(node);
//            vector_t state1 = prev_traj_.GetStates().at(node);
//            if (node == 0) {
//                state1 = state;
//            }
//
//            for (int ee = 0; ee < num_ee_; ee++) {
//                model_.GetFKLinearization(state1, state, prev_traj_.GetInputs(), ee, G, g);
//                if (node >= num_ineq_fk_) {
//                    data_.constraint_projections_.constraint_mat_.block(idx_eq,
//                                                                        node * num_states_ +
//                                                                        CentroidalModel::MOMENTUM_OFFSET, G.rows(),
//                                                                        G.cols()) = G;
//
//                    data_.constraint_projections_.b_.segment(idx_eq, POS_VARS) = -g;
//                } else {
//                    // We only need the margin in the z direction, in x and y we can treat it as an equality
//                    data_.constraint_mat_.SetMatrix(G.bottomRows<1>(),
//                                                    constraint_idx_ + idx_ineq,
//                                                    node * num_states_ + CentroidalModel::MOMENTUM_OFFSET);
//
//                    data_.fk_lb_(idx_ineq) = -g(2) - ((num_ineq_fk_ - node) * margin);
//                    data_.fk_ub_(idx_ineq) = -g(2) + ((num_ineq_fk_ - node) * margin);
//
//
//                    data_.constraint_mat_.SetMatrix(G.topRows<2>(),
//                                                    constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_,
//                                                    node * num_states_ + CentroidalModel::MOMENTUM_OFFSET);
//
//                    data_.fk_constants_.segment(idx_eq, 2) = -g.head<2>();
//                }
//
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    if (prev_traj_.IsSplineMutable(ee, coord)) {
//                        int vars_index, vars_affecting;
//                        std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
//
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            data_.constraint_projections_.constraint_mat_.block(
//                                    idx_eq, spline_offset + vars_index - vars_affecting,
//                                    1, vars_affecting) = -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(
//                                    time).transpose();
//                            idx_eq++;
//                        } else {
//                            data_.constraint_mat_.SetMatrix(
//                                    -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose(),
//                                    constraint_idx_ + idx_ineq,
//                                    spline_offset + vars_index - vars_affecting);
//                            idx_ineq++;
//                        }
//                    } else {
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            data_.constraint_projections_.b_(idx_eq) +=
//                                    prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
//                            idx_eq++;
//                        } else {
//                            data_.fk_lb_(idx_ineq) +=
//                                    prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
//                            data_.fk_ub_(idx_ineq) +=
//                                    prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
//                            idx_ineq++;
//                        }
//                    }
//                }
//            }
//        }
//        constraint_idx_ = constraint_idx_ + data_.num_fk_ineq_constraints_;
//    }
//
//    void MPC::AddFrictionConeConstraints() {
//        // -------------------- Friction pyramid -------------------- //
//        int force_offset = GetForceSplineStartIdx();
//        int idx = 0;
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            double time = GetTime(node);
//            for (int ee = 0; ee < num_ee_; ee++) {
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    if (prev_traj_.GetInputs().IsForceMutable(ee, coord, time)) {
//                        vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);
//
//                        int vars_index, vars_affecting;
//                        std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time,
//                                                                                                          coord);
//
//                        for (int fric_con = 0; fric_con < 4; fric_con++) {
//                            // All the friction constraints are effected by all 3 coordinates of the end effectors
//                            data_.constraint_mat_.SetMatrix(friction_pyramid_(fric_con, coord) * vars_lin.transpose(),
//                                                            constraint_idx_ + idx + fric_con,
//                                                            force_offset + vars_index - vars_affecting);
//                            data_.friction_cone_ub_(idx + fric_con) = 0;
//                            data_.friction_cone_lb_(idx + fric_con) = -qp_solver->GetInfinity(1)(0);
//                        }
//                    }
//                }
//                idx += 4;
//            }
//        }
//        constraint_idx_ += idx;
//    }
//
//    void MPC::AddBoxConstraints() {
//        int idx = 0;
//        for (int node = 0; node < info_.num_nodes; node++) {
//            // Velocity bounds
//            data_.box_ub_.segment(idx, num_joints_) = info_.vel_bounds;
//            data_.box_lb_.segment(idx, num_joints_) = -info_.vel_bounds;
//            data_.constraint_mat_.SetDiagonalMatrix(1, constraint_idx_ + idx, GetVelocityIndex(node), num_joints_);
//            idx += num_joints_;
//
//            // Configuration bounds
//            data_.box_ub_.segment(idx, num_joints_) = info_.joint_bounds_ub;
//            data_.box_lb_.segment(idx, num_joints_) = info_.joint_bounds_lb;
//            data_.constraint_mat_.SetDiagonalMatrix(1, constraint_idx_ + idx, GetJointIndex(node), num_joints_);
//            idx += num_joints_;
//        }
//        constraint_idx_ += idx;
//    }
//
//    void MPC::AddForceBoxConstraints() {
//        int force_idx = GetForceSplineStartIdx();
//        int row_idx = 0;
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            double time = GetTime(node);
//            for (int ee = 0; ee < num_ee_; ee++) {
//                const int coord = 2;
////                for (int coord = 0; coord < POS_VARS; coord++) {
//                if (prev_traj_.GetInputs().IsForceMutable(ee, coord, time)) {
//                    int vars_index, vars_affecting;
//                    std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time,
//                                                                                                      coord);
//                    vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);
//
//                    data_.constraint_mat_.SetMatrix(vars_lin.transpose(),
//                                                    constraint_idx_ + row_idx,
//                                                    force_idx + vars_index - vars_affecting);
//                    if (coord == 2) {
//                        data_.force_box_lb_(row_idx) = 0;
//                    } else {
//                        data_.force_box_lb_(row_idx) = -info_.force_bound;
//                    }
//                    data_.force_box_ub_(row_idx) = info_.force_bound;
//                    row_idx++;
//                }
////                }
//            }
//        }
//        constraint_idx_ += row_idx;
//    }
//
//    void MPC::AddQuadraticTrackingCost(const vector_t& state_des, const matrix_t& Q) {
//        if (Q.rows() != num_states_ || Q.cols() != num_states_) {
//            throw std::runtime_error("Supplied quadratic cost term is the wrong size.");
//        }
//
//        Q_ = Q;
//        w_ = -1 * Q * state_des;
//    }
//
//    void MPC::AddHessianApproxCost() {
//        for (int node = 0; node < info_.num_nodes; node++) {
//            data_.cost_mat_.SetMatrix(Q_, node * num_states_, node * num_states_);
//        }
//        if (Q_forces_.size() > 0) {
////            const Inputs& input = prev_traj_.GetInputs();
////            matrix_t force_map(info_.num_nodes+1, prev_traj_.GetInputs().GetTotalForceSplineVars());
////            for (int node = 0; node < info_.num_nodes + 1; node++) {
////                for (int ee = 0; ee < num_ee_; ee++) {
////                    for (int coord = 0; coord < POS_VARS; coord++) {
////                        if (input.IsForceMutable(ee, coord, GetTime(node))) {
////                            int vars_idx, vars_affecting;
////                            std::tie(vars_idx, vars_affecting) = input.GetForceSplineIndex(ee,
////                                                                                                            GetTime(node),
////                                                                                                            coord);
////
////                            force_map.block(node, vars_idx - vars_affecting, 1,
////                                            vars_affecting) =
////                                    input.GetForces().at(ee).at(coord).GetPolyVarsLin(GetTime(node)).transpose();
////                        }
////                    }
////                }
////            }
////            matrix_t temp = force_map.transpose() * Q_forces_ * force_map;
//            data_.cost_mat_.SetMatrix(Q_forces_, GetForceSplineStartIdx(), GetForceSplineStartIdx());
//        }
//    }
//
//    void MPC::AddGradientCost() {
//        for (int node = 0; node < info_.num_nodes; node++) {
//            data_.cost_linear.segment(node * num_states_, num_states_) = w_;
//        }
//    }
//
//    void MPC::AddFinalCost() {
//        data_.cost_mat_.SetMatrix(Phi_, info_.num_nodes * num_states_, info_.num_nodes * num_states_);
//
//        data_.cost_linear.segment(info_.num_nodes * num_states_, num_states_) = Phi_w_;
//    }
//
//    std::vector<std::vector<double>> MPC::CreateDefaultSwitchingTimes(int num_switches, int num_ee, double horizon) {
//        std::vector<std::vector<double>> switching_times;
//        for (int i = 0; i < num_ee; i++) {
//            std::vector<double> times;
//            for (int j = 0; j < num_switches; j++) {
//                times.push_back((j + 1) * horizon / num_switches);
//            }
//            switching_times.push_back(times);
//        }
//
//        return switching_times;
//    }
//
//    int MPC::GetForceSplineStartIdx() const {
//        return num_states_ * (1 + info_.num_nodes);     // initial condition and states for each node
//    }
//
//    int MPC::GetVelocityIndex(int node) const {
//        return GetForceSplineStartIdx() + prev_traj_.GetInputs().GetTotalForceSplineVars() + node * num_joints_;
//    }
//
//    int MPC::GetJointIndex(int node) const {
//        return (node + 1) * num_states_ + CentroidalModel::FLOATING_VEL_OFFSET + CentroidalModel::MOMENTUM_OFFSET;
//    }
//
//    int MPC::GetPosSplineStartIdx() const {
//        return GetVelocityIndex(info_.num_nodes);
//    }
//
//    Trajectory MPC::ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const {
//        // Start by copying the current trajectory to keep the switching times and what not
//        Trajectory traj(prev_traj_);
//
//        // Assign all the spline information
//        int force_idx = GetForceSplineStartIdx();
//        int pos_idx = GetPosSplineStartIdx();
//        for (int ee = 0; ee < num_ee_; ee++) {
//            for (int coord = 0; coord < POS_VARS; coord++) {
//                int vars_in_spline = traj.GetInputs().GetForces().at(ee).at(coord).GetTotalPolyVars();
//                traj.UpdateForceSpline(ee, coord, qp_sol.segment(force_idx, vars_in_spline));
//                force_idx += vars_in_spline;
//
//                if (traj.IsSplineMutable(ee, coord)) {
//                    vars_in_spline = traj.GetPositions().at(ee).at(coord).GetTotalPolyVars();
//                    traj.UpdatePositionSpline(ee, coord, qp_sol.segment(pos_idx, vars_in_spline));
//                    pos_idx += vars_in_spline;
//                }
//            }
//        }
//
//        // Set states and joint vels
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            // Need to convert the state back to the manifold valued state
//            vector_t man_state = CentroidalModel::ConvertAlgebraStateToManifoldState(
//                    qp_sol.segment(node * num_states_, num_states_), init_state);
//            // need to normalize quaternion returned by MPC
//            Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(man_state.segment(9, 4)));
//            // Note the warning on the pinocchio function!
//            pinocchio::quaternion::firstOrderNormalize(quat);
//
//            man_state(9) = quat.x();
//            man_state(10) = quat.y();
//            man_state(11) = quat.z();
//            man_state(12) = quat.w();
//
//            traj.SetState(node, man_state);
//        }
//
//        for (int node = 0; node < info_.num_nodes; node++) {
//            traj.SetInputVels(node, qp_sol.segment(GetVelocityIndex(node), num_joints_));
//        }
//
//        return traj;
//    }
//
//    void MPC::SetInitQPSizes() {
//        data_.num_decision_vars = (info_.num_nodes + 1) * num_states_ + info_.num_nodes * num_joints_ +
//                                  prev_traj_.GetInputs().GetTotalForceSplineVars() + prev_traj_.GetTotalPosSplineVars();
//        data_.num_dynamics_constraints = (info_.num_nodes + 1) * num_states_;
//
//        data_.num_cone_constraints_ = (info_.num_nodes + 1) * 4 * num_ee_;
//        data_.num_box_constraints_ = info_.num_nodes * 2 * num_joints_;
//        data_.num_force_box_constraints_ = GetNodeIntersectMutableForces();
//
//        data_.num_fk_ineq_constraints_ = num_ee_ * num_ineq_fk_;
//        data_.num_fk_constraints_ = POS_VARS * num_ee_ * (info_.num_nodes + 1) - data_.num_fk_ineq_constraints_;
//    }
//
//    void MPC::UpdateQPSizes() {
//        data_.num_decision_vars = (info_.num_nodes + 1) * num_states_ + info_.num_nodes * num_joints_ +
//                                  prev_traj_.GetInputs().GetTotalForceSplineVars() + prev_traj_.GetTotalPosSplineVars();
//
//        data_.num_force_box_constraints_ = GetNodeIntersectMutableForces();
//    }
//
//    void MPC::SetDefaultGaitTrajectory(Gaits gait, int num_polys, const std::array<std::array<double, 3>, 4>& ee_pos) {
//        if (num_ee_ != 4) {
//            throw std::runtime_error("Default gaits have only been implemented for 4 end effectors!");
//        }
//
//        std::vector<std::vector<double>> switching_times;
//
//        switch (gait) {
//            case Trot: {
//                std::vector<double> times;
//                times.push_back(0.2);
//                times.push_back(0.4);
//                times.push_back(0.6);
//                times.push_back(0.8);
//
//                Spline position1(2, times, true, Spline::Normal);
//                Spline position2(2, times, false, Spline::Normal);
//
//                Spline force1(num_polys, times, false, Spline::Force);
//                Spline force2(num_polys, times, true, Spline::Force);
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
//
//                prev_traj_.PrintTrajectoryToFile("trot_test.txt");
//                break;
//            }
//            case Amble: {
//                throw std::runtime_error("Amble not implemented yet!");
//                break;
//            }
//            case Static_Walk: {
//                throw std::runtime_error("Static Walk not implemented yet!");
//                break;
//            }
//            default:
//                throw std::runtime_error("Unsupported gait.");
//        }
//
//    }
//
//    void MPC::SetStateTrajectoryWarmStart(const std::vector<vector_t>& states) {
//        assert(states.size() == info_.num_nodes + 1 && states.at(0).size() == prev_traj_.GetStates().at(0).size());
//
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            prev_traj_.SetState(node, states.at(node));
//        }
//    }
//
//    vector_t MPC::GetTargetConfig(double time) const {
//        int node = floor((time - init_time_) / info_.integrator_dt);
//        return prev_traj_.GetStates().at(node).tail(num_states_ - 5);
//    }
//
//    vector_t MPC::GetNextTargetConfig() const {
//        return prev_traj_.GetState(1).tail(num_states_ - 5);
//    }
//
//    vector_t MPC::GetFullTargetState(double time) const {
//        int node = floor((time - init_time_) / info_.integrator_dt);
//        return prev_traj_.GetStates().at(node);
//    }
//
//    vector_t MPC::GetTargetVelocity(double time) const {
//        int node = floor((time - init_time_) / info_.integrator_dt);
//        return model_.ComputeBaseVelocities(prev_traj_.GetStates().at(node),
//                                            prev_traj_.GetInputs().GetVels(time));
//    }
//
//    vector_t MPC::GetForceTarget(double time) const {
//        controller::Contact contacts = GetDesiredContacts(time);
//        vector_t forces(contacts.GetNumContacts() * 3);
//        int idx = 0;
//        for (int ee = 0; ee < num_ee_; ee++) {
//            if (contacts.in_contact_.at(ee)) {
//                forces.segment(idx, POS_VARS) = prev_traj_.GetInputs().GetForce(ee, time);
//                idx += 3;
//            }
//        }
//        return forces;
//    }
//
//    // TODO: Check/think about this
//    vector_t MPC::GetTargetAcc(double time) const {
//        int node = floor((time - init_time_) / info_.integrator_dt);
//        return (model_.ComputeBaseVelocities(prev_traj_.GetStates().at(node),
//                                             prev_traj_.GetInputs().GetVels(time)) -
//                model_.ComputeBaseVelocities(prev_traj_.GetStates().at(node + 1),
//                                             prev_traj_.GetInputs().GetVels(time + info_.integrator_dt))) /
//               info_.integrator_dt;
//    }
//
//    const CentroidalModel& MPC::GetModel() const {
//        return model_;
//    }
//
//    double MPC::LineSearch(const vector_t& direction, const vector_t& init_state) {
//        double alpha = 1;
//
//        // Note: for now just using the equality constraints on the merit function
//        double merit = GetMeritValue(prev_qp_sol, mu_, init_state);
//        double merit_step = GetMeritValue(alpha * direction + prev_qp_sol, mu_, init_state);
//        double merit_directional = GetMeritGradient(prev_qp_sol, direction, mu_, init_state);
//
//        int i = 0;
//        while ((merit - merit_step) < -0.00001 * alpha * merit_directional && i < 5) {
//            alpha *= 0.5;
//            vector_t temp = (alpha * direction) + prev_qp_sol;
//            merit_step = GetMeritValue(temp, mu_, init_state);
//            i++;
//        }
//
//        return alpha;
//    }
//
//    double MPC::GetMeritValue(const vector_t& x, double mu, const vector_t& init_state) {
//        const Trajectory temp_traj = ConvertQPSolToTrajectory(x, init_state);
//
//        return mu * GetEqualityConstraintValues(temp_traj, init_state).lpNorm<1>() + GetCostValue(x);
//    }
//
//    double MPC::GetMeritValue(const Trajectory& traj, double mu, const mpc::vector_t& init_state) {
//        return mu * GetEqualityConstraintValues(traj, init_state).lpNorm<1>() + GetCostValue(traj.ConvertToQPVector());
//    }
//
//    double MPC::GetCostValue(const vector_t& x) const {
//        return 0.5 * x.dot(data_.sparse_cost_ * (x)) + data_.cost_linear.dot(x);
//    }
//
//    // TODO: Optimize this and make it so it is called less
//    vector_t MPC::GetEqualityConstraintValues(const Trajectory& traj, const vector_t& init_state) {
//        vector_t eq_constraints = vector_t::Zero(
//                data_.num_dynamics_constraints - num_states_); // data_.num_fk_constraints_
//
//        for (int node = 0; node < info_.num_nodes; node++) {
//            eq_constraints.segment(node * num_states_, num_states_) =
//                    CentroidalModel::ConvertManifoldStateToAlgebraState(traj.GetState(node + 1), init_state)
//                    - model_.GetDiscreteDynamics(
//                            CentroidalModel::ConvertManifoldStateToAlgebraState(traj.GetState(node), init_state),
//                            traj.GetInputs(), GetTime(node), init_state);
////            for (int ee = 0; ee < num_ee_; ee++) {
////                eq_constraints.segment(info_.num_nodes*num_states_ + node*3*num_ee_ + 3*ee, 3) =
////                        model_.GetEndEffectorLocationCOMFrame(traj.GetState(node), model_.GetEndEffectorFrame(ee))
////                        - traj.GetPosition(ee,GetTime(node));
////            }
//        }
//        return eq_constraints;
//    }
//
//    double MPC::GetTime(int node) const {
//        return node * info_.integrator_dt + init_time_;
//    }
//
//    double MPC::GetMeritGradient(const vector_t& x, const vector_t& p, double mu, const vector_t& init_state) {
//        const Trajectory temp_traj = ConvertQPSolToTrajectory(x, init_state);
//
//        return (data_.sparse_cost_ * x + data_.cost_linear).dot(p)
//               - mu * GetEqualityConstraintValues(temp_traj, init_state).lpNorm<1>();
//    }
//
//    // TODO: Make this weight only on force values, so we don't hit the derivatives
//    void MPC::AddForceCost(double weight) {
//        const int num_forces = prev_traj_.GetInputs().GetTotalForceSplineVars();
//        Q_forces_.resize(num_forces, num_forces);
//        for (int i = 0; i < num_forces; i++) {
//            Q_forces_(i, i) = weight * ((i + 1) % 2);
//        }
//    }
//
//    void MPC::RecordStats(double alpha, const vector_t& direction, const std::string& solve_type,
//                          const vector_t& ref_state, double solve_time) {
//        equality_constraint_violations_.push_back(GetEqualityConstraintValues(prev_traj_, ref_state).lpNorm<1>());
//        step_norm_.push_back(direction.norm());
//        alpha_.push_back(alpha);
//        cost_result_.push_back(GetCostValue(prev_qp_sol));
//        merit_result_.push_back(GetMeritValue(prev_traj_, mu_, ref_state));
//        merit_directional_deriv_.push_back(
//                GetMeritGradient(prev_qp_sol - alpha * direction, direction, mu_, ref_state));
//        solve_type_.push_back(solve_type);
//        ref_state_.push_back(ref_state);
//        solve_time_.push_back(solve_time);
//    }
//
//    void MPC::PrintStats() {
//        using std::setw;
//        using std::setfill;
//
//        const int col_width = 15;
//        const int table_width = 9 * col_width;
//
//        std::cout << setfill('-') << setw(table_width) << "" << std::endl;
//        std::cout << std::left << setfill(' ') << setw(table_width / 2 - 7) << "" << "MPC Statistics" << std::endl;
//        std::cout << setfill('-') << setw(table_width) << "" << std::endl;
//
//        std::cout << setfill(' ');
//        std::cout << setw(col_width) << "Solve #"
//                  << setw(col_width) << "Time (ms)"
//                  << setw(col_width) << "Constraints"
//                  << setw(col_width) << "Step Norm"
//                  << setw(col_width) << "Alpha"
//                  << setw(col_width) << "Cost"
//                  << setw(col_width) << "Merit"
//                  << setw(col_width) << "Merit dd"
//                  << setw(col_width) << "Solve Type" << std::endl;
//        std::cout << std::setfill('-') << setw(table_width) << "" << std::endl;
//        std::cout << setfill(' ');
//        for (int i = 0; i < alpha_.size(); i++) {
//            std::cout << setw(col_width) << i
//                      << setw(col_width) << solve_time_.at(i)
//                      << setw(col_width) << equality_constraint_violations_.at(i)
//                      << setw(col_width) << step_norm_.at(i)
//                      << setw(col_width) << alpha_.at(i)
//                      << setw(col_width) << cost_result_.at(i)
//                      << setw(col_width) << merit_result_.at(i)
//                      << setw(col_width) << merit_directional_deriv_.at(i)
//                      << setw(col_width) << solve_type_.at(i) << std::endl;
//        }
//        std::cout << std::endl;
//    }
//
//    controller::Contact MPC::GetDesiredContacts(double time) const {
//        std::vector<bool> in_contact = prev_traj_.GetContacts(time);
//        std::vector<int> contact_frames = model_.GetContactFrames();
//        controller::Contact contact;
//        contact.in_contact_ = in_contact;
////        contact.in_contact_.at(3) = !contact.in_contact_.at(3); // RR
////        contact.in_contact_.at(1) = !contact.in_contact_.at(1); // FR
////        contact.in_contact_.at(2) = !contact.in_contact_.at(2); // RL
////        contact.in_contact_.at(0) = !contact.in_contact_.at(0); // FL
//        contact.contact_frames_ = contact_frames;
//        return contact;
//    }
//
//    int MPC::GetNodeIntersectMutableForces() const {
//        int count = 0;
//        for (int node = 0; node < info_.num_nodes + 1; node++) {
//            for (int ee = 0; ee < num_ee_; ee++) {
//                const int coord = 2;
////                for (int coord = 0; coord < POS_VARS; coord++) {
//                if (prev_traj_.GetInputs().IsForceMutable(ee, coord, GetTime(node))) {
//                    count++;
//                }
////                }
//            }
//        }
//
//        return count;
//    }
//
//    Trajectory MPC::GetTrajectory() const {
//        return prev_traj_;
//    }
//
//    int MPC::GetNode(double time) const {
//        return std::ceil((time - init_time_) / info_.integrator_dt);
//    }
//
//    int MPC::GetNumDecisionVars() const {
//        return data_.num_decision_vars;
//    }
//
//    int MPC::GetNumConstraints() const {
//        return data_.GetTotalNumConstraints();
//    }
//
//    vector_t MPC::Getdx() {
//        if (solve_type_.at(solve_type_.size() - 1) == "Solved") { // TODO: Change to not be text based
//            return qp_solver->Getdx();
//        } else {
//            return vector_t::Zero(data_.num_decision_vars); // TODO: Change this return
//        }
//    }
//
//    bool MPC::ComputeDerivativeTerms() {
//        if (solve_type_.at(solve_type_.size() - 1) == "Solved") { // TODO: Change to not be text based
//            qp_solver->Computedx(data_.sparse_cost_, data_.cost_linear, prev_qp_sol);
//            vector_t dx = qp_solver->Getdx();
//            vector_t dy_lu = vector_t::Zero(data_.GetTotalNumConstraints());
//            qp_solver->SetupDerivativeCalcs(dx, dy_lu, dy_lu);
//            return true;
//        } else {
//            return false;
//        }
//    }
//
//    bool MPC::GetQPPartials(QPPartials& partials) const {
//        if (solve_type_.at(solve_type_.size() - 1) == "Solved") {
//            partials.dP.setZero();
//            partials.dA.setZero();
//            partials.dq.setZero();
//            partials.dl.setZero();
//            partials.du.setZero();
//
//            qp_solver->CalcDerivativeWrtMats(partials.dP, partials.dA);
//            qp_solver->CalcDerivativeWrtVecs(partials.dq, partials.dl, partials.du);
//
//            return true;
//        } else {
//            return false;
//        }
//    }
//
//    bool MPC::ComputeParamPartials(const Trajectory& traj, QPPartials& partials, int ee, int idx) {
//        if (solve_type_.at(solve_type_.size() - 1) == "Solved") {
//            QPData partial_data = data_;
//            partial_data.InitQPMats();
//
//            // For now the cost terms are not effected by the contact times
//
//            partials.dP.resize(data_.num_decision_vars, data_.num_decision_vars);
//            partials.dq.resize(data_.num_decision_vars);
//            partials.dA.resize(data_.GetTotalNumConstraints(), data_.num_decision_vars);
//            partials.du.resize(data_.GetTotalNumConstraints());
//            partials.dl.resize(data_.GetTotalNumConstraints());
//
//            // TODO: Make this more general
//            partials.dP *= 0;
//            partials.dq *= 0;
//
//
//            // Box constraints have no effect
//            // Force box constraints have no effect
//            // Cone constraints have on effect
//
//            // Dynamics constraints have an effect
////            utils::SparseMatrixBuilder A;
//            partial_data.constraint_mat_.Reserve(data_.sparse_constraint_.nonZeros());
////        const int contact_nodes = prev_traj_.GetNumContactNodes(ee);
//            matrix_t Adyn_deriv, Bdyn_deriv;
//            vector_t Cdyn_deriv;
//
//            for (int node = 0; node < info_.num_nodes; node++) {
//                matrix_t A, B;
//                vector_t C;
//                model_.GetLinearDiscreteDynamics(traj.GetState(node), traj.GetState(0), traj.GetInputs(),
//                                                 GetTime(node), A, B, C);
//
//                model_.ComputeLinearizationPartialWrtContactTimes(Adyn_deriv, Bdyn_deriv, Cdyn_deriv,
//                                                                  traj.GetState(node),
//                                                                  traj.GetInputs(), GetTime(node),
//                                                                  ee, idx);
//                // TODO: Is this sparsity pattern legit?
////                std::cout << "A: \n" << A << std::endl;
////                std::cout << "Adyn_deriv: \n" << Adyn_deriv << std::endl;
////                std::cout << "Bdyn_deriv: \n" << Bdyn_deriv << std::endl;
//
//                partial_data.constraint_mat_.SetMatrix(Adyn_deriv, (node + 1) * num_states_, node * num_states_);
//
//                const matrix_t B_spline = Bdyn_deriv.leftCols(traj.GetInputs().GetTotalForceSplineVars());
//                const matrix_t B_vels = Bdyn_deriv.rightCols(num_joints_);
//
//                partial_data.constraint_mat_.SetMatrix(B_vels, (node + 1) * num_states_, GetVelocityIndex(node));
//                partial_data.constraint_mat_.SetMatrix(B_spline, (node + 1) * num_states_, GetForceSplineStartIdx());
//
//                partial_data.dynamics_constants.segment((node + 1) * num_states_, num_states_) = -Cdyn_deriv;
//                partial_data.dynamics_constants.segment((node + 1) * num_states_, num_states_) = -Cdyn_deriv;
//            }
//
//
//            // FK constraints have an effect
//            const int spline_offset = GetPosSplineStartIdx();
//            int idx_eq = 0;
//            int idx_ineq = 0;
//
//            for (int node = 0; node < info_.num_nodes + 1; node++) {
//                const int time = GetTime(node);
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    if (traj.IsSplineMutable(ee, coord)) {
//                        int vars_index, vars_affecting;
//                        std::tie(vars_index, vars_affecting) = traj.GetPositionSplineIndex(ee, time, coord);
//
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            partial_data.constraint_mat_.SetMatrix(
//                                    -traj.GetPositions().at(ee).at(coord).ComputeCoefPartialWrtTime(time,
//                                                                                                    idx).transpose(),
//                                    constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_,
//                                    spline_offset + vars_index - vars_affecting);
//                            idx_eq++;
//                        } else {
//                            partial_data.constraint_mat_.SetMatrix(
//                                    -traj.GetPositions().at(ee).at(coord).ComputeCoefPartialWrtTime(time,
//                                                                                                    idx).transpose(),
//                                    constraint_idx_ + idx_ineq,
//                                    spline_offset + vars_index - vars_affecting);
//                            idx_ineq++;
//                        }
//                    } else {
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            partial_data.fk_constants_(idx_eq) +=
//                                    traj.GetPositions().at(ee).at(coord).ComputePartialWrtTime(time, idx);
//                            idx_eq++;
//                        } else {
//                            partial_data.fk_lb_(idx_ineq) +=
//                                    traj.GetPositions().at(ee).at(coord).ComputePartialWrtTime(time, idx);
//                            partial_data.fk_ub_(idx_ineq) +=
//                                    traj.GetPositions().at(ee).at(coord).ComputePartialWrtTime(time, idx);
//                            idx_ineq++;
//                        }
//                    }
//                }
//            }
//
//            partial_data.ConstructSparseMats();
//            partial_data.ConstructVectors();
//
//            partials.dA = partial_data.sparse_constraint_;
//            partials.dl = partial_data.lb_;
//            partials.du = partial_data.ub_;
//
//            return true;
//        } else {
//            return false;
//        }
//    }
//
//    vector_t MPC::GetQPSolution() const {
//        return prev_qp_sol;
//    }

} // mpc