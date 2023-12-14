//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <array>

#include "mpc.h"
#include "osqp_interface.h"

namespace mpc {

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
    }

    MPC::MPC(const MPCInfo& info, const std::string& robot_urdf) :
        info_(info),
        model_(robot_urdf, info.ee_frames, info.discretization_steps, info.integrator_dt),
        num_joints_(model_.GetNumJoints()),
        num_states_(model_.GetPinocchioNumConfig() + 5),
        num_ee_(model_.GetNumEndEffectors()),
        prev_traj_(info.num_nodes+1, num_states_+1, num_joints_,
               CreateDefaultSwitchingTimes(info.num_switches, num_ee_, info.integrator_dt*(info.num_nodes)),
               info.integrator_dt, info.swing_height) {

        assert(info_.ee_frames.size() == num_ee_);

        if (info_.vel_bounds.size() != num_joints_ || info_.joint_bounds_lb.size() != num_joints_
        || info_.joint_bounds_ub.size() != num_joints_) {
            throw std::runtime_error("Velocity or joint bounds do not match the number of joints on the robot.");
        }

        // --------------------------------- //
        // k = # of nodes
        // QP Decision vector:
        // [s_{n1}, ..., s_{nk}, u_{splines}, u_{j1}, ..., u_{jk}, pos_{splines}]
        // ^ in words: the state at each node, then the splines (forces), then the joint velocities at each node,
        // then the end effector position spline variables.
        // --------------------------------- //

        // state vector: [lcom, kcom, qb, q0, ..., qnj]
        // input vector: [f1, ..., fnee, v0, ..., vnj]. Note each f is a set of 3 splines.

        UpdateQPSizes();

        qp_solver = std::make_unique<OSQPInterface>(data_, false);

        SetFrictionPyramid();

        Phi_ = matrix_t::Zero(num_states_, num_states_);
        Phi_w_ = vector_t::Zero(num_states_);

        prev_qp_sol = vector_t::Zero(data_.num_decision_vars);

        init_time_ = 0;

        line_search_res_ = vector_t::Zero(data_.num_decision_vars);
        mu_ = 500;
        run_num_ = 0;
    }

    Trajectory MPC::Solve(const vector_t &centroidal_state, double init_time) {
        assert(centroidal_state.size() == num_states_ + 1);

        init_time_ = init_time;

        prev_traj_.SetInitTime(init_time);

        prev_traj_.AddPolys(info_.integrator_dt*info_.num_nodes + init_time);

        prev_traj_.RemoveUnusedPolys(init_time);

        UpdateQPSizes();

        // Reset to 0
        ResetQPMats();

        prev_traj_.SetState(0, centroidal_state);
        prev_qp_sol = prev_traj_.ConvertToQPVector();

        assert(prev_qp_sol.size() == data_.num_decision_vars);

        // initial condition
        data_.dynamics_constraints.topLeftCorner(num_states_, num_states_) = -matrix_t::Identity(num_states_, num_states_);
        data_.dynamics_constants.head(num_states_) =
                -CentroidalModel::ConvertManifoldStateToAlgebraState(centroidal_state, centroidal_state);

        // Constraints on all the nodes BUT the first
        for (int i = 0; i < info_.num_nodes; i++) {
            double time = GetTime(i);
            // ------------------------ Dynamics ------------------------ //
            AddDynamicsConstraints(centroidal_state, time, i);

            // ------------------------ Inequality Constraints ------------------------ //
//            AddBoxConstraints(centroidal_state, time, i);
        }

        // ----------------------- Costs ------------------------- //
        AddHessianApproxCost();
        AddGradientCost();
        AddFinalCost();

        // Constraints on ALL the nodes
        for (int i = 0; i < info_.num_nodes+1; i++) {
            double time = GetTime(i);
            // ------------------------ FK Constraint ------------------------ //
//            AddFKConstraints(centroidal_state, time, i);

            // ------------------------ Friction Cone Constraints ------------------------ //
//            AddFrictionConeConstraints(centroidal_state, time, i);

            // Positivity of z direction force
//            AddPositivityConstraints(time, i);

//            AddForceBoxConstraints(time, i);
        }

        // TODO: instead of adding a constraint consider removing them as decision variables
        // - None of these constraints will need to be added
        // - Above constraints that are effected by these will need to be changed
        // - Converting a trajectory to and from a QP vector will take some effort
        // - All the conversion functions for inputs and trajectories will need to be looked at
        // - These values will need to be set at some point
        // - When a new point in the spline is added it needs to be given the correct value
        // - Maybe denote a spline as mutable or not. For mutable splines denote certain polys as mutable
        // - Size of decision variables is different
//        AddSwingFootConstraints();  // All end effector z trajectories are set
//        AddForceConstraints();        // All force constants are set to 0
//        AddGroundIntersectConstraints();    // All position Z constants are set to 0


        qp_solver->SetupQP(data_, prev_qp_sol);
        vector_t sol = qp_solver->Solve(data_);
//        std::cout << "Dynamics lin violation on after solve: "
//                  << (data_.dynamics_constraints*sol - data_.dynamics_constants).lpNorm<1>() << std::endl;
//        std::cout << std::endl;


//        std::cout << "dynamics linearization violations norm: " << (data_.dynamics_constraints*sol - data_.dynamics_constants).lpNorm<1>() << std::endl;
//        Trajectory sol_traj = ConvertQPSolToTrajectory(sol, centroidal_state);
//        const double constraint_violation = GetEqualityConstraintValues(sol_traj, centroidal_state).lpNorm<1>();
//        std::cout << "dynamics violation norm: " << constraint_violation << std::endl;
//        std::cout << "state at node 31: \n" << sol_traj.GetState(31) << std::endl;
//        std::cout << "joint vel 0: \n" << sol_traj.GetInputs().GetVels(0.0) << std::endl;
//        std::cout << "Discrete dynamics for node 31: \n" <<
//            model_.GetDiscreteDynamics(CentroidalModel::ConvertManifoldStateToAlgebraState(sol_traj.GetState(31),
//                                   centroidal_state), sol_traj.GetInputs(), 31*info_.integrator_dt, centroidal_state) << std::endl;
//        std::cout << "Linear discrete dynamics for node 1: \n" << (data_.dynamics_constraints*sol - data_.dynamics_constants).segment(num_states_, num_states_) + sol.segment(2*num_states_, num_states_) << std::endl;
//        std::cout << "Generated state at node 32: \n" << sol.segment(32*num_states_, num_states_) << std::endl;

//        std::cout << "Force box violation (lb): " << data_.force_box_constraints_*sol - data_.force_box_lb_ << std::endl;

        vector_t p = sol - prev_qp_sol;
        double alpha = 0;
        if (sol.size() == prev_qp_sol.size()) {
            double max = 0;
            int max_index = 0;
            for (int i = 0; i < sol.size(); i++) {
                if (std::abs(sol(i) - prev_qp_sol(i)) > max) {
                    max = std::abs(sol(i) - prev_qp_sol(i));
                    max_index = i;
                }
            }
//            std::cout << "Max step: " << max << std::endl;
//            std::cout << "Max step index: " << max_index << std::endl;

            alpha = LineSearch(p, centroidal_state);
        }

//        if (run_num_ > 0){
//            std::cout << (alpha*(qp_solver->GetDualSolution() - prev_dual_sol_) + prev_dual_sol_).lpNorm<Eigen::Infinity>() << std::endl;
//        }

        prev_dual_sol_ = qp_solver->GetDualSolution();

        prev_qp_sol = alpha * p + prev_qp_sol;

        prev_traj_ = ConvertQPSolToTrajectory(prev_qp_sol, centroidal_state);

        prev_traj_.PrintTrajectoryToFile("prev_traj_after_solve.txt");

        RecordStats(alpha, p, qp_solver->GetSolveQuality(), centroidal_state);

        run_num_++;

        return prev_traj_;
    }

    void MPC::SetWarmStartTrajectory(const mpc::Trajectory &trajectory) {
        prev_traj_ = trajectory;

        prev_qp_sol = prev_traj_.ConvertToQPVector();
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

    void MPC::ResetQPMats() {
        // TODO: Only reset the ones that change
        data_.dynamics_constraints = matrix_t::Zero(data_.num_dynamics_constraints, data_.num_decision_vars);
        data_.dynamics_constants = vector_t::Zero(data_.num_dynamics_constraints);

        data_.cost_quadratic = matrix_t::Zero(data_.num_decision_vars, data_.num_decision_vars);
        data_.cost_linear = vector_t::Zero(data_.num_decision_vars);


        data_.fk_constraints_ = matrix_t::Zero(data_.num_fk_constraints_, data_.num_decision_vars);
        data_.fk_constants_ = vector_t::Zero(data_.num_fk_constraints_);

        data_.swing_force_constraints_ = matrix_t::Zero(data_.num_swing_foot_constraints_, data_.num_decision_vars);
        data_.swing_force_constants_ = vector_t::Zero(data_.num_swing_foot_constraints_);

        data_.foot_on_ground_constraints_ = matrix_t::Zero(data_.num_foot_on_ground_constraints_, data_.num_decision_vars);
        data_.foot_on_ground_constants_ = vector_t::Zero(data_.num_foot_on_ground_constraints_);

        data_.friction_cone_constraints_ = matrix_t::Zero(data_.num_cone_constraints_, data_.num_decision_vars);
        data_.friction_cone_lb_ = vector_t::Zero(data_.num_cone_constraints_);
        data_.friction_cone_ub_ = vector_t::Zero(data_.num_cone_constraints_);

//        data_.swing_trajectory_constraints_ = matrix_t::Zero(data_.num_swing_foot_pos_constraints_, data_.num_decision_vars);
//        data_.swing_trajectory_constants_ = vector_t::Zero(data_.num_swing_foot_pos_constraints_);

//        data_.positive_force_constraints_ = matrix_t::Zero(data_.num_positive_force_constraints_, data_.num_decision_vars);
//        data_.positive_force_lb_ = vector_t::Zero(data_.num_positive_force_constraints_);
//        data_.positive_force_ub_ = vector_t::Zero(data_.num_positive_force_constraints_);

        data_.box_constraints_ = matrix_t::Zero(data_.num_box_constraints_, data_.num_decision_vars);
        data_.box_lb_ = vector_t::Zero(data_.num_box_constraints_);
        data_.box_ub_ = vector_t::Zero(data_.num_box_constraints_);

        data_.foot_ground_inter_constraints_ = matrix_t::Zero(data_.num_foot_ground_inter_constraints_, data_.num_decision_vars);
        data_.foot_ground_inter_lb_ = vector_t::Zero(data_.num_foot_ground_inter_constraints_);
        data_.foot_ground_inter_ub_ = vector_t::Zero(data_.num_foot_ground_inter_constraints_);

        data_.force_box_constraints_ = matrix_t::Zero(data_.num_force_box_constraints_, data_.num_decision_vars);
        data_.force_box_lb_ = vector_t::Zero(data_.num_force_box_constraints_);
        data_.force_box_ub_ = vector_t::Zero(data_.num_force_box_constraints_);
    }

    void MPC::AddDynamicsConstraints(const vector_t& state, double time, int node) {
        matrix_t A, B;
        vector_t C;

        int force_spline_vars = prev_traj_.GetInputs().GetTotalForceSplineVars();

        // A can go in normally, B needs to be split
        model_.GetLinearDiscreteDynamics(prev_traj_.GetStates().at(node),
                                         state, prev_traj_.GetInputs(), time, A, B, C);
        data_.dynamics_constraints.block((node+1)*num_states_,
                                         node*num_states_,
                                         num_states_, num_states_) = A;
        data_.dynamics_constraints.block((node + 1) * num_states_,
                                         (node + 1) * num_states_,
                                         num_states_, num_states_) = -matrix_t::Identity(num_states_,
                                                                                         num_states_);

        matrix_t B_spline = B.leftCols(force_spline_vars);
        matrix_t B_vels = B.rightCols(num_joints_);

//        std::cout << "node: " << node << ", Bpsline: \n" << B_spline << std::endl;

        data_.dynamics_constraints.block((node+1)*num_states_,
                                         GetVelocityIndex(node),
                                         num_states_, num_joints_) = B_vels;
        data_.dynamics_constraints.block((node+1)*num_states_,
                                         GetForceSplineStartIdx(),
                                         num_states_, force_spline_vars) = B_spline;
        data_.dynamics_constants.segment((node+1)*num_states_, num_states_) = -C;

    }


    // TODO: When the nodes increase this start to look more sus
    // Need to confirm which legs should be fixed, etc...
    // TODO: Maybe just let the first constant of the z direction be mutable and all the other vars fixed
    void MPC::AddFKConstraints(const vector_t& state, double time, int node) {
        matrix_t G;
        vector_t g;

        const int spline_offset = GetPosSplineStartIdx();

        vector_t state1 = prev_traj_.GetStates().at(node);
        if (node == 0) {
            state1 = state;
        }

        // TODO: Remove this node restriction, then the initial condition will require we are at zero
        // TODO: I wonder if the FK would be better in the cost (would do with auto-diff)?
        // Or maybe this is better using the null space projection strategy?
        // This doesn't work well when I have a capped joint velocity limit
        if (node > 0) {
            for (int ee = 0; ee < num_ee_; ee++) {
                model_.GetFKLinearization(state1, state, prev_traj_.GetInputs(), ee, G, g);
                data_.fk_constraints_.block(node * POS_VARS * num_ee_ + ee * POS_VARS,
                                            node * num_states_ + CentroidalModel::MOMENTUM_OFFSET,
                                            POS_VARS, num_joints_ + CentroidalModel::FLOATING_VEL_OFFSET) = G;

                data_.fk_constants_.segment(node * POS_VARS * num_ee_ + ee * POS_VARS, POS_VARS) = -g;

                for (int coord = 0; coord < POS_VARS; coord++) {
                    if (prev_traj_.IsSplineMutable(ee, coord)) {
                        int vars_index, vars_affecting;
                        std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
                        data_.fk_constraints_.block(node * POS_VARS * num_ee_ + ee * POS_VARS + coord,
                                                    spline_offset + vars_index - vars_affecting,
                                                    1,
                                                    vars_affecting) =
                                -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose();
                    } else {
                        data_.fk_constants_(node * POS_VARS * num_ee_ + ee * POS_VARS + coord) +=
                                prev_traj_.GetPositions().at(ee).at(coord).ValueAt(time);
                    }
                }
            }
        }

    }

    void MPC::AddForceConstraints() {
        // TODO: Might be able to enforce only on Z then let the friction cone enforce it on the other coords
        // -------------------- Zero force when in swing ----------------- //
        int idx = 0; //(info_.num_nodes+1) * POS_VARS * num_ee_;
        int force_offset = GetForceSplineStartIdx();
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                const auto& poly_vars = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVars();
                const auto& poly_times = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyTimes();
                for (int i = 0; i < poly_times.size(); i++) {
                    if ((poly_vars.at(i).size() == 1 && i < poly_times.size()-1 && poly_vars.at(i+1).size() == 1)
                    || (i == 0 && poly_vars.at(i).size() == 1)) {
                        int vars_index, vars_affecting;
                        std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, poly_times.at(i), coord);
                        data_.swing_force_constraints_(idx, force_offset + vars_index - 1) = 1;
                        data_.swing_force_constants_(idx) = 0;
                        idx++;
                        i++;
                    }
                }
            }
        }

        // -------------------- Z force is positive when in contact (repeat I believe) -------------------- //
        // Need to enforce that fz is only a positive force in z for each end effector
//        int idx_in = data_.num_inequality_constraints - 1;
//        const int coord = 2;
//        for (int ee = 0; ee < num_ee_; ee++) {
//            const auto& poly_vars = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVars();
//            const auto& poly_times = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyTimes();
//            for (int i = 0; i < poly_times.size(); i++) {
//                if (poly_vars.at(i).size() != 1) {
//                    int vars_index, vars_affecting;
//                    std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, poly_times.at(i), coord);
//                    data_.inequality_constraints(idx_in, force_offset + vars_index - 2) = 1;   // TODO: will only work with POLY_ORDER=3
//                    data_.inequality_constants_lb(idx_in) = 0;
//                    data_.inequality_constants_ub(idx_in) = 1000; // TODO: make not hard coded
//                    idx_in--;
//                }
//            }
//        }
//        // Check that we get everything we expect
//        assert(data_.num_inequality_constraints - idx_in - 1 == prev_traj_.GetInputs().GetNumForceValsZ());
    }

    void MPC::AddGroundIntersectConstraints() {
        int idx = 0; //foot_on_ground_start_;

        const int coord = 2;
        int pos_offset = GetPosSplineStartIdx();

        for (int ee = 0; ee < num_ee_; ee++) {
            const auto& poly_vars = prev_traj_.GetPositions().at(ee).at(coord).GetPolyVars();
            const auto& poly_times = prev_traj_.GetPositions().at(ee).at(coord).GetPolyTimes();
            for (int i = 1; i < poly_times.size(); i++) { // TODO: Currently ignoring the first time
                // TODO: Change this to a function call to verify if its a constant
                if ((poly_vars.at(i).size() == 1 && i < poly_times.size() - 1 && poly_vars.at(i + 1).size() == 1)
                    || (i == 0 && poly_vars.at(i).size() == 1) ||
                    (poly_vars.at(i).size() == 1 && i == poly_times.size() - 1)) {
                    int vars_index, vars_affecting;
                    std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, poly_times.at(i),
                                                                                             coord);

                    data_.foot_on_ground_constraints_(idx, pos_offset + vars_index - 1) = 1;
                    data_.foot_on_ground_constants_(idx) = 0.0;

                    idx++;
                    i++;
                }
            }
        }
    }

    void MPC::AddFrictionConeConstraints(const vector_t& state, double time, int node) {
        // -------------------- Friction pyramid -------------------- //
        int force_offset = GetForceSplineStartIdx();
        int idx = 4*num_ee_*node;
        for (int ee = 0; ee < num_ee_; ee++) {
            std::vector<std::array<Spline, 3>> forces = prev_traj_.GetInputs().GetForces();
            for (int coord = 0; coord < POS_VARS; coord++) {
                vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);

                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time, coord);

                for (int fric_con = 0; fric_con < 4; fric_con++) {
                    // All the friction constraints are effected by all 3 coordinates of the end effectors
                    data_.friction_cone_constraints_.block(idx + fric_con,
                                                       force_offset + vars_index - vars_affecting,
                                                       1, vars_affecting) =
                                                               friction_pyramid_(fric_con, coord)*vars_lin.transpose();
                }
            }
            data_.friction_cone_ub_.segment(idx, 4) = vector_t::Zero(4);
            data_.friction_cone_lb_.segment(idx, 4) = -qp_solver->GetInfinity(4);
            idx += 4;
        }

        if (node == info_.num_nodes) {
            assert(idx == box_constraint_start_);
        }
    }

    // TODO: Currently assuming the end effector z parameterization only has 2 polynomials
    // Constrains thw swing foot trajectory
//    void MPC::AddSwingFootConstraints() {
//        // Assign the constant height and zero derivative (at the top) to the swing foot trajectory
//        int idx = 0;
//        const int coord = 2;
//        const int pos_offset = GetPosSplineStartIdx();
//        const auto& positions = prev_traj_.GetPositions();
//        for (int ee = 0; ee < num_ee_; ee++) {
//            const auto& poly_vars = positions.at(ee).at(coord).GetPolyVars();
//            const auto& poly_times = positions.at(ee).at(coord).GetPolyTimes();
//            for (int i = 0; i < poly_vars.size(); i++) {
//                if (poly_vars.at(i).size() == 2) {
//                    int vars_index, vars_affecting;
//                    std::tie(vars_index, vars_affecting) =
//                            prev_traj_.GetPositionSplineIndex(ee, poly_times.at(i), coord);
//                    data_.swing_trajectory_constraints_(idx, pos_offset + vars_index - 2) = 1;
//                    data_.swing_trajectory_constants_(idx) = info_.swing_height;
//                    data_.swing_trajectory_constraints_(idx+1, pos_offset + vars_index - 1) = 1;
//                    data_.swing_trajectory_constants_(idx+1) = 0;
//                    idx += 2;
//                }
//            }
//        }
//    }

    void MPC::AddPositivityConstraints(double time, int node) {
        // -------------------- Z force is non-negative at all time -------------------- //
//        int idx = node*num_ee_;
//        int force_offset = GetForceSplineStartIdx();
//        const int coord = 2;
//        for (int ee = 0; ee < num_ee_; ee++) {
//            int vars_index, vars_affecting;
//            std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time, coord);
//            vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);
//
//            data_.positive_force_constraints_.block(idx,
//                                                    force_offset + vars_index - vars_affecting,
//                                                    1, vars_affecting) = vars_lin.transpose();
//            data_.positive_force_lb_(idx) = 0;
//            data_.positive_force_ub_(idx) = qp_solver->GetInfinity(1)(0); // TODO: make not hard coded
//            idx++;
//        }

//        if (node == info_.num_nodes) {
//            assert(idx == ground_positive_start_);
//        }

        // -------------------- Z positions cannot intersect the ground anywhere -------------------- //
//        idx = node*num_ee_; // ground_positive_start_ +
//        int pos_offset = GetPosSplineStartIdx();
//        for (int ee = 0; ee < num_ee_; ee++) {
//            int vars_index, vars_affecting;
//            std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
//            vector_t vars_lin = prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time);
//
//            data_.foot_ground_inter_constraints_.block(idx, pos_offset + vars_index - vars_affecting,
//                                               1, vars_affecting) = vars_lin.transpose();
//            // TODO: Make this based on the decision variables.
//            data_.foot_ground_inter_lb_(idx) = -model_.GetCOMPosition(prev_traj_.GetStates().at(node))(2) - 0.05;
//            data_.foot_ground_inter_ub_(idx) = 1000; // TODO: make not hard coded
//            idx++;
//        }
//        if (node == info_.num_nodes) {
//            assert(idx == foot_on_ground_start_);
//        }
    }

    void MPC::AddBoxConstraints(const vector_t& state, double time, int node) {
        int idx = node*2*num_joints_; // box_constraint_start_ +

        // Velocity bounds
        data_.box_ub_.segment(idx,num_joints_) = info_.vel_bounds;
        data_.box_lb_.segment(idx,num_joints_) = -info_.vel_bounds;
        data_.box_constraints_.block(idx,
                                           GetVelocityIndex(node),
                                           num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);

        // Configuration bounds
        data_.box_ub_.segment(idx + num_joints_,num_joints_) = info_.joint_bounds_ub;
        data_.box_lb_.segment(idx + num_joints_,num_joints_) = info_.joint_bounds_lb;
        data_.box_constraints_.block(idx + num_joints_,
                                           GetJointIndex(node),
                                           num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);
    }

    void MPC::AddForceBoxConstraints(double time, int node) {
        int force_idx = GetForceSplineStartIdx();

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time, coord);
                vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);

                data_.force_box_constraints_.block(POS_VARS*num_ee_*node + ee*POS_VARS + coord,
                                                   force_idx + vars_index - vars_affecting,
                                                   1, vars_affecting) = vars_lin.transpose();
                if (coord == 2) {
                    data_.force_box_lb_(POS_VARS*num_ee_*node + ee*POS_VARS + coord) = 0;
                } else {
                    data_.force_box_lb_(POS_VARS*num_ee_*node + ee*POS_VARS + coord) = -info_.force_bound;
                }
                data_.force_box_ub_(POS_VARS*num_ee_*node + ee*POS_VARS + coord) = info_.force_bound;
            }
        }
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
            data_.cost_quadratic.block(node*num_states_,
                                   node*num_states_,
                                   num_states_, num_states_) = Q_;
        }
        if (Q_forces_.size() > 0) {
            data_.cost_quadratic.block(GetForceSplineStartIdx(), GetForceSplineStartIdx(),
                                       Q_forces_.rows(), Q_forces_.cols()) = Q_forces_;
        }
    }

    void MPC::AddGradientCost() {
        for (int node = 0; node < info_.num_nodes; node++) {
            data_.cost_linear.segment(node * num_states_, num_states_) = w_;
        }
    }

    void MPC::AddFinalCost() {
        data_.cost_quadratic.block(info_.num_nodes*num_states_, info_.num_nodes*num_states_,
                                   num_states_, num_states_) = Phi_;

        data_.cost_linear.segment(info_.num_nodes*num_states_, num_states_) = Phi_w_;
    }

    std::vector<std::vector<double>> MPC::CreateDefaultSwitchingTimes(int num_switches, int num_ee, double horizon) {
        std::vector<std::vector<double>> switching_times;
        for (int i = 0; i < num_ee; i++) {
            std::vector<double> times;
            for (int j = 0; j < num_switches; j++) {
                times.push_back((j+1)*horizon/num_switches);
            }
            switching_times.push_back(times);
        }

        return switching_times;
    }

    int MPC::GetForceSplineStartIdx() const {
        return num_states_*(1+info_.num_nodes);     // initial condition and states for each node
    }

    int MPC::GetVelocityIndex(int node) const {
        return GetForceSplineStartIdx() + prev_traj_.GetInputs().GetTotalForceSplineVars() + node*num_joints_;
    }

    int MPC::GetJointIndex(int node) const {
        return (node+1) * num_states_ + CentroidalModel::FLOATING_VEL_OFFSET + CentroidalModel::MOMENTUM_OFFSET;
    }

    int MPC::GetPosSplineStartIdx() const {
        return GetVelocityIndex(info_.num_nodes);
    }

    Trajectory MPC::ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const {
        // Start by copying the current trajectory to keep the switching times and what not
        Trajectory traj(prev_traj_);

        // Assign all the spline information
        int force_idx = GetForceSplineStartIdx();
        int pos_idx = GetPosSplineStartIdx(); //force_idx + info_.num_nodes*num_joints_ + traj.GetInputs().GetTotalForceSplineVars() - 1;
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_in_spline = traj.GetInputs().GetForces().at(ee).at(coord).GetTotalPolyVars();
//                int num_constant = traj.GetInputs().GetForces().at(ee).at(coord).GetNumConstant();
                traj.UpdateForceSpline(ee, coord, qp_sol.segment(force_idx, vars_in_spline));
                force_idx += vars_in_spline;

                vars_in_spline = traj.GetPositions().at(ee).at(coord).GetTotalPolyVars();
                traj.UpdatePositionSpline(ee, coord, qp_sol.segment(pos_idx, vars_in_spline));
                pos_idx += vars_in_spline;
            }
        }

        // Set states and joint vels
        for (int node = 0; node < info_.num_nodes + 1; node++) {
            // Need to convert the state back to the manifold valued state
            vector_t man_state = CentroidalModel::ConvertAlgebraStateToManifoldState(
                    qp_sol.segment(node*num_states_, num_states_), init_state);
            // need to normalize quaternion returned by MPC
            Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(man_state.segment(9, 4)));
            // Note the warning on the pinocchio function!
            pinocchio::quaternion::firstOrderNormalize(quat);
//            std::cout << "quat: \n" << quat << "\nnorm: " << quat.norm() << ", is normalized: " <<
//            pinocchio::quaternion::isNormalized(quat) << std::endl;
            man_state(9) = quat.x();
            man_state(10) = quat.y();
            man_state(11) = quat.z();
            man_state(12) = quat.w();

            traj.SetState(node, man_state);
        }

        for (int node = 0; node < info_.num_nodes; node++) {
            traj.SetInputVels(node, qp_sol.segment(GetVelocityIndex(node), num_joints_));
        }

//        traj.PrintTrajectoryToFile("internal_traj.txt");

        return traj;
    }

    void MPC::UpdateQPSizes() {
        data_.num_decision_vars = (info_.num_nodes+1)*num_states_ + info_.num_nodes*num_joints_ +
                                  prev_traj_.GetInputs().GetTotalForceSplineVars() + prev_traj_.GetTotalPosSplineVars();
        data_.num_dynamics_constraints = (info_.num_nodes+1)*num_states_;

        data_.num_cone_constraints_ = (info_.num_nodes+1)*4*num_ee_;
        data_.num_box_constraints_ = info_.num_nodes*2*num_joints_;
//        data_.num_positive_force_constraints_ = (info_.num_nodes+1)*num_ee_;
//        data_.num_swing_foot_pos_constraints_ = 2*prev_traj_.GetTotalPosNonConstantZ();
        data_.num_foot_ground_inter_constraints_ = (info_.num_nodes+1)*num_ee_;
        data_.num_foot_on_ground_constraints_ = prev_traj_.GetTotalPosConstantsZ();
        data_.num_force_box_constraints_ = POS_VARS*num_ee_*(info_.num_nodes+1);

        data_.num_fk_constraints_ = POS_VARS*num_ee_*(info_.num_nodes+1);
        data_.num_swing_foot_constraints_ = prev_traj_.GetInputs().GetTotalForceConstants();

        cone_constraint_start_ = 0;
        box_constraint_start_ = (info_.num_nodes+1)*4*num_ee_;
        positive_force_start_ = box_constraint_start_ + info_.num_nodes*2*num_joints_;  // TODO: need num_nodes+1?
        ground_positive_start_ = positive_force_start_ + (info_.num_nodes+1)*num_ee_;
        foot_on_ground_start_ = ground_positive_start_ + (info_.num_nodes+1)*num_ee_;

        data_.num_inequality_constraints = foot_on_ground_start_ + prev_traj_.GetTotalPosConstantsZ();
//                (info_.num_nodes+1)*(4*num_ee_) + info_.num_nodes*(2*num_joints_) + 4*num_ee_ +
//                                           prev_traj_.GetInputs().GetNumForceValsZ() + 2*(info_.num_nodes+1)*num_ee_ +
//                                           prev_traj_.GetTotalPosConstantsZ();

    }

    void MPC::SetDefaultGaitTrajectory(Gaits gait, int num_polys, const std::array<std::array<double, 3>, 4>& ee_pos) {
        if (num_ee_ != 4) {
            throw std::runtime_error("Default gaits have only been implemented for 4 end effectors!");
        }

//        assert(prev_traj_.GetTotalTime() >= 0.75);  // So we have time to make the gait

        std::vector<std::vector<double>> switching_times;

        switch (gait) {
            case Trot: {
                std::vector<double> times;
                // TODO: Change
//                times.push_back(prev_traj_.GetTotalTime()/2);
//                times.push_back(prev_traj_.GetTotalTime());
                times.push_back(0.35);
                times.push_back(0.75);


                Spline position1(num_polys, times, true, Spline::Normal);
                Spline position2(num_polys, times, false, Spline::Normal);

                Spline force1(num_polys, times, false, Spline::Force);
                Spline force2(num_polys, times, true, Spline::Force);

                for (int ee = 0; ee < num_ee_; ee++) {
                    if (ee == 0 || ee == 3) {
                        prev_traj_.SetEndEffectorSplines(ee, force1, position1);
                    } else {
                        prev_traj_.SetEndEffectorSplines(ee, force2, position2);
                    }
                    for (int coord = 0; coord < POS_VARS; coord++) {
                        prev_traj_.SetPositionsForAllTime(ee, ee_pos.at(ee));
                    }
                }

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

    vector_t MPC::GetFullTargetState(double time) const {
        int node = floor((time - init_time_)/info_.integrator_dt);
        return prev_traj_.GetStates().at(node);
    }

    vector_t MPC::GetTargetVelocity(double time) const {
        int node = floor((time - init_time_)/info_.integrator_dt);
        return model_.ComputeBaseVelocities(prev_traj_.GetStates().at(node),
                                                   prev_traj_.GetInputs().GetVels(time));
    }

    // TODO: Check/think about this
    vector_t MPC::GetTargetAcc(double time) const {
        int node = floor((time - init_time_)/info_.integrator_dt);
        return (model_.ComputeBaseVelocities(prev_traj_.GetStates().at(node),
                                            prev_traj_.GetInputs().GetVels(time)) -
                model_.ComputeBaseVelocities(prev_traj_.GetStates().at(node+1),
                                             prev_traj_.GetInputs().GetVels(time + info_.integrator_dt)))/info_.integrator_dt;
    }

    const CentroidalModel& MPC::GetModel() const {
        return model_;
    }

    // TODO: Consider removing the init state somehow
    double MPC::LineSearch(const vector_t& direction, const vector_t& init_state) {
        double alpha = 1;

       // Note: for now just using the equality constraints on the merit function
        double merit = GetMeritValue(prev_qp_sol, mu_, init_state);
        double merit_step = GetMeritValue(alpha*direction + prev_qp_sol, mu_, init_state);
        double merit_directional = GetMeritGradient(prev_qp_sol, direction, mu_, init_state);

//        std::cout << "merit directional derivative averaged per size: " << merit_directional / data_.num_decision_vars << std::endl;
//        std::cout << "merit value: " << merit << std::endl;
//        std::cout << "starting merit with the step: " << merit_step << std::endl;

//        if (merit_directional < -0.75) {

//        if (run_num_ >= 2) {
            int i = 0;
            while ((merit - merit_step) < -0.0001 * alpha * merit_directional && i < 5) {
                alpha *= 0.5;
                merit_step = GetMeritValue(alpha * direction + prev_qp_sol, mu_, init_state);
                i++;
            }

//        }
            if (i == 5) {
                alpha = 1;
            }
//        } else {
//            alpha = 0;
//        }

        // TODO: I would expect alpha = 1 to be the most common value

//        std::cout << "final merit decrease: " << merit - merit_step << std::endl;
//        std::cout << "final alpha: " << alpha << std::endl;
//        std::cout << "new merit value: " << merit_step << std::endl;

//        std::cout << std::endl;

        return alpha;
    }

    double MPC::GetMeritValue(const vector_t& x, double mu, const vector_t& init_state) const {
        // TODO: Do more efficiently
        const Trajectory temp_traj = ConvertQPSolToTrajectory(x, init_state);
        temp_traj.PrintTrajectoryToFile("temp_traj.txt");

        return GetCostValue(x) + mu*GetEqualityConstraintValues(temp_traj, init_state).lpNorm<1>();
    }

    double MPC::GetCostValue(const vector_t& x) const {
        return 0.5*x.dot(data_.cost_quadratic*(x)) + data_.cost_linear.dot(x);
    }

    // TODO: Something is very weird here. My QP solutions don't seem to satisfy this very well.
    // Given I am solving to completion, I would expect this to to be smaller on the QP solutions
    vector_t MPC::GetEqualityConstraintValues(const Trajectory& traj, const vector_t& init_state) const {
        vector_t eq_constraints = vector_t::Zero(data_.num_fk_constraints_ + data_.num_dynamics_constraints - num_states_);

        for (int node = 0; node < info_.num_nodes; node++) {
            eq_constraints.segment(node*num_states_, num_states_) =
                    CentroidalModel::ConvertManifoldStateToAlgebraState(traj.GetState(node+1), init_state)
                    - model_.GetDiscreteDynamics(CentroidalModel::ConvertManifoldStateToAlgebraState(traj.GetState(node), init_state),
                            traj.GetInputs(), GetTime(node), init_state);
//            eq_constraints.segment(node*num_states_ + 3, 3) = Eigen::Vector3d::Zero();

//            std::cout << "dynamics: \n" << model_.GetDiscreteDynamics(CentroidalModel::ConvertManifoldStateToAlgebraState(traj.GetState(node), init_state),
//                                                                      traj.GetInputs(), GetTime(node)) << std::endl;
//            std::cout << "trajectory: \n" << CentroidalModel::ConvertManifoldStateToAlgebraState(traj.GetState(node+1), init_state) << std::endl;

//            for (int ee = 0; ee < num_ee_; ee++) {
//                eq_constraints.segment(info_.num_nodes*num_states_ + node*3*num_ee_ + 3*ee, 3) =
//                        model_.GetEndEffectorLocationCOMFrame(traj.GetState(node), model_.GetEndEffectorFrame(ee))
//                        - traj.GetPosition(ee,GetTime(node));
//            }
        }

//        std::cout << "equality constraint violation: \n";
//        std::cout << "dynamics: " << eq_constraints.head(data_.num_dynamics_constraints - num_states_).lpNorm<1>() << std::endl;
//        std::cout << "FK: " << eq_constraints.tail(data_.num_fk_constraints_).lpNorm<1>() << std::endl;
//        std::cout << std::endl;

        return eq_constraints;
    }

    double MPC::GetTime(int node) const {
        return node * info_.integrator_dt + init_time_;
    }

    double MPC::GetMeritGradient(const vector_t& x, const vector_t& p, double mu, const vector_t& init_state) const {
        // TODO: Do more efficiently
        const Trajectory temp_traj = ConvertQPSolToTrajectory(x, init_state);

        return (data_.cost_quadratic*x + data_.cost_linear).dot(p)
        - mu * GetEqualityConstraintValues(temp_traj, init_state).lpNorm<1>();
    }

    void MPC::AddForceCost(double weight) {
        int num_forces = prev_traj_.GetInputs().GetTotalForceSplineVars();
        Q_forces_ = weight*matrix_t::Identity(num_forces, num_forces);
    }

    void MPC::RecordStats(double alpha, const vector_t& direction, const std::string& solve_type,
                          const vector_t& ref_state) {
        equality_constraint_violations_.push_back(GetEqualityConstraintValues(prev_traj_, ref_state).lpNorm<1>());
        step_norm_.push_back(direction.norm());
        alpha_.push_back(alpha);
        cost_result_.push_back(GetCostValue(prev_qp_sol));
        merit_result_.push_back(GetMeritValue(prev_qp_sol, mu_, ref_state));
        merit_directional_deriv_.push_back(GetMeritGradient(prev_qp_sol - alpha*direction, direction, mu_, ref_state));
        solve_type_.push_back(solve_type);
        ref_state_.push_back(ref_state);
    }

    void MPC::PrintStats() {
        using std::setw;
        using std::setfill;

        const int col_width = 15;
        const int table_width = 8*col_width;

        std::cout << setfill('-') << setw(table_width) << "" << std::endl;
        std::cout << std::left << setfill(' ') << setw(table_width/2 - 7) << "" << "MPC Statistics" << std::endl;
        std::cout << setfill('-') << setw(table_width) << "" << std::endl;

        std::cout << setfill(' ');
        std::cout << setw(col_width) << "Solve #"
                << setw(col_width) << "Constraints"
                << setw(col_width) << "Step Norm"
                << setw(col_width) << "Alpha"
                << setw(col_width) << "Cost"
                << setw(col_width) << "Merit"
                << setw(col_width) << "Merit dd"
                << setw(col_width) << "Solve Type" << std::endl;
        std::cout << std::setfill('-') << setw(table_width) << "" << std::endl;
        std::cout << setfill(' ');
        for (int i = 0; i < alpha_.size(); i++) {
            std::cout << setw(col_width) << i
                      << setw(col_width) << equality_constraint_violations_.at(i)
                      << setw(col_width) << step_norm_.at(i)
                      << setw(col_width) << alpha_.at(i)
                      << setw(col_width) << cost_result_.at(i)
                      << setw(col_width) << merit_result_.at(i)
                      << setw(col_width) << merit_directional_deriv_.at(i)
                      << setw(col_width) << solve_type_.at(i) << std::endl;
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

} // mpc