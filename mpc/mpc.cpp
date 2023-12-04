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
        time_horizon = info.time_horizon;
        num_qp_iterations = info.num_qp_iterations;
        num_contacts = info.num_contacts;
        friction_coef = info.friction_coef;
        vel_bounds = info.vel_bounds;
        joint_bounds = info.joint_bounds;
        ee_frames = info.ee_frames;
        discretization_steps = info.discretization_steps;
        num_switches = info.num_switches;
        integrator_dt = info.integrator_dt;
    }

    MPC::MPC(const MPCInfo& info, const std::string& robot_urdf) :
        info_(info),
        model_(robot_urdf, info.ee_frames, info.discretization_steps, info.integrator_dt),
        num_joints_(model_.GetNumJoints()),
        num_states_(model_.GetPinocchioNumConfig() + 5),
        num_ee_(model_.GetNumEndEffectors()),
        prev_traj_(info.num_nodes+1, num_states_+1, num_joints_,
               CreateDefaultSwitchingTimes(info.num_switches, num_ee_, info.time_horizon),
               info.time_horizon/info.num_nodes) {

        assert(info_.ee_frames.size() == num_ee_);

        if (info_.vel_bounds.size() != num_joints_ || info_.joint_bounds.size() != num_joints_) {
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

        qp_solver = std::make_unique<OSQPInterface>(data_, true);

        SetFrictionPyramid();

        Phi_ = matrix_t::Zero(num_states_, num_states_);
        Phi_w_ = vector_t::Zero(num_states_);
    }

    Trajectory MPC::Solve(const vector_t &centroidal_state, double init_time) {
        assert(centroidal_state.size() == num_states_ + 1);


        prev_traj_.SetInitTime(init_time);

        prev_traj_.AddPolys(info_.integrator_dt*info_.num_nodes + init_time);

        prev_traj_.RemoveUnusedPolys(init_time);

        UpdateQPSizes();

        // Reset to 0
        ResetQPMats();

        // initial condition
        data_.dynamics_constraints.topLeftCorner(num_states_, num_states_) = -matrix_t::Identity(num_states_, num_states_);
        data_.dynamics_constants.head(num_states_) =
                -CentroidalModel::ConvertManifoldStateToAlgebraState(centroidal_state, centroidal_state);


        // TODO: Current problems (12/1/23):
        // Even with NO FK constraints ee#1 x position is non-zero -- Fixed.
        // FK constraints seem to cause weird changes in joint angles even when we start on the desired position
        // Something is weird in the dynamics - z force is near zero (slightly negative), but position is constant.
        //      it should be falling with near zero force.
        // Trajectory printing switching times is only for positions or forces, not both -- Fixed
        // Positivity constraints on z forces do not seem to be working -- Fixed
        // Position variables not in use are non-zero -- Fixed
        // Last element of the states in the generated trajectory looks sus
        // Ground intersection constraint causing weird things with the forces, maybe fixed, issues with positivity + ground interaction
        // Friction cone constraints causing issues -- Fixed
        // Something is weird in the force dynamics, might just be integrator chains
        // As of 12/1 at 5:50pm box, dynamics, cone, positive, force constraints seem to work together (for only 20-30 nodes, more than that is weird).
        // As of 12/1 at 5:57pm dynamics, box, FK, ground intersect seem to work together for at least 50 nodes
        // I suspect there is a dynamics issue causing the weird forces and there is potentially some weird
        // indexing issue causing some of the weirdness when multiple constraints are enabled.

        // TODO: Need to re-do the indexing on the constraints to its less of a mess (row indexing)

        matrix_t G;
        vector_t g;
        // Constraints on all the nodes BUT the first
        for (int i = 0; i < info_.num_nodes; i++) {
            double time = i * info_.integrator_dt + init_time;
            // ------------------------ Dynamics ------------------------ //
            AddDynamicsConstraints(centroidal_state, time, i);

            // ------------------------ Inequality Constraints ------------------------ //
            AddBoxConstraints(centroidal_state, time, i);

            // ----------------------- Costs ------------------------- //
            AddHessianApproxCost(centroidal_state, time, i);
            AddGradientCost(centroidal_state, time, i);
        }

        AddFinalCost();

        // Constraints on ALL the nodes
        for (int i = 0; i < info_.num_nodes+1; i++) {
            double time = i * info_.integrator_dt + init_time;
            // ------------------------ FK Constraint ------------------------ //
            AddFKConstraints(centroidal_state, time, i);

            // ------------------------ Friction Cone Constraints ------------------------ //
//            AddFrictionConeConstraints(centroidal_state, time, i);

            // Positivity of z direction force
//            AddPositivityConstraints(time, i);
        }

        // TODO: instead of adding a constraint consider removing them as decision variables
//        AddForceConstraints();
        AddGroundIntersectConstraints();


//        PrintEqualityConstraints();
//        PrintInequalityConstraints();
//        PrintDynamicsConstraints();

        qp_solver->SetupQP(data_);
        vector_t sol = qp_solver->Solve(data_);

//        std::cout << "qp_sol: \n" << sol << std::endl;

        prev_traj_ = ConvertQPSolToTrajectory(sol, centroidal_state);

        prev_traj_.PrintTrajectoryToFile("prev_traj_after_solve.txt");

        return prev_traj_;
    }

    void MPC::SetWarmStartTrajectory(const mpc::Trajectory &trajectory) {
        prev_traj_ = trajectory;
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
        data_.dynamics_constraints = matrix_t::Zero(data_.num_dynamics_constraints, data_.num_decision_vars);
        data_.dynamics_constants = vector_t::Zero(data_.num_dynamics_constraints);

        data_.equality_constraints = matrix_t::Zero(data_.num_equality_constraints, data_.num_decision_vars);
        data_.equality_constants = vector_t::Zero(data_.num_equality_constraints);

        data_.inequality_constraints = matrix_t::Zero(data_.num_inequality_constraints, data_.num_decision_vars);
        data_.inequality_constants_ub = vector_t::Zero(data_.num_inequality_constraints);
        data_.inequality_constants_lb = vector_t::Zero(data_.num_inequality_constraints);

        data_.cost_quadratic = matrix_t::Zero(data_.num_decision_vars, data_.num_decision_vars);
        data_.cost_linear = vector_t::Zero(data_.num_decision_vars);
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

// TODO: update indexing
        data_.dynamics_constraints.block((node+1)*num_states_,
                                         node*num_joints_ + GetForceSplineStartIdx() + force_spline_vars - 1,
                                         num_states_, num_joints_) = B_vels;
        data_.dynamics_constraints.block((node+1)*num_states_,
                                         GetForceSplineStartIdx(),
                                         num_states_, force_spline_vars) = B_spline;
        data_.dynamics_constants.segment((node+1)*num_states_, num_states_) = -C;

    }


    void MPC::AddFKConstraints(const vector_t& state, double time, int node) {
        matrix_t G;
        vector_t g;

        int spline_offset = GetPosSplineStartIdx();

        vector_t state1 = prev_traj_.GetStates().at(node);
        if (node == 0) {
            state1 = state;
        }

        for (int ee = 0; ee < num_ee_; ee++) {
            model_.GetFKLinearization(state1, state, prev_traj_.GetInputs(), ee, G, g);
            data_.equality_constraints.block(node * POS_VARS * num_ee_ + ee * POS_VARS,
                                             node * num_states_ + CentroidalModel::MOMENTUM_OFFSET,
                                             POS_VARS, num_joints_ + CentroidalModel::FLOATING_VEL_OFFSET) = G;

            data_.equality_constants.segment(node * POS_VARS * num_ee_ + ee * POS_VARS, POS_VARS) = -g;

//            std::cout << "G for ee #" << ee << ": \n" << G << std::endl;
//            std::cout << "g for ee #" << ee << ": \n" << g << std::endl;

            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
                data_.equality_constraints.block(node * POS_VARS * num_ee_ + ee * POS_VARS + coord,
                                                 spline_offset + vars_index - vars_affecting,
                                                 1,
                                                 vars_affecting) =
                                                         -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose();
            }
        }

    }

    void MPC::AddForceConstraints() {
        // -------------------- Zero force when in swing ----------------- //
        int idx = (info_.num_nodes+1) * POS_VARS * num_ee_;
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
                        data_.equality_constraints(idx, force_offset + vars_index - 1) = 1;
                        data_.equality_constants(idx) = 0;
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

    // TODO: Consider how this works as the COM moves during traj opt
    // TODO: might want to consider using decision variables to get the COM at that point so this isn't pre-specified
    void MPC::AddGroundIntersectConstraints() {
        int idx = foot_on_ground_start_;

        const int coord = 2;
        int pos_offset = GetPosSplineStartIdx();
        for (int ee = 0; ee < num_ee_; ee++) {
            const auto& poly_vars = prev_traj_.GetPositions().at(ee).at(coord).GetPolyVars();
            const auto& poly_times = prev_traj_.GetPositions().at(ee).at(coord).GetPolyTimes();
            for (int i = 0; i < poly_times.size(); i++) {
                int node = floor(poly_times.at(i)/info_.integrator_dt); // Gets the node for the associated time
                if (node <= info_.num_nodes) {
                    if ((poly_vars.at(i).size() == 1 && i < poly_times.size() - 1 && poly_vars.at(i + 1).size() == 1)
                        || (i == 0 && poly_vars.at(i).size() == 1) ||
                        (poly_vars.at(i).size() == 1 && i == poly_times.size() - 1)) {
                        int vars_index, vars_affecting;
                        std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, poly_times.at(i),
                                                                                                 coord);
                        for (int j = 0; j < data_.inequality_constraints.row(idx).size(); j++) {
                            if (data_.inequality_constraints(idx, j) != 0) {
                                std::cout << "j: " << j << std::endl;
//                                throw std::runtime_error("indexing error!");
                            }
                        }

                        data_.inequality_constraints(idx, pos_offset + vars_index - 1) = 1;

                        data_.inequality_constants_ub(idx) =
                                -model_.GetCOMPosition(prev_traj_.GetStates().at(node))(2) + 0.05;
                        data_.inequality_constants_lb(idx) =
                                -model_.GetCOMPosition(prev_traj_.GetStates().at(node))(2) - 0.05;

                        idx++;
                        i++;
                    }
//                    std::cout << "COM pos at Node #" << node << ": \n" << model_.GetCOMPosition(prev_traj_.GetStates().at(node)) << std::endl;
                }
            }
        }
        assert(idx <= data_.num_inequality_constraints);
    }

    void MPC::AddFrictionConeConstraints(const vector_t& state, double time, int node) {
        // -------------------- Friction pyramid -------------------- //
        int force_offset = GetForceSplineStartIdx();
        int idx = cone_constraint_start_ + 4*num_ee_*node;
        for (int ee = 0; ee < num_ee_; ee++) {
            std::vector<std::array<Spline, 3>> forces = prev_traj_.GetInputs().GetForces();
            for (int coord = 0; coord < POS_VARS; coord++) {
                vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);

                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time, coord);

                for (int fric_con = 0; fric_con < 4; fric_con++) {
                    // All the friction constraints are effected by all 3 coordinates of the end effectors
                    data_.inequality_constraints.block(idx + fric_con,
                                                       force_offset + vars_index - vars_affecting,
                                                       1, vars_affecting) =
                                                               friction_pyramid_(fric_con, coord)*vars_lin.transpose();
                }
            }
            data_.inequality_constants_ub.segment(idx, 4) = vector_t::Zero(4);
            data_.inequality_constants_lb.segment(idx, 4) = -qp_solver->GetInfinity(4);
            idx += 4;
        }

        if (node == info_.num_nodes) {
            assert(idx == box_constraint_start_);
        }
    }

    void MPC::AddPositivityConstraints(double time, int node) {
        // -------------------- Z force is non-negative at all time -------------------- //
        int idx = positive_force_start_ + node*num_ee_;
        int force_offset = GetForceSplineStartIdx();
        const int coord = 2;
        for (int ee = 0; ee < num_ee_; ee++) {
            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time, coord);
            vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);

            data_.inequality_constraints.block(idx, force_offset + vars_index - vars_affecting, 1, vars_affecting) = vars_lin.transpose();
            data_.inequality_constants_lb(idx) = 0;
            data_.inequality_constants_ub(idx) = qp_solver->GetInfinity(1)(0); // TODO: make not hard coded
            idx++;
        }

        if (node == info_.num_nodes) {
            assert(idx == ground_positive_start_);
        }

        // -------------------- Z positions cannot intersect the ground anywhere -------------------- //
        idx = ground_positive_start_ + node*num_ee_;
        int pos_offset = GetPosSplineStartIdx();
        for (int ee = 0; ee < num_ee_; ee++) {
            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
            vector_t vars_lin = prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time);

            data_.inequality_constraints.block(idx, pos_offset + vars_index - vars_affecting,
                                               1, vars_affecting) = vars_lin.transpose();
            // TODO: Make this based on the decision variables.
            data_.inequality_constants_lb(idx) = -model_.GetCOMPosition(prev_traj_.GetStates().at(node))(2) - 0.05;
            data_.inequality_constants_ub(idx) = 1000; // TODO: make not hard coded
            idx++;
        }
        if (node == info_.num_nodes) {
            assert(idx == foot_on_ground_start_);
        }
    }

    void MPC::AddBoxConstraints(const vector_t& state, double time, int node) {
        int idx = box_constraint_start_ + node*2*num_joints_;

        // Velocity bounds
        data_.inequality_constants_ub.segment(idx,
                                              num_joints_) = info_.vel_bounds;
        data_.inequality_constants_lb.segment(idx,
                                              num_joints_) = -info_.vel_bounds;
        data_.inequality_constraints.block(idx,
                                           GetVelocityIndex(node),
                                           num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);

        // Configuration bounds
        data_.inequality_constants_ub.segment(idx,
                                              num_joints_) = info_.joint_bounds;
        data_.inequality_constants_lb.segment(idx,
                                              num_joints_) = -info_.joint_bounds;
        data_.inequality_constraints.block(idx,
                                           GetJointIndex(node),
                                           num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);

        idx += 2*num_joints_;
        if (node == info_.num_nodes - 1) {
            assert(idx == positive_force_start_);
        }
    }

    void MPC::AddQuadraticTrackingCost(const vector_t& state_des, const matrix_t& Q) {
        if (Q.rows() != num_states_ || Q.cols() != num_states_) {
            throw std::runtime_error("Supplied quadratic cost term is the wrong size.");
        }

        Q_ = Q;
        w_ = -1*Q*state_des;
    }

    void MPC::AddHessianApproxCost(const mpc::vector_t &state, double time, int node) {
        data_.cost_quadratic.block(node*num_states_,
                                   node*num_states_,
                                   num_states_, num_states_) = Q_;
    }

    void MPC::AddGradientCost(const mpc::vector_t &state, double time, int node) {
        data_.cost_linear.segment(node*num_states_, num_states_) = w_;
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

        traj.PrintTrajectoryToFile("internal_traj.txt");

        return traj;
    }

    void MPC::PrintDynamicsConstraints() const {
        std::cout << "dynamics constraints: \n" << data_.dynamics_constraints << std::endl;
        std::cout << "dynamics constants: \n" << data_.dynamics_constants << std::endl;
    }

    void MPC::PrintEqualityConstraints() const {
        std::cout << "equality constraints: \n" << data_.equality_constraints << std::endl;
        std::cout << "equality constants: \n" << data_.equality_constants << std::endl;
    }

    void MPC::PrintInequalityConstraints() const {
        std::cout << "inequality constraints: \n" << data_.inequality_constraints << std::endl;
        std::cout << "inequality ub: \n" << data_.inequality_constants_ub << std::endl;
        std::cout << "inequality lb: \n" << data_.inequality_constants_lb << std::endl;
    }

    void MPC::UpdateQPSizes() {
        data_.num_decision_vars = (info_.num_nodes+1)*num_states_ + info_.num_nodes*num_joints_ +
                                  prev_traj_.GetInputs().GetTotalForceSplineVars() + prev_traj_.GetTotalPosSplineVars();
        data_.num_dynamics_constraints = (info_.num_nodes+1)*num_states_;
        data_.num_equality_constraints = (info_.num_nodes+1)*POS_VARS*num_ee_ + prev_traj_.GetInputs().GetTotalForceConstants();

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

        assert(prev_traj_.GetTotalTime() >= 0.75);  // So we have time to make the gait

        std::vector<std::vector<double>> switching_times;

        switch (gait) {
            case Trot: {
                std::vector<double> times;
                times.push_back(0.35);
                times.push_back(0.75);

                Spline type1 = Spline(num_polys, times, true);
                Spline type2 = Spline(num_polys, times, false);

                for (int ee = 0; ee < num_ee_; ee++) {
                    if (ee == 0 || ee == 3) {
                        prev_traj_.SetEndEffectorSplines(ee, type1, type2);
                    } else {
                        prev_traj_.SetEndEffectorSplines(ee, type2, type1);
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

    vector_t MPC::GetTargetConfig() const {
        return prev_traj_.GetStates().at(1).tail(num_states_ - 5);
    }

    vector_t MPC::GetTargetVelocity() const {
        return model_.ComputeBaseVelocities(prev_traj_.GetStates().at(1),
                                                   prev_traj_.GetInputs().GetVel(0));
    }

    const CentroidalModel& MPC::GetModel() const {
        return model_;
    }


} // mpc