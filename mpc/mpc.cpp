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
        prev_traj_(info.num_nodes, num_states_, num_joints_,
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

        force_spline_vars_ = num_ee_ * 3 * prev_traj_.GetInputs().GetForces().at(0).at(0).GetTotalPolyVars();
        pos_spline_vars_ = num_ee_*3*prev_traj_.GetPositions().at(0).at(0).GetTotalPolyVars();
        data_.num_decision_vars = info_.num_nodes*(num_states_ + num_joints_) + force_spline_vars_ + pos_spline_vars_;
        data_.num_dynamics_constraints = (info_.num_nodes+1)*num_states_;
        data_.num_equality_constraints = info_.num_nodes*POS_VARS*num_ee_;
        data_.num_inequality_constraints = info_.num_nodes*(4*num_ee_ + 2*num_joints_);

        qp_solver = std::make_unique<OSQPInterface>();

        SetFrictionPyramid();
    }

    // TODO: See below
    // TODO: Switching times for all the spline will NOT be the same!
    // TODO: see above!
    // TODO: see above!
    Trajectory MPC::Solve(const vector_t &centroidal_state) {
        assert(centroidal_state.size() == num_states_ + 1);

        // Reset to 0
        ResetQPMats();

        // initial condition
        data_.dynamics_constraints.topLeftCorner(num_states_, num_states_) = -matrix_t::Identity(num_states_, num_states_);
        data_.dynamics_constants.head(num_states_) =
                -CentroidalModel::ConvertManifoldStateToAlgebraState(centroidal_state, centroidal_state);


        // TODO: Enforce no foot slide and no force at a distance.
        // Two options:
        //  1. Enforce equality constraints for certain variables. Probably easiest.
        //  2. Remove certain variables from the decision vector. Harder, but def the correct way in general.

        matrix_t G;
        vector_t g;
        for (int i = 0; i < info_.num_nodes; i++) {
            double time = i * info_.integrator_dt;
            // ------------------------ Dynamics ------------------------ //
            AddDynamicsConstraints(centroidal_state, time, i);

            // ------------------------ FK Constraint ------------------------ //
            AddFKConstraints(centroidal_state, time, i);

            // ------------------------ Inequality Constraints ------------------------ //
            AddInequalityConstraints(centroidal_state, time, i);

            // ----------------------- Costs ------------------------- //
            AddHessianApproxCost(centroidal_state, time, i);
            AddGradientCost(centroidal_state, time, i);
        }

        // TODO: look into doing this better
        EnforceFootSlipAndForceAtADistance();

        qp_solver->SetupQP(data_);
        vector_t sol = qp_solver->Solve();

        // TODO: Did I just compute a dx or a new x?
        // I believe I calculated x
        // So if I want to take a smaller step, I can recover dx. Should probably check this.

        prev_traj_ = ConvertQPSolToTrajectory(sol, centroidal_state);

//        std::cout << "Solution: \n" << sol << std::endl;

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
            throw std::runtime_error("Wrong sized final cost.");
        }

        Phi_ = Phi;
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

        matrix_t B_spline = B.leftCols(force_spline_vars_);
        matrix_t B_vels = B.rightCols(num_joints_);

        data_.dynamics_constraints.block((node+1)*num_states_,
                                         node*num_joints_ + num_states_*info_.num_nodes + force_spline_vars_,
                                         num_states_, num_joints_) = B_vels;
        data_.dynamics_constraints.block((node+1)*num_states_,
                                         num_states_*info_.num_nodes,
                                         num_states_, force_spline_vars_) = B_spline;
        data_.dynamics_constants.segment((node+1)*num_states_, num_states_) = -C;
    }

    void MPC::AddFKConstraints(const vector_t& state, double time, int node) {
        matrix_t G;
        vector_t g;
        for (int ee = 0; ee < num_ee_; ee++) {
            model_.GetFKLinearization(state, prev_traj_.GetInputs(), ee, G, g);
            data_.equality_constraints.block(node * POS_VARS * num_ee_ + ee * POS_VARS,
                                             node * num_states_ + CentroidalModel::MOMENTUM_OFFSET,
                                             POS_VARS, num_joints_ + CentroidalModel::FLOATING_VEL_OFFSET) = G;

            data_.equality_constants.segment(node * POS_VARS * num_ee_ + ee * POS_VARS, POS_VARS) = -g;

            for (int coord = 0; coord < POS_VARS; coord++) {
                data_.equality_constraints.block(node * POS_VARS * num_ee_ + ee * POS_VARS + coord,
                                                 GetPositionSplineIndex(ee, time, coord),
                                                 1,
                                                 Spline::POLY_ORDER) =
                                                         -prev_traj_.GetPositionsPolyVarsLin(ee, time).row(coord);
            }
        }
    }

    void MPC::AddInequalityConstraints(const vector_t& state, double time, int node) {

        // Friction pyramid
        for (int ee = 0; ee < num_ee_; ee++) {
            std::vector<std::array<Spline, 3>> forces = prev_traj_.GetInputs().GetForces();
            for (int coord = 0; coord < POS_VARS; coord++) {
                Eigen::Vector4d temp = friction_pyramid_.col(coord).transpose().cwiseProduct(
                        prev_traj_.GetInputs().GetForcePolyVarsLin(ee, time).row(coord));

                data_.inequality_constraints.block(node * (4 * num_ee_ + 2 * num_joints_),
                                                   GetForceSplineIndex(ee, time, coord),
                                                   4, 1) = temp;
            }
        }
        data_.inequality_constants_ub.segment(node*(4*num_ee_ + 2*num_joints_), 4*num_ee_) = vector_t::Zero(4*num_ee_);
        data_.inequality_constants_lb.segment(node*(4*num_ee_ + 2*num_joints_), 4*num_ee_) = -qp_solver->GetInfinity(4*num_ee_);


        // Velocity bounds
        data_.inequality_constants_ub.segment(node*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
                                              num_joints_) = info_.vel_bounds;
        data_.inequality_constants_lb.segment(node*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
                                              num_joints_) = -info_.vel_bounds;
        data_.inequality_constraints.block(node*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
                                           GetVelocityIndex(node),
                                           num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);

        // Configuration bounds
        data_.inequality_constants_ub.segment(node*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
                                             num_joints_) = info_.joint_bounds;
        data_.inequality_constants_lb.segment(node*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
                                              num_joints_) = -info_.joint_bounds;
        data_.inequality_constraints.block(node*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
                                           GetJointIndex(node),
                                           num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);

    }

    void MPC::AddHessianApproxCost(const mpc::vector_t &state, double time, int node) {
        data_.cost_quadratic.block(node*num_states_,
                                   node*num_states_,
                                   num_states_, num_states_) = Q_;
    }

    void MPC::AddGradientCost(const mpc::vector_t &state, double time, int node) {
        data_.cost_linear.segment(node*num_states_, num_states_) = w_;
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

    int MPC::GetForceSplineIndex(int end_effector, double time, int coord) const {
        int num_spline_vars_before = end_effector*POS_VARS*
                prev_traj_.GetInputs().GetForces().at(end_effector).at(coord).GetTotalPolyVars();
        int idx_into_spline_vars = coord*prev_traj_.GetInputs().GetForces().at(end_effector).at(coord).GetTotalPolyVars();
        int idx_into_poly = prev_traj_.GetInputs().GetForces().at(end_effector).at(coord).GetPolyIdx(time);

        return info_.num_nodes * num_states_ + num_spline_vars_before + idx_into_spline_vars +
                    idx_into_poly*Spline::POLY_ORDER;
    }

    int MPC::GetForceSplineIndexNoTime(int end_effector, int idx, int coord) const {
        int num_spline_vars_before = end_effector*POS_VARS*
                                     prev_traj_.GetInputs().GetForces().at(end_effector).at(coord).GetTotalPolyVars();
        int idx_into_spline_vars = coord*prev_traj_.GetInputs().GetForces().at(end_effector).at(coord).GetTotalPolyVars();

        return info_.num_nodes * num_states_ + num_spline_vars_before + idx_into_spline_vars +
               idx*Spline::POLY_ORDER;
    }

    int MPC::GetPositionSplineIndex(int end_effector, double time, int coord) const {
        int num_spline_vars_before = end_effector*POS_VARS*
                                     prev_traj_.GetPositions().at(end_effector).at(coord).GetTotalPolyVars();
        int idx_into_spline_vars = coord*prev_traj_.GetPositions().at(end_effector).at(coord).GetTotalPolyVars();
        int idx_into_poly = prev_traj_.GetPositions().at(end_effector).at(coord).GetPolyIdx(time);

        return info_.num_nodes * num_states_ + force_spline_vars_ + info_.num_nodes*num_joints_
                + num_spline_vars_before + idx_into_spline_vars + idx_into_poly*Spline::POLY_ORDER;
    }

    int MPC::GetPositionSplineIndexNoTime(int end_effector, int idx, int coord) const {
        int num_spline_vars_before = end_effector*POS_VARS*
                                     prev_traj_.GetPositions().at(end_effector).at(coord).GetTotalPolyVars();
        int idx_into_spline_vars = coord*prev_traj_.GetPositions().at(end_effector).at(coord).GetTotalPolyVars();

        return info_.num_nodes * num_states_ + force_spline_vars_ + info_.num_nodes*num_joints_
               + num_spline_vars_before + idx_into_spline_vars + idx*Spline::POLY_ORDER;
    }

    int MPC::GetVelocityIndex(int node) const {
        return info_.num_nodes * num_states_ + force_spline_vars_ + node*num_joints_;
    }

    int MPC::GetJointIndex(int node) const {
        return node * num_states_ + CentroidalModel::FLOATING_VEL_OFFSET;
    }

    Trajectory MPC::ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const {
        // Start by copying the current trajectory to keep the switching times and what not
        Trajectory traj(prev_traj_);

        // Assign all the spline information
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                std::vector<std::array<double, Spline::POLY_ORDER>> pos_vars;
                std::vector<std::array<double, Spline::POLY_ORDER>> force_vars;

                std::array<double, Spline::POLY_ORDER> pos_poly{};
                std::array<double, Spline::POLY_ORDER> force_poly{};

                // Form position spline
                for (int sp_idx = 0; sp_idx < prev_traj_.GetPositions().at(ee).at(coord).GetTotalPoly(); sp_idx++) {
                    // TODO: probably a better way to do this
                    for (int poly = 0; poly < Spline::POLY_ORDER; poly++) {
                        pos_poly.at(poly) = qp_sol.segment(GetPositionSplineIndexNoTime(ee, sp_idx, coord),
                                                           Spline::POLY_ORDER)(poly);
                    }
                    pos_vars.push_back(pos_poly);
                }

                // Form force spline
                for (int sp_idx = 0; sp_idx < prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetTotalPoly(); sp_idx++) {
                    // TODO: probably a better way to do this
                    for (int poly = 0; poly < Spline::POLY_ORDER; poly++) {
                        force_poly.at(poly) = qp_sol.segment(GetForceSplineIndexNoTime(ee, sp_idx, coord),
                                                             Spline::POLY_ORDER)(poly);
                    }
                    force_vars.push_back(force_poly);
                }

                traj.UpdatePosition(ee, coord, pos_vars);
                traj.UpdateForce(ee, coord, force_vars);
            }
        }

        // Set states and joint vels
        for (int node = 0; node < info_.num_nodes; node++) {
            // Need to convert the state back to the manifold valued state
            vector_t man_state = CentroidalModel::ConvertAlgebraStateToManifoldState(
                    qp_sol.segment(node*num_states_, num_states_), init_state);
            traj.SetState(node, man_state);
            traj.SetInputVels(node, qp_sol.segment(GetVelocityIndex(node), num_joints_));
        }

        return traj;
    }

    void MPC::EnforceFootSlipAndForceAtADistance() {
        // Go through each end effector
        // Grab the polynomials that correspond to constraints
        // Remove from the decision variables and constraints

        // For now, consider just enforcing equality constraints
        // i.e. when the ee is not in contact: f_ee(t) = {0,0,0,0}
        // when in contact: p_ee(t) = {k,k,0,0}

        // Note:
        // - Position polynomials during contact can be collapsed to a single decision variable (i.e. pos in each coord)
        // - z position is restricted to 0 in contact
        // - The surrounding polynomials need to have their derivative and end points match the next ones. This is true
        // for ALL polynomials. Can remove a lot of variables.
        // - May want to consider enforcing a swing height constraints
        // - During flight phase we should be able to completely remove its effect TBD on details
    }

} // mpc