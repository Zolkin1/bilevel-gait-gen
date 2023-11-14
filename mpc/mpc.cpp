//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "include/mpc.h"

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
        num_states_(model_.GetPinocchioNumConfig() + 6),
        num_ee_(model_.GetNumEndEffectors()),
        prev_traj_(info.num_nodes, num_states_, num_joints_,
               CreateDefaultSwitchingTimes(info.num_switches, num_ee_, info.time_horizon),
               info.time_horizon/info.num_nodes) {

        // TODO: deal with the parameterization. Is this really the size of the inputs?

        assert(info_.ee_frames.size() == num_ee_);

        if (info_.vel_bounds.size() != num_joints_ || info_.joint_bounds.size() != num_joints_) {
            throw std::runtime_error("Velocity or joint bounds do not match the number of joints on the robot.");
        }

        // state vector: [lcom, kcom, qb, q0, ..., qnj]
        com_start_idx_ = 0;
        config_start_idx_ = 6;

        // input vector: [f1, ..., fnee, pos1, ..., posee, v0, ..., vnj]
        force_start_idx_ = 0;
        pos_start_idx_ = num_ee_*POS_VARS;
        vel_start_idx = pos_start_idx_ + num_ee_*POS_VARS;
        num_inputs_ = num_joints_ + 6*num_ee_;

        assert(num_inputs_ == vel_start_idx + num_joints_);

        decision_vars_ = info_.num_nodes*(num_states_ + num_inputs_);
        dynamics_constraints_ = info_.num_nodes*num_states_;
        equality_constraints_ = info_.num_nodes*POS_VARS*num_ee_;
        inequality_constraints_ = info_.num_nodes*(4*num_ee_ + 2*num_joints_);

        SetFrictionPyramid();
    }

    Trajectory MPC::Solve(const mpc::vector_t &centroidal_state) {
        assert(centroidal_state.size() == num_states_);

        // Reset to 0
        data_.dynamics_constraints = matrix_t::Zero(dynamics_constraints_, dynamics_constraints_);
        data_.dynamics_constants = vector_t::Zero(dynamics_constraints_);

        data_.equality_constraints = matrix_t::Zero(equality_constraints_, equality_constraints_);
        data_.equality_constants = vector_t::Zero(equality_constraints_);

        data_.inequality_constraints = matrix_t::Zero(inequality_constraints_, inequality_constraints_);
        data_.inequality_constants_ub = vector_t::Zero(inequality_constraints_);
        data_.inequality_constants_lb = vector_t::Zero(inequality_constraints_);

        data_.cost_quadratic = matrix_t::Zero(decision_vars_, decision_vars_);
        data_.cost_linear = vector_t::Zero(decision_vars_);

        // initial condition
        data_.dynamics_constraints.topLeftCorner(num_states_, num_states_) = -matrix_t::Identity(num_states_, num_states_);
        data_.dynamics_constants.head(num_states_) = -centroidal_state;

        matrix_t A, B;
        vector_t C;
        model_.GetLinearDiscreteDynamics(centroidal_state, centroidal_state, prev_traj_.GetInputs(), 0, 0.1,
                                              A, B, C);


        std::cout << "A: \n" << A << std::endl;
        std::cout << "B: \n" << B << std::endl;
        std::cout << "C: \n" << C << std::endl;

        // xbar, ubar = prev_traj

        // modify the state
        vector_t state = prev_traj_.GetStates().at(0);
        vector_t xbar = CentroidalModel::ConvertManifoldStateToAlgebraState(state, centroidal_state);
        for (int i = 0; i < state.size(); i++) {
            if (i < 6) {
                state(i) += 0.1;
            } else if (i > 11) {
                state(i) += 0.01;
            }
        }

        // Modify the inputs
        auto switching_times = mpc::MPC::CreateDefaultSwitchingTimes(info_.num_switches, 4, info_.time_horizon);
        Inputs input_new(prev_traj_.GetInputs());
        std::array<mpc::Spline, 3> forces = {mpc::Spline(2, switching_times.at(0), true), mpc::Spline(2, switching_times.at(0), true),
                                             mpc::Spline(2, switching_times.at(0), true)};
        std::array<double, 4> vars = {1.1, 0, 1.1, 0};
        for (int coord = 0; coord < 3; coord++) {
            for (int poly = 0; poly < input_new.GetForces().at(0).at(0).GetTotalPoly(); poly++) {
                forces.at(coord).SetPolyVars(poly, vars);
            }
        }

        for (int i = 0; i < num_ee_; i++) {
            input_new.SetEndEffectorForce(i, forces);
        }

        vector_t joint_vels(num_joints_);
        joint_vels << .1,0,0, 0,.2,.3, 0,.1,.1, .4,.5,.1;
        input_new.SetJointVels(joint_vels, 0.1);

        vector_t state_alg = CentroidalModel::ConvertManifoldStateToAlgebraState(state, centroidal_state);

        vector_t xdot_true = model_.CalcDynamics(state, input_new, 0.1);
        vector_t xdot_approx = A*(state_alg - xbar) +
                B*(input_new.GetInputVector(0.1) - prev_traj_.GetInputs().GetInputVector(0.1)) + C;

        std::cout << "xdot true: \n" << xdot_true << std::endl;
        std::cout << "xdot approx: \n" << xdot_approx << std::endl;


//        for (int i = 0; i < info_->num_nodes; i++) {
//            // Get linear dynamics at each time node
//            matrix_t A = model_.GetLinearDiscreteDynamicsState(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
//            matrix_t B = model_.GetLinearDiscreteDynamicsInput(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
//            vector_t C = model_.GetConstantDiscreteDynamics(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
//
//            // Update dynamics equality constraint
//            data_.dynamics_constraints.block(i*num_states_ + num_states_, i*(num_states_+num_inputs_),
//                                             num_states_, num_states_) = A;
//            data_.dynamics_constraints.block(i*num_states_ + num_states_, i*(num_states_+num_inputs_) + num_states_,
//                                             num_states_, num_inputs_) = B;
//            data_.dynamics_constraints.block(i*num_states_ + num_states_, num_states_ + num_inputs_,
//                                             num_states_, num_states_) = -matrix_t::Identity(num_states_, num_states_);
//            data_.dynamics_constants.segment(i*num_states_ + num_states_, num_states_) = C;
//
//
//            // Need to add for each end effector
//            for (int j = 0; j < num_ee_; j++) {
//                // Get linearization of FK at each time node
//                matrix_t G = model_.GetFKLinearization(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
//                data_.equality_constraints.block(i * POS_VARS * num_ee_ + j * POS_VARS, i * (num_states_ + num_inputs_),
//                                                 POS_VARS, num_states_) = G;
//
//                vector_t g = model_.GetFK(prev_traj_.states_.at(i));
//                data_.equality_constants.segment(i * POS_VARS * num_ee_ + j * POS_VARS, POS_VARS) = g;
//            }
//            data_.equality_constraints.block(i * POS_VARS * num_ee_,
//                                             i * (num_states_ + num_inputs_) + num_states_ + pos_start_idx_,
//                                             POS_VARS*num_ee_,
//                                             num_inputs_) = matrix_t::Identity(POS_VARS*num_ee_, POS_VARS*num_ee_);
//
//            // Inequality constraints
//            // TODO: For now these are all constant, so I don't need to re-calc each time
//            // Friction cone constraints can be enforced at all times (because even 0 force is in the cone). Might
//            // want to be careful with the numerical instabilities doing that.
//            data_.inequality_constraints.block(i*(4*num_ee_ + 2*num_joints_),
//                                               i*(num_inputs_ + num_states_) + num_states_ + force_start_idx_,
//                                               4*num_ee_, num_ee_) = friction_pyramid_;
//            data_.inequality_constants_ub.segment(i*(4*num_ee_ + 2*num_joints_), 4*num_ee_) = vector_t::Zero(4*num_ee_);
//
//            // Velocity bounds
//            data_.inequality_constants_ub.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
//                                                  num_joints_) = info_->vel_bounds;
//            data_.inequality_constants_lb.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
//                                                  num_joints_) = -info_->vel_bounds;
//            data_.inequality_constraints.block(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
//                                               i*(num_inputs_ + num_states_) + num_states_ + vel_start_idx,
//                                               num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);
//
//            // Configuration bounds
//            data_.inequality_constants_ub.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
//                                                 num_joints_) = info_->joint_bounds;
//            data_.inequality_constants_lb.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
//                                                  num_joints_) = -info_->joint_bounds;
//            data_.inequality_constraints.block(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
//                                               i*(num_inputs_ + num_states_) + pos_start_idx_,
//                                               num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);
//
//
//            // Get hessian approx of discrete cost function (generalized gauss newton)
//            matrix_t J = GetCostHessianApprox(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
//            data_.cost_quadratic.block(i*(num_states_ + num_inputs_), i*(num_states_ + num_inputs_),
//                                       num_states_ + num_inputs_, num_states_ + num_inputs_) = J;
//
//
//            // Get gradient of discrete cost function
//            vector_t w = GetCostGradient(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
//            data_.cost_linear.segment(i*(num_states_ + num_inputs_), num_states_ + num_inputs_) = w;
//        }
//
//        // TODO: potentially search only null space of other equality constraint
//
//        // Setup QP
//        qp_solver.SetupQP(data_);
//
//        // solve qp
//        prev_traj_ = qp_solver.SolveQP();

        return prev_traj_;

        // --------------------------------------------------------------- //
        // ---------------------- At each time node: --------------------- //
        // --------------------------------------------------------------- //

        // ---------------- Dynamics equality constraints ---------------- //
        // Get discrete centroidal dynamics

        // ------------------ Other equality constraints ----------------- //
        // foot position input = forward kinematics of the joints at each point

        // -------------------- Inequality constraints ------------------- //
        // Joint velocity constraints (not sure if I need these)
        // Joint angle constraints
        // Friction cone constraints

        // ---------------------- Cost function ------------------------- //
        // Get discrete cost function

        // ---------------- Form lagrangian ---------------- //

        // ---------------- Gradient and Hessian approx of lagrangian ---------------- //
        // Note: these need to happen at each node independetly (I believe), so can potentially compute in parallel
        // Get linearization of discrete centroidal dynamics
        // Linearization of other equality constraints
        // Linearization of inequality constraints

        // ---------------- QP Solve ---------------- //
        // Update QP solver
        // Solve

        // line search?


        // Concise:
        // Initial conditions
        // Get linear dynamics at each time node
        // Get linearization of FK at each time node
        // Velocity constraints assumed polytopic and constant over time
        // Joint constraints assumed polytopic and constant over time
        // Friction pyramids are used
        // Get hessian approx of discrete cost function (generalized gauss newton)
        // Get gradient of discrete cost function
        // Contribution of constraints to the hessian may need to be ignored
        // Potentially search only in the null space of the fk constraint
        // Setup QP
        // Solve QP (to completion?)


    }

    void MPC::SetWarmStartTrajectory(const mpc::Trajectory &trajectory) {
        prev_traj_ = trajectory;
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
    };

} // mpc