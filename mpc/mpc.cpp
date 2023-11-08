//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "include/mpc.h"

namespace mpc {

    MPC::MPC(std::unique_ptr<MPCInfo> info, const std::string& robot_urdf) :
    info_(std::move(info)), model_(robot_urdf), num_joints_(model_.GetNumJoints()),
    num_states_(model_.GetNumConfig() + 6), num_ee_(model_.GetNumEndEffectors()), num_inputs_(num_joints_ + 6*num_ee_),
    prev_traj_(info_->num_nodes, num_states_, num_inputs_) {

        // TODO: deal with the parameterization. Is this really the size of the inputs?


        if (info_->vel_bounds.size() != num_joints_ || info_->joint_bounds.size() != num_joints_) {
            throw std::runtime_error("Velocity or joint bounds do not match the number of joints on the robot.");
        }

        // state vector: [lcom, kcom, qb, q0, ..., qnj]
        com_start_idx_ = 0;
        config_start_idx_ = 6;

        // input vector: [f1, ..., fnee, pos1, ..., posee, v0, ..., vnj]
        force_start_idx_ = 0;
        pos_start_idx_ = num_ee_*POS_VARS;
        vel_start_idx = pos_start_idx_ + num_ee_*POS_VARS;

        assert(num_inputs_ == vel_start_idx + num_joints_ + 1);

        decision_vars_ = info_->num_nodes*(num_states_ + num_inputs_);
        dynamics_constraints_ = info_->num_nodes*num_states_;
        equality_constraints_ = info_->num_nodes*POS_VARS*num_ee_;
        inequality_constraints_ = info_->num_nodes*(4*num_ee_ + 2*num_joints_);

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

        for (int i = 0; i < info_->num_nodes; i++) {
            // Get linear dynamics at each time node
            matrix_t A = model_.GetLinearDiscreteDynamicsState(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
            matrix_t B = model_.GetLinearDiscreteDynamicsInput(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
            vector_t C = model_.GetConstantDiscreteDynamics(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));

            // Update dynamics equality constraint
            data_.dynamics_constraints.block(i*num_states_ + num_states_, i*(num_states_+num_inputs_),
                                             num_states_, num_states_) = A;
            data_.dynamics_constraints.block(i*num_states_ + num_states_, i*(num_states_+num_inputs_) + num_states_,
                                             num_states_, num_inputs_) = B;
            data_.dynamics_constraints.block(i*num_states_ + num_states_, num_states_ + num_inputs_,
                                             num_states_, num_states_) = -matrix_t::Identity(num_states_, num_states_);
            data_.dynamics_constants.segment(i*num_states_ + num_states_, num_states_) = C;


            // Need to add for each end effector
            for (int j = 0; j < num_ee_; j++) {
                // Get linearization of FK at each time node
                // TODO: Specify which end effector.
                matrix_t G = model_.GetFKLinearization(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
                data_.equality_constraints.block(i * POS_VARS * num_ee_ + j * POS_VARS, i * (num_states_ + num_inputs_),
                                                 POS_VARS, num_states_) = G;

                vector_t g = model_.GetFK(prev_traj_.states_.at(i));
                data_.equality_constants.segment(i * POS_VARS * num_ee_ + j * POS_VARS, POS_VARS) = g;
            }
            data_.equality_constraints.block(i * POS_VARS * num_ee_,
                                             i * (num_states_ + num_inputs_) + num_states_ + pos_start_idx_,
                                             POS_VARS*num_ee_,
                                             num_inputs_) = matrix_t::Identity(POS_VARS*num_ee_, POS_VARS*num_ee_);

            // Inequality constraints
            // TODO: For now these are all constant, so I don't need to re-calc each time
            // Friction cone constraints can be enforced at all times (because even 0 force is in the cone). Might
            // want to be careful with the numerical instabilities doing that.
            data_.inequality_constraints.block(i*(4*num_ee_ + 2*num_joints_),
                                               i*(num_inputs_ + num_states_) + num_states_ + force_start_idx_,
                                               4*num_ee_, num_ee_) = friction_pyramid_;
            data_.inequality_constants_ub.segment(i*(4*num_ee_ + 2*num_joints_), 4*num_ee_) = vector_t::Zero(4*num_ee_);

            // Velocity bounds
            data_.inequality_constants_ub.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
                                                  num_joints_) = info_->vel_bounds;
            data_.inequality_constants_lb.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
                                                  num_joints_) = -info_->vel_bounds;
            data_.inequality_constraints.block(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_,
                                               i*(num_inputs_ + num_states_) + num_states_ + vel_start_idx,
                                               num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);

            // Configuration bounds
            data_.inequality_constants_ub.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
                                                 num_joints_) = info_->joint_bounds;
            data_.inequality_constants_lb.segment(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
                                                  num_joints_) = -info_->joint_bounds;
            data_.inequality_constraints.block(i*(4*num_ee_ + 2*num_joints_) + 4*num_ee_ + num_joints_,
                                               i*(num_inputs_ + num_states_) + pos_start_idx_,
                                               num_joints_, num_joints_) = matrix_t::Identity(num_joints_, num_joints_);


            // Get hessian approx of discrete cost function (generalized gauss newton)
            matrix_t J = GetCostHessianApprox(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
            data_.cost_quadratic.block(i*(num_states_ + num_inputs_), i*(num_states_ + num_inputs_),
                                       num_states_ + num_inputs_, num_states_ + num_inputs_) = J;


            // Get gradient of discrete cost function
            vector_t w = GetCostGradient(prev_traj_.states_.at(i), prev_traj_.inputs_.at(i));
            data_.cost_linear.segment(i*(num_states_ + num_inputs_), num_states_ + num_inputs_) = w;
        }

        // TODO: potentially search only null space of other equality constraint

        // Setup QP
        qp_solver.SetupQP(data_);

        // solve qp
        prev_traj_ = qp_solver.SolveQP();

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

    void MPC::SetFrictionPyramid() const {
        // Assumes flat ground
        Eigen::Vector3d h = {1, 0, 0};
        Eigen::Vector3d l = {0, 1, 0};
        Eigen::Vector3d n = {0, 0, 1};

        friction_pyramid_ << (h - n*info_->friction_coef).transpose(),
                            -(h + n*info_->friction_coef).transpose(),
                            (l - n*info_->friction_coef).transpose(),
                            -(l + n*info_->friction_coef).transpose();
    }

} // mpc