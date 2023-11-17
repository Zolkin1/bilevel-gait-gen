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

        data_.num_decision_vars = info_.num_nodes*(num_states_ + num_joints_) +
                prev_traj_.GetInputs().GetTotalForceSplineVars() + prev_traj_.GetTotalPosSplineVars();
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

        data_.dynamics_constraints.block((node+1)*num_states_,
                                         node*num_joints_ + num_states_*info_.num_nodes + force_spline_vars,
                                         num_states_, num_joints_) = B_vels;
        data_.dynamics_constraints.block((node+1)*num_states_,
                                         num_states_*info_.num_nodes,
                                         num_states_, force_spline_vars) = B_spline;
        data_.dynamics_constants.segment((node+1)*num_states_, num_states_) = -C;

        // ----------------------------- Testing --------------------------- //
        std::cout << "A: \n" << A << std::endl;
        std::cout << "B: \n" << B << std::endl;
        std::cout << "C: \n" << C << std::endl;

        // Modify state:
        vector_t state_new = prev_traj_.GetStates().at(0);
        vector_t xbar = CentroidalModel::ConvertManifoldStateToAlgebraState(state_new, state);
        for (int i = 0; i < state_new.size(); i++) {
            if (i < 6) {
                state_new(i) += 0.1;
            } else if (i > 11) {
                state_new(i) += 0.01;
            }
        }

        // Modify Inputs
        auto switching_times = mpc::MPC::CreateDefaultSwitchingTimes(info_.num_switches, 4, info_.time_horizon);
        Inputs input_new(prev_traj_.GetInputs());
        std::array<mpc::Spline, 3> forces = {mpc::Spline(2, switching_times.at(0), true), mpc::Spline(2, switching_times.at(0), true),
                                             mpc::Spline(2, switching_times.at(0), true)};

        for (int coord = 0; coord < 3; coord++) {
            for (int poly = 0; poly < input_new.GetForces().at(0).at(coord).GetNumPolyTimes(); poly++) {
                std::vector<double> vars;
                vars.push_back(1.1);
                if (input_new.GetForces().at(0).at(coord).GetPolyVars().at(poly).size() == 2) {
                    vars.push_back(0);
                }

                forces.at(coord).SetPolyVars(poly, vars);
            }
        }

        for (int i = 0; i < num_ee_; i++) {
            input_new.SetEndEffectorForce(i, forces);
        }

        vector_t joint_vels(num_joints_);
        joint_vels << .1,0,0, 0,.2,.3, 0,.1,.1, .4,.5,.1;
        input_new.SetJointVels(joint_vels, 0.1);
        vector_t state_alg = CentroidalModel::ConvertManifoldStateToAlgebraState(state_new, state);
        vector_t xdot_true = model_.CalcDynamics(state_new, input_new, 0.1);
        vector_t xdot_approx = A*(state_alg - xbar) +
                               B*(input_new.GetInputVector(0.1) - prev_traj_.GetInputs().GetInputVector(0.1)) + C;
        std::cout << "xdot true: \n" << xdot_true << std::endl;
        std::cout << "xdot approx: \n" << xdot_approx << std::endl;
    }


    void MPC::AddFKConstraints(const vector_t& state, double time, int node) {
        matrix_t G;
        vector_t g;

        int spline_offset = info_.num_nodes*num_states_ + prev_traj_.GetInputs().GetTotalForceSplineVars();

        for (int ee = 0; ee < num_ee_; ee++) {
            model_.GetFKLinearization(state, prev_traj_.GetInputs(), ee, G, g);
            data_.equality_constraints.block(node * POS_VARS * num_ee_ + ee * POS_VARS,
                                             node * num_states_ + CentroidalModel::MOMENTUM_OFFSET,
                                             POS_VARS, num_joints_ + CentroidalModel::FLOATING_VEL_OFFSET) = G;

            data_.equality_constants.segment(node * POS_VARS * num_ee_ + ee * POS_VARS, POS_VARS) = -g;

            std::cout << "G for ee #" << ee << ": \n" << G << std::endl;
            std::cout << "g for ee #" << ee << ": \n" << g << std::endl;

            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(ee, time, coord);
                std::cout << "vars index: " << vars_index << std::endl;

                data_.equality_constraints.block(node * POS_VARS * num_ee_ + ee * POS_VARS + coord,
                                                 spline_offset + vars_index - vars_affecting,
                                                 1,
                                                 vars_affecting) =
                                                         -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose();
            std::cout << -prev_traj_.GetPositions().at(ee).at(coord).GetPolyVarsLin(time).transpose() << std::endl;
            }
        }
    }

    void MPC::AddInequalityConstraints(const vector_t& state, double time, int node) {
        // Friction pyramid
        for (int ee = 0; ee < num_ee_; ee++) {
            std::vector<std::array<Spline, 3>> forces = prev_traj_.GetInputs().GetForces();
            for (int coord = 0; coord < POS_VARS; coord++) {
                vector_t vars_lin = prev_traj_.GetInputs().GetForces().at(ee).at(coord).GetPolyVarsLin(time);
                Eigen::Vector4d temp;

                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = prev_traj_.GetInputs().GetForceSplineIndex(ee, time, coord);

                for (int i = 0; i < vars_lin.size(); i++) {
                    temp = friction_pyramid_.col(coord)*vars_lin(i);

                    data_.inequality_constraints.block(node * (4 * num_ee_ + 2 * num_joints_),
                                                       info_.num_nodes*num_states_ + vars_index - vars_affecting,
                                                       4, 1) = temp;
                }
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

    int MPC::GetVelocityIndex(int node) const {
        return info_.num_nodes * num_states_ + prev_traj_.GetInputs().GetTotalForceSplineVars() + node*num_joints_;
    }

    int MPC::GetJointIndex(int node) const {
        return node * num_states_ + CentroidalModel::FLOATING_VEL_OFFSET;
    }

    Trajectory MPC::ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const {
        // Start by copying the current trajectory to keep the switching times and what not
        Trajectory traj(prev_traj_);

        // Assign all the spline information
        int force_idx = info_.num_nodes*num_states_;
        int pos_idx = info_.num_nodes*(num_states_ + num_joints_) + traj.GetInputs().GetTotalForceSplineVars();
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_in_spline = traj.GetInputs().GetForces().at(ee).at(coord).GetTotalPolyVars();
                traj.UpdateForceSpline(ee, coord, qp_sol.segment(force_idx, vars_in_spline));
                force_idx += vars_in_spline;

                vars_in_spline = traj.GetPositions().at(ee).at(coord).GetTotalPolyVars();
                traj.UpdatePositionSpline(ee, coord, qp_sol.segment(pos_idx, vars_in_spline));
                pos_idx += vars_in_spline;
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