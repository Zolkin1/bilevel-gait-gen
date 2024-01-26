//
// Created by zolkin on 1/19/24.
//

#include "mpc_single_rigid_body.h"

namespace mpc {

    /*
     * TODO:
     * - Debug orientation
     * - Determine why the end effector constraint causes the QP to be much harder to solve (bug?)
     */

    MPCSingleRigidBody::MPCSingleRigidBody(const mpc::MPCInfo& info, const std::string& robot_urdf) :
    MPC(info, robot_urdf) {
        InitalizeQPData();
        qp_solver = std::make_unique<OSQPInterface>(data_, false);
        num_run_ = 0;
    }

    Trajectory MPCSingleRigidBody::Solve(const mpc::vector_t& state, double init_time,
                                         const std::vector<vector_3t>& ee_start_locations) {
        utils::Timer solve_timer("mpc solve");
        solve_timer.StartTimer();

        assert(state.size() == num_states_ + 1);

        init_time_ = init_time;

        utils::Timer poly_update_timer("spline update");
        poly_update_timer.StartTimer();
        prev_traj_.SetInitTime(init_time);
        prev_traj_.AddPolys(info_.integrator_dt*info_.num_nodes + init_time);
        prev_traj_.RemoveUnusedPolys(init_time);
        poly_update_timer.StopTimer();
        UpdateNumInputs();

        utils::Timer data_update_timer("data update");
        data_update_timer.StartTimer();
        UpdateQPSizes();
        data_.InitQPMats();

        prev_traj_.SetState(0, state);
        prev_qp_sol = ConvertTrajToQPVec(prev_traj_);
        data_update_timer.StopTimer();

        assert(prev_qp_sol.size() == data_.num_decision_vars);

        utils::Timer constraint_costs_timer("constraints and costs");
        constraint_costs_timer.StartTimer();
        // ----------------------- Costs ------------------------- //
        AddHessianApproxCost();
        AddGradientCost();
        AddFinalCost();

        // -------------------- Constraints ---------------------- //
        // TODO: Go through the constraint enum
        constraint_idx_ = 0;
        utils::Timer dynamics_timer("dynamics constraints");
        for (const auto& constraint : data_.constraints_) {
            switch (constraint) {
                case Constraints::Dynamics:
                    dynamics_timer.StartTimer();
                    AddDynamicsConstraints(state);
                    dynamics_timer.StopTimer();
                    break;
                case Constraints::ForceBox:
                    AddForceBoxConstraints();
                    break;
                case Constraints::FrictionCone:
                    AddFrictionConeConstraints();
                    break;
                case Constraints::EndEffectorLocation:
                    AddEELocationConstraints(ee_start_locations);
                    break;
                default:
                    throw std::runtime_error("No such constraint exists.");
            }
        }
        constraint_costs_timer.StopTimer();

        data_.ConstructSparseMats();
        data_.ConstructVectors();
        if (constraint_projection_) {
            data_.ApplyProjection();
        }

        // ----------------------- Solve ------------------------- //
        qp_solver->SetupQP(data_, prev_qp_sol);
        utils::Timer qp_solve_timer("QP solve");
        qp_solve_timer.StartTimer();
        // TODO: DMA
        const vector_t sol = qp_solver->Solve(data_);
        qp_solve_timer.StopTimer();

        if (qp_solver->GetSolveQuality() != "Solved Inaccurate" && qp_solver->GetSolveQuality() != "Solved"
            && qp_solver->GetSolveQuality() != "Max Iter Reached") {
            std::cerr << "Warning: " << qp_solver->GetSolveQuality() << std::endl;
//            throw std::runtime_error("Bad solve.");
        }

        std::cout << "Solve type: " << qp_solver->GetSolveQuality() << std::endl;

        // TODO: DMA
        const vector_t p = sol - prev_qp_sol;
        double max = 0;
        int max_idx = 0;
        for (int i = 0; i < p.size(); i++) {
            if (std::abs(sol(i) - prev_qp_sol(i)) > max) {
                max = std::abs(sol(i) - prev_qp_sol(i));
                max_idx = i;
            }
        }

//        std::cout << "max difference is: " << max << " and occurs at index: " << max_idx << std::endl;
//        std::cout << "dynamics variables: " << num_states_*(info_.num_nodes+1) << std::endl;
//        std::cout << "force variables: " << GetVelocityIndex(0) << std::endl;

        utils::Timer line_search_timer("line search");
        double alpha = 1;
        if (num_run_ >= 0 && sol.size() == prev_qp_sol.size()) {
            line_search_timer.StartTimer();
            alpha = LineSearch(p, state);
            line_search_timer.StopTimer();
        }

        num_run_++;
//        if (run_num_ > 0){
//            std::cout << (alpha*(qp_solver->GetDualSolution() - prev_dual_sol_) + prev_dual_sol_).lpNorm<Eigen::Infinity>() << std::endl;
//        }

        prev_dual_sol_ = qp_solver->GetDualSolution();

        prev_qp_sol = ((alpha * p) + prev_qp_sol).eval();

        prev_traj_ = ConvertQPSolToTrajectory(prev_qp_sol, state);
        // TODO: Turn into unit test
        vector_t temp = ConvertTrajToQPVec(prev_traj_);
        assert(temp.size() == prev_qp_sol.size());
        for (int i = 0; i < temp.size(); i++) {
            if (std::abs(temp(i) - prev_qp_sol(i)) > 1e-4) {
                std::cout << "i: " << i << " error: " << std::abs(temp(i) - prev_qp_sol(i)) << std::endl;
                std::cout << "traj->vec: " << temp(i) << std::endl;
                std::cout << "vec: " << prev_qp_sol(i) << std::endl;
                std::cout << "----------" << std::endl;
            }
        }

        solve_timer.StopTimer();

        utils::Timer stats_timer("recording stats");
        stats_timer.StartTimer();
        RecordStats(alpha, p, qp_solver->GetSolveQuality(), state,
                    solve_timer.GetElapsedTimeMilliseconds());
        stats_timer.StopTimer();

        run_num_++;

        constraint_costs_timer.PrintElapsedTime();
        line_search_timer.PrintElapsedTime();
        poly_update_timer.PrintElapsedTime();
        data_update_timer.PrintElapsedTime();
        qp_solve_timer.PrintElapsedTime();
//        stats_timer.PrintElapsedTime();
        solve_timer.PrintElapsedTime();
        std::cout << "-----------" << std::endl;

//        for (int i = 0; i < info_.num_nodes; i++) {
//            vector_3t net_force;
//            net_force << 0, 0, -9.81* model_.GetMass();
//            for (int ee = 0; ee < num_ee_; ee++) {
//                net_force += prev_traj_.GetForce(ee, GetTime(i));
//            }
//
//            std::cout << "Node: " << i << ", net force: " << net_force.transpose() << std::endl;
//        }

        prev_traj_.PrintTrajectoryToFile("mpc_demo_traj.txt");

        return prev_traj_;
    }

    void MPCSingleRigidBody::AddDynamicsConstraints(const vector_t& state) {
        // initial condition
        data_.constraint_mat_.SetDiagonalMatrix(-1, 0, 0, num_states_);
        data_.dynamics_constants.head(num_states_) =
                -model_.ConvertManifoldStateToTangentState(state, state);

        const int force_spline_vars = prev_traj_.GetTotalForceSplineVars();
        const int position_spline_vars = prev_traj_.GetTotalPosSplineVars();

        A_.resize(num_states_, num_states_);
        B_.resize(num_states_, num_inputs_);
        C_.resize(num_states_);
        C2_.resize(num_states_);

        utils::Timer get_dyn_timer("linear discrete dynamics");

        for (int node = 0; node < info_.num_nodes; node++) {
            double time = GetTime(node);
            // A can go in normally, B needs to be split
            get_dyn_timer.StartTimer();
            model_.GetLinearDynamics(prev_traj_.GetState(node),
                                     state, prev_traj_,
                                     info_.integrator_dt,
                                     time, A_, B_, C_, C2_);


//            integrator_.DiscretizeLinearDynamics(A_, B_, C_, C2_);

            A_ = matrix_t::Identity(num_states_, num_states_) + integrator_.GetDt()*A_;
            B_ = integrator_.GetDt()*B_;
            C_ = integrator_.GetDt()*C_;

//            std::cout << "A: \n" << A_ << std::endl;
//            std::cout << "A*state: \n" << A_ * model_.ConvertManifoldStateToTangentState(prev_traj_.GetState(node), state) << std::endl;
//            std::cout << "B: \n" << B_ << std::endl;
//            std::cout << "C: \n" << C_ << std::endl;

//            B_.topRows(3).setZero();

            assert(B_.cols() == num_inputs_);
            get_dyn_timer.StopTimer();
//            get_dyn_timer.PrintElapsedTime();
            data_.constraint_mat_.SetMatrix(A_, constraint_idx_ + (node + 1) * num_states_, node * num_states_);

            data_.constraint_mat_.SetDiagonalMatrix(-1, constraint_idx_ + (node + 1) * num_states_,
                                                    (node + 1) * num_states_, num_states_);

//            const matrix_t B_forces = B_.leftCols(force_spline_vars);
//            const matrix_t B_positions = B_.rightCols(position_spline_vars);

            assert(B_.cols() == force_spline_vars + position_spline_vars);
            data_.constraint_mat_.SetMatrix(B_, constraint_idx_ + (node+1)*num_states_,
                                            GetForceSplineStartIdx());

//            data_.constraint_mat_.SetMatrix(B_forces, constraint_idx_ + (node + 1) * num_states_, GetForceSplineStartIdx());
//            data_.constraint_mat_.SetMatrix(B_positions, constraint_idx_ + (node + 1) * num_states_, GetPosSplineStartIdx());
            data_.dynamics_constants.segment(constraint_idx_ + (node + 1) * num_states_, num_states_) = -C_;
        }
        constraint_idx_ = (info_.num_nodes+1)*num_states_;
    }

    int MPCSingleRigidBody::GetForceSplineStartIdx() const {
        return num_states_*(1+info_.num_nodes);     // initial condition and states for each node
    }

    int MPCSingleRigidBody::GetPosSplineStartIdx() const {
        return GetForceSplineStartIdx() + prev_traj_.GetTotalForceSplineVars();
    }

    Trajectory MPC::ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const {
        // Start by copying the current trajectory to keep the switching times and what not
        Trajectory traj(prev_traj_);

        // Assign all the spline information
        int force_idx = GetForceSplineStartIdx();
        int pos_idx = GetPosSplineStartIdx();
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < POS_VARS; coord++) {
                int vars_in_spline = traj.GetTotalPolyVars(Trajectory::SplineTypes::Force, ee, coord);
                traj.UpdateForceSpline(ee, coord, qp_sol.segment(force_idx, vars_in_spline));
                force_idx += vars_in_spline;

                if (traj.IsSplineMutable(ee, coord)) {
                    vars_in_spline = traj.GetPositions().at(ee).at(coord).GetTotalPolyVars();
                    traj.UpdatePositionSpline(ee, coord, qp_sol.segment(pos_idx, vars_in_spline));
                    pos_idx += vars_in_spline;
                }
            }
        }

        // Set states and joint vels
        for (int node = 0; node < info_.num_nodes + 1; node++) {
            // Need to convert the state back to the manifold valued state
            // TODO: DMA
            vector_t man_state = model_.ConvertTangentStateToManifoldState(
                    qp_sol.segment(node*num_states_, num_states_), init_state);
            // need to normalize quaternion returned by MPC
            Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(
                    man_state.segment(SingleRigidBodyModel::QUAT_START, SingleRigidBodyModel::QUAT_SIZE)));

            // Note the warning on the pinocchio function!
            pinocchio::quaternion::firstOrderNormalize(quat);

            man_state(SingleRigidBodyModel::QUAT_START) = quat.x();
            man_state(SingleRigidBodyModel::QUAT_START + 1) = quat.y();
            man_state(SingleRigidBodyModel::QUAT_START + 2) = quat.z();
            man_state(SingleRigidBodyModel::QUAT_START + 3) = quat.w();

            traj.SetState(node, man_state);
        }

        return traj;
    }

    void MPCSingleRigidBody::SetInitQPSizes() {
        data_.num_decision_vars = (info_.num_nodes+1)*num_states_ + num_inputs_;
        data_.num_dynamics_constraints = (info_.num_nodes+1)*num_states_;

        data_.num_cone_constraints_ = (info_.num_nodes+1)*4*num_ee_;
        data_.num_force_box_constraints_ = GetNodeIntersectMutableForces();

        data_.num_ee_location_constraints_ = (info_.num_nodes+1)*2*num_ee_;
        data_.num_start_ee_constraints_ = 2*num_ee_;
    }

    vector_t MPCSingleRigidBody::ConvertTrajToQPVec(const Trajectory& traj) const {
        // TODO: DMA
        vector_t qp_vec = vector_t::Zero(traj.GetTotalVariables());

        const int num_states = traj.GetState(0).size();

        for (int i = 0; i < info_.num_nodes+1; i++) {
            qp_vec.segment(i*(num_states-1), num_states-1) =
                    model_.ConvertManifoldStateToTangentState(traj.GetState(i), traj.GetState(0));
        }

        qp_vec.tail(num_inputs_) = traj.SplinesAsVec();

        return qp_vec;
    }

    std::vector<std::vector<Eigen::Vector3d>> MPCSingleRigidBody::CreateVizData() {
        // TODO: DMA
        std::vector<std::vector<Eigen::Vector3d>> fk_traj(num_ee_+1);
        for (int ee = 0; ee < 5; ee++) {
            for (int node = 0; node < info_.num_nodes+1; node++) {
                if (ee == 4) {
                    fk_traj.at(ee).push_back(model_.GetCOMPosition(prev_traj_.GetState(node)));
                } else {
                    fk_traj.at(ee).push_back(prev_traj_.GetEndEffectorLocation(ee, GetTime(node)));
                }
            }
        }

        return fk_traj;
    }

    void MPCSingleRigidBody::InitalizeQPData() {
        SetInitQPSizes();
        data_.InitQPMats();
    }

    // TODO: Make this also constrain the z axis - maybe
    void MPCSingleRigidBody::AddEELocationConstraints(const std::vector<vector_3t>& ee_start_locations) {
        // TODO: Do I need a rotation matrix?
        int idx = 0;
        const int pos_start_idx = GetPosSplineStartIdx();

        // -------------- End Effector Box Constraints -------------- //
        // Note: Bounds are given for the 2D space centered under the hip
        constexpr int CONSTRAINT_COORDS = 2;
        const Eigen::Vector<double, CONSTRAINT_COORDS> bounds = {0.1, 0.2};//, 0.2};

        // TODO: DMA
        matrix_t A(data_.num_ee_location_constraints_, data_.num_decision_vars);
        A.setZero();


        for (int node = 0; node < info_.num_nodes+1; node++) {
            for (int ee = 0; ee < num_ee_; ee++) {
                switch (ee) {
                    case 0:
                        data_.ee_location_ub_.segment<CONSTRAINT_COORDS>(idx) << 0.15, 0.15;
                        data_.ee_location_lb_.segment<CONSTRAINT_COORDS>(idx) << 0.05, 0.05;
                        break;
                    case 1:
                        data_.ee_location_ub_.segment<CONSTRAINT_COORDS>(idx) << 0.15, -0.05;
                        data_.ee_location_lb_.segment<CONSTRAINT_COORDS>(idx) << 0.05, -0.15;
                        break;
                    case 2:
                        data_.ee_location_ub_.segment<CONSTRAINT_COORDS>(idx) << -0.1, 0.15;
                        data_.ee_location_lb_.segment<CONSTRAINT_COORDS>(idx) << -0.25, 0.05;
                        break;
                    case 3:
                        data_.ee_location_ub_.segment<CONSTRAINT_COORDS>(idx) << -0.1, -0.05;
                        data_.ee_location_lb_.segment<CONSTRAINT_COORDS>(idx) << -0.25, -0.15;
                        break;
                }
//                data_.ee_location_ub_.segment<CONSTRAINT_COORDS>(idx) = bounds + model_.GetCOMToHip(ee).head<CONSTRAINT_COORDS>();
//                data_.ee_location_lb_.segment<CONSTRAINT_COORDS>(idx) = -bounds + model_.GetCOMToHip(ee).head<CONSTRAINT_COORDS>();

                for (int coord = 0; coord < CONSTRAINT_COORDS; coord++) {
                    A(idx, node*num_states_ + coord) = -1;

                    if (coord < 2) {
                        int vars_idx, vars_affecting;
                        std::tie(vars_idx, vars_affecting) =
                                prev_traj_.GetPositionSplineIndex(ee, GetTime(node), coord);

                        // TODO: DMA
                        vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
                                                                    ee, coord, GetTime(node));

                        A.block(idx,
                                pos_start_idx + vars_idx - vars_affecting,
                                1, vars_affecting) = vars_lin.transpose();
                    } else {
                        data_.ee_location_ub_(idx) += prev_traj_.GetEndEffectorLocation(ee, GetTime(node))(2) - 0.45;
                        data_.ee_location_lb_(idx) += prev_traj_.GetEndEffectorLocation(ee, GetTime(node))(2) - 0.45 + 0.0;
                    }
                    idx++;
                }
            }
        }

        data_.constraint_mat_.SetMatrix(A, constraint_idx_, 0);
        constraint_idx_ += idx;
        assert(idx == data_.num_ee_location_constraints_);

        // -------------- End Effector Start Constraints -------------- //
        idx = 0;
        // TODO: DMA
        matrix_t M(data_.num_start_ee_constraints_, prev_traj_.GetTotalPosSplineVars());
        M.setZero();
        for (int ee = 0; ee < num_ee_; ee++) {
            data_.start_ee_constants_.segment<2>(idx) = ee_start_locations.at(ee).head<2>();

            for (int coord = 0; coord < 2; coord++) {
                int vars_idx, vars_affecting;
                std::tie(vars_idx, vars_affecting) =
                        prev_traj_.GetPositionSplineIndex(ee, GetTime(0), coord);

                // TODO: DMA
                vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
                                                            ee, coord, GetTime(0));

                assert(vars_lin.size() == vars_affecting);

//                vars_idx = (ee*2*3) + coord*3 + vars_affecting;
                M.block(idx, vars_idx - vars_affecting, 1, vars_affecting) =
                        vars_lin.transpose();
                idx++;
            }
        }

        data_.constraint_mat_.SetMatrix(M, constraint_idx_, GetPosSplineStartIdx());

        assert(idx == data_.start_ee_constants_.size());
    }

    vector_t MPCSingleRigidBody::GetFullTargetState(double time, const vector_t& prev_state) {
        int node = GetNode(time);

        vector_t state(model_.GetFullModelConfigSpace());

        // Craft the vector and return
        man_state_t mpc_state = prev_traj_.GetState(node);
        state.head<7>() << mpc_state.head<POS_VARS>(), mpc_state.segment<4>(6);

        int idx = 7;
        for (int ee = 0; ee < num_ee_; ee++) {
            // Grab state and end effector locations at that time
            vector_3t ee_location = prev_traj_.GetEndEffectorLocation(ee, time);

            // Compute IK to get the joint angles
            // TODO: Merge this into one IK call for all the joints
            vector_t ik_joints = model_.InverseKinematics(prev_traj_.GetState(node),
                                                             ee_location, ee, prev_state);
//            vector_t ik_joints = vector_t::Zero(3);

            state.segment(idx, ik_joints.size()) = ik_joints;
            idx += ik_joints.size();
        }

        return state;
    }

//    bool MPCSingleRigidBody::ComputeParamPartials(const Trajectory& traj, QPPartials& partials, int ee, int idx) {
//        if (solve_type_.at(solve_type_.size()-1) == "Solved") {
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
//                model_.GetLinearDiscreteDynamics(traj.GetState(node), traj.GetState(0), traj,
//                                                 GetTime(node), A, B, C);
//
//                model_.ComputeLinearizationPartialWrtContactTimes(Adyn_deriv, Bdyn_deriv, Cdyn_deriv,
//                                                                  traj.GetState(node),
//                                                                  traj, GetTime(node),
//                                                                  ee, idx);
//                // TODO: Is this sparsity pattern legit?
////                std::cout << "A: \n" << A << std::endl;
////                std::cout << "Adyn_deriv: \n" << Adyn_deriv << std::endl;
////                std::cout << "Bdyn_deriv: \n" << Bdyn_deriv << std::endl;
//
//                partial_data.constraint_mat_.SetMatrix(Adyn_deriv, (node + 1) * num_states_, node * num_states_);
//
//                const matrix_t B_spline = Bdyn_deriv.leftCols(traj.GetTotalForceSplineVars());
//                const matrix_t B_vels = Bdyn_deriv.rightCols(traj.GetTotalPosSplineVars());
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
//            for (int node = 0; node < info_.num_nodes+1; node++) {
//                const int time = GetTime(node);
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    if (traj.IsSplineMutable(ee, coord)) {
//                        int vars_index, vars_affecting;
//                        std::tie(vars_index, vars_affecting) = traj.GetPositionSplineIndex(ee, time, coord);
//
//                        if (node >= num_ineq_fk_ || coord != 2) {
//                            partial_data.constraint_mat_.SetMatrix(
//                                    -traj.GetPositions().at(ee).at(coord).ComputeCoefPartialWrtTime(time, idx).transpose(),
//                                    constraint_idx_ + idx_eq + data_.num_fk_ineq_constraints_,
//                                    spline_offset + vars_index - vars_affecting);
//                            idx_eq++;
//                        } else {
//                            partial_data.constraint_mat_.SetMatrix(
//                                    -traj.GetPositions().at(ee).at(coord).ComputeCoefPartialWrtTime(time, idx).transpose(),
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
} // mpc