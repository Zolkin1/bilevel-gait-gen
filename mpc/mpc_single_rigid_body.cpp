//
// Created by zolkin on 1/19/24.
//

#include "mpc_single_rigid_body.h"

namespace mpc {

    MPCSingleRigidBody::MPCSingleRigidBody(const mpc::MPCInfo& info, const std::string& robot_urdf) :
    MPC(info, robot_urdf) {
        InitalizeQPData();
        if (info_.verbose == Optimization || info_.verbose == All) {
            qp_solver = std::make_unique<ClarabelInterface>(data_, true);
        } else {
            qp_solver = std::make_unique<ClarabelInterface>(data_, false);
        }
        num_run_ = 0;

        prev_qp_sol = vector_t::Zero(data_.num_decision_vars);
        line_search_res_ = vector_t::Zero(data_.num_decision_vars);

        ee_bounds_ = info_.ee_box_size;
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

        if (info_.verbose == All) {
            std::cout << "run num: " << num_run_ << ", force vars: " << prev_traj_.GetTotalForceSplineVars() <<
                      ", position vars: " << prev_traj_.GetTotalPosSplineVars() << std::endl;
        }

        utils::Timer data_update_timer("data update");
        data_update_timer.StartTimer();
        UpdateQPSizes();
        data_.InitQPMats();

        prev_traj_.SetState(0, state); // if i change to 1 I still get errors, removing this entirely is also less error
        prev_qp_sol = ConvertTrajToQPVec(prev_traj_);
        data_update_timer.StopTimer();

        assert(prev_qp_sol.size() == data_.num_decision_vars);

        utils::Timer constraint_costs_timer("constraints and costs");
        constraint_costs_timer.StartTimer();
        // ----------------------- Costs ------------------------- //
        AddForceCost(info_.force_cost);
        AddHessianApproxCost();
        AddGradientCost();
        AddFinalCost();
        AddDiagonalCost();

//        AddEEPosCost();


//        data_.cost_mat_.SetDiagonalMatrix(1, 0, 0, data_.num_decision_vars);

        // -------------------- Constraints ---------------------- //
        constraint_idx_ = 0;
        utils::Timer dynamics_timer("dynamics constraints");
        for (const auto& constraint : data_.constraints_) {
            switch (constraint) {
                case Constraints::Dynamics:
                    dynamics_timer.StartTimer();
                    AddDynamicsConstraints(prev_traj_.GetState(0));
                    dynamics_timer.StopTimer();
                    break;
                case Constraints::ForceBox:
                    AddForceBoxConstraints();
                    break;
                case Constraints::FrictionCone:
                    AddFrictionConeConstraints();
                    break;
                case Constraints::EndEffectorLocation:
                    AddEELocationConstraints();
                    break;
                case Constraints::TDPosition:
                    AddTDPositionConstraints();
                    break;
                case Constraints::Raibert:
                    AddRaibertHeuristic();
                    break;
                case Constraints::EndEffectorStart:
                    AddEEStartConstraints(ee_start_locations);
                    break;
                default:
                    throw std::runtime_error("No such constraint exists.");
            }
        }
        constraint_costs_timer.StopTimer();

        data_.ConstructSparseMats();
        data_.ConstructVectors();

        // ----------------------- Solve ------------------------- //
        qp_solver->SetupQP(data_, prev_qp_sol);
        utils::Timer qp_solve_timer("QP solve");
        qp_solve_timer.StartTimer();
        // TODO: DMA
        vector_t sol;
        try {
            sol = qp_solver->Solve(data_);
        } catch (const std::string& error) {
            std::cerr << error << " Using previous solution instead." << std::endl;
            sol = prev_qp_sol;

//            std::cout << "Force box constraint: \n" <<
//                    data_.sparse_constraint_.middleRows(data_.num_dynamics_constraints,
//                                                        data_.num_force_box_constraints_) << std::endl;

//            std::cout << "Force box ub: " << data_.ub_.segment(data_.num_dynamics_constraints, data_.num_force_box_constraints_) .transpose()<< std::endl;
//            std::cout << "Force box lb: " << data_.lb_.segment(data_.num_dynamics_constraints, data_.num_force_box_constraints_).transpose() << std::endl;

//            throw std::runtime_error("Primal infeasible.");
        }

        qp_solve_timer.StopTimer();

//        std::cout << "friction cone constraint dual: " << qp_solver->GetDualSolution().transpose().segment(
//                data_.num_dynamics_constraints + data_.num_force_box_constraints_, data_.num_cone_constraints_) << std::endl;

        if (qp_solver->GetSolveQuality() != SolvedInacc && qp_solver->GetSolveQuality() != Solved
            && qp_solver->GetSolveQuality() != MaxIter) {
            std::cerr << "Warning: " << qp_solver->GetSolveQualityAsString() << std::endl;
            IncreaseEEBox();

//            throw std::runtime_error("Bad solve.");
        } else {
            DecreaseEEBox();
        }

        if (info_.verbose == All) {
            std::cout << "Solve type: " << qp_solver->GetSolveQualityAsString() << std::endl;
        }

        // TODO: DMA
        const vector_t p = sol - prev_qp_sol;

//        double max = 0;
//        int max_idx = 0;
//        for (int i = 0; i < p.size(); i++) {
//            if (std::abs(sol(i) - prev_qp_sol(i)) > max) {
//                max = std::abs(sol(i) - prev_qp_sol(i));
//                max_idx = i;
//            }
//        }
//        std::cout << "max difference is: " << max << " and occurs at index: " << max_idx << std::endl;
//        std::cout << "force variables start: " << num_states_*(info_.num_nodes+1) << std::endl;
//        std::cout << "position variables start: " << GetPosSplineStartIdx() << std::endl;

        utils::Timer line_search_timer("line search");
        double alpha = 1;
        if (num_run_ >= 0 && sol.size() == prev_qp_sol.size()) {
            line_search_timer.StartTimer();
            alpha = LineSearch(p, prev_traj_.GetState(0));
            line_search_timer.StopTimer();
        }

        prev_dual_sol_ = qp_solver->GetDualSolution();

        prev_qp_sol = ((alpha * p) + prev_qp_sol).eval(); // TODO: Remove eval?

        prev_traj_ = ConvertQPSolToTrajectory(prev_qp_sol, prev_traj_.GetState(0));

        solve_timer.StopTimer();

        utils::Timer stats_timer("recording stats");
        stats_timer.StartTimer();
        RecordStats(alpha, p, qp_solver->GetSolveQuality(), prev_traj_.GetState(0),
                    solve_timer.GetElapsedTimeMilliseconds(), GetCostValue(prev_qp_sol)); //GetCostValue(sol));
        stats_timer.StopTimer();

        run_num_++;

//        std::cout << "Cost: " << std::setprecision(17) << GetCost() << std::endl;

        if (info_.verbose == Timing || info_.verbose == All) {
            constraint_costs_timer.PrintElapsedTime();
            line_search_timer.PrintElapsedTime();
            poly_update_timer.PrintElapsedTime();
            data_update_timer.PrintElapsedTime();
            qp_solve_timer.PrintElapsedTime();
            solve_timer.PrintElapsedTime();
            std::cout << "-----------" << std::endl;
        }

//        for (int i = 0; i < info_.num_nodes; i++) {
//            vector_3t net_force;
//            net_force << 0, 0, -9.81* model_.GetMass();
//            for (int ee = 0; ee < num_ee_; ee++) {
//                net_force += prev_traj_.GetForce(ee, GetTime(i));
//            }
//
//            std::cout << "Node: " << i << ", net force: " << net_force.transpose() << std::endl;
//        }

//        prev_traj_.PrintTrajectoryToFile("mpc_demo_traj.txt");

        num_run_++;

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

            assert(B_.cols() == num_inputs_);
            get_dyn_timer.StopTimer();

            data_.constraint_mat_.SetMatrix(A_, constraint_idx_ + (node + 1) * num_states_, node * num_states_);

            data_.constraint_mat_.SetDiagonalMatrix(-1, constraint_idx_ + (node + 1) * num_states_,
                                                    (node + 1) * num_states_, num_states_);

            assert(B_.cols() == force_spline_vars + position_spline_vars);
            data_.constraint_mat_.SetMatrix(B_, constraint_idx_ + (node+1)*num_states_,
                                            GetForceSplineStartIdx());

            data_.dynamics_constants.segment(constraint_idx_ + (node + 1) * num_states_, num_states_) = -C_;
        }
        constraint_idx_ += (info_.num_nodes+1)*num_states_; // note the += rather than =
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

                if (coord < 2) {
                    vars_in_spline = traj.GetTotalPolyVars(Trajectory::SplineTypes::Position, ee ,coord); //traj.GetPositions().at(ee).at(coord).GetTotalPolyVars();
                    traj.UpdatePositionSpline(ee, coord, qp_sol.segment(pos_idx, vars_in_spline));
                    pos_idx += vars_in_spline;
                }
            }
        }

        assert(pos_idx == qp_sol.size());
        assert(force_idx == GetPosSplineStartIdx());

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

        data_.num_cone_constraints_ = GetNumFricConeConstraints();
//        data_.num_cone_constraints_ = (info_.num_nodes+1)*4*num_ee_;

        if (using_clarabel_) {
            data_.num_force_box_constraints_ = GetNumForceBoxConstraints(); //2*GetNodeIntersectMutableForces();
            data_.num_ee_location_constraints_ = 2*(info_.num_nodes-(EE_NODE_START-1))*2*num_ee_; // 2*(info_.num_nodes+1)*2*num_ee_
            data_.num_td_pos_constraints_ = GetNumTDConstraints();
            data_.num_raibert_constraints_ = GetNumRaibertConstraints();
        } else {
            data_.num_force_box_constraints_ = GetNodeIntersectMutableForces(); // TODO: Change
            data_.num_ee_location_constraints_ = (info_.num_nodes+1)*2*num_ee_;
        }

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


    void MPCSingleRigidBody::AddEELocationConstraints() {
        // TODO: Do I need a rotation matrix?
        const int pos_start_idx = GetPosSplineStartIdx();

        // -------------- End Effector Box Constraints -------------- //
        // Note: Bounds are given for the 2D space centered under the hip
        constexpr int CONSTRAINT_COORDS = 2;
        const Eigen::Vector<double, CONSTRAINT_COORDS> bounds = info_.ee_box_size/2;

        // TODO: DMA
        matrix_t A(data_.num_ee_location_constraints_, data_.num_decision_vars);
        A.setZero();

        int extra_runs = 1;
        if (using_clarabel_) {
            extra_runs = 2;
        }
        int idx = 0;
        for (int i = 0; i < extra_runs; i++) {
            for (int node = EE_NODE_START; node < info_.num_nodes + 1; node++) {
                for (int ee = 0; ee < num_ee_; ee++) {
                    if (i == 0) {
                        data_.ee_location_ub_.segment<CONSTRAINT_COORDS>(idx) =
                                bounds + model_.GetCOMToHip(ee).head<CONSTRAINT_COORDS>();
                        data_.ee_location_lb_.segment<CONSTRAINT_COORDS>(idx) =
                                -bounds + model_.GetCOMToHip(ee).head<CONSTRAINT_COORDS>();
                    }

                    for (int coord = 0; coord < CONSTRAINT_COORDS; coord++) {
                        if (i == 0) {
                            A(idx, (node) * num_states_ + coord) = -1;
                        } else {
                            A(idx, (node) * num_states_ + coord) = 1;
                        }

                        int vars_idx, vars_affecting;
                        std::tie(vars_idx, vars_affecting) =
                                prev_traj_.GetPositionSplineIndex(ee, GetTime(node), coord);

                        // TODO: DMA
                        vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
                                                                    ee, coord, GetTime(node));

                        if (i == 0) {
                            A.block(idx,
                                    pos_start_idx + vars_idx,
                                    1, vars_affecting) = vars_lin.transpose();
                        } else {
                            A.block(idx,
                                    pos_start_idx + vars_idx,
                                    1, vars_affecting) = -vars_lin.transpose();
                        }

                        idx++;
                    }
                }
            }
        }

        data_.constraint_mat_.SetMatrix(A, constraint_idx_, 0);
        assert(idx == data_.num_ee_location_constraints_);
        constraint_idx_ += idx;
    }

    void MPCSingleRigidBody::AddEEStartConstraints(const std::vector<vector_3t>& ee_start_locations) {
        // -------------- End Effector Start Constraints -------------- //
        int idx = 0;
//         TODO: DMA
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

                M.block(idx, vars_idx, 1, vars_affecting) =
                        vars_lin.transpose();
                idx++;
            }
        }

        data_.constraint_mat_.SetMatrix(M, constraint_idx_, GetPosSplineStartIdx());

        assert(idx == data_.start_ee_constants_.size());
        constraint_idx_ += idx;
    }

    vector_t MPCSingleRigidBody::GetFullTargetState(double time, const vector_t& prev_state) {
        int node = GetNode(time);

//        vector_t state(model_.GetFullModelConfigSpace());

        // Craft the vector and return
//        man_state_t mpc_state = prev_traj_.GetState(node);
//        state.head<7>() << mpc_state.head<POS_VARS>(), mpc_state.segment<4>(6);

//        int idx = 7;
        std::vector<vector_3t> ee_locations(num_ee_);
        for (int ee = 0; ee < num_ee_; ee++) {
            // Grab state and end effector locations at that time
            ee_locations.at(ee) = prev_traj_.GetEndEffectorLocation(ee, time);
        }

        // Compute IK to get the joint angles
        return model_.InverseKinematics(prev_traj_.GetState(node),
                                                      ee_locations, prev_state,
                                                      info_.joint_bounds_ub, info_.joint_bounds_lb);
//        state = ik_state;
//
//        return state;
    }

    std::vector<Eigen::Vector2d> MPCSingleRigidBody::GetEEBoxCenter() {
        std::vector<Eigen::Vector2d> box_centers(num_ee_);
        for (int ee = 0; ee < num_ee_; ee++) {
            box_centers.at(ee) = model_.GetCOMToHip(ee).head<2>();
        }

        return box_centers;
    }

    void MPCSingleRigidBody::UpdateWholeBodyTrajectory(mpc::Trajectory& traj) {
        // TODO: Figure out a better initial guess
        vector_t prev_state = vector_t::Zero(model_.GetNumManifoldStates());
        for (int i = 0; i < 20; i++) {
            std::vector<vector_3t> ee_locations(num_ee_);
            for (int ee = 0; ee < num_ee_; ee++) {
                // Grab state and end effector locations at that time
                ee_locations.at(ee) = prev_traj_.GetEndEffectorLocation(ee, traj.GetTime(i));
            }
            traj.UpdateFullConfig(i, model_.InverseKinematics(traj.GetState(i),
                                                              ee_locations, prev_state,
                                                              info_.joint_bounds_ub,
                                                              info_.joint_bounds_lb));

            // TODO: Do the velocity IK
            traj.UpdateFullVelocity(i, vector_t::Zero(model_.GetFullModelConfigSpace()-1));
        }
    }

    // TODO: Make this cleaner
    bool MPCSingleRigidBody::ComputeParamPartialsOSQP(const Trajectory& traj, QPPartials& partials, int ee, int contact_time_idx) {
        if (solve_type_.at(solve_type_.size()-1) == Solved) {
            QPData partial_data = data_;
            partial_data.InitQPMats();

            // For now the cost terms are not effected by the contact times

            partials.dP.resize(data_.num_decision_vars, data_.num_decision_vars);
            partials.dq.resize(data_.num_decision_vars);
            partials.dA.resize(data_.GetTotalNumConstraints(), data_.num_decision_vars);
            partials.du.resize(data_.GetTotalNumConstraints());
            partials.dl.resize(data_.GetTotalNumConstraints());

            // TODO: Make this more general
            partials.dP.setZero();
            partials.dq.setZero();

            // Force box constraints have no effect
            // Cone constraints have on effect

            // -------------- Dynamics constraints have an effect -------------- //
            partial_data.constraint_mat_.Reserve(data_.sparse_constraint_.nonZeros());

            matrix_t Adyn_deriv, Bdyn_deriv;
            vector_t Cdyn_deriv;

            for (int node = 0; node < info_.num_nodes; node++) {
                model_.ComputeLinearizationPartialWrtContactTimes(Adyn_deriv, Bdyn_deriv, Cdyn_deriv,
                                                                  traj.GetState(node),
                                                                  traj, GetTime(node),
                                                                  ee, contact_time_idx);
                Adyn_deriv = integrator_.GetDt()*Adyn_deriv;
                Bdyn_deriv = integrator_.GetDt()*Bdyn_deriv;
                Cdyn_deriv = integrator_.GetDt()*Cdyn_deriv;


                partial_data.constraint_mat_.SetMatrix(Adyn_deriv, (node + 1) * num_states_, node * num_states_);

                partial_data.constraint_mat_.SetMatrix(Bdyn_deriv, (node+1)*num_states_,
                                                GetForceSplineStartIdx());

                partial_data.dynamics_constants.segment((node + 1) * num_states_, num_states_) = -Cdyn_deriv;
            }


            // -------------- EE Box constraints have an effect -------------- //
            const int spline_offset = GetPosSplineStartIdx();
            int ee_constraint_start = data_.num_dynamics_constraints + data_.num_force_box_constraints_
                                            + data_.num_cone_constraints_;
            int idx = ee*(data_.num_ee_location_constraints_/num_ee_);

            for (int node = 0; node < info_.num_nodes+1; node++) {
                const double time = GetTime(node);
                for (int coord = 0; coord < 2; coord++) {
                    int vars_index, vars_affecting;
                    std::tie(vars_index, vars_affecting) = traj.GetPositionSplineIndex(ee, time, coord);

                    // position coef dependence
                    const vector_t pos_coef_partials = traj.GetPositionCoefPartialsWrtContactTime(ee, coord, time,
                                                                                                  contact_time_idx);

                    partial_data.constraint_mat_.SetMatrix(pos_coef_partials.transpose(),
                                                           ee_constraint_start + idx,
                                                           spline_offset + vars_index - vars_affecting);
                    idx++;
                }
            }

            assert(idx - ee*(data_.num_ee_location_constraints_/num_ee_) == data_.num_ee_location_constraints_/num_ee_);



            ee_constraint_start += data_.num_ee_location_constraints_;
            idx = 2*ee;

//         TODO: DMA
            matrix_t M(data_.num_start_ee_constraints_, traj.GetTotalPosSplineVars());
            M.setZero();
            idx = 0;
            for (int coord = 0; coord < 2; coord++) {
                int vars_idx, vars_affecting;
                std::tie(vars_idx, vars_affecting) =
                        prev_traj_.GetPositionSplineIndex(ee, GetTime(0), coord);

                // Position coef dependence

                // TODO: DMA
                vector_t pos_coef_partials = traj.GetPositionCoefPartialsWrtContactTime(ee, coord, GetTime(0), contact_time_idx);

                M.block(idx, vars_idx, 1, vars_affecting) =
                        pos_coef_partials.transpose();
                idx++;
            }

//            assert(idx == data_.num_start_ee_constraints_);

//            data_.constraint_mat_.SetMatrix(M, constraint_idx_, GetPosSplineStartIdx());

            partial_data.ConstructSparseMats();
            partial_data.ConstructVectors();

            partials.dA = partial_data.sparse_constraint_;
            partials.dl = partial_data.lb_;
            partials.du = partial_data.ub_;

            return true;
        } else {
            return false;
        }
    }

    bool MPCSingleRigidBody::ComputeParamPartialsClarabel(const mpc::Trajectory& traj, mpc::QPPartials& partials,
                                                          int ee, int contact_time_idx) {
        if (solve_type_.at(solve_type_.size()-1) == Solved) {
            QPData partial_data = data_;
            partial_data.InitQPMats();

            // For now the cost terms are not effected by the contact times

            // This is zero, so skip it for speed for now
//            partials.dP.resize(data_.num_decision_vars, data_.num_decision_vars);

            partials.dq.resize(data_.num_decision_vars);
            partials.dA.resize(data_.num_equality_, data_.num_decision_vars);
            partials.dG.resize(data_.num_inequality_, data_.num_decision_vars);
            partials.db.resize(data_.num_equality_);
            partials.dh.resize(data_.num_inequality_); //- data_.num_cone_constraints_);

            partials.SetZero();

            // Force box constraints have no effect
            // Cone constraints have on effect

            // -------------- Dynamics constraints have an effect -------------- //

            utils::SparseMatrixBuilder A_builder;
            A_builder.Reserve(data_.sparse_constraint_.nonZeros()); // TODO: Reserves too much

            utils::SparseMatrixBuilder G_builder;
            G_builder.Reserve(data_.sparse_constraint_.nonZeros()); // TODO: Reserves too much


            int eq_idx = 0;
            int ineq_idx = 0;
            for (int i = 0; i < data_.constraints_.size(); i++) {
                if (data_.constraints_.at(i) == Dynamics) {
                    utils::Timer dyn_partials_timer("param partials dynamics");
                    dyn_partials_timer.StartTimer();
                    matrix_t Adyn_deriv, Bdyn_deriv;
                    vector_t Cdyn_deriv;

                    for (int node = 0; node < info_.num_nodes; node++) {
                        model_.ComputeLinearizationPartialWrtContactTimes(Adyn_deriv, Bdyn_deriv, Cdyn_deriv,
                                                                          traj.GetState(node),
                                                                          traj, GetTime(node),
                                                                          ee, contact_time_idx);
                        Adyn_deriv = integrator_.GetDt() * Adyn_deriv;
                        Bdyn_deriv = integrator_.GetDt() * Bdyn_deriv;
                        Cdyn_deriv = integrator_.GetDt() * Cdyn_deriv;


                        A_builder.SetMatrix(Adyn_deriv, eq_idx + (node + 1) * num_states_, node * num_states_);

                        A_builder.SetMatrix(Bdyn_deriv, eq_idx + (node + 1) * num_states_,
                                            GetForceSplineStartIdx());

                        partials.db.segment(eq_idx + (node + 1) * num_states_, num_states_) = -Cdyn_deriv;
                    }
                    eq_idx += data_.num_dynamics_constraints;
                    dyn_partials_timer.StopTimer();
//                    dyn_partials_timer.PrintElapsedTime();
                } else if (data_.constraints_.at(i) == JointForwardKinematics) {
                    throw std::runtime_error("Joint Forward Kinematics not implemented with Carabel yet.");

                } else if (data_.constraints_.at(i) == EndEffectorLocation) {

                    // ----------- Inequality contribution ----------- //
                    const int spline_offset = GetPosSplineStartIdx();
                    int idx = 2*ee; //* (data_.num_ee_location_constraints_ / (2*num_ee_));
                    for (int node = EE_NODE_START; node < info_.num_nodes + 1; node++) {
                        const double time = GetTime(node);
                        for (int coord = 0; coord < 2; coord++) {
                            int vars_index, vars_affecting;
                            std::tie(vars_index, vars_affecting) = traj.GetPositionSplineIndex(ee, time, coord);

                            // position coef dependence
                            const vector_t pos_coef_partials = traj.GetPositionCoefPartialsWrtContactTime(ee, coord,
                                                                                                          time,
                                                                                                          contact_time_idx);

                            G_builder.SetMatrix(pos_coef_partials.transpose(),
                                                ineq_idx + idx,
                                                spline_offset + vars_index);

                            // Need to account for the derivative in the negative term
                            G_builder.SetMatrix(-pos_coef_partials.transpose(),
                                                ineq_idx + idx + data_.num_ee_location_constraints_/2,
                                                spline_offset + vars_index);
                            idx++;
                        }
                        idx += 2*(num_ee_-1);
                    }

//                    assert(idx  == data_.num_ee_location_constraints_);

                } else if (data_.constraints_.at(i) == EndEffectorStart) {
                    // ----------- Equality contribution ----------- //
                    matrix_t M(data_.num_start_ee_constraints_, traj.GetTotalPosSplineVars()); // TODO: DMA
                    M.setZero();
                    int idx = 0;
                    for (int coord = 0; coord < 2; coord++) {
                        int vars_idx, vars_affecting;
                        std::tie(vars_idx, vars_affecting) =
                                prev_traj_.GetPositionSplineIndex(ee, GetTime(0), coord);

                        // Position coef dependence
                        vector_t pos_coef_partials = traj.GetPositionCoefPartialsWrtContactTime(ee, coord,
                                                                                                GetTime(0),
                                                                                                contact_time_idx); // TODO: DMA

                        M.block(idx, vars_idx, 1, vars_affecting) =
                                pos_coef_partials.transpose();
                        idx++;
                    }

                    assert(idx == data_.num_start_ee_constraints_/num_ee_);

                    A_builder.SetMatrix(M, eq_idx, GetPosSplineStartIdx());
                    eq_idx += data_.num_start_ee_constraints_;
                    ineq_idx += data_.num_ee_location_constraints_;
                } else if (data_.constraints_.at(i) == ForceBox) {
                    AddForceBoxConstraintPartials(G_builder, contact_time_idx, ineq_idx, ee);
                    ineq_idx += data_.num_force_box_constraints_;
                } else if (data_.constraints_.at(i) == JointBox) {
                    throw std::runtime_error("Joint box not implemented with Clarabel yet.");
                } else if (data_.constraints_.at(i) == FrictionCone) {
                    AddFrictionConeConstraintPartials(G_builder, contact_time_idx, ineq_idx, ee);
                    ineq_idx += data_.num_cone_constraints_;
                } else if (data_.constraints_.at(i) == TDPosition) {
                    AddTDPositionConstraintPartial(A_builder, partials.db, contact_time_idx, eq_idx, ee);
                    eq_idx += data_.num_td_pos_constraints_;
                } else if (data_.constraints_.at(i) == Raibert) {
                    AddRaibertPartials(A_builder, contact_time_idx, eq_idx, ee);
                    eq_idx += data_.num_raibert_constraints_;
                }
            }

            // TODO: Investigate double free

            QPPartials partials2;
            partials2.dG.resize(data_.num_inequality_, data_.num_decision_vars);

            partials2.dG.setFromTriplets(G_builder.GetTriplet().begin(), G_builder.GetTriplet().end());
            partials.dA.setFromTriplets(A_builder.GetTriplet().begin(), A_builder.GetTriplet().end());

            partials.dG = partials2.dG;

            return true;
        } else {
            return false;
        }
    }

    double MPCSingleRigidBody::GetCost() const {
        return GetCostValue(prev_qp_sol);
    }

    double MPCSingleRigidBody::GetModifiedCost(int num_nodes) const {
        vector_t temp = prev_qp_sol;
//        temp.tail(temp.size() - num_states_*num_nodes).setZero();
        return GetCostValue(temp);
    }

    MPCSingleRigidBody::MPCSingleRigidBody(const mpc::MPCSingleRigidBody& other) : MPC(other) {
        *this = other;
        InitalizeQPData();
    }

    void MPCSingleRigidBody::AddEEPosCost() {
        const double weight = .1;

        const int pos_start_idx = GetPosSplineStartIdx();
        constexpr int CONSTRAINT_COORDS = 2;

        for (int node = EE_NODE_START; node < info_.num_nodes + 1; node++) {
            for (int ee = 0; ee < num_ee_; ee++) {
                model_.GetCOMToHip(ee); // need this to populate the vector

                data_.cost_mat_.SetDiagonalMatrix(weight, node*num_states_, node*num_states_, CONSTRAINT_COORDS);

                vector_2t k = model_.GetCOMHipOffset(ee); // + info_.ee_box_size/2;


                vector_2t q_c;
                q_c = (weight*k); //weight*prev_traj_.GetState(node).head<CONSTRAINT_COORDS>() -
                data_.cost_linear.segment<2>(node*num_states_) += q_c;

                for (int coord = 0; coord < CONSTRAINT_COORDS; coord++) {
                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) =
                            prev_traj_.GetPositionSplineIndex(ee, GetTime(node), coord);

                    vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
                                                                ee, coord, GetTime(node));

                    data_.cost_linear.segment(pos_start_idx + vars_idx, vars_affecting) += -q_c(coord)*vars_lin;

                    data_.cost_mat_.SetMatrix(weight*vars_lin, node * num_states_, pos_start_idx + vars_idx);
//                    data_.cost_mat_.SetDiagonalMatrix(weight, pos_start_idx + vars_idx, node * num_states_, vars_affecting);
                    data_.cost_mat_.SetMatrix(weight*vars_lin, pos_start_idx + vars_idx, node * num_states_);

                    data_.cost_mat_.SetMatrix(weight*vars_lin, pos_start_idx + vars_idx, pos_start_idx + vars_idx);
                }
            }
        }
    }

    // Note that if we have a late touch down then due to the start ee location constraint, the desired stance position will change to the current end effector location
    void MPCSingleRigidBody::AddTDPositionConstraints() {
        const int start_pos_idx = GetPosSplineStartIdx();
        int row_idx = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            if (prev_traj_.GetNextContactTime(ee, init_time_) - init_time_ < td_fraction_*prev_traj_.GetCurrentSwingTime(ee)) {
                const double td_time = prev_traj_.GetNextContactTime(ee, init_time_);
                data_.td_pos_constants_.segment<2>(row_idx) = prev_traj_.GetEndEffectorLocation(ee, td_time).head<2>();
                for (int coord = 0; coord < 2; coord++) {
                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) =
                            prev_traj_.GetPositionSplineIndex(ee, td_time, coord);

                    vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
                                                                ee, coord, td_time);
                    data_.constraint_mat_.SetMatrix(vars_lin.transpose(), constraint_idx_ + row_idx,
                                                    start_pos_idx + vars_idx);
                    row_idx++;
                }
            }
//            else if (prev_traj_.GetFirstTDTime(ee) < init_time_) { // there is a touch down in past time
//                const double td_time = prev_traj_.GetFirstTDTime(ee);
//                std::cout << "First TD for ee " << ee << " is at: " << td_time << std::endl;
//                data_.td_pos_constants_.segment<2>(row_idx) = prev_traj_.GetEndEffectorLocation(ee, td_time).head<2>();
//                for (int coord = 0; coord < 2; coord++) {
//                    int vars_idx, vars_affecting;
//                    std::tie(vars_idx, vars_affecting) =
//                            prev_traj_.GetPositionSplineIndex(ee, td_time, coord);
//
//                    vector_t vars_lin = prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
//                                                                ee, coord, td_time);
//                    data_.constraint_mat_.SetMatrix(vars_lin.transpose(), constraint_idx_ + row_idx,
//                                                    start_pos_idx + vars_idx);
//                    row_idx++;
//                }
//            }
        }
        assert(row_idx == data_.num_td_pos_constraints_);
        constraint_idx_ += row_idx;
    }

    void MPCSingleRigidBody::AddTDPositionConstraintPartial(utils::SparseMatrixBuilder& builder, vector_t& b, int contact_idx,
                                                            int eq_idx, int ee) {
        const int start_pos_idx = GetPosSplineStartIdx();
        int row_idx = 0; // need to start at the correct row
        for (int i = 0; i < ee; i++) {
            if (prev_traj_.GetNextContactTime(i, init_time_) - init_time_ < td_fraction_*prev_traj_.GetCurrentSwingTime(i)) {
                row_idx+=2;
            }
        }

        if (prev_traj_.GetNextContactTime(ee, init_time_) - init_time_ < prev_traj_.GetCurrentSwingTime(ee)/2) {
            const double td_time = prev_traj_.GetNextContactTime(ee, init_time_);
            b.segment<2>( eq_idx +row_idx) = prev_traj_.GetPositionPartialWrtContactTime(ee, td_time, contact_idx).head<2>();
            for (int coord = 0; coord < 2; coord++) {
                int vars_idx, vars_affecting;
                std::tie(vars_idx, vars_affecting) =
                        prev_traj_.GetPositionSplineIndex(ee, td_time, coord);

                vector_t vars_lin = prev_traj_.GetPositionCoefPartialsWrtContactTime(ee, coord, td_time, contact_idx);
                builder.SetMatrix(vars_lin.transpose(), eq_idx + row_idx,
                                                start_pos_idx + vars_idx);
                row_idx++;
            }
        }
//        else if (prev_traj_.GetFirstTDTime(ee) < init_time_) { // there is a touch down in past time
//            const double td_time = prev_traj_.GetFirstTDTime(ee);
//            b.segment<2>( eq_idx +row_idx) = prev_traj_.GetPositionPartialWrtContactTime(ee, td_time, contact_idx).head<2>();
//            for (int coord = 0; coord < 2; coord++) {
//                int vars_idx, vars_affecting;
//                std::tie(vars_idx, vars_affecting) =
//                        prev_traj_.GetPositionSplineIndex(ee, td_time, coord);
//
//                vector_t vars_lin = prev_traj_.GetPositionCoefPartialsWrtContactTime(ee, coord, td_time, contact_idx);
//                builder.SetMatrix(vars_lin.transpose(), eq_idx + row_idx,
//                                  start_pos_idx + vars_idx);
//                row_idx++;
//            }
//        }
    }

    void MPCSingleRigidBody::IncreaseEEBox() {
        info_.ee_box_size(0) += 0.05;
        info_.ee_box_size(1) += 0.05;
    }

    void MPCSingleRigidBody::DecreaseEEBox() {
        info_.ee_box_size(0) = std::max(info_.ee_box_size(0) - 0.05, ee_bounds_(0));
        info_.ee_box_size(1) = std::max(info_.ee_box_size(1) - 0.05, ee_bounds_(1));
    }

    void MPCSingleRigidBody::AddRaibertHeuristic() {
        // For each constant segment assign that foot position as the raibert heuristic (with no feedback term)
        // The raibert heuristic can be used for the flight phase, but those dynamics don't enter the MPC (the foot flight phase dynamics)
        // So we take the velocity measurements of the node closest to the beginning of the stance phase and assign that value
        // r_des = r_ref + v_com*deltaT/2
        // r_ref = COM + offset
        // Note that in this scenario the effect of the ground contact times is directly in the equation and thus no spline parameterization is necessary

        int row_idx = 0;
        int pos_start_idx = GetPosSplineStartIdx();

        const auto contact_times = prev_traj_.GetContactTimes();

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int j = 0; j < contact_times.at(ee).size(); j++) {
                if (contact_times.at(ee).at(j).GetType() == TouchDown
                    && contact_times.at(ee).at(j).GetTime() > init_time_ && contact_times.at(ee).at(j).GetTime() < info_.num_nodes*info_.integrator_dt + init_time_) {
                    for (int coord = 0; coord < 2; coord++) { // Only need the heuristic in the x-y plane
                        const double td_time = contact_times.at(ee).at(j).GetTime(); //prev_traj_.GetNextContactTime(ee, contact_times.at(ee).at(j).GetTime());

//                        std::cout << "touch down time: "<< td_time << std::endl;

                        int vars_affecting, vars_idx;
                        std::tie(vars_idx, vars_affecting) =
                                prev_traj_.GetPositionSplineIndex(ee, td_time, coord);

                        assert(vars_affecting == 1);

                        vector_t vars_lin =  prev_traj_.GetSplineLin(Trajectory::SplineTypes::Position,
                                                                     ee, coord, td_time);

                        // Foot position
                        data_.constraint_mat_.SetMatrix(vars_lin.transpose(), constraint_idx_ + row_idx,
                                                        pos_start_idx + vars_idx);

//                                                    A(idx, (node) * num_states_ + coord) = 1;
                        const int node = prev_traj_.GetNode(td_time - 0.01);

//                        if (row_idx == 0) {
//                            std::cout << "node: " << node << std::endl;
//                        }

                        assert(node < info_.num_nodes + 1);

                        // COM in that coord
                        data_.constraint_mat_.SetDiagonalMatrix(
                                -1.0, constraint_idx_ + row_idx, node*num_states_ + coord, 1);

                        int i = j + 1;
                        while (i < contact_times.at(ee).size() && contact_times.at(ee).at(i).GetType() != LiftOff) {
                            i++;
                        }

                        double contact_time;
                        if (i == contact_times.at(ee).size()) {
                                contact_time = 1.0; // Setting to a reasonable value
                        } else {
                            contact_time = contact_times.at(ee).at(i).GetTime() -
                                                        contact_times.at(ee).at(j).GetTime();
                        }

                        // Velocity COM
                        data_.constraint_mat_.SetDiagonalMatrix(
                                -0.00*contact_time/(model_.GetMass()*2), constraint_idx_ + row_idx, node*num_states_ + 3 + coord, 1);

                        // COM to hip offset
//                        std::cout << model_.GetCOMToHip(ee)(coord) << std::endl;
                        data_.raibert_constants_(row_idx) = model_.GetCOMToHip(ee)(coord);

                        row_idx++;
                    }
                }
            }
        }

        constraint_idx_ += row_idx;

        assert(row_idx == data_.num_raibert_constraints_);
    }

    void MPCSingleRigidBody::AddRaibertPartials(utils::SparseMatrixBuilder& builder, int contact_idx,
                                                int eq_idx, int ee) {
        int row_idx = 0;
        int pos_start_idx = GetPosSplineStartIdx();

        const auto contact_times = prev_traj_.GetContactTimes();

        if (contact_times.at(ee).at(contact_idx).GetType() == TouchDown
            && contact_times.at(ee).at(contact_idx).GetTime() > init_time_
            && contact_times.at(ee).at(contact_idx).GetTime() < info_.num_nodes*info_.integrator_dt + init_time_) {
            for (int coord = 0; coord < 2; coord++) { // Only need the heuristic in the x-y plane
                const double td_time = contact_times.at(ee).at(contact_idx).GetTime(); //prev_traj_.GetNextContactTime(ee, contact_times.at(ee).at(j).GetTime());


                const int node = prev_traj_.GetNode(td_time);

               assert(node < info_.num_nodes + 1);

                int i = contact_idx + 1;
                while (i < contact_times.at(ee).size() && contact_times.at(ee).at(i).GetType() != LiftOff) {
                    i++;
                }

                double contact_time;
                if (i == contact_times.at(ee).size()) {
                    contact_time = 1.0; // Setting to a reasonable value
                } else {
                    contact_time = contact_times.at(ee).at(i).GetTime() -
                                   contact_times.at(ee).at(contact_idx).GetTime();
                }

                // TODO: Fix
                // Velocity COM
                builder.SetDiagonalMatrix(
                        -0.00*contact_time/(model_.GetMass()*2), eq_idx + row_idx, node*num_states_ + 3 + coord, 1);

                row_idx++;
            }
        }
    }

    SingleRigidBodyModel MPCSingleRigidBody::GetModelCopy() const {
        return model_;
    }

} // mpc