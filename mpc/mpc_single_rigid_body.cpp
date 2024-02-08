//
// Created by zolkin on 1/19/24.
//

#include "mpc_single_rigid_body.h"

namespace mpc {

    /*
     * TODO:
     * - Determine why the QP is easier to solve when we have 24 position variables and 96 force variables
     *      this is the pattern of full solves vs inaccurate solves. We can also see jumps in the trajectory
     *      at certain places. I suspect this is realted to this problem.
     */

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

//        data_.cost_mat_.SetDiagonalMatrix(1, 0, 0, data_.num_decision_vars);

        // -------------------- Constraints ---------------------- //
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
//        if (constraint_projection_) {
//            data_.ApplyProjection();
//        }

        // ----------------------- Solve ------------------------- //
        qp_solver->SetupQP(data_, prev_qp_sol);
        utils::Timer qp_solve_timer("QP solve");
        qp_solve_timer.StartTimer();
        // TODO: DMA
        const vector_t sol = qp_solver->Solve(data_);
        qp_solve_timer.StopTimer();

        if (qp_solver->GetSolveQuality() != SolvedInacc && qp_solver->GetSolveQuality() != Solved
            && qp_solver->GetSolveQuality() != MaxIter) {
            std::cerr << "Warning: " << qp_solver->GetSolveQualityAsString() << std::endl;
            throw std::runtime_error("Bad solve.");
        }

        if (info_.verbose == All) {
            std::cout << "Solve type: " << qp_solver->GetSolveQualityAsString() << std::endl;
        }

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
//        std::cout << "force variables start: " << num_states_*(info_.num_nodes+1) << std::endl;
//        std::cout << "position variables start: " << GetPosSplineStartIdx() << std::endl;

        utils::Timer line_search_timer("line search");
        double alpha = 1;
        if (num_run_ >= 0 && sol.size() == prev_qp_sol.size()) {
            line_search_timer.StartTimer();
            alpha = LineSearch(p, state);
            line_search_timer.StopTimer();
        }

//        if (run_num_ > 0){
//            std::cout << (alpha*(qp_solver->GetDualSolution() - prev_dual_sol_) + prev_dual_sol_).lpNorm<Eigen::Infinity>() << std::endl;
//        }

        prev_dual_sol_ = qp_solver->GetDualSolution();

        prev_qp_sol = ((alpha * p) + prev_qp_sol).eval();

        prev_traj_ = ConvertQPSolToTrajectory(prev_qp_sol, state);
//        if (run_num_ == 50) {
//            std::cout << "ee location: " << prev_traj_.GetEndEffectorLocation(1, init_time).transpose() << std::endl;
//            std::cout << "lin vars: " << prev_traj_.GetSplineLin(Trajectory::Position, 1, 0, init_time).transpose() << std::endl;
//            int vars_index, vars_affecting;
//            std::tie(vars_index, vars_affecting) = prev_traj_.GetPositionSplineIndex(1, init_time, 0);
//            std::cout << "vars index, affecting: " << vars_index << ", " << vars_affecting << std::endl;
//            std::cout << "qp solve vars: " << sol.segment(GetPosSplineStartIdx() + vars_index - vars_affecting, vars_affecting).transpose() << std::endl;
//            std::cout << "product: " << prev_traj_.GetSplineLin(Trajectory::Position, 1, 0, init_time).dot(
//                    sol.segment(GetPosSplineStartIdx() + vars_index - vars_affecting, vars_affecting)) << std::endl;
//            std::cout << "ee start constraints: \n" << data_.sparse_constraint_.bottomRows<8>()*sol - data_.start_ee_constants_ << std::endl;
//            std::cout << "constraint mat: " << data_.sparse_constraint_.block(data_.GetTotalNumConstraints() - data_.num_start_ee_constraints_,
//                                                                              GetPosSplineStartIdx(), data_.num_start_ee_constraints_, prev_traj_.GetTotalPosSplineVars());
//            std::cout << "solution val: " << sol(GetPosSplineStartIdx() + vars_index - vars_affecting) << std::endl;
//        }

//        std::cout << "ee location ub constraints: \n" << data_.sparse_constraint_.middleRows(
//                data_.GetTotalNumConstraints() - data_.num_ee_location_constraints_,
//                data_.num_ee_location_constraints_)*sol - data_.ee_location_ub_ << std::endl;
//
//        std::cout << "ee location lb constraints: \n" << -data_.sparse_constraint_.middleRows(
//                data_.GetTotalNumConstraints() - data_.num_ee_location_constraints_,
//                data_.num_ee_location_constraints_)*sol + data_.ee_location_lb_ << std::endl;

//        std::cout << "constraint mat: \n" << data_.sparse_constraint_.middleRows(
//                data_.GetTotalNumConstraints() - (data_.num_ee_location_constraints_ + data_.num_start_ee_constraints_),
//                data_.num_ee_location_constraints_) << std::endl;

//        vector_t part_con = data_.sparse_constraint_.middleRows(
//                data_.GetTotalNumConstraints() - (data_.num_ee_location_constraints_ + data_.num_start_ee_constraints_),
//                data_.num_ee_location_constraints_)*sol;
//        double max_ub = 0, max_lb = 0;
//        int max_idx_ub = -1, max_idx_lb = -1;
//        for (int i = 0; i < data_.num_ee_location_constraints_; i++) {
//            if (part_con(i) - data_.ee_location_ub_(i) >= max_ub) {
//                max_ub = part_con(i);
//                max_idx_ub = i;
//            }
//
//            if (-part_con(i) + data_.ee_location_lb_(i) >= max_lb) {
//                max_lb = part_con(i);
//                max_idx_lb = i;
//            }
//        }
//        std::cout << "max upper bound violation: " << max_ub << ", at index " << max_idx_ub <<
//        "\n max lower bound violation: " << max_lb << ", at index " << max_idx_lb << std::endl;


        // TODO: Turn into unit test
//        vector_t temp = ConvertTrajToQPVec(prev_traj_);
//        assert(temp.size() == prev_qp_sol.size());
//        for (int i = 0; i < temp.size(); i++) {
//            if (std::abs(temp(i) - prev_qp_sol(i)) > 1e-4) {
//                std::cout << "i: " << i << " error: " << std::abs(temp(i) - prev_qp_sol(i)) << std::endl;
//                std::cout << "traj->vec: " << temp(i) << std::endl;
//                std::cout << "vec: " << prev_qp_sol(i) << std::endl;
//                std::cout << "----------" << std::endl;
//            }
//        }

        solve_timer.StopTimer();

        utils::Timer stats_timer("recording stats");
        stats_timer.StartTimer();
        RecordStats(alpha, p, qp_solver->GetSolveQuality(), state,
                    solve_timer.GetElapsedTimeMilliseconds(), GetCostValue(prev_qp_sol)); //GetCostValue(sol));
        stats_timer.StopTimer();

        run_num_++;


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

        prev_traj_.PrintTrajectoryToFile("mpc_demo_traj.txt");

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

        data_.num_cone_constraints_ = (info_.num_nodes+1)*4*num_ee_;
        if (using_clarabel_) {
            data_.num_force_box_constraints_ = 2*GetNodeIntersectMutableForces();
            data_.num_ee_location_constraints_ = 2*(info_.num_nodes)*2*num_ee_; // 2*(info_.num_nodes+1)*2*num_ee_
        } else {
            data_.num_force_box_constraints_ = GetNodeIntersectMutableForces();
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


    // TODO: Investigate more
    void MPCSingleRigidBody::AddEELocationConstraints(const std::vector<vector_3t>& ee_start_locations) {
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
            for (int node = 1; node < info_.num_nodes + 1; node++) {
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

        // -------------- End Effector Start Constraints -------------- //
        idx = 0;
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

    std::vector<Eigen::Vector2d> MPCSingleRigidBody::GetEEBoxCenter() const {
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

            data_.constraint_mat_.SetMatrix(M, constraint_idx_, GetPosSplineStartIdx());

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

            partials.dP.resize(data_.num_decision_vars, data_.num_decision_vars);
            partials.dq.resize(data_.num_decision_vars);
            partials.dA.resize(data_.num_equality_, data_.num_decision_vars);
            partials.dG.resize(data_.num_inequality_, data_.num_decision_vars);
            partials.db.resize(data_.num_equality_);
            partials.dh.resize(data_.num_inequality_);

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
                } else if (data_.constraints_.at(i) == JointForwardKinematics) {
                    throw std::runtime_error("Joint Forward Kinematics not implemented with Carabela yet.");

                } else if (data_.constraints_.at(i) == EndEffectorLocation) {

                    // ----------- Inequality contribution ----------- //
                    const int spline_offset = GetPosSplineStartIdx();
                    int idx = ee * (data_.num_ee_location_constraints_ / (2*num_ee_));
                    for (int node = 0; node < info_.num_nodes + 1; node++) {
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
                                                spline_offset + vars_index - vars_affecting);

                            // Need to account for the derivative in the negative term
                            G_builder.SetMatrix(-pos_coef_partials.transpose(),
                                                ineq_idx + idx + data_.num_ee_location_constraints_/2,
                                                spline_offset + vars_index - vars_affecting);
                            idx++;
                        }
                    }

                    assert(idx - ee * (data_.num_ee_location_constraints_ / (2*num_ee_)) ==
                           data_.num_ee_location_constraints_ / (2*num_ee_));


                    // ----------- Equality contribution ----------- //
                    matrix_t M(data_.num_start_ee_constraints_, traj.GetTotalPosSplineVars()); // TODO: DMA
                    M.setZero();
                    idx = 0;
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

                    assert(idx == data_.num_start_ee_constraints_);

                    A_builder.SetMatrix(M, eq_idx, GetPosSplineStartIdx());

                    partials.dG.setFromTriplets(G_builder.GetTriplet().begin(), G_builder.GetTriplet().end());
                    partials.dA.setFromTriplets(A_builder.GetTriplet().begin(), A_builder.GetTriplet().end());
                } else if (data_.constraints_.at(i) == ForceBox) {
                    ineq_idx += data_.num_force_box_constraints_;
                } else if (data_.constraints_.at(i) == JointBox) {
                    throw std::runtime_error("Joint box not implemented with Carabela yet.");

                } else if (data_.constraints_.at(i) == FrictionCone) {
                    ineq_idx += data_.num_cone_constraints_;
                }
            }

            return true;
        } else {
            return false;
        }
    }

    double MPCSingleRigidBody::GetCost() const {
        return GetCostValue(prev_qp_sol);
    }
} // mpc