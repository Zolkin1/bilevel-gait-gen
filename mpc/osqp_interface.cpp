//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//
#include "osqp_interface.h"

// TODO: Consider pre-conditioning

namespace mpc {
    OSQPInterface::OSQPInterface(QPData data, bool verbose)
        : QPInterface(data.num_decision_vars), verbose_(verbose), osqp_interface_timer_("osqp setup"),
          sparse_conversion_timer_("sparse conversion"),
          matrix_gather_timer_("matrix formation"),
          A_(1000,1000){

        // Set solver settings
        qp_solver_.settings()->setVerbosity(false);
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setPrimalInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setDualInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setAbsoluteTolerance(1e-4);
        qp_solver_.settings()->setRelativeTolerance(1e-4);
        qp_solver_.settings()->setScaledTerimination(false);
        qp_solver_.settings()->setMaxIteration(1000);
        qp_solver_.settings()->setRho(.01);
//        qp_solver_.settings()->setAlpha(1.6);
        qp_solver_.settings()->setWarmStart(true);     // TODO: figure out warm starting with changing sizes
        qp_solver_.settings()->setScaling(10);
        qp_solver_.settings()->setLinearSystemSolver(1);

        prev_dual_sol_ = vector_t::Zero(data.GetTotalNumConstraints());

        run = 0;
        prev_num_constraints_ = 0;
        prev_num_decision_ = 0;
    }

    void OSQPInterface::SetupQP(const mpc::QPData &data, const vector_t& warm_start) {
        osqp_interface_timer_.StartTimer();

        ConvertDataToOSQPConstraints(data);
        ConvertDataToOSQPCost(data);

        if (prev_num_constraints_ != data.GetTotalNumConstraints() || prev_num_decision_ != data.num_decision_vars) {
            qp_solver_.data()->setNumberOfVariables(data.num_decision_vars);
            qp_solver_.data()->setNumberOfConstraints(data.GetTotalNumConstraints());
            prev_num_constraints_ = data.GetTotalNumConstraints();
            prev_num_decision_ = data.num_decision_vars;
            prev_dual_sol_ = vector_t::Zero(data.GetTotalNumConstraints());

            qp_solver_.data()->clearLinearConstraintsMatrix();
            qp_solver_.data()->clearHessianMatrix();
            qp_solver_.clearSolver();

            if (!(qp_solver_.data()->setLinearConstraintsMatrix(A_) &&
                  qp_solver_.data()->setBounds(lb_, ub_))) {
                throw std::runtime_error("Unable to add the constraints to the QP solver.");
            }

            // Set solver costs
            if (!(qp_solver_.data()->setHessianMatrix(P_) && qp_solver_.data()->setGradient(w_))) {
                throw std::runtime_error("Unable to add the costs to the QP solver.");
            }

            // Re-init
            if (!qp_solver_.initSolver()) {
                throw std::runtime_error("Unable to initialize the solver.");
            }
            qp_solver_.setWarmStart(warm_start, prev_dual_sol_);
            std::cerr << "Inefficient start." << std::endl;
        } else {
            qp_solver_.updateHessianMatrix(P_);
            qp_solver_.updateGradient(w_);
            qp_solver_.updateLinearConstraintsMatrix(A_);
            qp_solver_.updateBounds(lb_, ub_);
            qp_solver_.setWarmStart(warm_start, prev_dual_sol_);
        }

//        for (int i = 0; i < A_.rows(); i++) {
//            for (int j = 0; j < A_.cols(); j++) {
//                assert(!std::isnan(A_(i,j)));
//                assert(std::abs(A_(i,j)) < 1e6);
//            }
//        }
//
//        for (int i = 0; i < P_.rows(); i++) {
//            for (int j = 0; j < P_.cols(); j++) {
//                assert(!std::isnan(P_(i,j)));
//                assert(std::abs(P_(i,j)) < 1e6);
//            }
//        }
//
//        for (int i = 0; i < lb_.size(); i++) {
//            assert(!std::isnan(lb_(i)));
//
//            assert(!std::isnan(ub_(i)));
//        }

        osqp_interface_timer_.StopTimer();
        std::cout << std::endl;
        osqp_interface_timer_.PrintElapsedTime();
//        sparse_conversion_timer_.PrintElapsedTime();
//        matrix_gather_timer_.PrintElapsedTime();
    }

    // TODO: remove data after debugging
    vector_t OSQPInterface::Solve(const QPData& data) {

        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            throw std::runtime_error("Could not solve QP.");
        }

        vector_t qp_sol = qp_solver_.getSolution();

        prev_qp_sol_ = qp_sol;

        prev_dual_sol_ = qp_solver_.getDualSolution();

        if (verbose_) {
            vector_t temp = A_ * qp_sol - ub_;
//            std::cout << "diff on upper bound: \n ---------- Dynamics ----------\n"
//                      << (temp).head(data.num_dynamics_constraints) <<
//                      "\n---------- Equality ----------\n"
//                      << (temp).segment(data.num_dynamics_constraints, data.num_equality_constraints) <<
//                      "\n---------- Inequality ----------\n" << (temp).tail(data.num_inequality_constraints)
//                      << std::endl;

            std::cout << "--- Solve Stats ---" << std::endl;
            std::cout << "max diff: " << temp.maxCoeff() << std::endl;
            int ind = 0;
            for (int i = 0; i < temp.size(); i++) {
                if (temp(i) == temp.maxCoeff()) {
                    ind = i;
                }
            }

            std::cout << "index of max: " << ind << std::endl;
            std::cout << "num dynamics constraints: " << data.num_dynamics_constraints << std::endl;
            std::cout << "num decision variables: " << data.num_decision_vars << std::endl;
            std::cout << std::endl;
        }

//        std::cout << "diff on uper bound: \n" << A_*qp_sol - ub_ << std::endl;

        run++;

//        qp_sol.head(data.num_dynamics_constraints) = qp_sol.head(data.num_dynamics_constraints)/10;
        return qp_sol;
    }

    void OSQPInterface::ConvertDataToOSQPConstraints(const mpc::QPData& data) {
        if (A_.rows() != data.GetTotalNumConstraints() || A_.cols() != data.num_decision_vars) {
            A_.resize(data.GetTotalNumConstraints(), data.num_decision_vars);
        }

        A_.setFromTriplets(data.constraint_mat_.GetTriplet().begin(), data.constraint_mat_.GetTriplet().end());

        lb_.resize(data.GetTotalNumConstraints());
        ub_.resize(data.GetTotalNumConstraints());

        lb_ << data.dynamics_constants,
                data.fk_lb_,
                data.fk_constants_,
                data.friction_cone_lb_,
                data.box_lb_,
                data.force_box_lb_;

        ub_ << data.dynamics_constants,
                data.fk_ub_,
                data.fk_constants_,
                data.friction_cone_ub_,
                data.box_ub_,
                data.force_box_ub_;
    }

    void OSQPInterface::ConvertDataToOSQPCost(const mpc::QPData& data) {
        if (P_.rows() != data.num_decision_vars || P_.cols() != data.num_decision_vars) {
            P_.resize(data.num_decision_vars, data.num_decision_vars);
        }

        P_.setFromTriplets(data.cost_mat_.GetTriplet().begin(), data.cost_mat_.GetTriplet().end());

        w_.resize(data.num_decision_vars);
        w_ << data.cost_linear;
    }

    vector_t OSQPInterface::GetInfinity(int size) const {
        vector_t infty = vector_t::Constant(size, OsqpEigen::INFTY);

        return infty;
    }

    std::string OSQPInterface::GetSolveQuality() const {
        switch (qp_solver_.getStatus()) {
            case OsqpEigen::Status::Solved:
                return "Solved";
                break;
            case OsqpEigen::Status::SolvedInaccurate:
                return "Solved Inaccurate";
                break;
            case OsqpEigen::Status::PrimalInfeasible:
                return "Primal Infeasible";
                break;
            case OsqpEigen::Status::DualInfeasible:
                return "Dual Infeasible";
                break;
            case OsqpEigen::Status::PrimalInfeasibleInaccurate:
                return "Primal Infeasible Inaccurate";
                break;
            case OsqpEigen::Status::DualInfeasibleInaccurate:
                return "Dual Infeasible Inaccurate";
                break;
            case OsqpEigen::Status::MaxIterReached:
                return "Max Iter Reached";
                break;
            case OsqpEigen::Status::Unsolved:
                return "Unsolved";
                break;
            default:
                return "Other";
                break;
        } 
    }

    vector_t OSQPInterface::GetDualSolution() const {
        return prev_dual_sol_;
    }

    void OSQPInterface::ConfigureForInitialRun() const {
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setAbsoluteTolerance(1e-5);
        qp_solver_.settings()->setRelativeTolerance(1e-5);
        qp_solver_.settings()->setMaxIteration(3000);
    }

    void OSQPInterface::ConfigureForRealTime() const {
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setAbsoluteTolerance(1e-4);
        qp_solver_.settings()->setRelativeTolerance(1e-4);
        qp_solver_.settings()->setMaxIteration(200);
    }

} // mpc