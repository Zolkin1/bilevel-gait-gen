//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//
#include "osqp_interface.h"

namespace mpc {
    OSQPInterface::OSQPInterface(QPData data, bool verbose)
        : sparse_conversion_timer_("sparse conversion"),
          matrix_gather_timer_("matrix formation"),
          osqp_interface_timer_("osqp setup") {

        prev_qp_sol_ = vector_t::Zero(data.num_decision_vars);
        dldx = vector_t::Zero(data.num_decision_vars);
        verbose_ = verbose;

        // Set solver settings
        qp_solver_.settings()->setVerbosity(true);
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setPrimalInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setDualInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setAbsoluteTolerance(1e-4);
        qp_solver_.settings()->setRelativeTolerance(1e-4);
        qp_solver_.settings()->setScaledTerimination(false);
        qp_solver_.settings()->setMaxIteration(1000);
        qp_solver_.settings()->setRho(.01);
//        qp_solver_.settings()->setAlpha(1.6);
        qp_solver_.settings()->setWarmStart(true);
        qp_solver_.settings()->setScaling(10); // TODO: This has a pretty big effect on the timing and the solution quality
        qp_solver_.settings()->setLinearSystemSolver(1);

        prev_dual_sol_ = vector_t::Zero(data.GetTotalNumConstraints());

        run = 0;
        prev_num_constraints_ = 0;
        prev_num_decision_ = 0;
    }

    void OSQPInterface::SetupQP(QPData& data, const vector_t& warm_start) {
        osqp_interface_timer_.StartTimer();

//        ConvertDataToOSQPConstraints(data);
//        ConvertDataToOSQPCost(data);

        // TODO: I will have an "inefficient start" in both branches if the sparsity pattern is NOT preserved.
        // TODO: Figure out a way to have the sparsity pattern preserved
//        if (prev_num_constraints_ != data.GetTotalNumConstraints() || prev_num_decision_ != data.num_decision_vars) {
            qp_solver_.data()->setNumberOfVariables(data.num_decision_vars);
            qp_solver_.data()->setNumberOfConstraints(data.GetTotalNumConstraints());
            prev_num_constraints_ = data.GetTotalNumConstraints();
            prev_num_decision_ = data.num_decision_vars;
            prev_dual_sol_ = vector_t::Zero(data.GetTotalNumConstraints());

            qp_solver_.data()->clearLinearConstraintsMatrix();
            qp_solver_.data()->clearHessianMatrix();
            qp_solver_.clearSolver();

            if (!(qp_solver_.data()->setLinearConstraintsMatrix(data.sparse_constraint_) &&
                  qp_solver_.data()->setBounds(data.lb_, data.ub_))) {
                throw std::runtime_error("Unable to add the constraints to the QP solver.");
            }

            // Set solver costs
            if (!(qp_solver_.data()->setHessianMatrix(data.sparse_cost_) &&
            qp_solver_.data()->setGradient(data.cost_linear))) {
                throw std::runtime_error("Unable to add the costs to the QP solver.");
            }

            // Re-init
            if (!qp_solver_.initSolver()) { // TODO: Dynamic memory is allocated here
                throw std::runtime_error("Unable to initialize the solver.");
            }
            qp_solver_.setWarmStart(warm_start, prev_dual_sol_);
//            std::cerr << "Inefficient start." << std::endl;
//        } else {
//            qp_solver_.updateHessianMatrix(data.sparse_cost_);
//            qp_solver_.updateGradient(data.cost_linear);
//            qp_solver_.updateLinearConstraintsMatrix(data.sparse_constraint_);
//            qp_solver_.updateBounds(data.lb_, data.ub_);
////            qp_solver_.setWarmStart(warm_start, prev_dual_sol_);
//        }

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
//        std::cout << std::endl;
        osqp_interface_timer_.PrintElapsedTime();
//        sparse_conversion_timer_.PrintElapsedTime();
//        matrix_gather_timer_.PrintElapsedTime();
    }

    // TODO: remove data after debugging
    vector_t OSQPInterface::Solve(const QPData& data) {

        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            throw std::runtime_error("Could not solve QP.");
        }

//        std::cout << "A constraint: \n" << data.sparse_constraint_.toDense().topLeftCorner<48+24,48+24>() << std::endl;

        vector_t qp_sol = qp_solver_.getSolution();

        prev_qp_sol_ = qp_sol;

        prev_dual_sol_ = qp_solver_.getDualSolution();

        if (verbose_) {
//            vector_t temp = data.sparse_constraint_ * qp_sol - data.ub_;
//
//            std::cout << "--- Solve Stats ---" << std::endl;
//            std::cout << "max diff: " << temp.maxCoeff() << std::endl;
//            int ind = 0;
//            for (int i = 0; i < temp.size(); i++) {
//                if (temp(i) == temp.maxCoeff()) {
//                    ind = i;
//                }
//            }
//
//            std::cout << "index of max: " << ind << std::endl;
//            std::cout << "num dynamics constraints: " << data.num_dynamics_constraints << std::endl;
//            std::cout << "num decision variables: " << data.num_decision_vars << std::endl;
//            std::cout << std::endl;
        }

        run++;

        return qp_sol;
    }

//    void OSQPInterface::ConvertDataToOSQPConstraints(const mpc::QPData& data) {
//        if (A_.rows() != data.GetTotalNumConstraints() || A_.cols() != data.num_decision_vars) {
//            A_.resize(data.GetTotalNumConstraints(), data.num_decision_vars);
//        }
//
//        A_.setFromTriplets(data.constraint_mat_.GetTriplet().begin(), data.constraint_mat_.GetTriplet().end());
//
////        Eigen::ColPivHouseholderQR<matrix_t> qr_mat(A_.middleRows(data.num_dynamics_constraints + data.num_fk_ineq_constraints_, data.num_fk_constraints_));
////        std::cout << "rank: " << qr_mat.rank() << ", rows: " << qr_mat.rows() << ", cols: " << qr_mat.cols() << std::endl;
//
//
//        lb_.resize(data.GetTotalNumConstraints());
//        ub_.resize(data.GetTotalNumConstraints());
//
//        lb_ << data.dynamics_constants,
//                data.fk_lb_,
//                data.fk_constants_,
//                data.friction_cone_lb_,
//                data.box_lb_,
//                data.force_box_lb_;
//
//        ub_ << data.dynamics_constants,
//                data.fk_ub_,
//                data.fk_constants_,
//                data.friction_cone_ub_,
//                data.box_ub_,
//                data.force_box_ub_;
//    }

//    void OSQPInterface::ConvertDataToOSQPCost(const mpc::QPData& data) {
//        if (P_.rows() != data.num_decision_vars || P_.cols() != data.num_decision_vars) {
//            P_.resize(data.num_decision_vars, data.num_decision_vars);
//        }
//
//        P_.setFromTriplets(data.cost_mat_.GetTriplet().begin(), data.cost_mat_.GetTriplet().end());
//
//        w_.resize(data.num_decision_vars);
//        w_ << data.cost_linear;
//    }

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

    void OSQPInterface::ConfigureForRealTime(double run_time_iters) const {
        qp_solver_.settings()->setPolish(false);
        qp_solver_.settings()->setAbsoluteTolerance(1e-5);
        qp_solver_.settings()->setRelativeTolerance(1e-6);
        qp_solver_.settings()->setMaxIteration(run_time_iters);
    }

    void OSQPInterface::SetupDerivativeCalcs(vector_t& dx, vector_t& dy_l, vector_t& dy_u) {
        // TODO: Deal with return errors
        qp_solver_.computeAdjointDerivative(dx, dy_l, dy_u);
    }

    void OSQPInterface::CalcDerivativeWrtMats(Eigen::SparseMatrix<double>& dP, Eigen::SparseMatrix<double>& dA) {
        // TODO: Deal with return errors
        OsqpEigen::ErrorExitFlag flag = qp_solver_.adjointDerivativeGetMat(dP, dA);
//        assert(!dA.toDense().array().isNaN());

        std::cout << "deriv mat result: ";

        switch (flag) {
            case OsqpEigen::ErrorExitFlag::DataValidationError :
                std::cout << "Data Validation error." << std::endl;
                break;
            case OsqpEigen::ErrorExitFlag::WorkspaceNotInitError :
                std::cout << "Workspace Not Initialized error." << std::endl;
                break;
            case OsqpEigen::ErrorExitFlag::NoError :
                std::cout << "No error." << std::endl;
                break;
            default:
                std::cout << "Other result." << std::endl;
                break;
        }
    }

    void OSQPInterface::CalcDerivativeWrtVecs(vector_t& dq, vector_t dl, vector_t du) {
        // TODO: Deal with return errors
        qp_solver_.adjointDerivativeGetVec(dq, dl, du);
    }

    void OSQPInterface::Computedx(const Eigen::SparseMatrix<double>& P, const vector_t& q,
                                  const vector_t& xstar) {
        dldx.noalias() = P*xstar + q;
//        std::cout << "reconstructed cost: " << 0.5*xstar.transpose()*P*xstar + q.dot(xstar) << std::endl;
    }

    vector_t OSQPInterface::Getdx() const {
        return dldx;
    }

} // mpc