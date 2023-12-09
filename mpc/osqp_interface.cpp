//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//
#include "osqp_interface.h"

// TODO: Consider pre-conditioning

namespace mpc {
    OSQPInterface::OSQPInterface(QPData data, bool verbose) : QPInterface(data.num_decision_vars), verbose_(verbose){

        // Set solver settings
        qp_solver_.settings()->setVerbosity(verbose_);
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setPrimalInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setDualInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setAbsoluteTolerance(1e-5);
        qp_solver_.settings()->setRelativeTolerance(1e-5);
        qp_solver_.settings()->setScaledTerimination(false);
        qp_solver_.settings()->setMaxIteration(5000);
//        qp_solver_.settings()->setTimeLimit(2e-3); -- Can't do this unless I somehow recompile osqp-eigen with PROFILING=1
        qp_solver_.settings()->setRho(.01);
//        qp_solver_.settings()->setAlpha(1.6);
        qp_solver_.settings()->setWarmStart(true);     // TODO: figure out warm starting with changing sizes
        qp_solver_.settings()->setScaling(10);

        prev_dual_sol_ = vector_t::Zero(data.GetTotalNumConstraints());

        run = 0;
    }

    void OSQPInterface::SetupQP(const mpc::QPData &data) {
        qp_solver_.data()->setNumberOfVariables(data.num_decision_vars);
        qp_solver_.data()->setNumberOfConstraints(data.GetTotalNumConstraints());

        qp_solver_.data()->clearLinearConstraintsMatrix();
        qp_solver_.data()->clearHessianMatrix();
        qp_solver_.clearSolver();

        ConvertDataToOSQPConstraints(data);
        ConvertDataToOSQPCost(data);

        // TODO: Remove
        for (int i = 0; i < A_.rows(); i++) {
            for (int j = 0; j < A_.cols(); j++) {
                assert(!std::isnan(A_(i,j)));
            }
        }

        for (int i = 0; i < P_.rows(); i++) {
            for (int j = 0; j < P_.cols(); j++) {
                assert(!std::isnan(P_(i,j)));
            }
        }

        Eigen::SparseMatrix<double> sparseA = A_.sparseView();
        Eigen::SparseMatrix<double> sparseP = P_.sparseView();

        if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
              qp_solver_.data()->setBounds(lb_, ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // Set solver costs
        if (!(qp_solver_.data()->setHessianMatrix(sparseP) && qp_solver_.data()->setGradient(w_))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }

        // Re-init
        if (!qp_solver_.initSolver()) {
            throw std::runtime_error("Unable to initialize the solver.");
        }
        qp_solver_.setWarmStart(prev_qp_sol_, prev_dual_sol_);
    }

    // TODO: remove data after debugging
    vector_t OSQPInterface::Solve(const QPData& data) {

        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            throw std::runtime_error("Could not solve QP.");
        }

        vector_t qp_sol = qp_solver_.getSolution();

        // TODO: Remove
        for (int i = 0; i < qp_sol.size(); i++) {
            assert(!std::isnan(qp_sol(i)));
        }

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
            std::cout << "num equality constraints: " << data.num_equality_constraints << std::endl;
            std::cout << "num inequality constraints: " << data.num_inequality_constraints << std::endl;
            std::cout << "num decision variables: " << data.num_decision_vars << std::endl;
            std::cout << std::endl;
        }

//        std::cout << "diff on uper bound: \n" << A_*qp_sol - ub_ << std::endl;

        run++;

//        qp_sol.head(data.num_dynamics_constraints) = qp_sol.head(data.num_dynamics_constraints)/10;
        return qp_sol;
    }

    void OSQPInterface::ConvertDataToOSQPConstraints(const mpc::QPData& data) {
        A_ = matrix_t::Zero(data.GetTotalNumConstraints(), data.num_decision_vars);

        A_ << data.dynamics_constraints, data.positive_force_constraints_, data.fk_constraints_,
                data.swing_force_constraints_, data.foot_on_ground_constraints_, data.friction_cone_constraints_,
                data.foot_ground_inter_constraints_, data.box_constraints_, data.force_box_constraints_;

        lb_ = vector_t::Zero(data.GetTotalNumConstraints());
        ub_ = lb_;

        lb_ << data.dynamics_constants, data.positive_force_lb_, data.fk_constants_, data.swing_force_constants_,
                data.foot_on_ground_lb_, data.friction_cone_lb_,
                data.foot_ground_inter_lb_, data.box_lb_, data.force_box_lb_;
        ub_ << data.dynamics_constants, data.positive_force_ub_, data.fk_constants_, data.swing_force_constants_,
                data.foot_on_ground_ub_, data.friction_cone_ub_,
                data.foot_ground_inter_ub_, data.box_ub_, data.force_box_ub_;
    }

    void OSQPInterface::ConvertDataToOSQPCost(const mpc::QPData& data) {
        P_ = matrix_t::Zero(data.num_decision_vars, data.num_decision_vars);
        P_ << data.cost_quadratic;


        w_ = vector_t::Zero(data.num_decision_vars);
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
            default:
                return "Other";
                break;
        } 
    }

    vector_t OSQPInterface::GetDualSolution() const {
        return prev_dual_sol_;
    }

} // mpc