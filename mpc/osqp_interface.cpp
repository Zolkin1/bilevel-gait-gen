//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//
#include "osqp_interface.h"

namespace mpc {
    OSQPInterface::OSQPInterface(QPData data) : QPInterface(data.num_decision_vars) {

        // Set solver settings
        qp_solver_.settings()->setVerbosity(true);
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setPrimalInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setDualInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setAbsoluteTolerance(1e-3);
        qp_solver_.settings()->setRelativeTolerance(1e-3);
        qp_solver_.settings()->setScaledTerimination(false);
        qp_solver_.settings()->setMaxIteration(4000);
//        qp_solver_.settings()->setTimeLimit(2e-3); -- Can't do this unless I somehow recompile osqp-eigen with PROFILING=1
        qp_solver_.settings()->setRho(.01);
        qp_solver_.settings()->setAlpha(1.6);
        qp_solver_.settings()->setWarmStart(true);

        prev_dual_sol_ = vector_t::Zero(data.num_dynamics_constraints + data.num_equality_constraints +
                data.num_inequality_constraints);

        run = 0;
    }

    void OSQPInterface::SetupQP(const mpc::QPData &data) {
        qp_solver_.data()->setNumberOfVariables(data.num_decision_vars);
        qp_solver_.data()->setNumberOfConstraints(data.num_dynamics_constraints + data.num_equality_constraints +
                                                    data.num_inequality_constraints);

        qp_solver_.data()->clearLinearConstraintsMatrix();
        qp_solver_.data()->clearHessianMatrix();
        qp_solver_.clearSolver();

        ConvertDataToOSQPConstraints(data);
        ConvertDataToOSQPCost(data);

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
        std::cout << "warm start result: " << qp_solver_.setWarmStart(prev_qp_sol_, prev_dual_sol_) << std::endl;
    }

    // TODO: remove data after debugging
    vector_t OSQPInterface::Solve(const QPData& data) {

        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            throw std::runtime_error("Could not solve QP.");
        }

        vector_t qp_sol = qp_solver_.getSolution();

        prev_qp_sol_ = qp_sol;
        prev_dual_sol_ = qp_solver_.getDualSolution();

        vector_t temp = A_*qp_sol - ub_;
        std::cout << "diff on upper bound: \n ---------- Dynamics ----------\n" << (temp).head(data.num_dynamics_constraints) <<
                  "\n---------- Equality ----------\n" << (temp).segment(data.num_dynamics_constraints, data.num_equality_constraints) <<
                  "\n---------- Inequality ----------\n" << (temp).tail(data.num_inequality_constraints)<< std::endl;

        std::cout << "max diff: " << temp.maxCoeff() << std::endl;
        int ind = 0;
        for (int i = 0; i < temp.size(); i++) {
            if (temp(i) == temp.maxCoeff()) {
                ind = i;
            }
        }

        std::cout << "index of max: " << ind << std::endl;

        std::cout << "row of constraint mat at max: " << A_.row(ind) << std::endl;

        std::cout << "num dynamics constraints: " << data.num_dynamics_constraints << std::endl;
        std::cout << "num equality constraints: " << data.num_equality_constraints << std::endl;
        std::cout << "num inequality constraints: " << data.num_inequality_constraints << std::endl;

//        std::cout << "diff on uper bound: \n" << A_*qp_sol - ub_ << std::endl;

        run++;

        return qp_sol;
    }

    void OSQPInterface::ConvertDataToOSQPConstraints(const mpc::QPData& data) {
        A_ = matrix_t::Zero(data.num_dynamics_constraints + data.num_equality_constraints + data.num_inequality_constraints,
                            data.num_decision_vars);
        A_ << data.dynamics_constraints, data.equality_constraints, data.inequality_constraints;

        lb_ = vector_t::Zero(data.num_dynamics_constraints + data.num_equality_constraints + data.num_inequality_constraints);
        ub_ = lb_;

        lb_ << data.dynamics_constants, data.equality_constants, data.inequality_constants_lb;
        ub_ << data.dynamics_constants, data.equality_constants, data.inequality_constants_ub;
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

} // mpc