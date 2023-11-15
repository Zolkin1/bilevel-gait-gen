//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "osqp_interface.h"

namespace mpc {
    OSQPInterface::OSQPInterface() : QPInterface() {
        // Set solver settings
        qp_solver_.settings()->setVerbosity(true);
        qp_solver_.settings()->setPolish(true);
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

        // Set solver constraints
        Eigen::SparseMatrix<double> sparseA = A_.sparseView();
        if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
              qp_solver_.data()->setBounds(lb_, ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // Set solver costs
        Eigen::SparseMatrix<double> sparseP = P_.sparseView();
        if (!(qp_solver_.data()->setHessianMatrix(sparseP) && qp_solver_.data()->setGradient(w_))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }

        // Re-init
        if (!qp_solver_.initSolver()) {
            throw std::runtime_error("Unable to initialize the solver.");
        }
    }

    vector_t OSQPInterface::Solve() {

        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            throw std::runtime_error("Could not solve QP.");
        }

        vector_t qp_sol = qp_solver_.getSolution();

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