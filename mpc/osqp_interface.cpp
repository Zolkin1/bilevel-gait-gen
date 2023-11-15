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
        qp_solver_.data()->clearLinearConstraintsMatrix();
        qp_solver_.data()->clearHessianMatrix();
        qp_solver_.clearSolver();

        matrix_t A = ConvertDataToOSQPConstraintMat(data);
        vector_t lb = ConvertDataToOSQPlb(data);
        vector_t ub = ConvertDataToOSQPub(data);

        matrix_t P = ConvertDataToOSQPCostMat(data);
        vector_t w = ConvertDataToOSQPCostVec(data);

        // Set solver constraints
        Eigen::SparseMatrix<double> sparseA = A.sparseView();
        if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
              qp_solver_.data()->setBounds(lb, ub))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // Set solver costs
        Eigen::SparseMatrix<double> sparseP = P.sparseView();
        if (!(qp_solver_.data()->setHessianMatrix(sparseP) && qp_solver_.data()->setGradient(w))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }

        // Re-init
        if (!qp_solver_.initSolver()) {
            throw std::runtime_error("Unable to initialize the solver.");
        }
    }

    Trajectory OSQPInterface::Solve() {
        vector_t qp_sol = qp_solver_.getSolution();

        std::cout << "qp sol: \n" << qp_sol << std::endl;

        // TODO: Parse into a trajectory
    }

} // mpc