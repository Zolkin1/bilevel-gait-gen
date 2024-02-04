//
// Created by zolkin on 1/15/24.
//

#include <iostream>
#include <iomanip>
#include <cassert>

#include "gait_optimizer.h"

namespace mpc {
    // TODO: Things to try:
    // Confirm parameter partials (write unit tests for spline partials, check parameter partials)
    // Check that the dHdth term is correct
    // Sanity check derivatives by adjusting the time and checking the cost change
    // Check LP math
    // Better costs

    GaitOptimizer::GaitOptimizer(int num_ee, int num_contact_nodes, int num_decision_vars, int num_constraints,
                                 double contact_time_ub, double min_time) :
    num_ee_(num_ee), num_decision_vars_(num_decision_vars),
    num_constraints_(num_constraints), contact_time_ub_(contact_time_ub), min_time_(min_time) {
        UpdateSizes(num_decision_vars, num_constraints);
        dHdth.resize(num_ee_*10);   // TODO: Don't hard code

        qp_solver_.settings()->setVerbosity(true);
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setPrimalInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setDualInfeasibilityTolerance(1e-6);
        qp_solver_.settings()->setAbsoluteTolerance(1e-5);
        qp_solver_.settings()->setRelativeTolerance(1e-5);
        qp_solver_.settings()->setScaledTerimination(false);
        qp_solver_.settings()->setMaxIteration(1000);
        qp_solver_.settings()->setRho(.01);
        qp_solver_.settings()->setWarmStart(true);
        qp_solver_.settings()->setScaling(10);
        qp_solver_.settings()->setLinearSystemSolver(1);

    }

    void GaitOptimizer::UpdateSizes(int num_decision_vars, int num_constraints) {
        num_decision_vars_ = num_decision_vars;
        num_constraints_ = num_constraints;

        qp_partials_.dA.resize(num_constraints_, num_decision_vars_);
        qp_partials_.dP.resize(num_decision_vars_, num_decision_vars_);
        qp_partials_.dq.resize(num_decision_vars_);
        qp_partials_.dl.resize(num_constraints_);
        qp_partials_.du.resize(num_constraints_);

        contact_times_.resize(num_ee_);
        contact_times_lb_.resize(num_ee_);
        contact_times_ub_.resize(num_ee_);

        for (int ee = 0; ee < num_ee_; ee++) {
//            std::vector<sp_matrix_t> dA_temp;
//            std::vector<sp_matrix_t> dP_temp;
//            std::vector<vector_t> dl_temp;
//            std::vector<vector_t> du_temp;
//            std::vector<vector_t> dq_temp;
//            for (int idx = 0; idx < num_contact_nodes_; idx++) {
//                dA_temp.emplace_back(num_constraints_, num_decision_vars_);
//                dP_temp.emplace_back(num_decision_vars_, num_decision_vars_);
//
//                dl_temp.emplace_back(num_constraints_);
//                du_temp.emplace_back(num_constraints_);
//                dq_temp.emplace_back(num_decision_vars_);
//            }

            // TODO: figure out a default QPPartials object
            param_partials_.emplace_back();
        }
    }

    void GaitOptimizer::SetCostFcnPartials(const QPPartials& partials) {
        qp_partials_ = partials;
    }

    void GaitOptimizer::Setdx(const vector_t& dx) {
        dldx = dx;
    }

    void GaitOptimizer::ComputeCostFcnDerivWrtContactTimes() {

        sp_matrix_t dldAth(num_constraints_, num_decision_vars_);
        sp_matrix_t dldPth(num_decision_vars_, num_decision_vars_);
        dHdth = vector_t::Zero(GetNumTimeNodes(num_ee_));


        for (int ee = 0; ee < num_ee_; ee++) {
            for (int idx = 0; idx < contact_times_.at(ee).size(); idx++) {
                const QPPartials& param_partial = param_partials_.at(ee).at(idx);
//                assert(!param_partial.dA.toDense().array().isNaN());

//                std::cout << "max param partial dA: " << param_partial.dA.coeffs().maxCoeff() << std::endl;

//                std::cout << "qp partial dA: \n" << qp_partials_.dA.toDense().topLeftCorner<48+24,48+24>() << std::endl;

                dldAth = qp_partials_.dA.cwiseProduct(param_partial.dA);
                dldPth = qp_partials_.dP.cwiseProduct(param_partial.dP);

                matrix_t A = dldAth.toDense();
                matrix_t P = param_partial.dP.toDense();

                for (int i = 0; i < param_partial.dl.size(); i++) {
                    assert(!isnanl(param_partial.dl(i)));
                    assert(!isnanl(param_partial.du(i)));
                    for (int j = 0; j < param_partial.dA.cols(); j++) {
                        if(isnanl(A(i, j))) {
                            std::cout << "qp partials: " << qp_partials_.dA.toDense()(i,j) << std::endl;
                            std::cout << "param partials: " << param_partial.dA.toDense()(i,j) << std::endl;
                            throw;
                        }
                    }
                }

                for (int i = 0; i < param_partial.dq.size(); i++) {
                    assert(!isnanl(param_partial.dq(i)));
                    for (int j = 0; j < param_partial.dP.cols(); j++) {
                        assert(!isnanl(P(i, j)));
                    }
                }

                // TODO: Fix du nans
                dHdth(GetNumTimeNodes(ee) + idx) = dldPth.sum() + dldAth.sum() +
                                                     qp_partials_.dl.dot(param_partial.dl) + qp_partials_.du.dot(param_partial.du) +
                                                     qp_partials_.dq.dot(param_partial.dq);
                assert(!isnanl(dHdth(GetNumTimeNodes(ee) + idx)));
            }
        }

//        std::cout << "max qp partial dA: " << qp_partials_.dA.coeffs().maxCoeff() << std::endl;
//        std::cout << "gradient term: \n" << dHdth << std::endl;

//        dHdth.normalize();
    }

    void GaitOptimizer::OptimizeContactTimes() {
        const int num_decision_vars = GetNumTimeNodes(num_ee_);
        const int num_constraints = num_decision_vars + num_decision_vars + num_ee_;

        A_builder_.Reserve(30);

        qp_solver_.data()->setNumberOfVariables(num_decision_vars);
        qp_solver_.data()->setNumberOfConstraints(num_constraints);

        qp_solver_.data()->clearLinearConstraintsMatrix();
        qp_solver_.data()->clearHessianMatrix();
        qp_solver_.clearSolver();

        Eigen::SparseMatrix<double> P(num_decision_vars, num_decision_vars);    // TODO: Try BFGS
        P.setZero();

        lb_.resize(num_constraints);
        lb_.setZero();
        ub_.resize(num_constraints);
        ub_.setZero();

        int next_row = CreatePolytopeConstraint(0);
        next_row = CreateStepBoundConstraint(next_row);
        CreateNextNodeConstraints(next_row);

        Eigen::SparseMatrix<double> A(num_constraints, num_decision_vars);
        A.setFromTriplets(A_builder_.GetTriplet().begin(), A_builder_.GetTriplet().end());

        // TODO: Print these all in the same rows so they can be more easily compared
//        std::cout << "A: \n" << A.toDense() << std::endl;
//        std::cout << "lb: \n" << lb_ << std::endl;
//        std::cout << "ub: \n" << ub_ << std::endl;

        PrintConstraints(A.toDense(), lb_, ub_);

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int i = 1; i < contact_times_.at(ee).size(); i++) {
                assert(contact_times_.at(ee).at(i-1).GetTime() <= contact_times_.at(ee).at(i).GetTime());
            }
        }

        if (!(qp_solver_.data()->setLinearConstraintsMatrix(A) &&
              qp_solver_.data()->setBounds(lb_, ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // Set solver costs
        if (!(qp_solver_.data()->setHessianMatrix(P) &&
              qp_solver_.data()->setGradient(dHdth))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }

        // Re-init
        if (!qp_solver_.initSolver()) { // TODO: Dynamic memory is allocated here
            throw std::runtime_error("Unable to initialize the solver.");
        }

        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            throw std::runtime_error("Could not solve QP.");
        }

        if (qp_solver_.getStatus() != OsqpEigen::Status::Solved) {
            std::cerr << "Could not solve the gait optimization problem." << std::endl;
            std::cerr << "Solve type: " << GetSolveQualityAsString() << std::endl;
            throw std::runtime_error("Bad gait optimization solve");
        }

        vector_t qp_sol = qp_solver_.getSolution();

        for (int i = 0; i < qp_sol.size(); i++) {
            assert(!isnanl(qp_sol(i)));
        }

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int idx = 0; idx < contact_times_.at(ee).size(); idx++) {
                contact_times_.at(ee).at(idx).SetTime(qp_sol(GetNumTimeNodes(ee) + idx));

                if (idx > 0) {
                    if (contact_times_.at(ee).at(idx-1).GetTime() - contact_times_.at(ee).at(idx).GetTime() <= 1e-3
                    && contact_times_.at(ee).at(idx-1).GetTime() - contact_times_.at(ee).at(idx).GetTime() > 0) {
                        contact_times_.at(ee).at(idx) = contact_times_.at(ee).at(idx-1);
                    }
                    assert(contact_times_.at(ee).at(idx-1).GetTime() <= contact_times_.at(ee).at(idx).GetTime());
                }
            }
        }
    }

    std::vector<time_v> GaitOptimizer::GetContactTimes() {
        return contact_times_;
    }

    QPPartials& GaitOptimizer::GetPartials() {
        return qp_partials_;
    }

    QPPartials& GaitOptimizer::GetParameterPartials(int ee, int idx) {
        return param_partials_.at(ee).at(idx);
    }

    void GaitOptimizer::SetNumContactTimes(int ee, int num_times) {
        contact_times_.at(ee).resize(num_times);
        contact_times_lb_.at(ee).resize(num_times);
        contact_times_ub_.at(ee).resize(num_times);

        param_partials_.at(ee).resize(num_times);
    }

    int GaitOptimizer::GetNumTimeNodes(int ee) const {
        int num_nodes = 0;
        for (int i = 0; i < ee; i++) {
            num_nodes += contact_times_.at(i).size();
        }

        return num_nodes;
    }

    void GaitOptimizer::SetContactTimes(std::vector<time_v> contact_times) {
        contact_times_ = contact_times;
    }

    int GaitOptimizer::CreatePolytopeConstraint(int start_row) {
        // Each node is constrained to be between the node before and after it. At the ends there are constant bounds
        int end_row = start_row;
        for (int ee = 0; ee < num_ee_; ee++) {
            const int nodes = contact_times_.at(ee).size();
            matrix_t A = matrix_t::Zero(nodes, nodes);
//            A(0,0) = 1;
            for (int i = 1; i < nodes; i++) {
                A(i-1,i-1) = 1;
                A(i-1,i) = -1;
            }
            A(nodes-1,nodes-1) = 1;

            vector_t lb = vector_t::Constant(nodes, -2);
//            lb(0) = contact_times_.at(ee).at(0);
            lb(nodes-1) = 0.5 + contact_times_.at(ee).at(0).GetTime();   // TODO: Make this not hard coded and time horizon - might not need this
            lb_.segment(start_row + GetNumTimeNodes(ee), nodes) = lb;
            vector_t ub = vector_t::Zero(nodes);
//            ub(0) = 0.2 + contact_times_.at(ee).at(contact_times_.at(ee).size()-1);
            // TODO: May want to make this a bound from the inital rather than from the final. when its from the final then it can run off to infinity over mutliple solves
            ub(nodes-1) = 0.2 + contact_times_.at(ee).at(contact_times_.at(ee).size()-1).GetTime();
            ub_.segment(start_row + GetNumTimeNodes(ee), nodes) = ub;

            if (ee > 0) {
                A_builder_.SetMatrix(A, start_row + GetNumTimeNodes(ee), GetNumTimeNodes(ee));
            } else {
                A_builder_.SetMatrix(A, start_row + GetNumTimeNodes(ee), GetNumTimeNodes(ee));
            }
            end_row += nodes;
        }
        return end_row;
    }

    int GaitOptimizer::CreateStepBoundConstraint(int start_row) {
        const double bound = 0.05; // TODO: Make not hard coded

        // Note: The first contact time node can never be changed.
        int idx = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int node = 0; node < contact_times_.at(ee).size(); node++) {
                if (node == 0) {
                    // TODO: This constraint ensures the first contact time never moves.
                    // TODO: The problem is that there will always be a node at 0 then.
                    ub_(start_row + idx) = contact_times_.at(ee).at(node).GetTime();
                    lb_(start_row + idx) = contact_times_.at(ee).at(node).GetTime();
                } else {
                    ub_(start_row + idx) = bound + contact_times_.at(ee).at(node).GetTime();
                    lb_(start_row + idx) = -bound + contact_times_.at(ee).at(node).GetTime();
                }
                idx++;
            }
        }
        A_builder_.SetDiagonalMatrix(1, start_row, 0, idx);
        return start_row + idx;
    }

    int GaitOptimizer::CreateNextNodeConstraints(int start_row) {
        int idx = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            if (contact_times_.at(ee).at(1).GetType() == TouchDown) {
                ub_(start_row + idx) = contact_times_.at(ee).at(1).GetTime();
                lb_(start_row + idx) = contact_times_.at(ee).at(1).GetTime();
                A_builder_.SetDiagonalMatrix(1, start_row + idx, GetNumTimeNodes(ee) + 1, 1);
                idx++;
            }
        }

        return start_row + idx;
    }

    void GaitOptimizer::ModifyQPPartials(const vector_t& xstar) {
        qp_partials_.dP = qp_partials_.dP + 0.5*xstar*xstar.transpose();
        qp_partials_.dq += xstar;
    }

    void GaitOptimizer::PrintConstraints(const matrix_t& A, const vector_t& lb, const vector_t& ub) {
        using std::setw;
        for (int row = 0; row < A.rows(); row++) {
            std::cout << setw(8) << std::setprecision(4) << lb(row) << " | ";
            for (int col = 0; col < A.cols(); col++) {
                std::cout << setw(3) << A(row, col) << " ";
            }
            std::cout << "| " << setw(8) << std::setprecision(4) << ub(row) << std::endl;
        }
    }

    std::string GaitOptimizer::GetSolveQualityAsString() const {
        switch (qp_solver_.getStatus()) {
            case OsqpEigen::Status::Solved:
                return "Solved";
            case OsqpEigen::Status::SolvedInaccurate:
                return "Solved Inaccurate";
            case OsqpEigen::Status::PrimalInfeasible:
                return "Primal Infeasible";
            case OsqpEigen::Status::DualInfeasible:
                return "Dual Infeasible";
            case OsqpEigen::Status::PrimalInfeasibleInaccurate:
                return "Primal Infeasible Inaccurate";
            case OsqpEigen::Status::DualInfeasibleInaccurate:
                return "Dual Infeasible Inaccurate";
            case OsqpEigen::Status::MaxIterReached:
                return "Max Iter Reached";
            case OsqpEigen::Status::Unsolved:
                return "Unsolved";
            default:
                return "Other";
        }
    }
}