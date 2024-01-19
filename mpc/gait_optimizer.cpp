//
// Created by zolkin on 1/15/24.
//

#include <iostream>
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
        qp_solver_.settings()->setAbsoluteTolerance(1e-4);
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

                std::cout << "max param partial dA: " << param_partial.dA.coeffs().maxCoeff() << std::endl;

//                std::cout << "qp partial dA: \n" << qp_partials_.dA.toDense().topLeftCorner<48+24,48+24>() << std::endl;

                dldAth = qp_partials_.dA.cwiseProduct(param_partial.dA);
                dldPth = qp_partials_.dP.cwiseProduct(param_partial.dP);
                // TODO: Fix du nans
                dHdth(GetNumTimeNodes(ee) + idx) = dldAth.sum() + dldPth.sum() +
                                                     qp_partials_.dl.dot(param_partial.dl) + qp_partials_.du.dot(param_partial.du) +
                                                     qp_partials_.dq.dot(param_partial.dq);
            }
        }
        std::cout << "max qp partial dA: " << qp_partials_.dA.coeffs().maxCoeff() << std::endl;


//        dHdth.normalize();
    }

    void GaitOptimizer::OptimizeContactTimes() {
        const int num_decision_vars = GetNumTimeNodes(num_ee_);
        const int num_constraints = num_ee_ + num_decision_vars;

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

        CreatePolytopeConstraint();

        Eigen::SparseMatrix<double> A(num_constraints, num_decision_vars);
        A.setFromTriplets(A_builder_.GetTriplet().begin(), A_builder_.GetTriplet().end());

//        std::cout << "A: \n" << A.toDense() << std::endl;

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

        vector_t qp_sol = qp_solver_.getSolution();

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int idx = 0; idx < contact_times_.at(ee).size(); idx++) {
                // TODO: Try normalizing the gradient step vector

                contact_times_.at(ee).at(idx) = qp_sol(GetNumTimeNodes(ee) + idx);

                // TODO: Min time not respected at the end of the time frame (i.e. at the upper bound)
                if (idx == 0) {
                    if (contact_times_.at(ee).at(idx) <= 0) {
                        contact_times_.at(ee).at(idx) = 0;
                    } else if (contact_times_.at(ee).at(idx) >= contact_time_ub_) {
                        contact_times_.at(ee).at(idx) = contact_time_ub_;
                    }
                } else {
                    if (contact_times_.at(ee).at(idx) <= contact_times_.at(ee).at(idx - 1)) {
                        contact_times_.at(ee).at(idx) = contact_times_.at(ee).at(idx-1) + min_time_;
                    } else if (contact_times_.at(ee).at(idx) >= contact_time_ub_) {
                        contact_times_.at(ee).at(idx) = contact_time_ub_;
                    }
                }

                // Project the times back onto the polytope
//                if (contact_times_.at(ee).at(idx) >= contact_times_ub_.at(ee).at(idx)) {
//                    contact_times_.at(ee).at(idx) = contact_times_ub_.at(ee).at(idx);
//                }
//
//                if (contact_times_.at(ee).at(idx) <= contact_times_lb_.at(ee).at(idx)) {
//                    contact_times_.at(ee).at(idx) = contact_times_lb_.at(ee).at(idx);
//                }
            }
        }
    }

    std::vector<std::vector<double>> GaitOptimizer::GetContactTimes() {
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

    void GaitOptimizer::SetContactTimes(std::vector<std::vector<double>> contact_times) {
        contact_times_ = contact_times;
    }

    void GaitOptimizer::CreatePolytopeConstraint() {
        // Each node is constrained to be between the node before and after it. At the ends there are constant bounds
        for (int ee = 0; ee < num_ee_; ee++) {
            const int nodes = contact_times_.at(ee).size();
            matrix_t A = matrix_t::Zero(nodes+1, nodes);
            A(0,0) = 1;
            for (int i = 1; i < nodes; i++) {
                A(i,i-1) = 1;
                A(i,i) = -1;
            }
            A(nodes,nodes-1) = 1;

            lb_.segment(GetNumTimeNodes(ee)+1, nodes+1).setZero();
            vector_t ub = vector_t::Zero(nodes+1);
            ub(0) = 1;
            ub(nodes) = 1;
            ub_.segment(GetNumTimeNodes(ee)+1, nodes+1) = ub;

            if (ee > 0) {
                A_builder_.SetMatrix(A, GetNumTimeNodes(ee)+ee, GetNumTimeNodes(ee));
            } else {
                A_builder_.SetMatrix(A, GetNumTimeNodes(ee), GetNumTimeNodes(ee));
            }
        }
    }

    void GaitOptimizer::ModifyQPPartials(const vector_t& xstar) {
        qp_partials_.dP = qp_partials_.dP + 0.5*xstar*xstar.transpose();
        qp_partials_.dq += xstar;
    }
}