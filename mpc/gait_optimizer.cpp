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
    num_constraints_(num_constraints), contact_time_ub_(contact_time_ub), min_time_(min_time),
    max_trust_region_(0.075) {
        UpdateSizes(num_decision_vars, num_constraints);
        dHdth.resize(num_ee_*10);   // TODO: Don't hard code

        qp_solver_.settings()->setVerbosity(false);
        qp_solver_.settings()->setPolish(true);
        qp_solver_.settings()->setPrimalInfeasibilityTolerance(1e-8);
        qp_solver_.settings()->setDualInfeasibilityTolerance(1e-8);
        qp_solver_.settings()->setAbsoluteTolerance(1e-10);
        qp_solver_.settings()->setRelativeTolerance(1e-10);
        qp_solver_.settings()->setScaledTerimination(false);
        qp_solver_.settings()->setMaxIteration(3000);
        qp_solver_.settings()->setRho(.01);
        qp_solver_.settings()->setWarmStart(true);
        qp_solver_.settings()->setScaling(10);
        qp_solver_.settings()->setLinearSystemSolver(1);

        gamma_ = 0.5;
        eta_ = 0.75;
        Delta_ = 1; // at 1e-8 I match predicted to actual perfectly!

        run_num_ = 0;
        past_decision_vars_ = 0;
    }

    void GaitOptimizer::UpdateSizes(int num_decision_vars, int num_constraints) {
        num_decision_vars_ = num_decision_vars;
        num_constraints_ = num_constraints;

//        qp_partials_.dA.resize(num_constraints_, num_decision_vars_);
//        qp_partials_.dP.resize(num_decision_vars_, num_decision_vars_);
//        qp_partials_.dq.resize(num_decision_vars_);
//        qp_partials_.dl.resize(num_constraints_);
//        qp_partials_.du.resize(num_constraints_);
//
//        contact_times_.resize(num_ee_);
//        contact_times_lb_.resize(num_ee_);
//        contact_times_ub_.resize(num_ee_);

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

        sp_matrix_t dldAth(qp_partials_.dA.rows(), qp_partials_.dA.cols());
        sp_matrix_t dldGth(qp_partials_.dG.rows(), qp_partials_.dG.cols());
        sp_matrix_t dldPth(num_decision_vars_, num_decision_vars_);
        dHdth = vector_t::Zero(GetNumTimeNodes(num_ee_));

        qp_partials_.dl.setZero();
        qp_partials_.du.setZero();

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int idx = 0; idx < contact_times_.at(ee).size(); idx++) {
                QPPartials& param_partial = param_partials_.at(ee).at(idx);
//                assert(!param_partial.dA.toDense().array().isNaN());

                if (param_partial.dA.size() == 0) {
                    param_partial.dA.resize(qp_partials_.dA.rows(), qp_partials_.dA.cols());
                    param_partial.dA.setZero();
                }

                if (param_partial.dG.size() == 0) {
                    param_partial.dG.resize(qp_partials_.dG.rows(), qp_partials_.dG.cols());
                    param_partial.dG.setZero();
                }

//                std::cout << "max param partial dA: " << param_partial.dA.toDense().maxCoeff() << std::endl;
//                std::cout << "max param partial dG: " << param_partial.dG.toDense().maxCoeff() << std::endl;

                //                std::cout << "qp partial dA: \n" << qp_partials_.dA.toDense().topLeftCorner<48+24,48+24>() << std::endl;

                dldAth = qp_partials_.dA.cwiseProduct(param_partial.dA);
                dldGth = qp_partials_.dG.cwiseProduct(param_partial.dG);
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
                            throw std::runtime_error("A partials have a NaN.");
                        }
                    }
                }

                for (int i = 0; i < param_partial.dq.size(); i++) {
                    assert(!isnanl(param_partial.dq(i)));
                    for (int j = 0; j < param_partial.dP.cols(); j++) {
                        assert(!isnanl(P(i, j)));
                    }
                }

                dHdth(GetNumTimeNodes(ee) + idx) = dldGth.sum( )+ dldAth.sum() + dldPth.sum() +
                                                     qp_partials_.dl.dot(param_partial.dl) + qp_partials_.du.dot(param_partial.du) +
                                                     qp_partials_.dq.dot(param_partial.dq) + qp_partials_.db.dot(param_partial.db) +
                                                     qp_partials_.dh.dot(param_partial.dh);
                assert(!isnanl(dHdth(GetNumTimeNodes(ee) + idx)));
            }
        }

//        std::cout << "max qp partial dA: " << qp_partials_.dA.coeffs().maxCoeff() << std::endl;
//        std::cout << "max qp partial dG: " << qp_partials_.dG.coeffs().maxCoeff() << std::endl;
        double max = -1e10;
        int max_row = -1;
        int max_col = -1;
        matrix_t dG = qp_partials_.dG.toDense();
        for (int i = 0; i < qp_partials_.dG.rows(); i++) {
            for (int j = 0; j < qp_partials_.dG.cols(); j++) {
                if (dG(i,j) > max) {
                    max = qp_partials_.dG.toDense()(i,j);
                    max_row = i;
                    max_col = j;
                }
            }
        }

//        std::cout << "dG max coeff row, col: " << max_row << ", " << max_col << std::endl;
        std::cout << "gradient term: " << dHdth.transpose() << std::endl;

//        dHdth.normalize();
    }

    void GaitOptimizer::OptimizeContactTimes(double time, double actual_red_cost) {
        OptimizeContactTimes(time, actual_red_cost, 1, true);
    }

    void GaitOptimizer::OptimizeContactTimes(double time, double actual_red_cost, double alpha, bool adapt_trust_region) {
        // ------------------- Setup ------------------- //
        int num_decision_vars = GetNumTimeNodes(num_ee_);
        int num_constraints = num_decision_vars + num_decision_vars + 3*num_ee_;

        A_builder_.Reserve(30);

        // ------------------- Trust Region Modification and Step Acceptance ------------------- //
        // Note: constraints are linear and thus always fully satisfied after the QP solve
        //       so acceptance and trust region sizes are based only on the cost function
        // Note: this happens here (rather than after the QP solve) because we can't get the new cost
        //       until after an MPC solve. So we update the trust region here.
        if (run_num_ > 0 && adapt_trust_region) {
            double rho = actual_red_cost / pred_red_cost_;
            std::cout << "rho: " << rho << std::endl;
            if (rho > eta_) {
//                IncreaseTrustRegion(rho);
                old_contact_times_ = contact_times_;
            } else {
                assert(old_contact_times_.size() == contact_times_.size());
//                DecreaseTrustRegion(step_);

// TODO: Put trust region changing back

//                contact_times_ = old_contact_times_;
//                num_decision_vars = GetNumTimeNodes(num_ee_);
//                num_constraints = num_decision_vars + num_decision_vars + 3*num_ee_;
            }
        }

        qp_solver_.data()->setNumberOfVariables(num_decision_vars);
        qp_solver_.data()->setNumberOfConstraints(num_constraints);

        qp_solver_.data()->clearLinearConstraintsMatrix();
        qp_solver_.data()->clearHessianMatrix();
        qp_solver_.clearSolver();

        // ------------------- Create Constraints ------------------- //
        lb_.resize(num_constraints);
        lb_.setZero();
        ub_.resize(num_constraints);
        ub_.setZero();

        int next_row = CreatePolytopeConstraint(0);
//        next_row = CreateStepBoundConstraint(next_row);
        next_row = CreateStartConstraint(next_row);
        next_row = CreateTrustRegionConstraint(next_row);
        CreateNextNodeConstraints(next_row, time);

        sp_matrix_t A(num_constraints, num_decision_vars);
        A.setFromTriplets(A_builder_.GetTriplet().begin(), A_builder_.GetTriplet().end());


        // ------------------- BFGS ------------------- //
        // TODO: Deal with changing size of decision variables
        if (dual_.size() != num_constraints) {
            dual_ = vector_t::Zero(num_constraints);
        }
        if (gradkp1_.size() != num_decision_vars) {
            gradkp1_ = vector_t::Zero(num_decision_vars);
        }

        if (xk_.size() != num_decision_vars) {
            xk_ = vector_t::Zero(num_decision_vars);
        }

//        if (run_num_ == 0) {
//            old_grad_ = dHdth;
//        }

        // TODO: Fix
//        if (dHdth.size() != num_decision_vars) {
//            dHdth = old_grad_;  // TODO: Will this stick me in a loop somehow?
//        }

//        UpdateLagrangianGradients(A);
//
////         BFGS does appear to help
//        if (past_decision_vars_ == num_decision_vars) {
//            AdjustBSize(num_decision_vars);
//
//            DampedBFGSUpdate();
//
//            const auto ldlt = Bk_.ldlt();
//            if (ldlt.info() == Eigen::NumericalIssue) {
//                std::cerr << "Bk is not PSD." << std::endl;
//            }
//        } else {
//            Bk_ = matrix_t::Identity(num_decision_vars, num_decision_vars);
//        }


        Bk_ = matrix_t::Zero(num_decision_vars, num_decision_vars);
        std::cout << "No BFGS" << std::endl;

        // TODO: Consider building a sparse matrix in the first place
        sp_matrix_t P = Bk_.sparseView();

        // ------------------- Run Solver ------------------- //
//        PrintConstraints(A.toDense(), lb_, ub_);

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

        if (qp_solver_.getStatus() != OsqpEigen::Status::Solved && qp_solver_.getStatus() != OsqpEigen::Status::SolvedInaccurate
            && qp_solver_.getStatus() != OsqpEigen::Status::MaxIterReached) {
            std::cerr << "Could not solve the gait optimization problem." << std::endl;
            std::cerr << "Solve type: " << GetSolveQualityAsString() << std::endl;
            throw std::runtime_error("Bad gait optimization solve");
        }

        if (qp_solver_.getStatus() == OsqpEigen::Status::MaxIterReached) {
            std::cerr << "Max iterations reached on the gait optimization." << std::endl;
        }

        step_ = qp_solver_.getSolution();
        dual_ = qp_solver_.getDualSolution();

        // ------------------- Numerical Conditioning ------------------- //
        for (int i = 0; i < step_.size(); i++) {
            assert(!isnanl(step_(i)));
        }


        // Note: the step is always accepted since I can't get the actual cost reduction without an MPC solve.
        //      so for now we should just start with a small trust region and expand out because when the region
        //      is small enough we guaruntee to have some reduction.

        xk_ = xkp1_;

        // TODO: Remove after finish debugging
        int max_idx = -1;
        double max = 1e10;
        for (int i = 0; i < dHdth.size(); i++) {
            if (dHdth(i) < max && (i%2) != 0) {
                max = dHdth(i);
                max_idx = i;
            }
        }

//        std::cout << "Max gradient occurs at: " << max_idx << " with value: " << dHdth(max_idx) << std::endl;
//        step_.setZero();
//        double step_len = 1e-3;
//        step_(max_idx) = step_len;

        // Finite differencing!
//        step_.setZero();
//        const double dt = 1e-8;     // TODO: Finite difference surprisingly sensitive to this (can I use it as ground truth?)
//        const int step_idx = 4; // some are worse than others (i.e. 4)
//        step_(step_idx) = dt;

        step_ = alpha*step_; // used for line-search/debugging
        xkp1_ = xk_ + step_;

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int idx = 0; idx < contact_times_.at(ee).size(); idx++) {
                contact_times_.at(ee).at(idx).SetTime(xkp1_(GetNumTimeNodes(ee) + idx));

                if (idx > 0) {
                    if (contact_times_.at(ee).at(idx-1).GetTime() - contact_times_.at(ee).at(idx).GetTime() <= 1e-3
                        && contact_times_.at(ee).at(idx-1).GetTime() - contact_times_.at(ee).at(idx).GetTime() > 0) {
                        contact_times_.at(ee).at(idx) = contact_times_.at(ee).at(idx-1);
                    }
                    assert(contact_times_.at(ee).at(idx-1).GetTime() <= contact_times_.at(ee).at(idx).GetTime());
                }
            }
        }

//        std::cout << "approx grad: ";
//        for (int i = 0; i < last_step_.size(); i++) {
//            std::cout << actual_red_cost/last_step_(i) << ", ";
//        }
//        std::cout << std::endl;

        last_step_ = step_;

        pred_red_cost_ = -dHdth.dot(step_) - 0.5*step_.transpose()*Bk_*step_;

        past_decision_vars_ = num_decision_vars;

//        old_grad_ = dHdth;

//        std::cout << "Actual cost reduction (previous step): " << actual_red_cost << std::endl;
//        if (actual_red_cost < 0) {
//            std::cerr << "Bad cost reduction!" << std::endl;
//        }
//        std::cout << "Predicted cost reduction: " << pred_red_cost_ << std::endl;
//        std::cout << "trust region size: " << Delta_ << std::endl;
//        std::cout << "gradient norm: " << dHdth.norm() << std::endl;
        std::cout << "step norm: " << step_.norm() << std::endl;

//        std::cout << "finite difference: " << -actual_red_cost/dt << std::endl;
//        std::cout << "step idx: " << step_idx << std::endl;
        std::cout << std::endl;

        run_num_++;
    }

    std::vector<time_v>& GaitOptimizer::GetContactTimes() {
        return contact_times_;
    }

    QPPartials& GaitOptimizer::GetQPPartials() {
        return qp_partials_;
    }

    QPPartials& GaitOptimizer::GetParameterPartials(int ee, int idx) {
        return param_partials_.at(ee).at(idx);
    }

    void GaitOptimizer::SetNumContactTimes(int ee, int num_times) {
        contact_times_.at(ee).resize(num_times);
//        contact_times_lb_.at(ee).resize(num_times);
//        contact_times_ub_.at(ee).resize(num_times);

        param_partials_.at(ee).resize(num_times);
    }

    int GaitOptimizer::GetNumTimeNodes(int ee) const {
        int num_nodes = 0;
        for (int i = 0; i < ee; i++) {
            num_nodes += contact_times_.at(i).size();
        }

        return num_nodes;
    }

    void GaitOptimizer::SetContactTimes(const std::vector<time_v>& contact_times) {
        if (old_contact_times_.size() != 0) {
            for (int ee = 0; ee < num_ee_; ee++) {
                for (int node = old_contact_times_.at(ee).size(); node < contact_times.at(ee).size(); node++) {
                    old_contact_times_.at(ee).push_back(contact_times.at(ee).at(node));
                }
            }
        } else {
            old_contact_times_ = contact_times;
        }

        contact_times_ = contact_times;
        xkp1_ = ContactTimesToQPVec();    // TODO: Is this the vector I want to assign to?
    }

    int GaitOptimizer::CreatePolytopeConstraint(int start_row) {
        // Each node is constrained to be between the node before and after it. At the ends there are constant bounds
        double constexpr MIN_TIME = 0.1;

        int end_row = start_row;
        for (int ee = 0; ee < num_ee_; ee++) {
            const int nodes = contact_times_.at(ee).size();
            matrix_t A = matrix_t::Zero(nodes, nodes);
            for (int i = 1; i < nodes; i++) {
                A(i-1,i-1) = 1;
                A(i-1,i) = -1;
                ub_(start_row + GetNumTimeNodes(ee) + i-1) = contact_times_.at(ee).at(i).GetTime()
                        - contact_times_.at(ee).at(i-1).GetTime() - MIN_TIME;
                lb_(start_row + GetNumTimeNodes(ee) + i-1) = -2;
            }
            A(nodes-1,nodes-1) = 1;

            lb_(start_row + GetNumTimeNodes(ee) + nodes - 1) = 0; //-0.2; // TODO: Check this
            ub_(start_row + GetNumTimeNodes(ee) + nodes - 1) = 0.2;

            assert(GetNumTimeNodes(ee) + nodes <= GetNumTimeNodes(num_ee_));

//            PrintConstraints(A, lb_, ub_);

            A_builder_.SetMatrix(A, start_row + GetNumTimeNodes(ee), GetNumTimeNodes(ee));
            end_row += nodes;
        }
        return end_row;
    }

    int GaitOptimizer::CreateStepBoundConstraint(int start_row) {
        // TODO: Make this bound update with each step
        // TODO: Consider changing this infinity norm bound to an L1 norm bound
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

    int GaitOptimizer::CreateStartConstraint(int start_row) {
        for (int ee = 0; ee < num_ee_; ee++) {
            lb_(start_row + ee) = 0; //contact_times_.at(ee).at(0).GetTime();
            ub_(start_row + ee) = 0; //contact_times_.at(ee).at(0).GetTime();
            A_builder_.SetDiagonalMatrix(1, start_row + ee, GetNumTimeNodes(ee), 1);
        }

        return start_row + num_ee_;
    }

    int GaitOptimizer::CreateTrustRegionConstraint(int start_row) {
        const int num_decision_vars = GetNumTimeNodes(num_ee_);
        // Note: this is an infinity norm trust region
        lb_.segment(start_row, num_decision_vars) = vector_t::Constant(num_decision_vars, -Delta_);
        ub_.segment(start_row, num_decision_vars) = vector_t::Constant(num_decision_vars, Delta_);
        A_builder_.SetDiagonalMatrix(1, start_row, 0, num_decision_vars);

        return start_row + num_decision_vars;
    }

    int GaitOptimizer::CreateNextNodeConstraints(int start_row, double time) {

        int idx = 0;
        for (int ee = 0; ee < num_ee_; ee++) {
            int next_node = -1;
            for (int i = 1; i < contact_times_.at(ee).size(); i++) {
                if (contact_times_.at(ee).at(i).GetTime() >= time) {
                    next_node = i;
                    break;
                }
            }

            if (contact_times_.at(ee).at(next_node).GetType() == TouchDown) {
                ub_(start_row + idx + 1) = 0; //contact_times_.at(ee).at(next_node).GetTime();
                lb_(start_row + idx + 1) = 0; //contact_times_.at(ee).at(next_node).GetTime();
                ub_(start_row + idx) = 0; //contact_times_.at(ee).at(next_node-1).GetTime();
                lb_(start_row + idx) = 0; //contact_times_.at(ee).at(next_node-1).GetTime();
                A_builder_.SetDiagonalMatrix(1, start_row + idx, GetNumTimeNodes(ee) + next_node - 1, 2);
                idx+=2;
            }
        }

        return start_row + idx;
    }

    void GaitOptimizer::ModifyQPPartials(const vector_t& xstar) {
        qp_partials_.dP = qp_partials_.dP + 0.5*xstar*xstar.transpose();
        qp_partials_.dq += xstar;
    }

    void GaitOptimizer::DampedBFGSUpdate() {
        // sk = xkp1 - xk
        const vector_t sk = xkp1_ - xk_;

        // yk = GLkp1 - GLk
        const vector_t yk = gradkp1_ - gradk_;

        // thetak = option
        double thetak;
        if (sk.dot(yk) >=  0.2*sk.transpose()*Bk_*sk) {
            thetak = 1;
        } else {
            thetak = static_cast<double>((0.8*sk.transpose()*Bk_*sk))/(sk.transpose()*Bk_*sk - sk.dot(yk));
        }

        // rk = thetak*yk + (1-thetak)*Bk*sk
        const vector_t rk = thetak*yk + (1 - thetak)*Bk_*sk;

        // Bkp1 = Bk - (Bk*sk*sk^T*Bk)/(sk^T*Bk*sk) + (rk*rk^T)/sk^T*rk;
        Bk_ = Bk_ - (Bk_*sk*sk.transpose()*Bk_)/static_cast<double>(sk.transpose()*Bk_*sk) + (rk*rk.transpose())/(sk.dot(rk));
    }

    void GaitOptimizer::UpdateLagrangianGradients(const sp_matrix_t& A) {
        gradk_ = gradkp1_;
        gradkp1_ = dHdth + A.transpose()*dual_;
    }

    void GaitOptimizer::IncreaseTrustRegion(double rho) {
        if (rho > 0.75 && step_.lpNorm<Eigen::Infinity>() >= Delta_) {
            Delta_ = std::min(2*Delta_, max_trust_region_);
        }
    }

    void GaitOptimizer::DecreaseTrustRegion(const vector_t& step) {
        Delta_ = gamma_*step.lpNorm<Eigen::Infinity>();
        if (Delta_ <= 1e-5) {
            Delta_ = 1e-5;
        }
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

    vector_t GaitOptimizer::ContactTimesToQPVec() const {
        vector_t x(GetNumTimeNodes(num_ee_));
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int node = 0; node < contact_times_.at(ee).size(); node++) {
                x(GetNumTimeNodes(ee) + node) = contact_times_.at(ee).at(node).GetTime();
            }
        }

        return x;
    }

    void GaitOptimizer::AdjustBSize(int num_decision_vars) {
        if (run_num_ == 1) {
            Bk_ = matrix_t::Identity(num_decision_vars, num_decision_vars);
        }

        // TODO: Do something smarter where we keep the curvature for the parts that we know
        if (Bk_.rows() != num_decision_vars || Bk_.cols() != num_decision_vars) {
            Bk_ = matrix_t::Identity(num_decision_vars, num_decision_vars);
        }
    }

    double GaitOptimizer::GetStepNorm() const {
        return step_.norm();
    }

    void GaitOptimizer::ResetBFGSCondition() {
        past_decision_vars_ = 0;
    }

    std::vector<time_v> GaitOptimizer::GetContactTimes(double alpha) const {
        vector_t new_times = xk_ + alpha*step_;

        std::vector<time_v> contacts = contact_times_;

        for (int ee = 0; ee < num_ee_; ee++) {
            for (int idx = 0; idx < contacts.at(ee).size(); idx++) {
                contacts.at(ee).at(idx).SetTime(new_times(GetNumTimeNodes(ee) + idx));

                if (idx > 0) {
                    if (contacts.at(ee).at(idx-1).GetTime() - contacts.at(ee).at(idx).GetTime() <= 1e-3
                        && contacts.at(ee).at(idx-1).GetTime() - contacts.at(ee).at(idx).GetTime() > 0) {
                        contacts.at(ee).at(idx) = contacts.at(ee).at(idx-1);
                    }
                    assert(contacts.at(ee).at(idx-1).GetTime() <= contacts.at(ee).at(idx).GetTime());
                }
            }
        }

        return contacts;
    }
}