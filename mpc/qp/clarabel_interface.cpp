//
// Created by zolkin on 2/7/24.
//

#include "qp/clarabel_interface.h"

namespace mpc {
    using namespace clarabel;

    // TODO: Time the constructor. Worried about data copy time
    ClarabelInterface::ClarabelInterface(const QPData& data, bool verbose) : QPInterface(data.num_decision_vars) {
        settings_ = DefaultSettings<double>::default_settings();
        settings_.verbose = verbose;
        settings_.tol_gap_rel = 1e-5;
        settings_.tol_gap_abs = 1e-5;
//        settings_.max_iter = 10;
    }

    void ClarabelInterface::SetupQP(mpc::QPData& data, const mpc::vector_t& warm_start) {
        cones_ = {};
        num_equality_constraints_ = 0;
        num_inequality_constraints_ = 0;
        for (int i = 0; i < data.constraints_.size(); i++) {
            switch (data.constraints_.at(i)) {
                case Constraints::Dynamics:
                    cones_.push_back(ZeroConeT<double>(data.num_dynamics_constraints));
                    num_equality_constraints_ += data.num_dynamics_constraints;
                    break;
                case Constraints::JointForwardKinematics:
                        throw std::runtime_error("not supported yet");
                    break;
                case Constraints::EndEffectorLocation:
                    cones_.push_back(NonnegativeConeT<double>(data.num_ee_location_constraints_));
                    num_inequality_constraints_ += data.num_ee_location_constraints_;

                    cones_.push_back(ZeroConeT<double>(data.num_start_ee_constraints_));
                    num_equality_constraints_ += data.num_start_ee_constraints_;
                    break;
                case Constraints::ForceBox:
                    cones_.push_back(NonnegativeConeT<double>(data.num_force_box_constraints_));
                    num_inequality_constraints_ += data.num_force_box_constraints_;
                    break;
                case Constraints::JointBox:
                    throw std::runtime_error("Joint box not implemented with clarabel yet");

                    // 2x to split the lb and ub
                    cones_.push_back(NonnegativeConeT<double>(data.num_box_constraints_));
                    break;
                case Constraints::FrictionCone:
                    cones_.push_back(NonnegativeConeT<double>(data.num_cone_constraints_));
                    num_inequality_constraints_ += data.num_cone_constraints_;
                    break;
            }
        }

//        data.sparse_cost_.makeCompressed();
//        data.sparse_constraint_.makeCompressed();
        solver_ = std::make_unique<DefaultSolver<double>>(data.sparse_cost_, data.cost_linear,
                data.sparse_constraint_,data.ub_, cones_, settings_);
    }

    vector_t ClarabelInterface::Solve(const QPData& data) {
        solver_->solve();

        DefaultSolution<double> sol = solver_->solution();

        assert(sol.x.size() == data.num_decision_vars);
        assert(sol.z.size() == data.GetTotalNumConstraints());

        switch (sol.status) {
            case clarabel::SolverStatus::Solved:
                solve_quality_ = Solved;
                break;
            case clarabel::SolverStatus::Unsolved:
                solve_quality_ = Unsolved;
                break;
            case clarabel::SolverStatus::PrimalInfeasible:
                solve_quality_ = PrimalInfeasible;
                break;
            case clarabel::SolverStatus::DualInfeasible:
                solve_quality_ = DualInfeasible;
                break;
            case clarabel::SolverStatus::MaxIterations:
                solve_quality_ = MaxIter;
                break;
            case clarabel::SolverStatus::AlmostSolved:
                solve_quality_ = SolvedInacc;
                break;
            case clarabel::SolverStatus::AlmostPrimalInfeasible:
                solve_quality_ = PrimalInfeasibleInacc;
                break;
            case clarabel::SolverStatus::AlmostDualInfeasible:
                solve_quality_ = DualInfeasibleInacc;
                break;
            default:
                solve_quality_ = Other;
        }

        dual_ = sol.z;
        primal_ = sol.x;
        return primal_;
    }

    SolveQuality ClarabelInterface::GetSolveQuality() const {
        return solve_quality_;
    }

    vector_t ClarabelInterface::GetDualSolution() const {
        return dual_;
    }

    void ClarabelInterface::ConfigureForInitialRun() const {}

    void ClarabelInterface::ConfigureForRealTime(double run_time_iters) const {}


    vector_t ClarabelInterface::Getdx() const {
        return dx_;
    }

    void ClarabelInterface::CalcDerivativeWrtVecs(vector_t& dq, vector_t& db, vector_t& dh) {
        dq = d_.head(primal_.size());

        dh.noalias() = lam_.asDiagonal()*-1*d_.segment(primal_.size(), num_inequality_constraints_);

        db = -d_.tail(num_equality_constraints_);
    }

    void ClarabelInterface::CalcDerivativeWrtMats(sp_matrix_t& dP, sp_matrix_t& dA, sp_matrix_t& dG) {

        const vector_t& dz = d_.head(primal_.size());
        const vector_t& dnu = d_.tail(num_equality_constraints_);
        const vector_t& dlam = d_.segment(primal_.size(), num_inequality_constraints_);

        dP = 0.5*(dz*primal_.transpose() + primal_*dz.transpose()).sparseView();    // TODO: Is sparse view efficient?

        dA = (dnu*primal_.transpose() + nu_*dz.transpose()).sparseView();

        dG = (lam_.asDiagonal()*dlam*primal_.transpose() + lam_*dz.transpose()).sparseView();
    }

    void ClarabelInterface::SetupDerivativeCalcs(mpc::vector_t& dx, mpc::vector_t& dy_l, mpc::vector_t& dy_u,
                                                 const QPData& data) {
        /*
         * Setup the system to solve for the differentials.
         *
         *   [Q   G^T*D(lam)   A^T]^(-1)  [dldz^T]
         * - [G   D(G*z - h)    0 ]       [  0  ]
         *   [A        0        0 ]       [  0  ]
         *
         *   Where Q is the hessian, A is the equality constraint matrix,
         *   G is the inequality constraint matrix, h is the inequality upper bounds,
         *   z is the optimal primal solution, lam is the optimal dual for the inequality constraints
         *
         *   The first row goes into dz (differential of the primal), the second row is dlam (differential
         *   of the inequality dual), and the third row is dnu (differential of the equality dual).
         */

        // Construct the matrix
        // - Need to split the constraint matrix into the equality and inequality constraints.
        //      can split them based on the cones

        matrix_t A(num_equality_constraints_, data.num_decision_vars);
        matrix_t G(num_inequality_constraints_, data.num_decision_vars);
        lam_.resize(num_inequality_constraints_);
        nu_.resize(num_equality_constraints_);
        vector_t h(num_inequality_constraints_);

        int eq_idx = 0;
        int ineq_idx = 0;
        int gen_idx = 0;
        for (int i = 0; i < data.constraints_.size(); i++) {
            switch (data.constraints_.at(i)) {
                case Constraints::Dynamics:
                    A.middleRows(eq_idx, data.num_dynamics_constraints) =
                            data.sparse_constraint_.middleRows(gen_idx, data.num_dynamics_constraints);
                    nu_.segment(eq_idx, data.num_dynamics_constraints) =
                            dual_.segment(gen_idx, data.num_dynamics_constraints);
                    eq_idx += data.num_dynamics_constraints;
                    gen_idx += data.num_dynamics_constraints;
                    break;
                case Constraints::JointForwardKinematics:
                    throw std::runtime_error("not supported yet");
                    break;
                case Constraints::EndEffectorLocation:
                    G.middleRows(ineq_idx, data.num_ee_location_constraints_) =
                            data.sparse_constraint_.middleRows(gen_idx, data.num_ee_location_constraints_);
                    lam_.segment(ineq_idx, data.num_ee_location_constraints_) =
                            dual_.segment(gen_idx, data.num_ee_location_constraints_);
                    h.segment(ineq_idx, data.num_ee_location_constraints_) =
                            data.ub_.segment(gen_idx, data.num_ee_location_constraints_);
                    ineq_idx += data.num_ee_location_constraints_;
                    gen_idx += data.num_ee_location_constraints_;

                    A.middleRows(eq_idx, data.num_start_ee_constraints_) =
                            data.sparse_constraint_.middleRows(gen_idx, data.num_start_ee_constraints_);
                    nu_.segment(eq_idx, data.num_start_ee_constraints_) =
                            dual_.segment(gen_idx, data.num_start_ee_constraints_);
                    eq_idx += data.num_start_ee_constraints_;
                    gen_idx += data.num_start_ee_constraints_;
                    break;
                case Constraints::ForceBox:
                    G.middleRows(ineq_idx, data.num_force_box_constraints_) =
                            data.sparse_constraint_.middleRows(gen_idx, data.num_force_box_constraints_);
                    lam_.segment(ineq_idx, data.num_force_box_constraints_) =
                            dual_.segment(gen_idx, data.num_force_box_constraints_);
                    h.segment(ineq_idx, data.num_force_box_constraints_) =
                            data.ub_.segment(gen_idx, data.num_force_box_constraints_);
                    ineq_idx += data.num_force_box_constraints_;
                    gen_idx += data.num_force_box_constraints_;
                    break;
                case Constraints::JointBox:
                    throw std::runtime_error("Joint box not implemented with clarabel yet");
                    break;
                case Constraints::FrictionCone:
                    G.middleRows(ineq_idx, data.num_cone_constraints_) =
                            data.sparse_constraint_.middleRows(gen_idx, data.num_cone_constraints_);
                    lam_.segment(ineq_idx, data.num_cone_constraints_) =
                            dual_.segment(gen_idx, data.num_cone_constraints_);
                    h.segment(ineq_idx, data.num_cone_constraints_) =
                            data.ub_.segment(gen_idx, data.num_cone_constraints_);
                    ineq_idx += data.num_cone_constraints_;
                    gen_idx += data.num_cone_constraints_;
                    break;
            }
        }

        matrix_t DiffMat (data.num_decision_vars + data.GetTotalNumConstraints(),
                          data.num_decision_vars + data.GetTotalNumConstraints());
        DiffMat.setZero();

        DiffMat.topLeftCorner(data.num_decision_vars, data.num_decision_vars) = -data.sparse_cost_;
        DiffMat.block(0, data.num_decision_vars,
                      data.num_decision_vars, num_inequality_constraints_) = -G.transpose() * lam_.asDiagonal();
        DiffMat.block(0, data.num_decision_vars + num_inequality_constraints_,
                      data.num_decision_vars, num_equality_constraints_) = -A.transpose();

        DiffMat.block(data.num_decision_vars, 0, G.rows(), G.cols()) = -G;
        DiffMat.block(data.num_decision_vars, G.cols(),
                      num_inequality_constraints_, num_inequality_constraints_) =
                              -1*(G*primal_ - h).asDiagonal();

        DiffMat.bottomLeftCorner(num_equality_constraints_, data.num_decision_vars) = -A;

        // Solve the system
        d_.resize(data.num_decision_vars + num_equality_constraints_ + num_inequality_constraints_);
        d_ << dx, vector_t::Zero(num_equality_constraints_ + num_inequality_constraints_);

//        DiffMat.llt().solveInPlace(d_);
        DiffMat.ldlt().solveInPlace(d_);
    }

    vector_t ClarabelInterface::Computedx(const Eigen::SparseMatrix<double>& P, const mpc::vector_t& q,
                                      const mpc::vector_t& xstar) {
        dx_.noalias() = P*xstar + q;
        return dx_;
    }

} // mpc