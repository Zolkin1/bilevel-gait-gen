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
//        settings_.tol_gap_rel = 5e-5;
//        settings_.tol_gap_abs = 5e-5;
//        settings_.max_iter = 10;
    }

    void ClarabelInterface::SetupQP(mpc::QPData& data, const mpc::vector_t& warm_start) {
        cones_ = {};

        for (int i = 0; i < data.constraints_.size(); i++) {
            switch (data.constraints_.at(i)) {
                case Constraints::Dynamics:
                    cones_.push_back(ZeroConeT<double>(data.num_dynamics_constraints));
                    break;
                case Constraints::JointForwardKinematics:
                        throw std::runtime_error("not supported yet");
                    break;
                case Constraints::EndEffectorLocation:
                    cones_.push_back(NonnegativeConeT<double>(data.num_ee_location_constraints_));

                    cones_.push_back(ZeroConeT<double>(data.num_start_ee_constraints_));

                    break;
                case Constraints::ForceBox:
                    cones_.push_back(NonnegativeConeT<double>(data.num_force_box_constraints_));
                    break;
                case Constraints::JointBox:
                    throw std::runtime_error("Joint box not implemented with clarabel yet");

                    // 2x to split the lb and ub
                    cones_.push_back(NonnegativeConeT<double>(data.num_box_constraints_));
                    break;
                case Constraints::FrictionCone:
                    cones_.push_back(NonnegativeConeT<double>(data.num_cone_constraints_));
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

        return sol.x;
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
        throw std::runtime_error("GetDx not impelemnted yet for Clarabel.");
    }

    void ClarabelInterface::CalcDerivativeWrtVecs(mpc::vector_t& dq, mpc::vector_t dl, mpc::vector_t du) {
        throw std::runtime_error("Not impelemnted yet for Clarabel.");
    }

    void ClarabelInterface::CalcDerivativeWrtMats(Eigen::SparseMatrix<double>& dP, Eigen::SparseMatrix<double>& dA) {
        throw std::runtime_error("Not impelemnted yet for Clarabel.");
    }

    void ClarabelInterface::SetupDerivativeCalcs(mpc::vector_t& dx, mpc::vector_t& dy_l, mpc::vector_t& dy_u) {
        throw std::runtime_error("Not impelemnted yet for Clarabel.");
    }

    void ClarabelInterface::Computedx(const Eigen::SparseMatrix<double>& P, const mpc::vector_t& q,
                                      const mpc::vector_t& xstar) {
        throw std::runtime_error("Not impelemnted yet for Clarabel.");
    }

} // mpc