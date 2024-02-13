//
// Created by zolkin on 2/7/24.
//

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>

#include "qp/clarabel_interface.h"

namespace mpc {
    using namespace clarabel;

    // TODO: Time the constructor. Worried about data copy time
    ClarabelInterface::ClarabelInterface(const QPData& data, bool verbose) : QPInterface(data.num_decision_vars) {
        settings_ = DefaultSettings<double>::default_settings();
        settings_.verbose = verbose;
        settings_.tol_gap_rel = 1e-8;
        settings_.tol_gap_abs = 1e-8;
//        settings_.linesearch_backtrack_step = 0.2;
//        settings_.min_switch_step_length = 1e-2;
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

        slacks_ = sol.s;
        dual_ = sol.z;
        primal_ = sol.x;

        if (solve_quality_ == PrimalInfeasible) {
            vector_t temp = data.sparse_constraint_*primal_ + slacks_ - data.ub_;
            vector_t primal_product = data.sparse_constraint_*primal_;

            const matrix_t A = data.sparse_constraint_.toDense();

//            std::cout << "A top rows: " << std::endl;
//            std::cout << A.topRows<6>() << std::endl;

            const int nodes = 6;

//            std::cout << "Constraint matrix: \n";
//            for (int j = 0; j < data.GetTotalNumConstraints(); j++) {
//                for (int i = 0; i < nodes*12; i++) {
//                    std::cout << std::setw(15) << A(j, i);
//                }
//
//                std::cout << " |||| ";
//
//                for (int i = nodes*12; i < A.cols(); i++) {
//                    std::cout << std::setw(15) << A(j, i);
//                }
//
//                std::cout << std::endl;
//            }

//            std::cout << "Constraint matrix: \n" << data.sparse_constraint_.toDense() << std::endl;

//            std::cout << "A*z^* | s | b" << std::endl;
//            for (int i = 0; i < temp.size(); i++) {
//                std::cout << std::setw(15) << primal_product(i) << " | " <<
//                    std::setw(15) << slacks_(i) << " | " <<
//                    std::setw(15) << data.ub_(i) << std::endl;
//            }
        }

        return primal_;
    }

    SolveQuality ClarabelInterface::GetSolveQuality() const {
        return solve_quality_;
    }

    vector_t ClarabelInterface::GetDualSolution() const {
        return dual_;
    }

    void ClarabelInterface::ConfigureForInitialRun() {
        settings_.tol_gap_rel = 1e-11;
        settings_.tol_gap_abs = 1e-11;
    }

    void ClarabelInterface::ConfigureForRealTime(double run_time_iters) {
         settings_.tol_gap_rel = 1e-8; //1e-11;
         settings_.tol_gap_abs = 1e-8; //1e-11;
    }


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

        dP.resize(primal_.size(), primal_.size());
        dA.resize(nu_.size(), primal_.size());
        dG.resize(lam_.size(), primal_.size());

//        std::cout << "dz*z^*: " << dz*primal_.transpose() << std::endl;

        utils::SparseMatrixBuilder P_builder;
        P_builder.Reserve(1000);    // TODO: Adjust this number
        P_builder.SetMatrix(0.5*(dz*primal_.transpose() + primal_*dz.transpose()), 0, 0);
        dP.setFromTriplets(P_builder.GetTriplet().begin(), P_builder.GetTriplet().end());

        utils::SparseMatrixBuilder A_builder;
        A_builder.Reserve(1000);
        A_builder.SetMatrix((dnu*primal_.transpose() + nu_*dz.transpose()), 0, 0);
        dA.setFromTriplets(A_builder.GetTriplet().begin(), A_builder.GetTriplet().end());

        utils::SparseMatrixBuilder G_builder;
        G_builder.Reserve(1000);
        G_builder.SetMatrix((lam_.asDiagonal()*dlam*primal_.transpose() + lam_*dz.transpose()), 0, 0);
        dG.setFromTriplets(G_builder.GetTriplet().begin(), G_builder.GetTriplet().end());
//        std::cout << dG.toDense() << std::endl;
        matrix_t dG_dense = dG.toDense();
//        std::cout << "dG max coef: " << dG_dense.maxCoeff() << std::endl;
//        std::cout << "dA max coef: " << dA.toDense().maxCoeff() << std::endl;
//        const double dG_max = dG_dense.maxCoeff();
//        for (int i = 0; i < dG.rows(); i++) {
//            for (int j = 0; j < dG.cols(); j++) {
//                if (dG_max == dG_dense(i,j)) {
//                    std::cout << "dG max coeff row: " << i << " col: " << j << std::endl;
//                }
//            }
//        }

//        std::cout << "lam_ max: " << lam_.maxCoeff() << ", note: large value represents high sensitivity" << std::endl;
//        std::cout << "dlam max: " << dlam.maxCoeff() << std::endl;
//        for (int i = 0; i < lam_.size(); i++) {
//            if (lam_.maxCoeff() == lam_(i)) {
//                std::cout << "lam_ max occurs at index: " << i << std::endl;
//            }
//        }

//        std::cout << "dnu max: " << dnu.maxCoeff() << std::endl;
//        std::cout << "nu_ max: " << nu_.maxCoeff() << std::endl;
//
//        std::cout << "dz max: " << dz.maxCoeff() << std::endl;
//        std::cout << "primal max: " << primal_.maxCoeff() << std::endl;
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

        utils::Timer timer("KKT mat solve");
        timer.StartTimer();

        // TODO: Speed up! The setup for the solve is about 20ms
        matrix_t A(data.num_equality_, data.num_decision_vars);
        matrix_t G(data.num_inequality_, data.num_decision_vars);
        lam_.resize(data.num_inequality_);
        nu_.resize(data.num_equality_);
        vector_t h(data.num_inequality_);
        vector_t s_ineq(data.num_inequality_);

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
                    s_ineq.segment(ineq_idx, data.num_ee_location_constraints_) =
                            slacks_.segment(gen_idx, data.num_ee_location_constraints_);
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
                    s_ineq.segment(ineq_idx, data.num_force_box_constraints_) =
                            slacks_.segment(gen_idx, data.num_force_box_constraints_);
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
                    s_ineq.segment(ineq_idx, data.num_cone_constraints_) =
                            slacks_.segment(gen_idx, data.num_cone_constraints_);
                    ineq_idx += data.num_cone_constraints_;
                    gen_idx += data.num_cone_constraints_;
                    break;
            }
        }

        // TODO: Check that lam or eq constraint is always 0!

        assert(gen_idx == data.GetTotalNumConstraints());

        G = -G; // TODO: What's up with this negative sign?

        utils::SparseMatrixBuilder mat_builder;
        mat_builder.Reserve(3*data.sparse_constraint_.nonZeros() +
                        3*data.sparse_cost_.nonZeros()); // TODO: might be too much

        mat_builder.SetMatrix(data.sparse_cost_.toDense(), 0, 0);

//        mat_builder.SetMatrix(G.transpose(), 0, data.num_decision_vars);
        mat_builder.SetMatrix(G.transpose()*lam_.asDiagonal(), 0, data.num_decision_vars);

        mat_builder.SetMatrix(A.transpose(),
                              0, data.num_decision_vars + data.num_inequality_);

//        mat_builder.SetMatrix((lam_.asDiagonal()*G), data.num_decision_vars, 0);
//        mat_builder.SetMatrix((G.array().colwise() * lam_.array()), data.num_decision_vars, 0);
        mat_builder.SetMatrix(G, data.num_decision_vars, 0);

        mat_builder.SetVectorDiagonally(s_ineq, data.num_decision_vars, data.num_decision_vars);
        mat_builder.SetMatrix(A, data.num_decision_vars + data.num_inequality_, 0);

        assert(data.num_inequality_ + data.num_equality_ == data.GetTotalNumConstraints());

//        mat_builder.SetDiagonalMatrix(1e-3, 0, 0, data.num_decision_vars + data.GetTotalNumConstraints());

        sp_matrix_t diff_mat(data.num_decision_vars + data.GetTotalNumConstraints(),
                            data.num_decision_vars + data.GetTotalNumConstraints());

        diff_mat.setFromTriplets(mat_builder.GetTriplet().begin(), mat_builder.GetTriplet().end());

        // Getting some rank info
//        Eigen::FullPivHouseholderQR<matrix_t> qr;
//        qr.compute(A.transpose());
//        std::cout << "A^T rank: " << qr.rank() << ", A^T rows: " << qr.rows() << ", A^T cols: " << qr.cols() << std::endl;
//
//
//        qr.compute(G.transpose());
//        std::cout << "G^T rank: " << qr.rank() << ", G^T rows: " << qr.rows() << ", G^T cols: " << qr.cols() << std::endl;

        int non_zero_lam = 0;
        int non_zero_slacks = 0;
        constexpr double ZERO_THRESHOLD = 1e-6;
        for (int i = 0; i < lam_.size(); i++) {
            if (lam_(i) > ZERO_THRESHOLD) {
                non_zero_lam++;
            }

            if (s_ineq(i) > ZERO_THRESHOLD) {
                non_zero_slacks++;
            }
        }
        std::cout << "Nonzero lambda values: " << non_zero_lam << ", lambda length: " << lam_.size() << std::endl;
        std::cout << "Nonzero inequality slack values: " << non_zero_slacks << ", slacks length: " << s_ineq.size() << std::endl;
        std::cout << "Sum of nonzeros: " << non_zero_lam + non_zero_slacks << std::endl;

//        Eigen::JacobiSVD<matrix_t> svd(A);
//        double cond = svd.singularValues()(0)
//                      / svd.singularValues()(svd.singularValues().size()-1);
//
//        std::cout << "Equality constraint mat max svd value: " << svd.singularValues()(0) << std::endl;
//        std::cout << "Equality constraint mat min svd value: " << svd.singularValues()(svd.singularValues().size()-1) << std::endl;
//        std::cout << "Equality constraint mat Condition number: " << cond << std::endl;
//
//        Eigen::JacobiSVD<matrix_t> svd2(G);
//        cond = svd2.singularValues()(0)
//                      / svd2.singularValues()(svd2.singularValues().size()-1);
//
//        std::cout << "Inequality constraint mat max svd value: " << svd2.singularValues()(0) << std::endl;
//        std::cout << "Inequality constraint mat min svd value: " << svd2.singularValues()(svd2.singularValues().size()-1) << std::endl;
//        std::cout << "Inequality constraint mat Condition number: " << cond << std::endl;
//
//        Eigen::JacobiSVD<matrix_t> svd3(data.sparse_constraint_);
//        cond = svd3.singularValues()(0)
//                      / svd3.singularValues()(svd3.singularValues().size()-1);
//
//        std::cout << "Constraint mat max svd value: " << svd3.singularValues()(0) << std::endl;
//        std::cout << "Constraint mat min svd value: " << svd3.singularValues()(svd3.singularValues().size()-1) << std::endl;
//        std::cout << "Constraint mat Condition number: " << cond << std::endl;
//
//        Eigen::JacobiSVD<matrix_t> svd4(data.sparse_cost_);
//        cond = svd4.singularValues()(0)
//                      / svd4.singularValues()(svd4.singularValues().size()-1);
//
//        std::cout << "Cost mat max svd value: " << svd4.singularValues()(0) << std::endl;
//        std::cout << "Cost mat min svd value: " << svd4.singularValues()(svd4.singularValues().size()-1) << std::endl;
//        std::cout << "Cost mat Condition number: " << cond << std::endl;

        double max = -1e10;
        int max_row = -1;
        int max_col = -1;
        matrix_t mat = diff_mat.toDense();
        for (int i = 0; i < diff_mat.rows(); i++) {
            for (int j = 0; j < diff_mat.cols(); j++) {
                if (mat(i,j) > max) {
                    max = mat(i,j);
                    max_row = i;
                    max_col = j;
                }
            }
        }

        std::cout << "KKT mat max val: " << max << ", row: " << max_row << ", col: " << max_col << std::endl;

        // Solve the system
        // The system solve is only about 1-2ms
        d_.resize(data.num_decision_vars + data.num_equality_ + data.num_inequality_);
        d_ << dx, vector_t::Zero(data.num_equality_ + data.num_inequality_);


        assert(diff_mat.rows() == diff_mat.cols());

        sp_matrix_t temp = diff_mat.transpose()*diff_mat;

        temp.makeCompressed();
//        temp.pruned(1e-6);

        // TODO: Might be bad as it wants the A^TA to be PD
        Eigen::LeastSquaresConjugateGradient<sp_matrix_t> solver;
//        Eigen::ConjugateGradient<sp_matrix_t, Eigen::Lower|Eigen::Upper> solver;
//        Eigen::SimplicialLDLT<sp_matrix_t> solver;

//        solver.compute(temp);

//        diff_mat.pruned(ZERO_THRESHOLD);

        // TODO: Check if I am row rank or column rank deficient and think about those implications
        solver.analyzePattern(diff_mat);
        solver.factorize(diff_mat);
        // TODO: Can't solve the high tolerance problem with 50 nodes
        solver.setTolerance(1e-7); // TODO: Higher tolerances make this much better! Ideally AT LEAST 1e-8!
        // TODO: Note it seems more related to the conditioning of the problem rather than the tolerance.
        solver.setMaxIterations(10*temp.rows());

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Could not factorize the differential matrix.");
        }

        vector_t temp2 = d_; //diff_mat.transpose()*d_;
        d_ = -solver.solve(temp2);
        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error()      << std::endl;
        std::cout << "max elemnt of d_: " << d_.maxCoeff() << std::endl;

//        std::cout << "d_: " << d_ << std::endl;
//        d_ = diff_mat.transpose()*d_;
//        std::cout << "again d_: " << d_ << std::endl;
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Could not solve the differential matrix");
        }
        timer.StopTimer();

        timer.PrintElapsedTime();

        num_equality_constraints_ = data.num_equality_;
        num_inequality_constraints_ = data.num_inequality_;

        std::cout << "Force box constraints: " << data.num_force_box_constraints_ << std::endl;
        std::cout << "Friction cone constraints: " << data.num_cone_constraints_ << std::endl;
        std::cout << "EE box constraints: " << data.num_ee_location_constraints_ << std::endl;
        std::cout << "EE start constraints: " << data.num_start_ee_constraints_ << std::endl;
    }

    vector_t ClarabelInterface::Computedx(const Eigen::SparseMatrix<double>& P, const mpc::vector_t& q,
                                      const mpc::vector_t& xstar) {
        dx_.noalias() = P*xstar + q;
        return dx_;
    }

} // mpc