//
// Created by zolkin on 2/7/24.
//

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>

#include "lsqr-cpp/lsqr.h"

#include "qp/clarabel_interface.h"

namespace mpc {
    using namespace clarabel;

    // TODO: Time the constructor. Worried about data copy time
    ClarabelInterface::ClarabelInterface(const QPData& data, bool verbose) : QPInterface(data.num_decision_vars) {
        settings_ = DefaultSettings<double>::default_settings();
        settings_.verbose = verbose;
        settings_.tol_gap_rel = 1e-8;
        settings_.tol_gap_abs = 1e-8;
        settings_.tol_feas = 1e-10;
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
                case Constraints::TDPosition:
                    cones_.push_back(ZeroConeT<double>(data.num_td_pos_constraints_));
                    break;
                case Constraints::EndEffectorStart:
                    cones_.push_back(ZeroConeT<double>(data.num_start_ee_constraints_));
                    break;
                case Constraints::Raibert:
                    cones_.push_back(ZeroConeT<double>(data.num_raibert_constraints_));
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
            const std::string& error = "Primal infeasible.";
            throw (error);

            // TODO: Is this the behavior I want?
            primal_ = vector_t::Constant(sol.x.size(), 1e10);

            vector_t temp = data.sparse_constraint_*primal_ + slacks_ - data.ub_;
            vector_t primal_product = data.sparse_constraint_*primal_;

//            const matrix_t A = data.sparse_constraint_.toDense();

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
        settings_.tol_gap_rel = 1e-15;
        settings_.tol_gap_abs = 1e-15;
    }

    void ClarabelInterface::ConfigureForRealTime(double run_time_iters) {
        // Note: the higher quality definetly shows
        // TODO: Tune tolerances
         settings_.tol_gap_rel = 1e-15; //1e-15;
         settings_.tol_gap_abs = 1e-15; //1e-15;
    }


    vector_t ClarabelInterface::Getdx() const {
        return dx_;
    }

    void ClarabelInterface::CalcDerivativeWrtVecs(vector_t& dq, vector_t& db, vector_t& dh) {
        dq = d_.head(primal_.size());

        dh.noalias() = lam_.asDiagonal()*-1*d_.segment(primal_.size(), num_inequality_constraints_ - 0);

        db = -d_.tail(num_equality_constraints_);

//        qp_file << "dq: \n" << dq << std::endl;
//        qp_file << std::endl;
//        qp_file << "dh: \n" << dh << std::endl;
//        qp_file << std::endl;
//        qp_file << "db: \n" << db << std::endl;
//        qp_file << std::endl;
//        qp_file.close();
    }

    void ClarabelInterface::CalcDerivativeWrtMats(matrix_t& dP, matrix_t& dA, matrix_t& dG) {

        const int num_cone_constraints = 0;

        const vector_t& dz = d_.head(primal_.size());
        const vector_t& dnu = d_.tail(num_equality_constraints_);
        const vector_t& dlam = d_.segment(primal_.size(), num_inequality_constraints_ - num_cone_constraints);

        // TODO: Resizing takes a non-trivial amount of time
//        dP.resize(primal_.size(), primal_.size());
        dA.resize(nu_.size(), primal_.size());
        dG.resize(lam_.size(), primal_.size());

        assert(dnu.size() == nu_.size());
        assert(dlam.size() == lam_.size());
        assert(dz.size() == primal_.size());

//        std::cout << "dz*z^*: " << dz*primal_.transpose() << std::endl;

// P is always zero for now, so ignore it for speed.
//        utils::SparseMatrixBuilder P_builder;
//        P_builder.Reserve(1000);
//        P_builder.SetMatrix(0.5*(dz*primal_.transpose() + primal_*dz.transpose()), 0, 0);
//        dP.setFromTriplets(P_builder.GetTriplet().begin(), P_builder.GetTriplet().end());

// Building the sparse matrix is very slow
//        utils::Timer Asparse_timer("dA sparse creation");
//        Asparse_timer.StartTimer();
//        utils::SparseMatrixBuilder A_builder;
//        A_builder.Reserve(500000); // ~471,200 elements actually in here, should this be a dense matrix?
//        A_builder.SetMatrix((dnu*primal_.transpose() + nu_*dz.transpose()), 0, 0);
//        dA.setFromTriplets(A_builder.GetTriplet().begin(), A_builder.GetTriplet().end());
//        Asparse_timer.StopTimer();
//        Asparse_timer.PrintElapsedTime();

        utils::Timer Adense_timer("dA dense creation");
        Adense_timer.StartTimer();
        dA.noalias() = dnu*primal_.transpose() + nu_*dz.transpose();
        Adense_timer.StopTimer();
//        Adense_timer.PrintElapsedTime();

//        const int row = 0; // 1
//        const int col = 2; // 750
//        std::cout << "dA at (" << row << "," << col << "): " << dA.toDense()(row, col) << std::endl;

//        utils::SparseMatrixBuilder G_builder;
//        G_builder.Reserve(500000);
//        G_builder.SetMatrix((lam_.asDiagonal()*dlam*primal_.transpose() + lam_*dz.transpose()), 0, 0);
//        dG.setFromTriplets(G_builder.GetTriplet().begin(), G_builder.GetTriplet().end());

        utils::Timer Gdense_timer("dG dense creation");
        Gdense_timer.StartTimer();
        dG.noalias() = lam_.asDiagonal()*dlam*primal_.transpose() + lam_*dz.transpose();
        Gdense_timer.StopTimer();
//        Gdense_timer.PrintElapsedTime();

//        qp_file << "dP: \n" << dP.toDense() << std::endl;
//        qp_file << std::endl;
//        qp_file << "dA: \n" << dA.toDense() << std::endl;
//        qp_file << std::endl;
//        qp_file << "dG: \n" << dG.toDense() << std::endl;
//        qp_file << std::endl;
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

        num_equality_constraints_ = data.num_equality_; // - data.num_td_pos_constraints_
        num_inequality_constraints_ = data.num_inequality_;

        // TODO: Speed up! The setup for the solve is about 20ms
        utils::Timer mat_build_timer("matrix building");
        mat_build_timer.StartTimer();
//        matrix_t A(data.num_equality_, data.num_decision_vars);
//        matrix_t G(data.num_inequality_, data.num_decision_vars);
        lam_.resize(data.num_inequality_);
        nu_.resize(num_equality_constraints_);
//        vector_t h(data.num_inequality_);
//        vector_t s_ineq(data.num_inequality_);
//        vector_t b(data.num_equality_);


        utils::SparseMatrixBuilder mat_builder;
        mat_builder.Reserve(3*data.sparse_constraint_.nonZeros() +
                            3*data.sparse_cost_.nonZeros()); // TODO: might be too much

        int eq_idx = 0;
        int ineq_idx = 0;
        int gen_idx = 0;
        // TODO: should be able to build most of this straight into the builder rather than into a dense matrix first.
        for (int i = 0; i < data.constraints_.size(); i++) {
            switch (data.constraints_.at(i)) {
                case Constraints::Dynamics:
//                    A.middleRows(eq_idx, data.num_dynamics_constraints) =
//                            data.sparse_constraint_.middleRows(gen_idx, data.num_dynamics_constraints);

                    // Put this in where all the equality constraints go
                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_dynamics_constraints).transpose(),
                                          0, data.num_decision_vars + data.num_inequality_ + eq_idx);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_dynamics_constraints),
                                          data.num_decision_vars + data.num_inequality_ + eq_idx,
                                          0);

                    nu_.segment(eq_idx, data.num_dynamics_constraints) =
                            dual_.segment(gen_idx, data.num_dynamics_constraints);

//                    b.segment(eq_idx, data.num_dynamics_constraints) =
//                            data.ub_.segment(gen_idx, data.num_dynamics_constraints);

                    eq_idx += data.num_dynamics_constraints;
                    gen_idx += data.num_dynamics_constraints;
                    break;
                case Constraints::JointForwardKinematics:
                    throw std::runtime_error("not supported yet");
                    break;
                case Constraints::EndEffectorLocation:
//                    G.middleRows(ineq_idx, data.num_ee_location_constraints_) =
//                            data.sparse_constraint_.middleRows(gen_idx, data.num_ee_location_constraints_);

                    lam_.segment(ineq_idx, data.num_ee_location_constraints_) =
                            dual_.segment(gen_idx, data.num_ee_location_constraints_);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_ee_location_constraints_),
                                          data.num_decision_vars + ineq_idx, 0);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_ee_location_constraints_).transpose()
                                            * dual_.segment(gen_idx, data.num_ee_location_constraints_).asDiagonal(),
                                          0, data.num_decision_vars + ineq_idx);

//                    h.segment(ineq_idx, data.num_ee_location_constraints_) =
//                            data.ub_.segment(gen_idx, data.num_ee_location_constraints_);

//                    s_ineq.segment(ineq_idx, data.num_ee_location_constraints_) =
//                            slacks_.segment(gen_idx, data.num_ee_location_constraints_);

                    mat_builder.SetVectorDiagonally(slacks_.segment(gen_idx, data.num_ee_location_constraints_),
                                                    data.num_decision_vars + ineq_idx, data.num_decision_vars + ineq_idx);

                    ineq_idx += data.num_ee_location_constraints_;
                    gen_idx += data.num_ee_location_constraints_;
                    break;
                case Constraints::EndEffectorStart:
                    //                    A.middleRows(eq_idx, data.num_start_ee_constraints_) =
//                            data.sparse_constraint_.middleRows(gen_idx, data.num_start_ee_constraints_);

                    // Put this in where all the equality constraints go
                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_start_ee_constraints_).transpose(),
                                          0, data.num_decision_vars + data.num_inequality_ + eq_idx);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_start_ee_constraints_),
                                          data.num_decision_vars + data.num_inequality_ + eq_idx,
                                          0);

                    nu_.segment(eq_idx, data.num_start_ee_constraints_) =
                            dual_.segment(gen_idx, data.num_start_ee_constraints_);
//                    b.segment(eq_idx, data.num_start_ee_constraints_) =
//                            data.ub_.segment(gen_idx, data.num_start_ee_constraints_);
                    eq_idx += data.num_start_ee_constraints_;
                    gen_idx += data.num_start_ee_constraints_;
                    break;
                case Constraints::ForceBox:
//                    G.middleRows(ineq_idx, data.num_force_box_constraints_) =
//                            data.sparse_constraint_.middleRows(gen_idx, data.num_force_box_constraints_);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_force_box_constraints_),
                                          data.num_decision_vars + ineq_idx, 0);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_force_box_constraints_).transpose()
                                          * dual_.segment(gen_idx, data.num_force_box_constraints_).asDiagonal(),
                                          0, data.num_decision_vars + ineq_idx);

                    lam_.segment(ineq_idx, data.num_force_box_constraints_) =
                            dual_.segment(gen_idx, data.num_force_box_constraints_);

//                    h.segment(ineq_idx, data.num_force_box_constraints_) =
//                            data.ub_.segment(gen_idx, data.num_force_box_constraints_);

//                    s_ineq.segment(ineq_idx, data.num_force_box_constraints_) =
//                            slacks_.segment(gen_idx, data.num_force_box_constraints_);

                    mat_builder.SetVectorDiagonally(slacks_.segment(gen_idx, data.num_force_box_constraints_),
                                                    data.num_decision_vars + ineq_idx, data.num_decision_vars + ineq_idx);

                    ineq_idx += data.num_force_box_constraints_;
                    gen_idx += data.num_force_box_constraints_;
                    break;
                case Constraints::JointBox:
                    throw std::runtime_error("Joint box not implemented with clarabel yet");
                    break;
                case Constraints::FrictionCone:
//                    G.middleRows(ineq_idx, data.num_cone_constraints_) =
//                            data.sparse_constraint_.middleRows(gen_idx, data.num_cone_constraints_);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_cone_constraints_),
                                          data.num_decision_vars + ineq_idx, 0);

                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_cone_constraints_).transpose()
                                          * dual_.segment(gen_idx, data.num_cone_constraints_).asDiagonal(),
                                          0, data.num_decision_vars + ineq_idx);

                    lam_.segment(ineq_idx, data.num_cone_constraints_) =
                            dual_.segment(gen_idx, data.num_cone_constraints_);

//                    h.segment(ineq_idx, data.num_cone_constraints_) =
//                            data.ub_.segment(gen_idx, data.num_cone_constraints_);

//                    s_ineq.segment(ineq_idx, data.num_cone_constraints_) =
//                            slacks_.segment(gen_idx, data.num_cone_constraints_);

                    mat_builder.SetVectorDiagonally(slacks_.segment(gen_idx, data.num_cone_constraints_),
                                                    data.num_decision_vars + ineq_idx, data.num_decision_vars + ineq_idx);

                    ineq_idx += data.num_cone_constraints_;
                    gen_idx += data.num_cone_constraints_;
                    break;
                case Constraints::TDPosition:
                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_td_pos_constraints_).transpose(),
                                        0, data.num_decision_vars + data.num_inequality_ + eq_idx);
                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_td_pos_constraints_),
                                          data.num_decision_vars + data.num_inequality_ + eq_idx,
                                          0);

                    nu_.segment(eq_idx, data.num_td_pos_constraints_) =
                            dual_.segment(gen_idx, data.num_td_pos_constraints_);

                    eq_idx += data.num_td_pos_constraints_;
                    gen_idx += data.num_td_pos_constraints_;
                    break;

                case Constraints::Raibert:
                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_raibert_constraints_).transpose(),
                                          0, data.num_decision_vars + data.num_inequality_ + eq_idx);
                    mat_builder.SetMatrix(data.sparse_constraint_.middleRows(gen_idx, data.num_raibert_constraints_),
                                          data.num_decision_vars + data.num_inequality_ + eq_idx,
                                          0);

                    nu_.segment(eq_idx, data.num_raibert_constraints_) =
                            dual_.segment(gen_idx, data.num_raibert_constraints_);

                    eq_idx += data.num_raibert_constraints_;
                    gen_idx += data.num_raibert_constraints_;
                    break;
            }
        }

        mat_builder.SetMatrix(data.sparse_cost_, 0, 0);

        sp_matrix_t diff_mat(data.num_decision_vars + data.GetTotalNumConstraints(),
                             data.num_decision_vars + data.GetTotalNumConstraints());

        diff_mat.setFromTriplets(mat_builder.GetTriplet().begin(), mat_builder.GetTriplet().end());

        // Print slack and lambda
//        std::cout << std::setw(15) << "slacks | " << std::setw(15) << "lambda" << std::endl;
//        for (int i = 0; i < lam_.size(); i++) {
//            std::cout << std::setw(15) << s_ineq(i) << " | " << std::setw(15) << lam_(i) << std::endl;
//        }

//        qp_file.open("clarabel_qp_data.txt");
//        qp_file << "A: \n" << A << std::endl;
//        qp_file << std::endl;
//        qp_file << "G: \n" << G << std::endl;
//        qp_file << std::endl;
//        qp_file << "b: \n" << b << std::endl;
//        qp_file << std::endl;
//        qp_file << "h: \n" << h << std::endl;
//        qp_file << std::endl;
//        qp_file << "P: \n" << data.sparse_cost_.toDense() << std::endl;
//        qp_file << std::endl;
//        qp_file << "q: \n" << data.cost_linear << std::endl;
//        qp_file << std::endl;
//        qp_file << "x: \n" << primal_ << std::endl;
//        qp_file << std::endl;
//        qp_file << "dual: \n" << dual_ << std::endl;
//        qp_file << std::endl;

        // TODO: Check that lam or eq constraint is always 0!

//        assert(gen_idx == data.GetTotalNumConstraints());
//
//        assert(G.rows() == data.num_inequality_);
//        assert(A.rows() == data.num_equality_);
//
//        mat_builder.Reserve(3*data.sparse_constraint_.nonZeros() +
//                            3*data.sparse_cost_.nonZeros());
//
//        mat_builder.SetMatrix(data.sparse_cost_.toDense(), 0, 0);
//
//        mat_builder.SetMatrix(G.transpose()*lam_.asDiagonal(), 0, data.num_decision_vars);
//
//        mat_builder.SetMatrix(A.transpose(),
//                              0, data.num_decision_vars + G.rows());
//
//        mat_builder.SetMatrix(G, data.num_decision_vars, 0);
//
//        mat_builder.SetVectorDiagonally(s_ineq, data.num_decision_vars, data.num_decision_vars);
//        mat_builder.SetMatrix(A, data.num_decision_vars + G.rows(), 0);
//
//        assert(data.num_inequality_ + data.num_equality_ == data.GetTotalNumConstraints());
//
//        sp_matrix_t diff_mat(data.num_decision_vars + data.GetTotalNumConstraints(),
//                            data.num_decision_vars + data.GetTotalNumConstraints());
//
//        diff_mat.setFromTriplets(mat_builder.GetTriplet().begin(), mat_builder.GetTriplet().end());

        mat_build_timer.StopTimer();
        mat_build_timer.PrintElapsedTime();

//        std::cout << "approx: " << diff_mat.isApprox(diff_mat1, 1e-12) << std::endl;

//        int non_zero_lam = 0;
//        int non_zero_slacks = 0;
//        constexpr double ZERO_THRESHOLD = 1e-6;
//        for (int i = 0; i < lam_.size(); i++) {
//            if (lam_(i) > ZERO_THRESHOLD) {
//                non_zero_lam++;
//            }
//
//            if (s_ineq(i) > ZERO_THRESHOLD) {
//                non_zero_slacks++;
//            }
//        }
//        std::cout << "Nonzero lambda values: " << non_zero_lam << ", lambda length: " << lam_.size() << std::endl;
//        std::cout << "Nonzero inequality slack values: " << non_zero_slacks << ", slacks length: " << s_ineq.size() << std::endl;
//        std::cout << "Sum of nonzeros: " << non_zero_lam + non_zero_slacks << std::endl;
//        std::cout << "A rows: " << A.rows() << std::endl;
//        std::cout << "Decision vars: " << data.num_decision_vars << std::endl;


        // Solve the system
        // The system solve is only about 1-2ms
        d_.resize(data.num_decision_vars + num_equality_constraints_ + data.num_inequality_);
        d_ << dx, vector_t::Zero(num_equality_constraints_ + data.num_inequality_);

//        diff_mat.pruned(1e-8); // doesn't make a big difference

        assert(diff_mat.rows() == diff_mat.cols());

        // Transpose method doesn't work well!
//        sp_matrix_t temp = diff_mat.transpose()*diff_mat;

//        temp.makeCompressed();
//        temp.pruned(1e-6);

        Eigen::SparseLU<sp_matrix_t> solver; // works really well!


        vector_t temp2 = d_;// diff_mat.transpose()*d_;

        utils::Timer lu_timer("sparse lu");

        lu_timer.StartTimer();
        solver.analyzePattern(diff_mat); // this is about 1-2ms of compute time

        solver.factorize(diff_mat); // diff_mat

        // TODO: Handle gracefully
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Could not factor the differential matrix.");
        }

        vector_t x = solver.solve(temp2);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Could not solve the equation.");
        }
        lu_timer.StopTimer();
        lu_timer.PrintElapsedTime();

//        std::cout << "Sparse LU error norm: " << (diff_mat*x - d_).norm() << std::endl;
//        std::cout << "Sparse LU norm: " << x.norm() << std::endl;

        d_ = -x;

        timer.StopTimer();

        timer.PrintElapsedTime();
    }

    vector_t ClarabelInterface::Computedx(const Eigen::SparseMatrix<double>& P, const mpc::vector_t& q,
                                      const mpc::vector_t& xstar) {
        dx_.noalias() = P*xstar + q;

        // TODO: This is making it only consider the state tracking costs
//        dx_.tail(dx_.size() - 12*45).setZero(); // not using all the states (only the first 45 nodes)

        return dx_;
    }

    ClarabelInterface& ClarabelInterface::operator=(const mpc::ClarabelInterface& other) {
        settings_ = other.settings_;
        cones_ = other.cones_;
        solve_quality_ = other.solve_quality_;
        num_equality_constraints_ = other.num_equality_constraints_;
        num_inequality_constraints_ = other.num_inequality_constraints_;
        dual_ = other.dual_;
        primal_ = other.primal_;
        slacks_ = other.slacks_;
        d_ = other.d_;
        dx_ = other.dx_;
        lam_ = other.lam_;
        nu_ = other.nu_;

        return *this;
    }

    ClarabelInterface::ClarabelInterface(const mpc::ClarabelInterface& other) : QPInterface(other.prev_qp_sol_.size()){
        *this = other;
    }

    void ClarabelInterface::SetVerbosity(bool verbose) {
        settings_.verbose = verbose;
    }

} // mpc