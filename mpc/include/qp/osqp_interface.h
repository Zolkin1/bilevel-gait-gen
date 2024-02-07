//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
#define BILEVEL_GAIT_GEN_OSQP_INTERFACE_H

#include <Eigen/SparseCore>

#include "OsqpEigen/OsqpEigen.h"

#include "qp_interface.h"
#include "timer.h"

namespace mpc {
    class OSQPInterface {
    public:
        enum SolveQuality {
            Solved = 0,
            SolvedInacc = 1,
            MaxIter = 2,
            PrimalInfeasible = 3,
            DualInfeasible = 4,
            PrimalInfeasibleInacc = 5,
            DualInfeasibleInacc = 6,
            Unsolved = 7,
            Other = 8
        };

        OSQPInterface(const QPData& data, bool verbose);

        void SetupQP(QPData& data, const vector_t& warm_start);

        // TODO: remove data
        vector_t Solve(const QPData& data);

        vector_t GetInfinity(int size) const;

        SolveQuality GetSolveQuality() const;

        std::string GetSolveQualityAsString() const;

        vector_t GetDualSolution() const;

        void ConfigureForInitialRun() const;

        void ConfigureForRealTime(double run_time_iters) const;

        void SetupDerivativeCalcs(vector_t& dx, vector_t& dy_l, vector_t& dy_u);

        void CalcDerivativeWrtMats(Eigen::SparseMatrix<double>& dP, Eigen::SparseMatrix<double>& dA);

        void CalcDerivativeWrtVecs(vector_t& dq, vector_t dl, vector_t du);

        void Computedx(const Eigen::SparseMatrix<double>& P, const vector_t& q, const vector_t& xstar);

        void SetSolveTolerances(double abs_tol, double rel_tol);

        /**
         * Calculates the partial of the quadratic cost fcn wrt the decision variables, x.
         * Cost fcn of the form (1/2)*x^T*P*x + q^T*x
         * @param P cost fcn hessian term
         * @param q cost fcn gradient term
         * @param xstar optimal solution to the QP
         */
        vector_t Getdx() const;

    protected:
    private:
//        void ConvertDataToOSQPConstraints(const QPData& data);

//        void ConvertDataToOSQPCost(const QPData& data);

        // QP Solver
        OsqpEigen::Solver qp_solver_;

        vector_t prev_dual_sol_;
        vector_t prev_qp_sol_;

        vector_t dldx;

        int run;

        bool verbose_;

        utils::Timer osqp_interface_timer_;
        utils::Timer sparse_conversion_timer_;
        utils::Timer matrix_gather_timer_;

        int prev_num_constraints_;
        int prev_num_decision_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
