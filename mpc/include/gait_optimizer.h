//
// Created by zolkin on 1/15/24.
//

#ifndef BILEVEL_GAIT_GEN_GAIT_OPTIMIZER_H
#define BILEVEL_GAIT_GEN_GAIT_OPTIMIZER_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "OsqpEigen/OsqpEigen.h"
#include "sparse_matrix_builder.h"
#include "spline/end_effector_splines.h"

namespace mpc {
    using matrix_t = Eigen::MatrixXd;
    using vector_t = Eigen::VectorXd;
    using sp_matrix_t = Eigen::SparseMatrix<double>;

    struct QPPartials {
        sp_matrix_t dA;
        sp_matrix_t dP;
        sp_matrix_t dG;
        vector_t dl;
        vector_t du;
        vector_t dq;
        vector_t db;
        vector_t dh;

        void SetZero() {
            dP.setZero();
            dA.setZero();
            dG.setZero();
            dq.setZero();
            dl.setZero();
            du.setZero();
            db.setZero();
            dh.setZero();
        }
    };

    class GaitOptimizer {
    public:
        GaitOptimizer(int num_ee, int num_contact_nodes, int num_decision_vars,
                      int num_constraints, double contact_time_ub, double min_time);

        void Setdx(const vector_t& dx);

        QPPartials& GetPartials();

        /**
         * Computes the partial derivative of the cost function wrt all the QP parameters:
         * dP, dq, dA, dlu, dub
         *
         * @note Setdx MUST be called first.
         */
        void SetCostFcnPartials(const QPPartials& partials);


        /**
         * Sets all the partial derivatives of the QP parameters wrt the contact times.
         * @param dA
         * @param dP
         * @param dlb
         * @param dub
         * @param dq
         * @param ee end effector
         * @param contact_idx the index of the contact time that the partial is wrt to
         */
//        void SetParameterPartials(const sp_matrix_t& dA, const sp_matrix_t& dP, const vector_t& dlb,
//                                  const vector_t& dub, const vector_t& dq, int ee, int contact_idx);

        QPPartials& GetParameterPartials(int ee, int idx);

        /**
         * Computes the derivatives of the cost function wrt the contact times.
         *
         * @note All partials MUST be set prior to calling this function.
         */
        void ComputeCostFcnDerivWrtContactTimes();

        /**
         * Performs the optimization on the contact times.
         * In the current version, this is a gradient descent step and a projection back onto the polytope.
         *
         * @note ComputeCostFcnDerivWrtContactTimes MUST be called prior to this function.
         */
        void OptimizeContactTimes(double time, double actual_red_cost);

        /**
         * @return the current computed contact times
         */
        std::vector<time_v> GetContactTimes();

        void SetNumContactTimes(int ee, int num_times);

        void UpdateSizes(int num_decision_vars, int num_constraints);

        void SetContactTimes(const std::vector<time_v>& contact_times);

        void ModifyQPPartials(const vector_t& xstar);

    protected:
    private:
        int GetNumTimeNodes(int ee) const;

        int CreatePolytopeConstraint(int start_row);

        int CreateStepBoundConstraint(int start_row);

        int CreateStartConstraint(int start_row);

        int CreateTrustRegionConstraint(int start_row);

        void UpdateLagrangianGradients(const sp_matrix_t& A);

        void AdjustBSize(int num_decision_vars);

        void IncreaseTrustRegion(double rho);

        void DecreaseTrustRegion(const vector_t& step);

        int CreateNextNodeConstraints(int start_row, double time);

        void DampedBFGSUpdate();

        vector_t ContactTimesToQPVec() const;

        void PrintConstraints(const matrix_t& A, const vector_t& lb, const vector_t& ub);

        std::string GetSolveQualityAsString() const;

        std::vector<time_v> contact_times_;
        std::vector<time_v> old_contact_times_;
        std::vector<std::vector<double>> contact_times_lb_, contact_times_ub_;

        vector_t dHdth;

        std::vector<std::vector<QPPartials>> param_partials_;

//        std::vector<std::vector<sp_matrix_t>> dAdth;         // partial deriv of A wrt contact times
//        std::vector<std::vector<sp_matrix_t>> dPdth;         // partial deriv of P wrt contact times
//        std::vector<std::vector<vector_t>> dlbdth;        // partial deriv of lb wrt contact times
//        std::vector<std::vector<vector_t>> dubdth;        // partial deriv of ub wrt contact times
//        std::vector<std::vector<vector_t>> dqdth;         // partial deriv of q wrt contact times


        QPPartials qp_partials_;

//        sp_matrix_t dldA;          // partial deriv of cost fcn wrt A
//        sp_matrix_t dldP;          // partial deriv of cost fcn wrt P
//        vector_t dldq;          // partial deriv of cost fcn wrt q
//        vector_t dldlb;         // partial deriv of cost fcn wrt lb
//        vector_t dldub;         // partial deriv of cost fcn wrt ub
        vector_t dldx;          // partial deriv of cost fcn wrt x

        const int num_ee_;
//        int num_contact_nodes_; // TODO: I want to make this const
        int num_decision_vars_; // TODO: I want to make this const
        int num_constraints_;   // TODO: I want to make this const

        int past_decision_vars_;
//        vector_t old_grad_;

        double contact_time_ub_;
        double min_time_;

        int run_num_;

        OsqpEigen::Solver qp_solver_;

        utils::SparseMatrixBuilder A_builder_;
        vector_t lb_, ub_;

        vector_t xk_, xkp1_;
        vector_t gradk_, gradkp1_;
        vector_t dual_;
        matrix_t Bk_;
        vector_t step_;

        double pred_red_cost_;
        double gamma_; // Trust region shrinking scalar
        double eta_;    // Reduction threshold
        double Delta_;  // Trust region size
        const double max_trust_region_;
    };
}


#endif //BILEVEL_GAIT_GEN_GAIT_OPTIMIZER_H
