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
#include "end_effector_splines.h"

namespace mpc {
    using matrix_t = Eigen::MatrixXd;
    using vector_t = Eigen::VectorXd;
    using sp_matrix_t = Eigen::SparseMatrix<double>;

    struct QPPartials {
        sp_matrix_t dA;
        sp_matrix_t dP;
        vector_t dl;
        vector_t du;
        vector_t dq;
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
        void OptimizeContactTimes();

        /**
         * @return the current computed contact times
         */
        std::vector<time_v> GetContactTimes();

        void SetNumContactTimes(int ee, int num_times);

        void UpdateSizes(int num_decision_vars, int num_constraints);

        void SetContactTimes(std::vector<time_v> contact_times);

        void ModifyQPPartials(const vector_t& xstar);

    protected:
    private:
        int GetNumTimeNodes(int ee) const;

        int CreatePolytopeConstraint(int start_row);

        int CreateStepBoundConstraint(int start_row);

        void PrintConstraints(const matrix_t& A, const vector_t& lb, const vector_t& ub);

        std::string GetSolveQualityAsString() const;

        std::vector<time_v> contact_times_;
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

        double contact_time_ub_;
        double min_time_;

        OsqpEigen::Solver qp_solver_;

        utils::SparseMatrixBuilder A_builder_;
        vector_t lb_, ub_;
    };
}


#endif //BILEVEL_GAIT_GEN_GAIT_OPTIMIZER_H
