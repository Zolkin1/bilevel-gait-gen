//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_QPDATA_H
#define BILEVEL_GAIT_GEN_QPDATA_H

#include <Eigen/Core>

#include "sparse_matrix_builder.h"

namespace mpc {

    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;

    enum Constraints {
        Dynamics,
        JointForwardKinematics,
        EndEffectorLocation,
        ForceBox,
        JointBox,
        FrictionCone
    };

    struct QPSettings {
        const bool constraint_projection_;
        const int constraint_mat_nnz_;
        const int cost_mat_nnz_;

        QPSettings(bool constraint_projection, int constraint_mat_nnz, int cost_mat_nnz);
    };

    struct QPConstraintProjection {
        matrix_t constraint_mat_;
        vector_t b_;

        matrix_t null_space_mat_;
        vector_t p_;

        QPConstraintProjection(int rows, int cols);

        void CalculateConstraintProjection();
    };

    struct QPData {
        QPSettings settings_;
        QPConstraintProjection constraint_projections_;
        const std::vector<Constraints> constraints_;

        // Eigen triplet for filling the sparse matrix
        utils::SparseMatrixBuilder constraint_mat_;
        utils::SparseMatrixBuilder cost_mat_;

        Eigen::SparseMatrix<double> sparse_constraint_;
        Eigen::SparseMatrix<double> sparse_cost_;

        vector_t lb_;
        vector_t ub_;

        vector_t dynamics_constants;

        vector_t fk_constants_;

        vector_t fk_lb_;
        vector_t fk_ub_;

        vector_t friction_cone_ub_;
        vector_t friction_cone_lb_;

        vector_t box_lb_;
        vector_t box_ub_;

        vector_t cost_linear;

        vector_t force_box_lb_;
        vector_t force_box_ub_;

        vector_t ee_location_lb_;
        vector_t ee_location_ub_;

        vector_t start_ee_constants_;

        int num_dynamics_constraints;
        int num_decision_vars;

        int num_cone_constraints_;
        int num_box_constraints_;
        int num_fk_constraints_;
        int num_force_box_constraints_;
        int num_fk_ineq_constraints_;
        int num_ee_location_constraints_;
        int num_start_ee_constraints_;

        QPData();
        QPData(bool constraint_projection, int constraint_mat_nnz, int cost_mat_nnz,
               std::vector<Constraints> constraints);
        QPData(const QPSettings settings, std::vector<Constraints> constraints);
        QPData(const QPSettings settings, int projection_rows, int projection_cols,
               std::vector<Constraints> constraints);

        int GetTotalNumConstraints() const;

        void InitQPMats();

        void ConstructSparseMats();

        void ConstructVectors();

        void ApplyProjection();
    };
} // mpc

#endif //BILEVEL_GAIT_GEN_QPDATA_H
