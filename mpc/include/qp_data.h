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


    // mpc data, stored as a struct to be cache friendly
    // - gradients and hessians
    // - constant terms
    // - warm start
    // TODO: is this really the struct I want?
    struct QPData {
        // NOTE: I ultimately want sparse matricies, so I will just store triplets here
        // and construct the sparse matrix later

        // Eigen triplet for filling the sparse matrix
        utils::SparseMatrixBuilder constraint_mat_;
        utils::SparseMatrixBuilder cost_mat_;

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

        int num_dynamics_constraints;
        int num_decision_vars;

        int num_cone_constraints_;
        int num_box_constraints_;
        int num_fk_constraints_;
        int num_force_box_constraints_;
        int num_fk_ineq_constraints_;

        int GetTotalNumConstraints() const {
            return num_cone_constraints_ + num_box_constraints_ +
            num_fk_constraints_ + num_force_box_constraints_ +
            num_dynamics_constraints + num_fk_ineq_constraints_;
        }

        void InitQPMats() {
            constraint_mat_.Reserve(18000);     // TODO: Make not hard coded
            cost_mat_.Reserve(2000);            // TODO: Make not hard coded

            dynamics_constants = vector_t::Zero(num_dynamics_constraints);

            cost_linear = vector_t::Zero(num_decision_vars);

            fk_constants_ = vector_t::Zero(num_fk_constraints_);

            fk_lb_ = vector_t::Zero(num_fk_ineq_constraints_);
            fk_ub_ = vector_t::Zero(num_fk_ineq_constraints_);

            friction_cone_lb_ = vector_t::Zero(num_cone_constraints_);
            friction_cone_ub_ = vector_t::Zero(num_cone_constraints_);

            box_lb_ = vector_t::Zero(num_box_constraints_);
            box_ub_ = vector_t::Zero(num_box_constraints_);

            force_box_lb_ = vector_t::Zero(num_force_box_constraints_);
            force_box_ub_ = vector_t::Zero(num_force_box_constraints_);
        }
    };
} // mpc

#endif //BILEVEL_GAIT_GEN_QPDATA_H
