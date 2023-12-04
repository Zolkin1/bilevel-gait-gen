//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_QPDATA_H
#define BILEVEL_GAIT_GEN_QPDATA_H

#include <Eigen/Core>

namespace mpc {

    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;


    // mpc data, stored as a struct to be cache friendly
    // - gradients and hessians
    // - constant terms
    // - warm start
    // TODO: is this really the struct I want?
    struct QPData {
        matrix_t dynamics_constraints;
        vector_t dynamics_constants;

//        matrix_t equality_constraints;
//        vector_t equality_constants;
//
//        matrix_t inequality_constraints;
//        vector_t inequality_constants_ub;
//        vector_t inequality_constants_lb;

        matrix_t fk_constraints_;
        vector_t fk_constants_;

        matrix_t swing_force_constraints_;
        vector_t swing_force_constants_;

        matrix_t foot_on_ground_constraints_;
        vector_t foot_on_ground_lb_;
        vector_t foot_on_ground_ub_;

        matrix_t friction_cone_constraints_;
        vector_t friction_cone_ub_;
        vector_t friction_cone_lb_;

        matrix_t positive_force_constraints_;
        vector_t positive_force_lb_;
        vector_t positive_force_ub_;

        matrix_t foot_ground_inter_constraints_;
        vector_t foot_ground_inter_lb_;
        vector_t foot_ground_inter_ub_;

        matrix_t box_constraints_;
        vector_t box_lb_;
        vector_t box_ub_;

        matrix_t cost_quadratic;
        vector_t cost_linear;

        int num_dynamics_constraints;
        int num_equality_constraints;
        int num_inequality_constraints;
        int num_decision_vars;

        int num_cone_constraints_;
        int num_box_constraints_;
        int num_positive_force_constraints_;
        int num_foot_ground_inter_constraints_;
        int num_foot_on_ground_constraints_;
        int num_fk_constraints_;
        int num_swing_foot_constraints_;
    };
} // mpc

#endif //BILEVEL_GAIT_GEN_QPDATA_H
