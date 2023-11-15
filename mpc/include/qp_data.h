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

        matrix_t equality_constraints;
        vector_t equality_constants;

        matrix_t inequality_constraints;
        vector_t inequality_constants_ub;
        vector_t inequality_constants_lb;

        matrix_t cost_quadratic;
        vector_t cost_linear;

        int num_dynamics_constraints;
        int num_equality_constraints;
        int num_inequality_constraints;
        int num_decision_vars;
    };
} // mpc

#endif //BILEVEL_GAIT_GEN_QPDATA_H
