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

        matrix_t fk_constraints_;
        vector_t fk_constants_;

        matrix_t fk_ineq_constraints_;
        vector_t fk_lb_;
        vector_t fk_ub_;

        matrix_t friction_cone_constraints_;
        vector_t friction_cone_ub_;
        vector_t friction_cone_lb_;

        matrix_t box_constraints_;
        vector_t box_lb_;
        vector_t box_ub_;

        matrix_t cost_quadratic;
        vector_t cost_linear;

        matrix_t force_box_constraints_;
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
            dynamics_constraints = matrix_t::Zero(num_dynamics_constraints, num_decision_vars);
            dynamics_constants = vector_t::Zero(num_dynamics_constraints);

            cost_quadratic = matrix_t::Zero(num_decision_vars, num_decision_vars);
            cost_linear = vector_t::Zero(num_decision_vars);

            fk_constraints_ = matrix_t::Zero(num_fk_constraints_, num_decision_vars);
            fk_constants_ = vector_t::Zero(num_fk_constraints_);

            fk_ineq_constraints_ = matrix_t::Zero(num_fk_ineq_constraints_, num_decision_vars);
            fk_lb_ = vector_t::Zero(num_fk_ineq_constraints_);
            fk_ub_ = vector_t::Zero(num_fk_ineq_constraints_);

            friction_cone_constraints_ = matrix_t::Zero(num_cone_constraints_, num_decision_vars);
            friction_cone_lb_ = vector_t::Zero(num_cone_constraints_);
            friction_cone_ub_ = vector_t::Zero(num_cone_constraints_);

            box_constraints_ = matrix_t::Zero(num_box_constraints_, num_decision_vars);
            box_lb_ = vector_t::Zero(num_box_constraints_);
            box_ub_ = vector_t::Zero(num_box_constraints_);

            force_box_constraints_ = matrix_t::Zero(num_force_box_constraints_, num_decision_vars);
            force_box_lb_ = vector_t::Zero(num_force_box_constraints_);
            force_box_ub_ = vector_t::Zero(num_force_box_constraints_);
        }

        void ResizeQPMats() {
            dynamics_constraints.resize(num_dynamics_constraints, num_decision_vars);
            dynamics_constants.resize(num_dynamics_constraints);

            cost_quadratic.resize(num_decision_vars, num_decision_vars);
            cost_linear.resize(num_decision_vars);

            fk_constraints_.resize(num_fk_constraints_, num_decision_vars);
            fk_constants_.resize(num_fk_constraints_);

            fk_ineq_constraints_.resize(num_fk_ineq_constraints_, num_decision_vars);
            fk_lb_.resize(num_fk_ineq_constraints_);
            fk_ub_.resize(num_fk_ineq_constraints_);

            friction_cone_constraints_.resize(num_cone_constraints_, num_decision_vars);
            friction_cone_lb_.resize(num_cone_constraints_);
            friction_cone_ub_.resize(num_cone_constraints_);

            box_constraints_.resize(num_box_constraints_, num_decision_vars);
            box_lb_.resize(num_box_constraints_);
            box_ub_.resize(num_box_constraints_);

            force_box_constraints_.resize(num_force_box_constraints_, num_decision_vars);
            force_box_lb_.resize(num_force_box_constraints_);
            force_box_ub_.resize(num_force_box_constraints_);
        }
    };
} // mpc

#endif //BILEVEL_GAIT_GEN_QPDATA_H
