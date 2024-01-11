//
// Created by zolkin on 1/11/24.
//

#include <Eigen/Core>
#include <Eigen/QR>
#include "qp_data.h"

namespace mpc {
    QPSettings::QPSettings(bool constraint_projection, int constraint_mat_nnz, int cost_mat_nnz)
    : constraint_projection_(constraint_projection),
    constraint_mat_nnz_(constraint_mat_nnz),
    cost_mat_nnz_(cost_mat_nnz) {}

    QPConstraintProjection::QPConstraintProjection(int num_constraints, int num_inputs) {
        constraint_mat_.resize(num_constraints, num_inputs);
        b_.resize(num_constraints);

        null_space_mat_.resize(num_inputs - num_constraints, num_inputs);
        p_.resize(num_inputs - num_constraints);
    }

    // TODO: Fix
    void QPConstraintProjection::CalculateConstraintProjection() {
        const int num_constraints = constraint_mat_.rows();
        const int num_inputs = constraint_mat_.cols();

        const Eigen::HouseholderQR<matrix_t> QR(constraint_mat_.transpose());
        const matrix_t Q = QR.householderQ();
        null_space_mat_ = Q.rightCols(num_inputs - num_constraints);
        const auto R = QR.matrixQR().topRows(num_constraints).triangularView<Eigen::Upper>();
        const matrix_t pinvA = R.solve(Q.leftCols(num_constraints).transpose());
        p_ = pinvA.transpose() * b_;
    }

    QPData::QPData() :
    settings_(false, 25000, 2000),
    constraint_projections_(0,0),
    lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000) {
    }

    QPData::QPData(bool constraint_projection, int constraint_mat_nnz, int cost_mat_nnz) :
    settings_(constraint_projection, constraint_mat_nnz, cost_mat_nnz), constraint_projections_(10, 100),
    lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000) {
    }

    QPData::QPData(const mpc::QPSettings settings) : settings_(settings),
    constraint_projections_(10, 100),
    lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000){}

    QPData::QPData(const QPSettings settings, int projection_rows, int projection_cols) :
            settings_(settings), constraint_projections_(projection_rows, projection_cols),
            lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000) {}

    int QPData::GetTotalNumConstraints() const {
        if (!settings_.constraint_projection_) {
            return num_cone_constraints_ + num_box_constraints_ +
                   num_force_box_constraints_ + num_dynamics_constraints +
                   num_fk_constraints_ + num_fk_ineq_constraints_;
        } else {
            return num_cone_constraints_ + num_box_constraints_ +
                   num_force_box_constraints_ + num_dynamics_constraints +
                   + num_fk_ineq_constraints_;
        }

    }

    void QPData::InitQPMats() {
        constraint_mat_.Reserve(settings_.constraint_mat_nnz_);
        cost_mat_.Reserve(settings_.cost_mat_nnz_);

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

        lb_.resize(GetTotalNumConstraints());
        ub_.resize(GetTotalNumConstraints());

        constraint_projections_.constraint_mat_.resize(num_fk_constraints_, num_decision_vars);
        constraint_projections_.null_space_mat_.resize(num_decision_vars - num_fk_constraints_, num_decision_vars);
        constraint_projections_.p_.resize(num_decision_vars);
        constraint_projections_.b_.resize(num_fk_constraints_);
    }

    void QPData::ConstructSparseMats() {
        if (sparse_constraint_.rows() != GetTotalNumConstraints() || sparse_constraint_.cols() != num_decision_vars) {
            sparse_constraint_.resize(GetTotalNumConstraints(), num_decision_vars);
        }
        if (sparse_cost_.rows() != num_decision_vars || sparse_cost_.cols() != num_decision_vars) {
            sparse_cost_.resize(num_decision_vars, num_decision_vars);
        }
        sparse_constraint_.setFromTriplets(constraint_mat_.GetTriplet().begin(), constraint_mat_.GetTriplet().end());
        sparse_cost_.setFromTriplets(cost_mat_.GetTriplet().begin(), cost_mat_.GetTriplet().end());
    }

    // TODO: Fix
    void QPData::ApplyProjection() {
        constraint_projections_.CalculateConstraintProjection();

        cost_linear = (cost_linear.transpose()*constraint_projections_.null_space_mat_).transpose();
                        //+ constraint_projections_.p_.transpose()*sparse_cost_*constraint_projections_.null_space_mat_).transpose();

        lb_ -= sparse_constraint_*constraint_projections_.p_;
        ub_ -= sparse_constraint_*constraint_projections_.p_;

        constraint_mat_.Reserve(settings_.constraint_mat_nnz_);
        constraint_mat_.SetMatrix(sparse_constraint_ * constraint_projections_.constraint_mat_.transpose(), 0, 0);

        cost_mat_.Reserve(settings_.cost_mat_nnz_);
        cost_mat_.SetMatrix(constraint_projections_.null_space_mat_.transpose() *
            sparse_cost_ * constraint_projections_.null_space_mat_, 0, 0);

        ConstructSparseMats();
    }

    void QPData::ConstructVectors() {
        if (!settings_.constraint_projection_) {
            lb_ << dynamics_constants,
                    fk_lb_,
                    fk_constants_,
                    friction_cone_lb_,
                    box_lb_,
                    force_box_lb_;

            ub_ << dynamics_constants,
                    fk_ub_,
                    fk_constants_,
                    friction_cone_ub_,
                    box_ub_,
                    force_box_ub_;
        } else {
            lb_ << dynamics_constants,
                    fk_lb_,
                    friction_cone_lb_,
                    box_lb_,
                    force_box_lb_;

            ub_ << dynamics_constants,
                    fk_ub_,
                    friction_cone_ub_,
                    box_ub_,
                    force_box_ub_;
        }
    }

}