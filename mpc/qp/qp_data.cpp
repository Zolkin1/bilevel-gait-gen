//
// Created by zolkin on 1/11/24.
//

#include <Eigen/Core>
#include <Eigen/QR>
#include "qp/qp_data.h"

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

    // Note: Uses all constraints if none are provided
    QPData::QPData() :
        settings_(false, 25000, 2000),
        constraint_projections_(0,0),
        lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000),
        constraints_({Constraints::Dynamics, Constraints::JointForwardKinematics, Constraints::EndEffectorLocation,
                      Constraints::ForceBox, Constraints::JointBox, Constraints::FrictionCone}),
                      using_clarabel_(false) {}

    QPData::QPData(bool constraint_projection, int constraint_mat_nnz, int cost_mat_nnz,
                   std::vector<Constraints> constraints) :
        settings_(constraint_projection, constraint_mat_nnz, cost_mat_nnz), constraint_projections_(10, 100),
        lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000),
        constraints_(constraints), using_clarabel_(false) {}

    QPData::QPData(const mpc::QPSettings settings, std::vector<Constraints> constraints) : settings_(settings),
        constraint_projections_(10, 100),
        lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000),
        constraints_(constraints), using_clarabel_(false) {}

    QPData::QPData(const QPSettings settings, int projection_rows, int projection_cols, std::vector<Constraints> constraints) :
        settings_(settings), constraint_projections_(projection_rows, projection_cols),
        lb_(), ub_(), sparse_constraint_(10000,10000), sparse_cost_(10000,10000),
        constraints_(constraints), using_clarabel_(false) {}

    int QPData::GetTotalNumConstraints() const {
        int num_constraints = 0;
        for (int i = 0; i < constraints_.size(); i++) {
            switch (constraints_.at(i)) {
                case Constraints::Dynamics:
                    num_constraints += num_dynamics_constraints;
                    break;
                case Constraints::JointForwardKinematics:
                    num_constraints += num_fk_constraints_ + num_fk_ineq_constraints_;
                    break;
                case Constraints::EndEffectorLocation:
                    num_constraints += num_ee_location_constraints_ + num_start_ee_constraints_;
                    break;
                case Constraints::ForceBox:
                    num_constraints += num_force_box_constraints_;
                    break;
                case Constraints::JointBox:
                    num_constraints += num_box_constraints_;
                    break;
                case Constraints::FrictionCone:
                    num_constraints += num_cone_constraints_;
                    break;
            }
        }

        return num_constraints;

    }

    void QPData::InitQPMats() {
        constraint_mat_.Reserve(settings_.constraint_mat_nnz_);
        cost_mat_.Reserve(settings_.cost_mat_nnz_);

        for (int i = 0; i < constraints_.size(); i++) {
            switch (constraints_.at(i)) {
                case Constraints::Dynamics:
                    dynamics_constants = vector_t::Zero(num_dynamics_constraints);
                    break;
                case Constraints::JointForwardKinematics:
                    if (!settings_.constraint_projection_) {
                        fk_constants_ = vector_t::Zero(num_fk_constraints_);
                        fk_lb_ = vector_t::Zero(num_fk_ineq_constraints_);
                        fk_ub_ = vector_t::Zero(num_fk_ineq_constraints_);
                    } else {
                        fk_lb_ = vector_t::Zero(num_fk_ineq_constraints_);
                        fk_ub_ = vector_t::Zero(num_fk_ineq_constraints_);
                        constraint_projections_.constraint_mat_.resize(num_fk_constraints_, num_decision_vars);
                        constraint_projections_.null_space_mat_.resize(num_decision_vars - num_fk_constraints_, num_decision_vars);
                        constraint_projections_.p_.resize(num_decision_vars);
                        constraint_projections_.b_.resize(num_fk_constraints_);
                    }
                    break;
                case Constraints::EndEffectorLocation:
                    if (!using_clarabel_) {
                        ee_location_lb_ = vector_t::Zero(num_ee_location_constraints_);
                        ee_location_ub_ = vector_t::Zero(num_ee_location_constraints_);
                    } else {
                        ee_location_lb_ = vector_t::Zero(num_ee_location_constraints_/2);
                        ee_location_ub_ = vector_t::Zero(num_ee_location_constraints_/2);
                    }
                    start_ee_constants_ = vector_t::Zero(num_start_ee_constraints_);
                    break;
                case Constraints::ForceBox:
                    if (!using_clarabel_) {
                        force_box_lb_ = vector_t::Zero(num_force_box_constraints_);
                        force_box_ub_ = vector_t::Zero(num_force_box_constraints_);
                    } else {
                        force_box_lb_ = vector_t::Zero(num_force_box_constraints_/2);
                        force_box_ub_ = vector_t::Zero(num_force_box_constraints_/2);
                    }
                    break;
                case Constraints::JointBox:
                    box_lb_ = vector_t::Zero(num_box_constraints_);
                    box_ub_ = vector_t::Zero(num_box_constraints_);
                    break;
                case Constraints::FrictionCone:
                    friction_cone_lb_ = vector_t::Zero(num_cone_constraints_);
                    friction_cone_ub_ = vector_t::Zero(num_cone_constraints_);
                    break;
            }
        }


        cost_linear = vector_t::Zero(num_decision_vars);

        lb_.resize(GetTotalNumConstraints());
        lb_.setZero(); // TODO: should not need this
        ub_.resize(GetTotalNumConstraints());
        ub_.setZero(); // TODO: should not need this
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
        int idx = 0;
        num_inequality_ = 0;
        num_equality_ = 0;
        for (int i = 0; i < constraints_.size(); i++) {
            switch (constraints_.at(i)) {
                case Constraints::Dynamics:
                    if (!using_clarabel_) {
                        lb_.segment(idx, num_dynamics_constraints) = dynamics_constants;
                        ub_.segment(idx, num_dynamics_constraints) = dynamics_constants;
                    } else {
                        ub_.segment(idx, num_dynamics_constraints) = dynamics_constants;
                    }
                    idx += num_dynamics_constraints;
                    num_equality_ += num_dynamics_constraints;
                    break;
                case Constraints::JointForwardKinematics:
                    if (!settings_.constraint_projection_) {
                        lb_.segment(idx, num_fk_ineq_constraints_) = fk_lb_;
                        ub_.segment(idx, num_fk_ineq_constraints_) = fk_ub_;
                        idx += num_fk_ineq_constraints_;
                        lb_.segment(idx, num_fk_constraints_) = fk_constants_;
                        ub_.segment(idx, num_fk_constraints_) = fk_constants_;
                        idx += num_fk_constraints_;
                    } else {
                        lb_.segment(idx, num_fk_ineq_constraints_) = fk_lb_;
                        ub_.segment(idx, num_fk_ineq_constraints_) = fk_ub_;
                        idx += num_fk_ineq_constraints_;
                    }
                    // TODO: inequality/equality
                    break;
                case Constraints::EndEffectorLocation:
                    if (!using_clarabel_) {
                        lb_.segment(idx, num_ee_location_constraints_) = ee_location_lb_;
                        ub_.segment(idx, num_ee_location_constraints_) = ee_location_ub_;
                        idx += num_ee_location_constraints_;
                        lb_.segment(idx, num_start_ee_constraints_) = start_ee_constants_;
                        ub_.segment(idx, num_start_ee_constraints_) = start_ee_constants_;
                        idx += num_start_ee_constraints_;
                    } else {
                        ub_.segment(idx, num_ee_location_constraints_) << ee_location_ub_, -ee_location_lb_;
                        idx += num_ee_location_constraints_;
                        ub_.segment(idx, num_start_ee_constraints_) = start_ee_constants_;
                        idx += num_start_ee_constraints_;
                    }
                    num_inequality_ += num_ee_location_constraints_;
                    num_equality_ += num_start_ee_constraints_;
                    break;
                case Constraints::ForceBox:
                    if (!using_clarabel_) {
                        lb_.segment(idx, num_force_box_constraints_) = force_box_lb_;
                        ub_.segment(idx, num_force_box_constraints_) = force_box_ub_;
                        idx += num_force_box_constraints_;
                    } else {
                        ub_.segment(idx, num_force_box_constraints_) << force_box_ub_, -force_box_lb_;
                        idx += num_force_box_constraints_;
                    }
                    num_inequality_ += num_force_box_constraints_;
                    break;
                case Constraints::JointBox:
                    lb_.segment(idx, num_box_constraints_) = box_lb_;
                    ub_.segment(idx, num_box_constraints_) = box_ub_;
                    idx += num_box_constraints_;
                    num_inequality_ += num_box_constraints_;
                    break;
                case Constraints::FrictionCone:
                    if (!using_clarabel_) {
                        lb_.segment(idx, num_cone_constraints_) = friction_cone_lb_;
                        ub_.segment(idx, num_cone_constraints_) = friction_cone_ub_;
                    } else {
                        ub_.segment(idx, num_cone_constraints_) = friction_cone_ub_;
                    }
                    idx += num_cone_constraints_;
                    num_inequality_ += num_cone_constraints_;
                    break;
            }
        }
    }

}