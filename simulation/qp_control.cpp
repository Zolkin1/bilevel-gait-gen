//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <Eigen/Sparse>

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/frames.hpp"

#include "qp_control.h"

namespace simulator {
    QPControl::QPControl(double control_rate, std::string robot_urdf, const std::string &foot_type, int nv,
                         const Eigen::VectorXd& torque_bounds) :
                         Controller(control_rate, robot_urdf, foot_type), num_vel_(nv) {

        torque_bounds_ = torque_bounds;
        // TODO: Adjust any solver settings
        //qp_solver_.settings()->setVerbosity(false);
    }

    void QPControl::InitSolver(const mjData* data) {
        // assume each joint is actuated
        num_actuators_ = num_vel_ - FLOATING_VEL_OFFSET;

        num_decision_vars_ = num_vel_ + CONSTRAINT_PER_FOOT*GetNumContacts();

        qp_solver_.data()->setNumberOfVariables(num_decision_vars_);

        num_constraints_ = num_vel_ + 7*GetNumContacts() + num_actuators_;
        qp_solver_.data()->setNumberOfConstraints(num_constraints_);

        A_ = Eigen::MatrixXd::Zero(num_constraints_, num_decision_vars_);
        lb_ = Eigen::VectorXd::Zero(num_constraints_);
        ub_ = lb_;

        P_ = Eigen::MatrixXd::Zero(num_decision_vars_, num_decision_vars_);
        w_ = Eigen::VectorXd::Zero(num_decision_vars_);

        const Eigen::VectorXd q = ConvertMujocoConfigToPinocchio(data);
        const Eigen::VectorXd v = ConvertMujocoVelToPinocchio(data);

        // Add the constraints to the matricies
        ComputeDynamicsTerms(q, v);
        AddDynamicsConstraints(v);
        AddContactMotionConstraints(v);
        AddTorqueConstraints(v);

        // TODO
//        AddFrictionConeConstraints();

        // Init solver constraints
        Eigen::SparseMatrix<double> sparseA = A_.sparseView();
        if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
              qp_solver_.data()->setLowerBound(lb_) &&
              qp_solver_.data()->setUpperBound(ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // TODO
        // Init the cost matricies
//        AddLegTrackingCost();
//        AddTorsoCost();
//        AddForceTrackingCost();

        // Init solver costs
        Eigen::SparseMatrix<double> sparseP = P_.sparseView();
        if (!(qp_solver_.data()->setHessianMatrix(sparseP) && qp_solver_.data()->setGradient(w_))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }

        // Init solver
        if (!qp_solver_.initSolver()) {
            throw std::runtime_error("Unable to initialize the solver.");
        }
    }

    std::vector<mjtNum> QPControl::ComputeControlAction(const mjModel* model, const mjData* data) {
        UpdateConstraintsAndCost(data);

        // Solve qp
        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            std::cerr << "Could not solve WBC QP. Returning 0 control action." << std::endl;
            std::vector<mjtNum> temp;
            for (int i = 0; i < model->njnt; i++) {
                temp.at(i) = 0;
            }
            return temp;
        }

        Eigen::VectorXd qp_sol = qp_solver_.getSolution();
        std::vector<mjtNum> control_action = RecoverControlInputs(qp_sol);

        return control_action;
    }

    void QPControl::ComputeDynamicsTerms(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);

        // Get M
        pinocchio::crba(pin_model_, *pin_data_,  q);
        // Make M symmetric
        pin_data_->M.triangularView<Eigen::StrictlyLower>() =
                pin_data_->M.transpose().triangularView<Eigen::StrictlyLower>();

        // Get C
        pinocchio::computeCoriolisMatrix(pin_model_, *pin_data_, q, v);

        // Get g
        pinocchio::computeGeneralizedGravity(pin_model_, *pin_data_, q);

        // Js
        Js_ = Eigen::MatrixXd::Zero(Js_.rows(), Js_.cols());
        Js_ = GetConstraintJacobian(q);

        // Jsdot
        Jsdot_ = Eigen::MatrixXd::Zero(Jsdot_.rows(), Jsdot_.cols());
        Jsdot_ = GetConstraintJacobianDerivative(q, v);

    }

    void QPControl::AddDynamicsConstraints(const Eigen::VectorXd& v) {
        // Add equality constraints to the solver params
        lb_.head(FLOATING_VEL_OFFSET) = -(pin_data_->g.head(FLOATING_VEL_OFFSET) +
                (pin_data_->C*v).head(FLOATING_VEL_OFFSET));
        ub_.head(FLOATING_VEL_OFFSET) = lb_.head(FLOATING_VEL_OFFSET);

        A_.topLeftCorner(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET) =
                pin_data_->M.topLeftCorner(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET);

        A_.topRightCorner(FLOATING_VEL_OFFSET, CONSTRAINT_PER_FOOT*GetNumContacts()) =
                -Js_.transpose().topLeftCorner(FLOATING_VEL_OFFSET, Js_.rows());

    }

    void QPControl::AddContactMotionConstraints(const Eigen::VectorXd& v) {
        // Assign to solver params
        if (Js_.size() > 0) {
            lb_.segment(FLOATING_VEL_OFFSET, Js_.rows()) = -Jsdot_ * v;
            ub_.segment(FLOATING_VEL_OFFSET, Js_.rows()) = lb_.segment(Js_.rows(), FLOATING_VEL_OFFSET);

            A_.block(FLOATING_VEL_OFFSET, 0, Js_.rows(), Js_.cols()) = Js_;
        }
    }

    void QPControl::AddTorqueConstraints(const Eigen::VectorXd& v) {
        // TODO: Double check. I believe I just use the actuated degrees here

        if (Js_.size() > 0) {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), 0, num_decision_vars_,
                     num_actuators_) <<
                                     pin_data_->M.block(FLOATING_VEL_OFFSET, 0, num_actuators_, pin_data_->M.cols()),
                    -Js_.transpose().block(FLOATING_VEL_OFFSET, 0, num_actuators_, Js_.rows());
        } else {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), 0, num_actuators_,
                     num_decision_vars_) << pin_data_->M.block(FLOATING_VEL_OFFSET, 0, num_actuators_, pin_data_->M.cols());
        }

        lb_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), num_actuators_) =
                (pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) - torque_bounds_;
        ub_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), num_actuators_) =
                -lb_.segment( FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), num_actuators_);
    }

    // TODO: Can I anything on the update? Maybe some of the costs...
    void QPControl::UpdateConstraintsAndCost(const mjData* data) {
        num_decision_vars_ = num_vel_ + CONSTRAINT_PER_FOOT*GetNumContacts();

        qp_solver_.data()->setNumberOfVariables(num_decision_vars_);

        num_constraints_ = num_vel_ + 7*GetNumContacts() + num_actuators_;
        qp_solver_.data()->setNumberOfConstraints(num_constraints_);

        A_ = Eigen::MatrixXd::Zero(num_constraints_, num_decision_vars_);
        lb_ = Eigen::VectorXd::Zero(num_constraints_);
        ub_ = lb_;

        P_ = Eigen::MatrixXd::Zero(num_decision_vars_, num_decision_vars_);
        w_ = Eigen::VectorXd::Zero(num_decision_vars_);

        const Eigen::VectorXd q = ConvertMujocoConfigToPinocchio(data);
        const Eigen::VectorXd v = ConvertMujocoVelToPinocchio(data);

        // Add the constraints to the matricies
        ComputeDynamicsTerms(q, v);
        AddDynamicsConstraints(v);
        AddContactMotionConstraints(v);
        AddTorqueConstraints(v);

//        AddFrictionConeConstraints();

        // Init solver constraints
        Eigen::SparseMatrix<double> sparseA = A_.sparseView();
        if (!(qp_solver_.updateLinearConstraintsMatrix(sparseA) &&
              qp_solver_.updateBounds(lb_, ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // Init the cost matricies
//        AddLegTrackingCost();
//        AddTorsoCost();
//        AddForceTrackingCost();

        // Init solver costs
        Eigen::SparseMatrix<double> sparseP = P_.sparseView();
        if (!(qp_solver_.updateHessianMatrix(sparseP) && qp_solver_.updateGradient(w_))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }
    }

    // TODO: Implement
    std::vector<double> QPControl::RecoverControlInputs(const Eigen::VectorXd& qp_sol) {
        std::vector<double> control_input;
        for (int i = 0; i < num_actuators_; i++) {
            control_input.push_back(0);
        }
    }

    Eigen::MatrixXd QPControl::GetConstraintJacobian(const Eigen::VectorXd& q) {
        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);

        // Js is the "support Jacobian". Js \in R^{constraint_per_foot*ncontacts x nv}
        // Js = [Jc^T ... Jc^T]^T, Jc \in R^{nv x constraint_per_foot}

        int num_contacts = GetNumContacts();
        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact_frames_.size(); i++) {
            if (in_contact_.at(i)) {
                int frame_id = pin_model_.getFrameId(pin_model_.frames.at(contact_frames_.at(i)).name,
                                                     pin_model_.frames.at(contact_frames_.at(i)).type);

                pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_id, pinocchio::LOCAL_WORLD_ALIGNED, J);
                Js.block(j * CONSTRAINT_PER_FOOT, 0, CONSTRAINT_PER_FOOT, J.cols()) =
                        J.topLeftCorner(CONSTRAINT_PER_FOOT, pin_model_.nv);
                j++;
                J = Eigen::MatrixXd::Zero(J.rows(), J.cols());
            }
        }

        return Js;
    }

    Eigen::MatrixXd QPControl::GetConstraintJacobianDerivative(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        pinocchio::computeJointJacobiansTimeVariation(pin_model_, *pin_data_, q, v);

        int num_contacts = GetNumContacts();
        Eigen::MatrixXd Jsdot = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd Jdot = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact_frames_.size(); i++) {
            if (in_contact_.at(i)) {
                int frame_id = pin_model_.getFrameId(pin_model_.frames.at(contact_frames_.at(i)).name,
                                                     pin_model_.frames.at(contact_frames_.at(i)).type);

                pinocchio::getFrameJacobianTimeVariation(pin_model_, *pin_data_,
                                                         frame_id, pinocchio::LOCAL_WORLD_ALIGNED, Jdot);
                Jsdot.block(j * CONSTRAINT_PER_FOOT, 0, CONSTRAINT_PER_FOOT, Jdot.cols()) =
                        Jdot.topLeftCorner(CONSTRAINT_PER_FOOT, pin_model_.nv);
                j++;
                Jdot = Eigen::MatrixXd::Zero(Jdot.rows(), Jdot.cols());
            }
        }

        return Jsdot;
    }

} // simulator