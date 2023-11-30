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

namespace controller {
    // TODO: Change the vector to be const refs
    // TODO: Need a way to update the desired forces and make sure it is using the correct number of contacts
    QPControl::QPControl(double control_rate, std::string robot_urdf, const std::string &foot_type, int nv,
                         const Eigen::VectorXd& torque_bounds, double friction_coef,
                         std::vector<double> base_pos_gains,
                         std::vector<double> base_ang_gains,
                         std::vector<double> joint_gains,
                         double leg_weight,
                         double torso_weight,
                         double force_weight) :
                         Controller(control_rate, robot_urdf, foot_type), num_vel_(nv), friction_coef_(friction_coef) {

        torque_bounds_ = torque_bounds;

        torso_tracking_weight_ = torso_weight;
        leg_tracking_weight_ = leg_weight;
        force_tracking_weight_ = force_weight;

        // Set gains
        SetBasePosGains(base_pos_gains.at(0), base_pos_gains.at(1));
        SetBaseAngleGains(base_ang_gains.at(0), base_ang_gains.at(1));
        SetJointGains(joint_gains.at(0), joint_gains.at(1));

        // Initialize these to 0
        config_target_ = Eigen::VectorXd::Zero(num_vel_ + 1);
        vel_target_ = Eigen::VectorXd::Zero(num_vel_);
        acc_target_ = Eigen::VectorXd::Zero(num_vel_);

        // assume each joint is actuated
        num_actuators_ = num_vel_ - FLOATING_VEL_OFFSET;

        // Set solver settings
        qp_solver_.settings()->setVerbosity(false);
        qp_solver_.settings()->setPolish(true);     // Makes a HUGE difference

        prev_num_contacts_ = -1;
    }

    Eigen::VectorXd QPControl::ComputeControlAction(const Eigen::VectorXd& q,
                                                    const Eigen::VectorXd& v,
                                                    const Eigen::VectorXd& a,
                                                    const Contact& contact) {
        UpdateConstraintsAndCost(q, v, a, contact);

        // Solve qp
        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            std::cerr << "Could not solve WBC QP. Returning 0 control action." << std::endl;
            Eigen::VectorXd temp = Eigen::VectorXd::Zero(num_inputs_);
            return temp;
        }

        Eigen::VectorXd qp_sol = qp_solver_.getSolution();

        Eigen::VectorXd control_action = Eigen::VectorXd::Zero(3*num_inputs_);

        RecoverControlInputs(qp_sol, v, control_action, contact);

        // since I am passing to a PID controller, I still need to pass the configuration and velocity
        AssignPositionControl(control_action);
        AssignVelocityControl(control_action);

        return control_action;
    }

    void QPControl::SetBasePosGains(double kv, double kp) {
        kv_pos_ = kv;
        kp_pos_ = kp;
    }

    void QPControl::SetBaseAngleGains(double kv, double kp) {
        kv_ang_ = kv;
        kp_ang_ = kp;
    }

    void QPControl::SetJointGains(double kv, double kp) {
        kv_joint_ = kv;
        kp_joint_ = kp;
    }

    void QPControl::ComputeDynamicsTerms(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Contact& contact) {
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
        Js_ = GetConstraintJacobian(q, contact);

        // Jsdot
        Jsdot_ = Eigen::MatrixXd::Zero(Jsdot_.rows(), Jsdot_.cols());
        Jsdot_ = GetConstraintJacobianDerivative(q, v, contact);

    }

    void QPControl::AddDynamicsConstraints(const Eigen::VectorXd& v, const Contact& contact) {
        // Add equality constraints to the solver params
        lb_.head(FLOATING_VEL_OFFSET) = -(pin_data_->g.head(FLOATING_VEL_OFFSET) +
                (pin_data_->C*v).head(FLOATING_VEL_OFFSET));
        ub_.head(FLOATING_VEL_OFFSET) = lb_.head(FLOATING_VEL_OFFSET);

        A_.topLeftCorner(FLOATING_VEL_OFFSET, num_vel_) =
                pin_data_->M.topLeftCorner(FLOATING_VEL_OFFSET, num_vel_);

        A_.topRightCorner(FLOATING_VEL_OFFSET, CONSTRAINT_PER_FOOT*contact.GetNumContacts()) =
                -Js_.transpose().topLeftCorner(FLOATING_VEL_OFFSET, Js_.rows());

    }

    void QPControl::AddContactMotionConstraints(const Eigen::VectorXd& v) {
        // Assign to solver params
        if (Js_.size() > 0) {
            lb_.segment(FLOATING_VEL_OFFSET, Js_.rows()) = -Jsdot_ * v;
            ub_.segment(FLOATING_VEL_OFFSET, Js_.rows()) = lb_.segment(FLOATING_VEL_OFFSET, Js_.rows());

            A_.block(FLOATING_VEL_OFFSET, 0, Js_.rows(), Js_.cols()) = Js_;
        }

    }

    void QPControl::AddTorqueConstraints(const Eigen::VectorXd& v, const Contact& contact) {
        int num_contacts = contact.GetNumContacts();
        if (Js_.size() > 0) {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * num_contacts, 0,
                     num_actuators_, num_decision_vars_) <<
                                     pin_data_->M.block(FLOATING_VEL_OFFSET, 0, num_actuators_, pin_data_->M.cols()),
                    -Js_.transpose().block(FLOATING_VEL_OFFSET, 0, num_actuators_, Js_.rows());
        } else {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * num_contacts, 0, num_actuators_,
                     num_decision_vars_) << pin_data_->M.block(FLOATING_VEL_OFFSET, 0, num_actuators_, pin_data_->M.cols());
        }

        lb_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * num_contacts, num_actuators_) =
               -(pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) - torque_bounds_;
        ub_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * num_contacts, num_actuators_) =
                -(pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) + torque_bounds_;
    }

    void QPControl::AddFrictionConeConstraints(const Contact& contact) {
        // Since we assume flat ground:
        Eigen::Vector3d h = {1, 0, 0};
        Eigen::Vector3d l = {0, 1, 0};
        Eigen::Vector3d n = {0, 0, 1};

        int num_contacts = contact.GetNumContacts();
        for (int i = 0; i < num_contacts; i++) {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i, num_vel_ + 3*i, 4, 3) <<
                    (h - n*friction_coef_).transpose(),
                    -(h + n*friction_coef_).transpose(),
                    (l - n*friction_coef_).transpose(),
                    -(l + n*friction_coef_).transpose();
            ub_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i, 4) =
                    Eigen::Vector4d::Zero();
            lb_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i, 4) =
                    Eigen::Vector4d::Ones() * -OsqpEigen::INFTY;
        }

    }

    // TODO: Examine 2's in the cost function
    // TODO: Make the leg tracking able to only track subset's of legs (i.e. swing legs)
    void QPControl::AddLegTrackingCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        // For now, just make all the joint accelerations go to 0.
        P_.block(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET, num_inputs_, num_inputs_) =
                P_.block(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET, num_inputs_, num_inputs_) + leg_tracking_weight_ *
                2*Eigen::MatrixXd::Identity(num_inputs_, num_inputs_);

        Eigen::VectorXd target = acc_target_.tail(num_inputs_) +
                                 kv_joint_*(vel_target_.tail(num_inputs_) - v.tail(num_inputs_)) +
                                 kp_joint_*(config_target_.tail(num_inputs_) - q.tail(num_inputs_));

        w_.segment(FLOATING_VEL_OFFSET, num_inputs_) = w_.segment(FLOATING_VEL_OFFSET, num_inputs_) +
                -2*target*leg_tracking_weight_;
    }

    void QPControl::AddTorsoCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        // Add the position costs
        P_.topLeftCorner(POS_VARS, POS_VARS) = P_.topLeftCorner(POS_VARS, POS_VARS) +
                torso_tracking_weight_*2*Eigen::MatrixXd::Identity(POS_VARS, POS_VARS);

        Eigen::VectorXd target = acc_target_.head(POS_VARS) + kv_pos_ * (vel_target_.head(POS_VARS) - v.head(POS_VARS)) +
                kp_pos_*(config_target_.head(POS_VARS) - q.head(POS_VARS));

        w_.head(POS_VARS) = w_.head(POS_VARS) -2*target*torso_tracking_weight_;

        // Add the orientation costs
        P_.block(POS_VARS, POS_VARS, 3, 3) = P_.block(POS_VARS, 0, 3, 3) + torso_tracking_weight_*2*Eigen::MatrixXd::Identity(3,3);

        Eigen::Quaternion<double> orientation(static_cast<Eigen::Vector4d>(q.segment(POS_VARS, 4)));
        Eigen::Quaternion<double> des_orientation(static_cast<Eigen::Vector4d>(config_target_.segment(POS_VARS, 4)));
        Eigen::VectorXd angle_target = kv_ang_*(vel_target_.segment(POS_VARS, 3) - v.segment(POS_VARS, 3)) +
                kp_ang_*pinocchio::quaternion::log3(orientation.inverse()*des_orientation);
        w_.segment(POS_VARS, 3) = w_.segment(POS_VARS, 3) - 2*angle_target*torso_tracking_weight_;
    }

    void QPControl::AddForceTrackingCost(const Contact& contact) {
        int num_contacts = contact.GetNumContacts();
        P_.block(FLOATING_VEL_OFFSET + num_inputs_, num_vel_,
                 CONSTRAINT_PER_FOOT*num_contacts, CONSTRAINT_PER_FOOT*num_contacts) =
                         2*Eigen::MatrixXd::Identity(CONSTRAINT_PER_FOOT*num_contacts, CONSTRAINT_PER_FOOT*num_contacts);

        Eigen::VectorXd target = force_target_;

        w_.segment(FLOATING_VEL_OFFSET + num_inputs_, CONSTRAINT_PER_FOOT*num_contacts) = -2*target*force_tracking_weight_;
    }

    void QPControl::UpdateConstraintsAndCost(const Eigen::VectorXd& q,
                                             const Eigen::VectorXd& v,
                                             const Eigen::VectorXd& a,
                                             const Contact& contact) {
        int num_contacts = contact.GetNumContacts();
        num_decision_vars_ = num_vel_ + CONSTRAINT_PER_FOOT*num_contacts;

        qp_solver_.data()->setNumberOfVariables(num_decision_vars_);

        num_constraints_ = num_vel_ + 7*num_contacts + num_actuators_;
        qp_solver_.data()->setNumberOfConstraints(num_constraints_);

        A_ = Eigen::MatrixXd::Zero(num_constraints_, num_decision_vars_);
        lb_ = Eigen::VectorXd::Zero(num_constraints_);
        ub_ = lb_;

        P_ = Eigen::MatrixXd::Zero(num_decision_vars_, num_decision_vars_);
        w_ = Eigen::VectorXd::Zero(num_decision_vars_);

        // TODO: Change. For now assuming we are always trying to apply 0 contact force
        force_target_ = Eigen::VectorXd::Zero(3*num_contacts);

        // Add the constraints to the matricies
        ComputeDynamicsTerms(q, v, contact);
        AddDynamicsConstraints(v, contact);
        AddContactMotionConstraints(v);
        AddTorqueConstraints(v, contact);
        AddFrictionConeConstraints(contact);

        // Init the cost matricies
        AddLegTrackingCost(q, v);
        AddTorsoCost(q, v);
        AddForceTrackingCost(contact);

        // Check if the sparsity pattern changed
        if (num_contacts != prev_num_contacts_) {
            qp_solver_.data()->clearLinearConstraintsMatrix();
            qp_solver_.data()->clearHessianMatrix();
            qp_solver_.clearSolver();

            // Set solver constraints
            Eigen::SparseMatrix<double> sparseA = A_.sparseView();
            if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
                  qp_solver_.data()->setBounds(lb_, ub_))) {
                throw std::runtime_error("Unable to add the constraints to the QP solver.");
            }

            // Set solver costs
            Eigen::SparseMatrix<double> sparseP = P_.sparseView();
            if (!(qp_solver_.data()->setHessianMatrix(sparseP) && qp_solver_.data()->setGradient(w_))) {
                throw std::runtime_error("Unable to add the costs to the QP solver.");
            }

            // Re-init
            if (!qp_solver_.initSolver()) {
                throw std::runtime_error("Unable to initialize the solver.");
            }
        } else {
            // Set solver constraints
            Eigen::SparseMatrix<double> sparseA = A_.sparseView();
            if (!(qp_solver_.updateLinearConstraintsMatrix(sparseA) &&
                  qp_solver_.updateBounds(lb_, ub_))) {
                throw std::runtime_error("Unable to add the constraints to the QP solver.");
            }

            // Set solver costs
            Eigen::SparseMatrix<double> sparseP = P_.sparseView();
            if (!(qp_solver_.updateHessianMatrix(sparseP) && qp_solver_.updateGradient(w_))) {
                throw std::runtime_error("Unable to add the costs to the QP solver.");
            }
        }

        prev_num_contacts_ = num_contacts;
    }

    void QPControl::RecoverControlInputs(const Eigen::VectorXd& qp_sol, const Eigen::VectorXd& v,
                                         Eigen::VectorXd& control, const Contact& contact) {
        // Recover torques using ID
        Eigen::VectorXd tau =
                (pin_data_->M*qp_sol.head(num_vel_) - Js_.transpose()*qp_sol.tail(CONSTRAINT_PER_FOOT*contact.GetNumContacts())
                                                      + pin_data_->C*v + pin_data_->g).tail(num_inputs_);

        for (int i = 2*num_inputs_; i < 3*num_inputs_; i++) {
            control(i) = tau(i - 2*num_inputs_);
        }
    }

    Eigen::MatrixXd QPControl::GetConstraintJacobian(const Eigen::VectorXd& q, const Contact& contact) {
        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);

        // Js is the "support Jacobian". Js \in R^{constraint_per_foot*ncontacts x nv}
        // Js = [Jc^T ... Jc^T]^T, Jc \in R^{nv x constraint_per_foot}

        int num_contacts = contact.GetNumContacts();
        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact.contact_frames_.size(); i++) {
            if (contact.in_contact_.at(i)) {
                int frame_id = static_cast<int>(pin_model_.getFrameId(pin_model_.frames.at(contact.contact_frames_.at(i)).name,
                                                     pin_model_.frames.at(contact.contact_frames_.at(i)).type));

                pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_id, pinocchio::LOCAL_WORLD_ALIGNED, J);
                Js.block(j * CONSTRAINT_PER_FOOT, 0, CONSTRAINT_PER_FOOT, J.cols()) =
                        J.topLeftCorner(CONSTRAINT_PER_FOOT, pin_model_.nv);
                j++;
                J = Eigen::MatrixXd::Zero(J.rows(), J.cols());
            }
        }

        return Js;
    }

    Eigen::MatrixXd QPControl::GetConstraintJacobianDerivative(const Eigen::VectorXd& q, const Eigen::VectorXd& v,
                                                               const Contact& contact) {
        pinocchio::computeJointJacobiansTimeVariation(pin_model_, *pin_data_, q, v);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);  // might not need this

        int num_contacts = contact.GetNumContacts();
        Eigen::MatrixXd Jsdot = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd Jdot = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact.contact_frames_.size(); i++) {
            if (contact.in_contact_.at(i)) {
                int frame_id = static_cast<int>(pin_model_.getFrameId(pin_model_.frames.at(contact.contact_frames_.at(i)).name,
                                                     pin_model_.frames.at(contact.contact_frames_.at(i)).type));

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

} // controller