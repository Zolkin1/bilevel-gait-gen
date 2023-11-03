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
                         const Eigen::VectorXd& torque_bounds, double friction_coef) :
                         Controller(control_rate, robot_urdf, foot_type), num_vel_(nv), friction_coef_(friction_coef) {

        torque_bounds_ = torque_bounds;

        // TODO: be able to adjust these
        torso_tracking_weight_ = 10;
        leg_tracking_weight_ = 10;
        force_tracking_weight_ = 1;

        // TODO: be able to adjust these
        kv_pos_ = 1;
        kp_pos_ = 2;
        kv_ang_ = 1;
        kp_ang_ = 2;

        // Initialize these to 0
        config_target_ = Eigen::VectorXd::Zero(num_vel_ + 1);
        vel_target_ = Eigen::VectorXd::Zero(num_vel_);
        acc_target_ = Eigen::VectorXd::Zero(num_vel_);

        // TODO: Adjust any solver settings
        qp_solver_.settings()->setVerbosity(false);
    }

    void QPControl::InitSolver(const mjModel* model, const mjData* data) {
        int num_contacts = UpdateContacts(model, data);

        // assume each joint is actuated
        num_actuators_ = num_vel_ - FLOATING_VEL_OFFSET;

        num_decision_vars_ = num_vel_ + CONSTRAINT_PER_FOOT*num_contacts;

        qp_solver_.data()->setNumberOfVariables(num_decision_vars_);

        num_constraints_ = num_vel_ + 7*num_contacts + num_actuators_;
        qp_solver_.data()->setNumberOfConstraints(num_constraints_);

        A_ = Eigen::MatrixXd::Zero(num_constraints_, num_decision_vars_);
        lb_ = Eigen::VectorXd::Zero(num_constraints_);
        ub_ = lb_;

        P_ = Eigen::MatrixXd::Zero(num_decision_vars_, num_decision_vars_);
        w_ = Eigen::VectorXd::Zero(num_decision_vars_);

        const Eigen::VectorXd q = ConvertMujocoConfigToPinocchio(data);
        const Eigen::VectorXd v = ConvertMujocoVelToPinocchio(data);
        const Eigen::VectorXd a = ConvertMujocoAccToPinocchio(data);

        // Add the constraints to the matricies
        ComputeDynamicsTerms(q, v, a);
        AddDynamicsConstraints(v);
        AddContactMotionConstraints(v);
        AddTorqueConstraints(v);
        AddFrictionConeConstraints(v);

        // Init solver constraints
        Eigen::SparseMatrix<double> sparseA = A_.sparseView();
        if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
              qp_solver_.data()->setLowerBound(lb_) &&
              qp_solver_.data()->setUpperBound(ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // TODO
        // Init the cost matricies
        AddLegTrackingCost(q,v);
        AddTorsoCost(q, v);
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
        int num_contacts = UpdateContacts(model, data);

        UpdateConstraintsAndCost(data);

        // Solve qp
        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            std::cerr << "Could not solve WBC QP. Returning 0 control action." << std::endl;
            std::vector<mjtNum> temp;
            for (int i = 0; i < model->njnt; i++) {
                temp.push_back(0);
            }
            return temp;
        }

        Eigen::VectorXd qp_sol = qp_solver_.getSolution();

        std::vector<mjtNum> control_action(3*num_inputs_);

        const Eigen::VectorXd v = ConvertMujocoVelToPinocchio(data);
        RecoverControlInputs(qp_sol, v, control_action);

        // since i am passing to a PID controller, I still need to pass the configuration and velocity
        AssignPositionControl(control_action);
        AssignVelocityControl(control_action);

//        std::cout << "feed forward: " << std::endl;
//        for (int i = 0; i < num_inputs_; i++) {
//            std::cout << control_action.at(i + 2*num_inputs_) << std::endl;
//        }

        return control_action;
    }

    void QPControl::ComputeDynamicsTerms(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Eigen::VectorXd& a) {
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
        JsR_ = GetConstraintRotationJacobian(q, v);

        // Contact acceleration
        contact_acc_ = GetContactAcc(q, v, a);

    }

    void QPControl::AddDynamicsConstraints(const Eigen::VectorXd& v) {
        // Add equality constraints to the solver params
        lb_.head(FLOATING_VEL_OFFSET) = -(pin_data_->g.head(FLOATING_VEL_OFFSET) +
                (pin_data_->C*v).head(FLOATING_VEL_OFFSET));
        ub_.head(FLOATING_VEL_OFFSET) = lb_.head(FLOATING_VEL_OFFSET);

        A_.topLeftCorner(FLOATING_VEL_OFFSET, num_vel_) =
                pin_data_->M.topLeftCorner(FLOATING_VEL_OFFSET, num_vel_);

        A_.topRightCorner(FLOATING_VEL_OFFSET, CONSTRAINT_PER_FOOT*GetNumContacts()) =
                -Js_.transpose().topLeftCorner(FLOATING_VEL_OFFSET, Js_.rows());

    }

    void QPControl::AddContactMotionConstraints(const Eigen::VectorXd& v) {
        // Assign to solver params
        if (Js_.size() > 0) {
            // TODO: which one to use?
            lb_.segment(FLOATING_VEL_OFFSET, Js_.rows()) = -Jsdot_ * v; // -contact_acc_;
            ub_.segment(FLOATING_VEL_OFFSET, Js_.rows()) = lb_.segment(FLOATING_VEL_OFFSET, Js_.rows());

            A_.block(FLOATING_VEL_OFFSET, 0, Js_.rows(), Js_.cols()) = Js_;

//            std::cout << "Jsdot*v: \n" << Jsdot_*v << std::endl;
//            std::cout << "contact acceleration: \n" << contact_acc_ << std::endl;
//            std::cout << "Js: \n" << Js_ << std::endl;
        }

    }

    void QPControl::AddTorqueConstraints(const Eigen::VectorXd& v) {
        if (Js_.size() > 0) {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), 0,
                     num_actuators_, num_decision_vars_) <<
                                     pin_data_->M.block(FLOATING_VEL_OFFSET, 0, num_actuators_, pin_data_->M.cols()),
                    -Js_.transpose().block(FLOATING_VEL_OFFSET, 0, num_actuators_, Js_.rows());
        } else {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), 0, num_actuators_,
                     num_decision_vars_) << pin_data_->M.block(FLOATING_VEL_OFFSET, 0, num_actuators_, pin_data_->M.cols());
        }

        lb_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), num_actuators_) =
               -(pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) - torque_bounds_;
        ub_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * GetNumContacts(), num_actuators_) =
                -(pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) + torque_bounds_;
    }

    void QPControl::AddFrictionConeConstraints(const Eigen::VectorXd& v) {
        // Since we assume flat ground:
        Eigen::Vector3d h = {1, 0, 0};
        Eigen::Vector3d l = {0, 1, 0};
        Eigen::Vector3d n = {0, 0, 1};

        int num_contacts = GetNumContacts();
        for (int i = 0; i < num_contacts; i++) {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i, num_vel_, 4, 3) <<
                    (h - n*friction_coef_).transpose(),
                    -(h + n*friction_coef_).transpose(),
                    (l - n*friction_coef_).transpose(),
                    -(l + n*friction_coef_).transpose();
            ub_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i, 4) = Eigen::Vector4d::Zero();
            lb_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i, 4) =
                    Eigen::Vector4d::Ones() * -OsqpEigen::INFTY;
        }

    }

    // TODO: Implement properly
    void QPControl::AddLegTrackingCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        // For now, just make all the joint accelerations go to 0.
        P_.block(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET, num_inputs_, num_inputs_) =
                P_.block(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET, num_inputs_, num_inputs_) + leg_tracking_weight_ *
                Eigen::MatrixXd::Identity(num_inputs_, num_inputs_);
    }

    void QPControl::AddTorsoCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        // Add the position costs
        P_.topLeftCorner(POS_VARS, POS_VARS) = P_.topLeftCorner(POS_VARS, POS_VARS) +
                torso_tracking_weight_*Eigen::MatrixXd::Identity(POS_VARS, POS_VARS);

        Eigen::VectorXd target = acc_target_.head(POS_VARS) + kv_pos_ * (vel_target_.head(POS_VARS) - v.head(POS_VARS)) +
                kp_pos_*(config_target_.head(POS_VARS) - q.head(POS_VARS));

        w_.head(POS_VARS) = w_.head(POS_VARS) -2*target*torso_tracking_weight_;

        // Add the orientation costs
        P_.block(POS_VARS, 0, 3, 3) = P_.block(POS_VARS, 0, 3, 3) + torso_tracking_weight_*Eigen::MatrixXd::Identity(3,3);

        Eigen::Vector4d temp = config_target_.segment(POS_VARS, 4) - q.segment(POS_VARS, 4);
        Eigen::Quaternion<double> quat_diff(temp);
        Eigen::VectorXd angle_target = -kv_ang_*v.segment(POS_VARS, 3) +
                kp_ang_*pinocchio::quaternion::log3(quat_diff);
        w_.segment(POS_VARS, 3) = w_.segment(POS_VARS, 3) - 2*angle_target*torso_tracking_weight_;
    }

    // TODO: Can I anything on the update? Maybe some of the costs...
    // TODO: if the number of contacts don't change then we can re-use the sparsity
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
        const Eigen::VectorXd a = ConvertMujocoAccToPinocchio(data);

        // Add the constraints to the matricies
        ComputeDynamicsTerms(q, v, a);
        AddDynamicsConstraints(v);
        AddContactMotionConstraints(v);
        AddTorqueConstraints(v);
        AddFrictionConeConstraints(v);

        // Init solver constraints
        Eigen::SparseMatrix<double> sparseA = A_.sparseView();
        // I think if the solve size changes I might need to clear then set
        qp_solver_.data()->clearLinearConstraintsMatrix();
        qp_solver_.data()->clearHessianMatrix();
        qp_solver_.clearSolver();

        if (!(qp_solver_.data()->setLinearConstraintsMatrix(sparseA) &&
              qp_solver_.data()->setBounds(lb_, ub_))) {
            throw std::runtime_error("Unable to add the constraints to the QP solver.");
        }

        // Init the cost matricies
        AddLegTrackingCost(q, v);
        AddTorsoCost(q, v);
//        AddForceTrackingCost();

        // Init solver costs
        Eigen::SparseMatrix<double> sparseP = P_.sparseView();
        if (!(qp_solver_.data()->setHessianMatrix(sparseP) && qp_solver_.data()->setGradient(w_))) {
            throw std::runtime_error("Unable to add the costs to the QP solver.");
        }

        // Maybe I need to re-init
        if (!qp_solver_.initSolver()) {
            throw std::runtime_error("Unable to initialize the solver.");
        }
    }

    void QPControl::RecoverControlInputs(const Eigen::VectorXd& qp_sol, const Eigen::VectorXd& v,
                                         std::vector<mjtNum>& control) {
        // Recover torques using ID
        Eigen::VectorXd tau =
                (pin_data_->M*qp_sol.head(num_vel_) - Js_.transpose()*qp_sol.tail(CONSTRAINT_PER_FOOT*GetNumContacts())
                                                      + pin_data_->C*v + pin_data_->g).tail(num_inputs_);

        Eigen::VectorXd mujoco_tau = ConvertPinocchioJointToMujoco(tau);
//        std::cout << "mujoco_tau: " << mujoco_tau << std::endl;
        for (int i = 2*num_inputs_; i < 3*num_inputs_; i++) {
            control.at(i) = mujoco_tau(i - 2*num_inputs_);
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

                // TODO: Can I use getFrameVelocity()?
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
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);  // might not need this

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

    Eigen::MatrixXd QPControl::GetConstraintRotationJacobian(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);  // might not need this

        int num_contacts = GetNumContacts();
        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact_frames_.size(); i++) {
            if (in_contact_.at(i)) {
                int frame_id = pin_model_.getFrameId(pin_model_.frames.at(contact_frames_.at(i)).name,
                                                     pin_model_.frames.at(contact_frames_.at(i)).type);

                pinocchio::getFrameJacobian(pin_model_, *pin_data_,
                                                         frame_id, pinocchio::LOCAL_WORLD_ALIGNED, J);
                Js.block(j * 3, 0, 3, J.cols()) = J.block(3, 0, 3, pin_model_.nv);
                j++;
                J = Eigen::MatrixXd::Zero(J.rows(), J.cols());
            }
        }

        return Js;
    }

    Eigen::VectorXd QPControl::GetContactAcc(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Eigen::VectorXd& a) {
        // second order fk:
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q, v, a);
        Eigen::VectorXd contact_acc = Eigen::VectorXd::Zero(CONSTRAINT_PER_FOOT*GetNumContacts());

        int j = 0;
        for (int i = 0; i < contact_frames_.size(); i++) {
            if (in_contact_.at(i)) {
                int frame_id = pin_model_.getFrameId(pin_model_.frames.at(contact_frames_.at(i)).name,
                                                     pin_model_.frames.at(contact_frames_.at(i)).type);

                pinocchio::Motion::Vector6 frame_acc =
                        pinocchio::getFrameAcceleration(pin_model_, *pin_data_, frame_id, pinocchio::LOCAL_WORLD_ALIGNED);
                contact_acc.segment(j * CONSTRAINT_PER_FOOT, CONSTRAINT_PER_FOOT) = frame_acc.head(CONSTRAINT_PER_FOOT);
                j++;
            }
        }

        return contact_acc;
    }

} // simulator