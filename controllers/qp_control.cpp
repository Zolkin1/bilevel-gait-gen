//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <Eigen/Sparse>

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/center-of-mass.hpp"

#include "qp_control.h"
#include "timer.h"

namespace controller {
    // TODO: Change the vector to be const refs
    // TODO: Need a way to update the desired forces and make sure it is using the correct number of contacts
    QPControl::QPControl(double control_rate, std::string robot_urdf, const std::string &foot_type, int nv,
                         const Eigen::VectorXd& torque_bounds, double friction_coef,
                         const std::vector<double>& base_pos_gains,
                         const std::vector<double>& base_ang_gains,
                         const vector_t& kp_joint_gains,
                         const vector_t& kd_joint_gains,
                         double leg_weight,
                         double torso_weight,
                         double force_weight,
                         int num_contacts) :
                         Controller(control_rate, robot_urdf, foot_type), num_vel_(nv),
                         friction_coef_(friction_coef), des_contact_(num_contacts) {

        torque_bounds_ = torque_bounds;

        torso_tracking_weight_ = torso_weight;
        leg_tracking_weight_ = leg_weight;
        force_tracking_weight_ = force_weight;

        // Set gains
        SetBasePosGains(base_pos_gains.at(0), base_pos_gains.at(1));
        SetBaseAngleGains(base_ang_gains.at(0), base_ang_gains.at(1));
        SetJointGains(kd_joint_gains, kp_joint_gains);

        // Initialize these to 0
        config_target_ = Eigen::VectorXd::Zero(num_vel_ + 1);
        vel_target_ = Eigen::VectorXd::Zero(num_vel_);
        acc_target_ = Eigen::VectorXd::Zero(num_vel_);

        // assume each joint is actuated
        num_actuators_ = num_vel_ - FLOATING_VEL_OFFSET;

        // Set solver settings
        qp_solver_.settings()->setRelativeTolerance(1e-10);
        qp_solver_.settings()->setAbsoluteTolerance(1e-10);
        qp_solver_.settings()->setVerbosity(false);
        qp_solver_.settings()->setPolish(true);     // Makes a HUGE difference

        prev_num_contacts_ = -1;

        log_file_.open("qp_control_log.txt");
        log_file_ << std::setw(15) << "time "
                << std::setw(15) << "q_des "
                << std::setw(15) << "v_des "
                << std::setw(15) << "q "
                << std::setw(15) << "v "
                << std::setw(15) << "a "
                << std::setw(15) << "grf"
                << std::endl;
    }

    Eigen::VectorXd QPControl::ComputeControlAction(const Eigen::VectorXd& q,
                                                    const Eigen::VectorXd& v,
                                                    const Eigen::VectorXd& a,
                                                    const Contact& contact,
                                                    double time) {
        utils::Timer timer("qp control");
        timer.StartTimer();
        Contact contact2(4);
//        contact2.in_contact_ = {true, true, true, true};
//        contact2.contact_frames_ = contact.contact_frames_;
        contact2 = des_contact_;
//        des_contact_ = contact2;

        UpdateConstraintsAndCost(q, v, a, contact2);

        // Solve qp
        if (qp_solver_.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
            std::cerr << "Could not solve WBC QP. Returning 0 control action." << std::endl;
            Eigen::VectorXd temp = Eigen::VectorXd::Zero(num_inputs_);
            return temp;
        }

        Eigen::VectorXd qp_sol = qp_solver_.getSolution();

        Eigen::VectorXd control_action = Eigen::VectorXd::Zero(3*num_inputs_);

        // Question: Wouldn't we expect there to be a non-zero accleration if we don't consider the effect of the ground reaction forces?
        // If there was no downward acceleration then we could never get a force.

        RecoverControlInputs(qp_sol, v, control_action, contact2);

        // since I am passing to a PID controller, I still need to pass the configuration and velocity
        AssignPositionControl(control_action);
        AssignVelocityControl(control_action);

        vector_t grf(12);
        int j = 0;
        for (int i = 0; i < contact.in_contact_.size(); i++) {
            if (contact.in_contact_.at(i)) {
                grf.segment(i*3, 3) << qp_sol.segment(qp_sol.size() - 3*(contact.GetNumContacts() - j), 3);
            } else {
                grf.segment(i*3, 3) << Eigen::Vector3d::Zero();
            }
        }

        // Note: Do I really want the COM?
        // Get the com
//        pinocchio::centerOfMass(pin_model_, *pin_data_, q, v, a, false);

// Currently using the floating base coordinates and NOT the COM
        Eigen::Vector3d com = q.head<3>(); //pin_data_->com[0];
        Eigen::Vector3d vcom = v.head<3>(); //pin_data_->vcom[0];
        Eigen::Vector3d acom = a.head<3>(); //pin_data_->acom[0];

        LogInfo(time, config_target_, vel_target_, com, q.tail(q.size() - 3), vcom, v.tail(v.size() - 3), acom, a.tail(a.size()-3), grf);

        timer.StopTimer();
        timer.PrintElapsedTime();

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

    void QPControl::SetJointGains(const vector_t & kv, const vector_t& kp) {
        kv_joint_ = kv;
        kp_joint_ = kp;
    }

    void QPControl::UpdateDesiredContacts(const Contact& contact) {
        des_contact_ = contact;
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
        // TODO: g really causes a lot of problems
        lb_.head(FLOATING_VEL_OFFSET) = -(pin_data_->g.head(FLOATING_VEL_OFFSET) +
                (pin_data_->C*v).head(FLOATING_VEL_OFFSET));
        ub_.head(FLOATING_VEL_OFFSET) = lb_.head(FLOATING_VEL_OFFSET);

        A_.topLeftCorner(FLOATING_VEL_OFFSET, num_vel_) =
                pin_data_->M.topLeftCorner(FLOATING_VEL_OFFSET, num_vel_);

//        A_.topRightCorner(FLOATING_VEL_OFFSET, CONSTRAINT_PER_FOOT*GetNumBothContacts(contact, des_contact_)) =
        A_.block(0, num_vel_, FLOATING_VEL_OFFSET,
               CONSTRAINT_PER_FOOT* GetNumBothContacts(contact, des_contact_)) =
            -Js_.transpose().topRows(FLOATING_VEL_OFFSET);
    }

    void QPControl::AddContactMotionConstraints(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Contact& contact) {
        // Assign to solver params
        // TODO: This appears to be wrong
        if (Js_.size() > 0) {
//            std::cout << "Jsdot: \n" << Jsdot_ << std::endl;
//            std::cout << "Js: \n" << Js_ << std::endl;

            pinocchio::forwardKinematics(pin_model_, *pin_data_, q, v, 0*v);
            for (int i = 0; i < des_contact_.in_contact_.size(); i++) {
                if (contact.in_contact_.at(i) && des_contact_.in_contact_.at(i)) {
                    int frame_id = static_cast<int>(pin_model_.getFrameId(pin_model_.frames.at(des_contact_.contact_frames_.at(i)).name,
                                                                          pin_model_.frames.at(des_contact_.contact_frames_.at(i)).type));
                    lb_.segment(FLOATING_VEL_OFFSET + i*3, 3) =
                            -pinocchio::getFrameClassicalAcceleration(pin_model_, *pin_data_,
                                                                     frame_id, pinocchio::LOCAL_WORLD_ALIGNED).linear(); //-Jsdot_ * v;
                }
            }
            const int num_contacts = GetNumBothContacts(des_contact_, contact);
            ub_.segment(FLOATING_VEL_OFFSET, num_contacts*3) =
                    lb_.segment(FLOATING_VEL_OFFSET, num_contacts*3);

            A_.block(FLOATING_VEL_OFFSET, 0, Js_.rows(), Js_.cols()) = Js_;
        }

    }

    void QPControl::AddTorqueConstraints(const Eigen::VectorXd& v, const Contact& contact) {
        int num_contacts = GetNumBothContacts(contact, des_contact_);
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
               -(pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) - torque_bounds_; // removed + pin_data_->g
        ub_.segment(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT * num_contacts, num_actuators_) =
                -(pin_data_->C*v + pin_data_->g).segment(FLOATING_VEL_OFFSET, num_actuators_) + torque_bounds_; // removed + pin_data_->g
    }

    void QPControl::AddFrictionConeConstraints(const Contact& contact) {
        // Since we assume flat ground:
        Eigen::Vector3d h = {1, 0, 0};
        Eigen::Vector3d l = {0, 1, 0};
        Eigen::Vector3d n = {0, 0, 1};

        int num_contacts = GetNumBothContacts(contact, des_contact_);
        for (int i = 0; i < num_contacts; i++) {
            A_.block(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*i,
                     num_vel_ + 3*i, 4, 3) <<
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

    void QPControl::AddPositiveGRFConstraints(const controller::Contact& contact) {
        int num_contacts = GetNumBothContacts(contact, des_contact_);
        for (int i = 0; i < num_contacts; i++) {
            A_(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*num_contacts + i,
                     num_vel_ + 2 + 3*i) = 1;
            ub_(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*num_contacts + i) = 200;
            lb_(FLOATING_VEL_OFFSET + CONSTRAINT_PER_FOOT*num_contacts + num_actuators_ + 4*num_contacts + i) = 0;
        }
    }

    // TODO: Make the leg tracking able to only track subset's of legs (i.e. swing legs)
    void QPControl::AddLegTrackingCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        // For now, just make all the joint accelerations go to 0.
        P_.block(FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET, num_inputs_, num_inputs_) =
                 leg_tracking_weight_ * 2*Eigen::MatrixXd::Identity(num_inputs_, num_inputs_);

        Eigen::VectorXd target = acc_target_.tail(num_inputs_) +
                                 kv_joint_.cwiseProduct(vel_target_.tail(num_inputs_) - v.tail(num_inputs_)) +
                                 kp_joint_.cwiseProduct(config_target_.tail(num_inputs_) - q.tail(num_inputs_));

        w_.segment(FLOATING_VEL_OFFSET, num_inputs_) = -2*target*leg_tracking_weight_;
    }

    void QPControl::AddTorsoCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
        // Add the position costs
        P_.topLeftCorner(POS_VARS, POS_VARS) =
                torso_tracking_weight_*2*Eigen::MatrixXd::Identity(POS_VARS, POS_VARS);

        // Calculate the com
        pinocchio::centerOfMass(pin_model_, *pin_data_, q, v, false);

        Eigen::VectorXd target = acc_target_.head(POS_VARS)
                + kv_pos_ * (vel_target_.head(POS_VARS) - v.head<3>()) //pin_data_->vcom[0])
                + kp_pos_*(config_target_.head(POS_VARS) - q.head<3>()); //pin_data_->com[0]);

        w_.head(POS_VARS) = -2*target*torso_tracking_weight_;

        // Add the orientation costs
        // TODO: Orientation costs break the system
        P_.block(POS_VARS, POS_VARS, 3, 3) =
                torso_tracking_weight_*2*Eigen::MatrixXd::Identity(3,3);

        Eigen::Quaternion<double> orientation(static_cast<Eigen::Vector4d>(q.segment(POS_VARS, 4)));
        orientation.normalize();
        Eigen::Quaternion<double> des_orientation(static_cast<Eigen::Vector4d>(config_target_.segment(POS_VARS, 4)));
        des_orientation.normalize();

        // Note: velocity orientations are in different frames. Need to convert.
        // Should be able to exp the vel back to the surface. Then invert to the local frame then log back to the tangent space
        Eigen::Quaterniond vel_quat;
        Eigen::Vector3d temp = vel_target_.segment<3>(POS_VARS);
        pinocchio::quaternion::exp3(temp, vel_quat);
        vel_quat = orientation.inverse()*vel_quat;
        Eigen::Vector3d vel_frame = pinocchio::quaternion::log3(vel_quat);

        Eigen::VectorXd angle_target = kv_ang_*(vel_frame - v.segment(POS_VARS, 3)) +
                kp_ang_*pinocchio::quaternion::log3(orientation.inverse()*des_orientation);

        w_.segment(POS_VARS, 3) = -2*angle_target*torso_tracking_weight_;
    }

    void QPControl::AddForceTrackingCost(const Contact& contact) {
        int num_contacts = GetNumBothContacts(contact, des_contact_);
        P_.block(FLOATING_VEL_OFFSET + num_inputs_, num_vel_,
                 CONSTRAINT_PER_FOOT*num_contacts, CONSTRAINT_PER_FOOT*num_contacts) =
                         force_tracking_weight_*2*Eigen::MatrixXd::Identity(CONSTRAINT_PER_FOOT*num_contacts, CONSTRAINT_PER_FOOT*num_contacts);

        Eigen::VectorXd target = force_target_;

        w_.segment(FLOATING_VEL_OFFSET + num_inputs_, CONSTRAINT_PER_FOOT*num_contacts) = -2*target*force_tracking_weight_;
    }

    void QPControl::UpdateForceTargets(const Eigen::VectorXd& force) {
        force_target_ = force;
    }

    void QPControl::UpdateConstraintsAndCost(const Eigen::VectorXd& q,
                                             const Eigen::VectorXd& v,
                                             const Eigen::VectorXd& a,
                                             const Contact& contact) {
        int num_contacts = GetNumBothContacts(contact, des_contact_);
        num_decision_vars_ = num_vel_ + CONSTRAINT_PER_FOOT*num_contacts;

        qp_solver_.data()->setNumberOfVariables(num_decision_vars_);

        num_constraints_ = FLOATING_VEL_OFFSET + 7*num_contacts + num_actuators_ + num_contacts;
        qp_solver_.data()->setNumberOfConstraints(num_constraints_);

        A_ = Eigen::MatrixXd::Zero(num_constraints_, num_decision_vars_);
        lb_ = Eigen::VectorXd::Zero(num_constraints_);
        ub_ = lb_;

        P_ = Eigen::MatrixXd::Zero(num_decision_vars_, num_decision_vars_);
        w_ = Eigen::VectorXd::Zero(num_decision_vars_);

        // Add the constraints to the matricies
        ComputeDynamicsTerms(q, v, contact);
        /// TODO: Keep an eye on the nonlinear terms
        acc_target_.setZero();  // pin_data_->M*a + pinocchio::nonLinearEffects(pin_model_, *pin_data_, q, v);
        acc_target_.head(FLOATING_VEL_OFFSET).setZero();

        AddDynamicsConstraints(v, contact);
        AddContactMotionConstraints(q, v, contact);
        AddTorqueConstraints(v, contact);
        AddFrictionConeConstraints(contact);
        AddPositiveGRFConstraints(contact);

        // Init the cost matricies
        AddLegTrackingCost(q, v);
        AddTorsoCost(q, v); // TODO: Investigate/fix -- need to check the scaling on this constraint and check that the dynamics are good
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
                (pin_data_->M*qp_sol.head(num_vel_) - Js_.transpose()*qp_sol.tail(CONSTRAINT_PER_FOOT*
                                                              GetNumBothContacts(contact, des_contact_))
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

        int num_contacts = GetNumBothContacts(contact, des_contact_);
        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact.contact_frames_.size(); i++) {
            if (contact.in_contact_.at(i) && des_contact_.in_contact_.at(i)) {
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

        int num_contacts = GetNumBothContacts(contact, des_contact_);
        Eigen::MatrixXd Jsdot = Eigen::MatrixXd::Zero(CONSTRAINT_PER_FOOT*num_contacts, num_vel_);
        Eigen::MatrixXd Jdot = Eigen::MatrixXd::Zero(FLOATING_VEL_OFFSET, num_vel_);

        int j = 0;
        for (int i = 0; i < contact.contact_frames_.size(); i++) {
            if (contact.in_contact_.at(i) && des_contact_.in_contact_.at(i)) {
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

    int QPControl::GetNumBothContacts(const Contact& contact1, const Contact& contact2) const {
        assert(contact1.in_contact_.size() == contact2.in_contact_.size());
        int num_contacts = 0;
        for (int i = 0; i < contact1.in_contact_.size(); i++) {
            if (contact1.in_contact_.at(i) && contact2.in_contact_.at(i)) {
                num_contacts++;
            }
        }
        return num_contacts;
    }

    void QPControl::LogInfo(double time,
                            const Eigen::VectorXd& q_des,
                            const Eigen::VectorXd& v_des,
                            const Eigen::Vector3d& com,
                            const Eigen::VectorXd& q,
                            const Eigen::Vector3d& vcom,
                            const Eigen::VectorXd& v,
                            const Eigen::Vector3d& acom,
                            const Eigen::VectorXd& a,
                            const Eigen::VectorXd& grf) {
        log_file_ << std::setw(15) << time
                << std::setw(15) << q_des.transpose()
                << std::setw(15) << v_des.transpose()
                << std::setw(15) << com.transpose()
                << std::setw(15) << q.transpose()
                << std::setw(15) << vcom.transpose()
                << std::setw(15) << v.transpose()
                << std::setw(15) << acom.transpose()
                << std::setw(15) << a.transpose()
                << std::setw(15) << grf.transpose() << std::endl;

        log_file_ << std::endl;

    }

} // controller