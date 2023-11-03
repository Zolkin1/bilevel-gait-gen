//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <Eigen/Core>
#include <Eigen/QR>
#include <utility>

#include "mujoco.h"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/frames.hpp"

#include "pd_grav_comp.h"


namespace simulator {
    PDGravComp::PDGravComp(double control_freq, std::string robot_urdf, const std::string& foot_type,
                           const Eigen::VectorXd& config_set_point, const Eigen::VectorXd& vel_set_point):
            Controller(control_freq, std::move(robot_urdf), foot_type) {
        acc_target_ = Eigen::VectorXd::Zero(num_inputs_);

        config_target_ = config_set_point;
        vel_target_ = vel_set_point;
    }

    std::vector<mjtNum> PDGravComp::ComputeControlAction(const mjModel* model, const mjData* data) {
        // Compute position and velocity set point

        // Update contacts
        UpdateContacts(model, data);

        // Compute feedforward torque
        ComputeFeedForwardValue(data);

        // Assign the values in order to the vector
        std::vector<mjtNum> control(3*num_inputs_);     // 3 gives position, velocity, feedforward
        AssignPositionControl(control);
        AssignVelocityControl(control);
        AssignFeedForward(control);

        return control;
    }

    void PDGravComp::ComputeFeedForwardValue(const mjData* data) {
        const Eigen::VectorXd q = ConvertMujocoConfigToPinocchio(data);
        const Eigen::VectorXd v = ConvertMujocoVelToPinocchio(data);

        // Update the forward kinematics
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);

        // S = input selection matrix = [0, I] (nv - FLOATING_VEL x nv)
        Eigen::MatrixXd tempZero =
                Eigen::MatrixXd::Zero(pin_model_.nv - FLOATING_VEL_OFFSET, FLOATING_VEL_OFFSET);
        Eigen::MatrixXd tempIdentity =
                Eigen::MatrixXd::Identity(pin_model_.nv - FLOATING_VEL_OFFSET, pin_model_.nv - FLOATING_VEL_OFFSET);

        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(pin_model_.nv - FLOATING_VEL_OFFSET, pin_model_.nv);

        S << tempZero, tempIdentity;

        // M = floating base inertia matrix
        pinocchio::crba(pin_model_, *pin_data_,  q);

        // Make M symmetric
        pin_data_->M.triangularView<Eigen::StrictlyLower>() = pin_data_->M.transpose().triangularView<Eigen::StrictlyLower>();

        // h = floating base centripetal, coriolis and gravity forces
        pinocchio::computeCoriolisMatrix(pin_model_, *pin_data_, q, v);
        pinocchio::computeGeneralizedGravity(pin_model_, *pin_data_, q);


        // Jc = Jacobian of k linearly indep. constraints (k X nv)
        int k = 0;
        for (auto && i : in_contact_) {
            if (i) {
                k += CONSTRAINT_PER_FOOT;
            }
        }
        Eigen::MatrixXd Jc = Eigen::MatrixXd::Zero(k, pin_model_.nv);

        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);

        Eigen::Matrix<double, 6, Eigen::Dynamic> J =
                Eigen::Matrix<double, 6, Eigen::Dynamic>::Zero(6, pin_model_.nv);

        // Only grab the frames for feet in contact
        int j = 0;
        for (int i = 0; i < contact_frames_.size(); i++) {
            if (in_contact_.at(i)) {
                int frame_id = pin_model_.getFrameId(pin_model_.frames.at(contact_frames_.at(i)).name, pin_model_.frames.at(contact_frames_.at(i)).type);

                pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_id, pinocchio::LOCAL_WORLD_ALIGNED, J);
                Jc.block(j * CONSTRAINT_PER_FOOT, 0, CONSTRAINT_PER_FOOT, J.cols()) =
                        J.topLeftCorner(CONSTRAINT_PER_FOOT, pin_model_.nv);
                j++;
                J = Eigen::MatrixXd::Zero(J.rows(), J.cols());
            }
        }

        // Su = [0 I]
        Eigen::MatrixXd Su = Eigen::MatrixXd::Zero(pin_model_.nv - k, pin_model_.nv);

        tempZero = Eigen::MatrixXd::Zero(Su.rows(), Jc.rows());
        tempIdentity = Eigen::MatrixXd::Identity(Su.rows(), Su.rows());
        Su << tempZero, tempIdentity;

        // Q = QR(J_C)
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(Jc.transpose());
        Eigen::MatrixXd Q = qr.householderQ();

        // tau = pinv(Su*Q^T*S^T)*Su*Q^T*(M*qddot + h)
        Eigen::MatrixXd temp = Su*Q.transpose()*S.transpose();

        // TODO: Take in non-zero acceleration
        Eigen::VectorXd des_qddot = Eigen::VectorXd::Zero(pin_model_.nv);       // For now, it is zero desired acceleration

        acc_target_ = temp.fullPivHouseholderQr().solve(Su * Q.transpose() * (pin_data_->M * des_qddot + pin_data_->C * v + pin_data_->g));

        acc_target_ = ConvertPinocchioJointToMujoco(acc_target_);

    }

    void PDGravComp::AssignFeedForward(std::vector<mjtNum> &control) {
        for (int i = 2*num_inputs_; i < 3*num_inputs_; i++) {
            control.at(i) = acc_target_(i - 2 * num_inputs_);
        }
    }
}