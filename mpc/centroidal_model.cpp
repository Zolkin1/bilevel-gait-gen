//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/centroidal-derivatives.hpp"
#include "pinocchio/parsers/urdf.hpp"

#include "centroidal_model.h"


namespace mpc {
    using matrix6x_t = Eigen::Matrix<double, 6, Eigen::Dynamic>;

    // TODO: Move to a xyz orientation representation for the configuration. Drop the quaternion.

    CentroidalModel::CentroidalModel(const std::string &robot_urdf, const std::vector<std::string>& frames,
                                     int discretization_steps) :
            GRAVITY(0., 0., -9.81) {
        // create the pinocchio model - always a free flyer
        pinocchio::urdf::buildModel(robot_urdf, pinocchio::JointModelFreeFlyer(), pin_model_, false);

        // create the pinocchio data
        pin_data_ = std::make_unique<pinocchio::Data>(pin_model_);

        num_joints_ = pin_model_.nq - FLOATING_BASE_OFFSET;
        num_total_states_ = MOMENTUM_OFFSET + pin_model_.nq - 1;

        robot_mass_ = pin_data_->mass[0];

        num_ee_ = frames.size();
        frames_= frames;

        CreateFrameMap(frames);

        if (discretization_steps != 1) {
            throw std::runtime_error("Only discretization step of 1 is currently supported.");
        }
        discretization_steps_ = discretization_steps;
    }

    matrix_t CentroidalModel::GetLinearDiscreteDynamicsState(const vector_t& state, const Inputs& input,
                                                             int node, double time) {
        matrix_t A = matrix_t::Zero(num_total_states_, num_total_states_);

        // TODO: double check all the frames

        // TODO: is this zero or the jacobian of the FK?
        vector_t q = ConvertMPCStateToPinocchioState(state);
        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);    // compute here so we don't need to recompute each time
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);
        Eigen::Matrix3Xd cross_sum = Eigen::Matrix3Xd::Zero(3, q.size());
        for (int i = 0; i < num_ee_; i++) {
            // Get the jacobian of the FK
            matrix_t J = GetFKJacobianForEndEffector(q, frames_.at(i), false);
            Eigen::Vector3d force = input.GetForce(i, time);
            for (int j = 0; j < q.size(); j++) {
                // Cross product with the corresponding input force, add these for all the end effectors
                cross_sum.col(j) += static_cast<Eigen::Vector3d>(-J.col(j)).cross(force);
            }
        }
        A.middleRows<3>(3) = cross_sum;

        // Deal with centroidal dynamics
        matrix6x_t CMM = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);  // apparently dhdotda = CMM according to pinocchio docs
        matrix6x_t dhdq = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
        matrix6x_t dhdotdq = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
        matrix6x_t dhdotdv = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);

        vector_t a = vector_t::Zero(pin_model_.nv);

        // Get all the terms I need
        pinocchio::computeCentroidalDynamicsDerivatives(pin_model_, *pin_data_, q, input.GetVels(time), a,
                                                        dhdq, dhdotdq, dhdotdv, CMM);

        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
        matrix_t CMMj = CMM.rightCols(num_joints_);

        matrix_t CMMbinv = CMMb.inverse();
        A.middleRows<MOMENTUM_OFFSET>(FLOATING_VEL_OFFSET) << CMMbinv, -CMMbinv * dhdq; // TODO: Check negative sign


        // TODO: For now I'm just assuming the CMM has no partial wrt q. Will need to come back and change this.
        // I need the partial w.r.t to each configuration variable
        // Note: (A^-1)' = (A^-1)*A'*(A^-1)
        // So lets get CMM'


        // Discretize with the integrator
        return integrator_->CalcDerivWrtStateSingleStep(state, A);

    }

    matrix_t CentroidalModel::GetLinearDiscreteDynamicsInput(const vector_t& state, const vector_t& input) {
        // Calculate the jacobian wrt the input and evaluate it at the input and state given
        // Remember to discretize!
    }

    vector_t CentroidalModel::GetConstantDiscreteDynamics(const vector_t& state, const vector_t& input) {
        // Evaluate the discrete dynamics at the state and input to get constant term
        // Remember to discretize!

    }

    matrix_t CentroidalModel::GetFKJacobianForEndEffector(const vector_t& q, const std::string& frame, bool compute_jac) {
        if (compute_jac) {
            pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
            pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);
        }

        matrix6x_t J(6,pin_model_.nv);
        pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_map_.at(frame), pinocchio::LOCAL_WORLD_ALIGNED, J);

        // Only return the position coordinates
        return J.topRows<3>();
    }

    Eigen::Vector3d CentroidalModel::GetEndEffectorLocationCOMFrame(const vector_t& state, const std::string& frame) const {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, state.tail(num_joints_));
        return pin_data_->oMi.at(frame_map_.at(frame)).translation();
    }

    // Note: The MPC state is NOT using a quaternion representation
    int CentroidalModel::GetNumConfig() const {
        return pin_model_.nq - 1;
    }

    int CentroidalModel::GetNumJoints() const {
        return num_joints_;
    }

    int CentroidalModel::GetNumEndEffectors() const {
        return num_ee_;
    }

    vector_t CentroidalModel::ConvertMPCStateToPinocchioState(const mpc::vector_t &state) const {
        Eigen::Vector4d quat = ConvertZYXRotToQuaternion(state.segment(3,3));
        vector_t q(state.size() + 1);
        q.head(POS_VARS) = state.head(POS_VARS);
        q.segment(POS_VARS, 4) = quat;
        q.tail(num_joints_) = state.tail(num_joints_);

        return quat;
    }

    // converts to a x,y,z,w quaternion
    Eigen::Vector4d CentroidalModel::ConvertZYXRotToQuaternion(const Eigen::Vector3d& zyx_rot) {
        Eigen::Quaterniond quat = static_cast<Eigen::Quaterniond>(Eigen::AngleAxisd(zyx_rot(0), Eigen::Matrix<double, 3, 1>::UnitZ()) *
                                           Eigen::AngleAxisd(zyx_rot(1), Eigen::Matrix<double, 3, 1>::UnitY()) *
                                           Eigen::AngleAxisd(zyx_rot(2), Eigen::Matrix<double, 3, 1>::UnitX()));
        Eigen::Vector4d quat_vec;
        quat_vec(0) = quat.x();
        quat_vec(1) = quat.y();
        quat_vec(2) = quat.z();
        quat_vec(3) = quat.w();
        return quat_vec;
    }

    Eigen::Vector3d CentroidalModel::ConvertQuaternionToZYXRot(const Eigen::Vector4d& quat) {
        Eigen::Vector3d rot;
        double q_w = quat(3);
        double q_x = quat(0);
        double q_y = quat(1);
        double q_z = quat(2);
        // x
        rot(2) = std::atan2(2*(q_w*q_x + q_y*q_z), 1 - 2*(q_x*q_x + q_y*q_y));
        // y
        rot(1) = -M_PI/2 + 2*std::atan2(std::sqrt(1 + 2*(q_w*q_y - q_x*q_z)), std::sqrt(1 - 2*(q_w*q_y - q_x*q_z)));
        // z
        rot(0) = std::atan2(2*(q_w*q_z + q_x*q_y), 1 - 2*(q_y*q_y + q_z*q_z));

        return rot;
    }

    void CentroidalModel::CreateFrameMap(const std::vector<std::string>& frames) {
        for (int i = 0; i < num_ee_; i++) {
            for (int j = 0; j < pin_model_.frames.size(); j++)
            if (pin_model_.frames.at(j).name == frames.at(i)) {
                frame_map_.insert(std::pair<std::string, int>(frames.at(i), pin_model_.getFrameId(frames.at(i))));
            }
        }
    }

    vector_t CentroidalModel::CalcDynamics(const vector_t& state, const Inputs& input, double time) const {
        vector_t xdot = vector_t::Zero(num_total_states_);

        xdot.topRows<3>() = robot_mass_*GRAVITY;

        for (int i = 0; i < num_ee_; i++) {
            xdot.topRows<3>() += input.GetForce(i, time);

            // Note that the input ee position can deviate from the state ee position.
            // So I use the state ee position.
            xdot.middleRows<3>(3) += GetEndEffectorLocationCOMFrame(state, frames_.at(i)).cross(input.GetForce(i, time));
        }

        matrix_t CMM = pinocchio::computeCentroidalMap(pin_model_, *pin_data_, state.tail(num_joints_));
        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
        matrix_t CMMj = CMM.rightCols(num_joints_);
        xdot.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) =
                CMMb.householderQr().solve(state.head(MOMENTUM_OFFSET) - CMMj*input.GetVels(time));

        xdot.bottomRows(num_joints_) = input.GetVels(time);

        return xdot;
    }

} // mpc