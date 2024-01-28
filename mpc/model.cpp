//
// Created by zolkin on 1/18/24.
//

#include "pinocchio/algorithm/center-of-mass.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"

#include "model.h"

namespace mpc {
    Model::Model(const std::string& robot_urdf, const std::vector<std::string>& frames, int discretization_steps,
                 double dt, bool uses_joints, const std::vector<Constraints>& constraints) :
            GRAVITY(0., 0., -9.81), num_ee_(frames.size()), uses_joints_(uses_joints),
            constraints_(constraints) {

//        integrator_ = std::make_unique<RKIntegrator>(dt);

        // create the pinocchio model - always a free flyer
        pinocchio::urdf::buildModel(robot_urdf, pinocchio::JointModelFreeFlyer(), pin_model_, false);

        // create the pinocchio data
        pin_data_ = std::make_unique<pinocchio::Data>(pin_model_);

        robot_mass_ = pinocchio::computeTotalMass(pin_model_);

        frames_= frames;

        CreateFrameMap(frames);

        if (discretization_steps != 1) {
            throw std::runtime_error("Only discretization step of 1 is currently supported.");
        }
        discretization_steps_ = discretization_steps;
    }

    void Model::CreateFrameMap(const std::vector<std::string>& frames) {
        for (int i = 0; i < num_ee_; i++) {
            for (int j = 0; j < pin_model_.frames.size(); j++) {
                if (pin_model_.frames.at(j).name == frames.at(i)) {
                    frame_map_.insert(std::pair<std::string, int>(frames.at(i), pin_model_.getFrameId(frames.at(i))));
                }
            }
        }
    }

    int Model::GetNumEndEffectors() const {
        return num_ee_;
    }

    std::vector<int> Model::GetContactFrames() const {
        std::vector<int> frames;
        for (const auto& frame : frames_) {
            frames.push_back(frame_map_.at(frame));
        }

        return frames;
    }

    double Model::GetMass() const {
        return robot_mass_;
    }

    bool Model::UsesJoints() const {
        return uses_joints_;
    }

    const std::vector<Constraints>& Model::GetApplicableConstraints() const {
        return constraints_;
    }

    int Model::GetNumJoints() const {
        return 0;
    }

    int Model::GetFullModelConfigSpace() const {
        return pin_model_.nq;
    }

    matrix_33t Model::GetEEJacobian(int end_effector, const mpc::vector_t& q) {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);
        pinocchio::updateFramePlacements(pin_model_, *pin_data_);

        pinocchio::Data::Matrix6x J(6, pin_model_.nv);
        J.setZero();

        const int id = frame_map_.at(frames_.at(end_effector));
        pinocchio::computeFrameJacobian(pin_model_, *pin_data_, q, id, J);

        return J.block<3,3>(0, 6 + 3*end_effector);    // TODO: Make not hard coded
    }

    matrix_33t Model::GetEEJacobianDeriv(int end_effector, const vector_t& q, const vector_t& v) {
        pinocchio::computeJointJacobiansTimeVariation(pin_model_, *pin_data_, q, v);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);

        pinocchio::Data::Matrix6x Jdot(6, pin_model_.nv);
        Jdot.setZero();

        const int id = frame_map_.at(frames_.at(end_effector));
        pinocchio::getFrameJacobianTimeVariation(pin_model_, *pin_data_, id, pinocchio::LOCAL_WORLD_ALIGNED, Jdot);

        return Jdot.block<3,3>(0, 6 + 3*end_effector);    // TODO: Make not hard coded
    }

    matrix_33t Model::GetBaseRotationMatrix(const vector_t& q) {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);
        return pin_data_->oMi[pin_model_.getJointId("root_joint")].rotation().matrix();
    }

    matrix_3t Model::GetCoriolisMat(const vector_t& q, const vector_t& v) {
        // Get C
        pinocchio::computeCoriolisMatrix(pin_model_, *pin_data_, q, v);
        return pin_data_->C;
    }

    vector_3t Model::GetGravityVec(const mpc::vector_t& q) {
        // Get g
        pinocchio::computeGeneralizedGravity(pin_model_, *pin_data_, q);

        return pin_data_->g;
    }

    matrix_33t Model::GetOperationalSpaceInertia(int end_effector, const mpc::vector_t& q) {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);

        // Get M
        pinocchio::crba(pin_model_, *pin_data_,  q);

        // Make M symmetric
        pin_data_->M.triangularView<Eigen::StrictlyLower>() =
        pin_data_->M.transpose().triangularView<Eigen::StrictlyLower>();

        matrix_33t J = GetEEJacobian(end_effector, q);

        return J*pin_data_->M.ldlt().solve(J.transpose());
    }
} // mpc