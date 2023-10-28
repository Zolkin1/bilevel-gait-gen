//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <string>
#include <utility>

#include "pinocchio/algorithm/kinematics.hpp"

#include "controller.h"

namespace simulator {
    Controller::Controller(double control_rate, std::string robot_urdf) : rate_(1.0/control_rate),
    robot_urdf_(std::move(robot_urdf)) {
        // create the pinocchio model - always a free flyer
        pinocchio::urdf::buildModel(robot_urdf_, pinocchio::JointModelFreeFlyer(), pin_model_, false);

        // create the pinocchio data
        pin_data_ = std::make_unique<pinocchio::Data>(pin_model_);

        // call forward kinematics so we can ensure it has been called
        Eigen::VectorXd q = Eigen::VectorXd::Zero(pin_model_.nq);
        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);

        // Assuming every joint (except the floating base) is actuated
        num_inputs_ = pin_model_.nq - FLOATING_BASE_OFFSET;
    }

    void Controller::DefineContacts(const std::vector<std::string>& frames, const std::vector<int>& mujoco_bodies) {
        if (frames.size() != mujoco_bodies.size()) {
            std::cerr << "Number of pinocchio contact frames and number of mujoco contact bodies do not match!"
            << "Not assigning contacts in the controller." << std::endl;
            return;
        }

        for (int i = 0; i < frames.size(); i++) {
            for (int j = 0; j < pin_model_.frames.size(); j++) {
                // Check if the names given match a frame in pinocchio then take that index
                if (pin_model_.frames.at(j).name == frames.at(i)) {
                    contact_frames_.push_back(j);
                    mujoco_bodies_.push_back(mujoco_bodies.at(i));
                    in_contact_.push_back(false);
                }
            }
        }

        // Check that we got everything
        if (contact_frames_.size() != frames.size()) {
            std::cerr << "Could not find " << frames.size() - contact_frames_.size()
            << " frames. Check provided frames." << std::endl;
        }
    }

    double Controller::GetRate() const {
        return rate_;
    }

    int Controller::GetNumInputs() const {
        return num_inputs_;
    }

    void Controller::PrintConfigNames() const {
        for (const auto & name : pin_model_.names) {
            std::cout << name << std::endl;
        }
    }

    // TODO: Make this not hard coded
    Eigen::VectorXd Controller::ConvertMujocoConfigToPinocchio(const mjData* data) const {
        Eigen::VectorXd q = Eigen::VectorXd::Zero(pin_model_.nq);

        // re-organize joints
        for (int i = 0; i < 3; i++) {
            q(i+7) = data->qpos[i+10];
            q(i+10) = data->qpos[i+7];
            q(i+13) = data->qpos[i+16];
            q(i+16) = data->qpos[i+13];
        }

        // floating base position
        for (int i = 0; i < 3; i++) {
            q(i) = data->qpos[i];
        }

        // floating base quaternion, note pinocchio uses (x,y,z,w) and mujoco uses (w,x,y,z)
        q(6) = data->qpos[3];
        q(3) = data->qpos[4];
        q(4) = data->qpos[5];
        q(5) = data->qpos[6];


        return q;
    }

    // TODO: Make this not hard coded
    Eigen::VectorXd Controller::ConvertMujocoVelToPinocchio(const mjData* data) const {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(pin_model_.nv);

        for (int i = 0; i < 3; i++) {
            v(i+6) = data->qvel[i+9];
            v(i+9) = data->qvel[i+6];
            v(i+12) = data->qvel[i+15];
            v(i+15) = data->qvel[i+12];
        }

        for (int i = 0; i < 6; i++) {
            v(i) = data->qvel[i];
        }

        return v;
    }

    void Controller::UpdateContacts(const mjModel* model, const mjData* data) {
        for (auto && i : in_contact_) {
            i = false;
        }

        for (int i = 0; i < data->ncon; i++) {
            for (int j = 0; j < mujoco_bodies_.size(); j++) {
                if (model->geom_bodyid[data->contact[i].geom2] == mujoco_bodies_.at(j)) {
                    in_contact_.at(j) = true;
                }
            }
        }
    }

    // TODO: Make this not hard coded
    Eigen::VectorXd Controller::ConvertPinocchioJointToMujoco(const Eigen::VectorXd& joints) {
        Eigen::VectorXd mujoco_joints(joints.size());
        for (int i = 0; i < 3; i++) {
            mujoco_joints(i) = joints(i+3);
            mujoco_joints(i+3) = joints(i);
            mujoco_joints(i+6) = joints(i+9);
            mujoco_joints(i+9) = joints(i+6);
        }

        return mujoco_joints;
    }
} // simulator