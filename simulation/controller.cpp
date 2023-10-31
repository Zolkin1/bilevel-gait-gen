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

    Eigen::VectorXd Controller::ConvertMujocoConfigToPinocchio(const mjData* data) const {
        Eigen::VectorXd q = Eigen::VectorXd::Zero(pin_model_.nq);

        // Floating base config
        // floating base position
        for (int i = 0; i < 3; i++) {
            q(i) = data->qpos[i];
        }

        // floating base quaternion, note pinocchio uses (x,y,z,w) and mujoco uses (w,x,y,z)
        q(6) = data->qpos[3];
        q(3) = data->qpos[4];
        q(4) = data->qpos[5];
        q(5) = data->qpos[6];


        // Joints
        for (int i = 0; i < pin_model_.nv - FLOATING_VEL_OFFSET; i++) {
            q(mujoco_to_pinocchio_joint_map_.at(i+1) - 2 + FLOATING_BASE_OFFSET) = data->qpos[i + FLOATING_BASE_OFFSET];
        }


        return q;
    }

    Eigen::VectorXd Controller::ConvertMujocoVelToPinocchio(const mjData* data) const {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(pin_model_.nv);

        // Floating base velocities
        for (int i = 0; i < FLOATING_VEL_OFFSET; i++) {
            v(i) = data->qvel[i];
        }

        // Joint velocities
        for (int i = 0; i < pin_model_.nv - FLOATING_VEL_OFFSET; i++) {
            v(mujoco_to_pinocchio_joint_map_.at(i+1) - 2 + FLOATING_VEL_OFFSET) = data->qvel[i + FLOATING_VEL_OFFSET];
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

    Eigen::VectorXd Controller::ConvertPinocchioJointToMujoco(const Eigen::VectorXd& joints) {
        Eigen::VectorXd mujoco_joints(joints.size());

        for (int i = 0; i < joints.size(); i++) {
            mujoco_joints(i) = joints(mujoco_to_pinocchio_joint_map_.at(i+1) - 2); // -2 removes the two extra pinocchio joints
        }

        return mujoco_joints;
    }

    void Controller::CreateJointMap(const mjModel* model) {
        for (int i = 0; i < model->njnt; i++) {
            for (int j = 0; j < pin_model_.njoints; j++) {
                if (model->names + model->name_jntadr[i] == pin_model_.names.at(j)) {
                    mujoco_to_pinocchio_joint_map_.insert(std::pair<int, int>(i,j));
                }
            }
        }
    }
} // simulator