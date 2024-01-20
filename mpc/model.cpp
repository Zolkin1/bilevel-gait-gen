//
// Created by zolkin on 1/18/24.
//

#include "pinocchio/algorithm/center-of-mass.hpp"

#include "model.h"

namespace mpc {
    Model::Model(const std::string& robot_urdf, const std::vector<std::string>& frames, int discretization_steps,
                 double dt, bool uses_joints, const std::vector<Constraints>& constraints) :
            GRAVITY(0., 0., -9.81), num_ee_(frames.size()), uses_joints_(uses_joints),
            constraints_(constraints) {

        integrator_ = std::make_unique<RKIntegrator>(dt);

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

} // mpc