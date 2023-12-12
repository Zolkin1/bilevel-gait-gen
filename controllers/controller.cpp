//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <string>
#include <utility>

#include "pinocchio/algorithm/kinematics.hpp"

#include "include/controller.h"

namespace controller {
    int Contact::GetNumContacts() const {
        int num_contacts = 0;
        for (const bool i : in_contact_) {
            if (i) {
                num_contacts++;
            }
        }

        return num_contacts;
    }

    Controller::Controller(double control_rate, std::string robot_urdf, const std::string& foot_type) : rate_(1.0/control_rate),
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

        if (foot_type == "FLAT_FEET") {
            CONSTRAINT_PER_FOOT = CONSTRAINT_PER_FLAT_FOOT;
        } else if (foot_type == "POINT_FEET") {
            CONSTRAINT_PER_FOOT = CONSTRAINT_PER_POINT_FOOT;
        } else {
            std::cerr << "Foot type not found. Using point feet." << std::endl;
        }
    }

    void Controller::InitSolver() {}

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

    const pinocchio::Model& Controller::GetPinocchioModel() const {
        return pin_model_;
    }

    void Controller::UpdateTargetConfig(const Eigen::VectorXd& q) {
        assert(q.size() == config_target_.size());
        config_target_ = q;
    }

    void Controller::UpdateTargetVel(const Eigen::VectorXd& v) {
        assert(v.size() == vel_target_.size());
        vel_target_ = v;
    }

    void Controller::UpdateTargetAcc(const Eigen::VectorXd& a) {
        assert(a.size() == acc_target_.size());
        acc_target_ = a;
    }

    void Controller::AssignPositionControl(Eigen::VectorXd& control) {
        for (int i = 0; i < num_inputs_; i++) {
            control(i) = config_target_(i + FLOATING_BASE_OFFSET);
        }
    }

    void Controller::AssignVelocityControl(Eigen::VectorXd& control) {
        for (int i = num_inputs_; i < 2*num_inputs_; i++) {
            control(i) = vel_target_(i - num_inputs_ + FLOATING_VEL_OFFSET);
        }
    }

    Contact::Contact() {}

    Contact::Contact(int num_contacts) : in_contact_(num_contacts, true), contact_frames_(num_contacts, 0) {}

} // controller