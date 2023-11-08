//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/parsers/urdf.hpp"

#include "centroidal_model.h"


namespace mpc {
    CentroidalModel::CentroidalModel(const std::string &robot_urdf) {
        // create the pinocchio model - always a free flyer
        pinocchio::urdf::buildModel(robot_urdf, pinocchio::JointModelFreeFlyer(), pin_model_, false);

        // create the pinocchio data
        pin_data_ = std::make_unique<pinocchio::Data>(pin_model_);

        num_joints_ = pin_model_.nq - FLOATING_BASE_OFFSET;


        // TODO: Get the number of end effectors and their data
        num_ee_ = 0;
    }

    // TODO: Add in frame
    vector_t CentroidalModel::GetFK(const vector_t& state) {
        // Forward kinematics

        // Get a specific frame
    }

    int CentroidalModel::GetNumConfig() const {
        return pin_model_.nq;
    }

    int CentroidalModel::GetNumJoints() const {
        return num_joints_;
    }

    int CentroidalModel::GetNumEndEffectors() const {
        return num_ee_;
    }

} // mpc