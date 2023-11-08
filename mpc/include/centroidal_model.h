//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
#define BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H

#include <memory>

#include <Eigen/Core>
#include "pinocchio/parsers/urdf.hpp"

#include "integrator.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;

    class CentroidalModel {
    public:
        explicit CentroidalModel(const std::string& robot_urdf);

        matrix_t GetLinearDiscreteDynamicsState(const vector_t& state, const vector_t& input);

        matrix_t GetLinearDiscreteDynamicsInput(const vector_t& state, const vector_t& input);

        vector_t GetConstantDiscreteDynamics(const vector_t& state, const vector_t& input);

        // See computeForwardKinematicsDerivatives
        // Needs to return with the frame transition built in
        matrix_t GetFKLinearization(const vector_t& state, const vector_t& input);

        vector_t GetFK(const vector_t& state);

        int GetNumConfig() const;

        int GetNumJoints() const;

        int GetNumEndEffectors() const;

        // conversion from full state to centroidal state
    protected:
    private:
        // pinocchio model
        pinocchio::Model pin_model_;

        // pinocchio data
        std::unique_ptr<pinocchio::Data> pin_data_;

        int num_joints_;
        int num_ee_;

        // integrator to be used

        // Constants
        static int constexpr FLOATING_BASE_OFFSET = 7;
        static int constexpr FLOATING_VEL_OFFSET = 6;
    };
}


#endif //BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
