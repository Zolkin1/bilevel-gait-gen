//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
#define BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H

#include <memory>

#include <Eigen/Core>
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"

#include "integrator.h"
#include "inputs.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;

    class CentroidalModel {
    public:
        explicit CentroidalModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                                 int discretization_steps);

        matrix_t GetLinearDiscreteDynamicsState(const vector_t& state, const Inputs& input, int node, double time);

        matrix_t GetLinearDiscreteDynamicsInput(const vector_t& state, const vector_t& input);

        vector_t GetConstantDiscreteDynamics(const vector_t& state, const vector_t& input);

        // See computeForwardKinematicsDerivatives
        // Needs to return with the frame transition built in
        matrix_t GetFKLinearization(const vector_t& state, const vector_t& input);

        vector_t GetFKJacobianForEndEffector(const vector_t& state, const std::string& frame, bool compute_jac);

        Eigen::Vector3d GetEndEffectorLocationCOMFrame(const vector_t& state, const std::string& frame) const;

        int GetNumConfig() const;

        int GetNumJoints() const;

        int GetNumEndEffectors() const;

        vector_t CalcDynamics(const vector_t& state, const Inputs& input, double time) const;

        // conversion from full state to centroidal state
    protected:
    private:
        void CreateFrameMap(const std::vector<std::string>& frames);

        // pinocchio model
        pinocchio::Model pin_model_;

        // pinocchio data
        std::unique_ptr<pinocchio::Data> pin_data_;

        double robot_mass_;

        int num_total_states_;     // [h, qb, qj]
        int num_joints_;
        int num_ee_;

        int discretization_steps_;

        vector_t prev_fk_q_;

        std::map<std::string, int> frame_map_;
        std::vector<std::string> frames_;

        // integrator to be used
        std::unique_ptr<Integrator> integrator_;

        // Constants
        static int constexpr FLOATING_BASE_OFFSET = 7;
        static int constexpr FLOATING_VEL_OFFSET = 6;
        static int constexpr MOMENTUM_OFFSET = 6;

        const Eigen::Vector3d GRAVITY;
    };
}


#endif //BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
