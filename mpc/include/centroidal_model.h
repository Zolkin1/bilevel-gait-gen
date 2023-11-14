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
    using vector6_t = Eigen::Vector<double, 6>;

    class CentroidalModel {
        // Note: I am using tait-bryan z-y-x angles
    public:
        explicit CentroidalModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                                 int discretization_steps, double dt);

        /**
         * The state given should be: [h, qb, qj]. h is the linear and angular momenetum at the COM
         * qb is the floating base. x,y,z then a quaternion in pinocchio form. qj is the joint configurations.
         * Then the pinocchio configuration vector is the same as the configuration vector here.
         * The MPC decision variables will be [h, qb_alg, qj]. qb_alg is the floating base in the lie algebra.
         * @param state
         * @param input
         * @param node
         * @param time
         * @return
         */
         // TODO: remove node argument
        void GetLinearDiscreteDynamics(const vector_t& state, const vector_t& ref_state, const Inputs& input, int node, double time,
                                           matrix_t& A, matrix_t& B, vector_t& C);

        // See computeForwardKinematicsDerivatives
        // Needs to return with the frame transition built in
        matrix_t GetFKLinearization(const vector_t& state, const vector_t& input);

        matrix_t GetFKJacobianForEndEffector(const vector_t& state, const std::string& frame, bool compute_jac);

        Eigen::Vector3d GetEndEffectorLocationCOMFrame(const vector_t& state, const std::string& frame) const;

        int GetPinocchioNumConfig() const;

        int GetNumJoints() const;

        int GetNumEndEffectors() const;

        vector_t CalcDynamics(const vector_t& state, const Inputs& input, double time) const;

        static Eigen::Vector4d ConvertZYXRotToQuaternion(const Eigen::Vector3d& zyx_rot);
        static Eigen::Vector3d ConvertQuaternionToZYXRot(const Eigen::Vector4d& quat);

        // Note: State expected to have pinocchio/eigen quaternion representation
        static vector_t ConvertManifoldStateToAlgebraState(const vector_t& state, const vector_t& ref_state);
        static vector_t ConvertAlgebraStateToManifoldState(const vector_t& state, const vector_t& ref_state);

        // conversion from full state to centroidal state
    protected:
    private:
        // Constants
        static int constexpr FLOATING_BASE_OFFSET = 7;      // 7 for quaternion, 6 for euler
        static int constexpr FLOATING_VEL_OFFSET = 6;
        static int constexpr MOMENTUM_OFFSET = 6;
        static int constexpr POS_VARS = 3;

        vector_t ComputePinocchioVelocities(const vector_t& state, const vector_t& joint_vels, const matrix_t& CMMbinv,
                                            const matrix_t& CMMj);

        vector_t ConvertMPCStateToPinocchioState(const vector_t& state) const;
        vector_t ConvertPinocchioStateToMPCState(const Eigen::Vector<double, MOMENTUM_OFFSET>& momentum,
                                                                  const vector_t& state) const;

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

        const Eigen::Vector3d GRAVITY;
    };
}


#endif //BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
