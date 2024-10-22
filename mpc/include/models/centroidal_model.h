//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
#define BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H

#include <memory>

#include <Eigen/Core>
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"

//#include "integrator.h"
#include "rk_integrator.h"
#include "timer.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;
    using vector6_t = Eigen::Vector<double, 6>;
    using matrix_3t = Eigen::Matrix3Xd;

    class CentroidalModel {
        // Note: I am using tait-bryan z-y-x angles
    public:
        // Constants
        static int constexpr FLOATING_BASE_OFFSET = 7;      // 7 for quaternion, 6 for euler
        static int constexpr FLOATING_VEL_OFFSET = 6;
        static int constexpr MOMENTUM_OFFSET = 6;
        static int constexpr POS_VARS = 3;

        explicit CentroidalModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                                 int discretization_steps, double dt);

        /**
         * @param state
         * @param input
         * @param node
         * @param time
         * @return
         */
//        void GetLinearDiscreteDynamics(const vector_t& state, const vector_t& ref_state, const Inputs& input, double time,
//                                           matrix_t& A, matrix_t& B, vector_t& C);

        // See computeForwardKinematicsDerivatives
        // Needs to return with the frame transition built in
//        void GetFKLinearization(const vector_t& state, const vector_t& ref_state, const Inputs& input, int end_effector,
//                                    matrix_t& A, vector_t& C);

        matrix_3t GetFKJacobianForEndEffector(const vector_t& state, const std::string& frame, bool compute_jac);

        Eigen::Vector3d GetEndEffectorLocationCOMFrame(const vector_t& state,
                                                       const std::string& frame) const;

        int GetPinocchioNumConfig() const;

        int GetNumJoints() const;

        int GetNumEndEffectors() const;

//        void SetDynamicsRefState(const vector_t& state);

//        vector_t CalcDynamics(const vector_t& state, const Inputs& input, double time,
//                              const vector_t& ref_state);

        static Eigen::Vector4d ConvertZYXRotToQuaternion(const Eigen::Vector3d& zyx_rot);
        static Eigen::Vector3d ConvertQuaternionToZYXRot(const Eigen::Vector4d& quat);

        // Note: State expected to have pinocchio/eigen quaternion representation
        static vector_t ConvertManifoldStateToAlgebraState(const vector_t& state, const vector_t& ref_state);
        static vector_t ConvertAlgebraStateToManifoldState(const vector_t& state, const vector_t& ref_state);

        // Take in MPC state
        Eigen::Vector3d GetCOMPosition(const vector_t& state) const;

        vector_t ComputeBaseVelocities(const vector_t& state, const vector_t& vel) const;

        double GetMass() const;

//        vector_t GetDiscreteDynamics(const vector_t& state, const Inputs& input, double time,
//                                     const vector_t& ref_state);

        const std::string& GetEndEffectorFrame(int ee) const;

        std::vector<int> GetContactFrames() const;

        /**
         * Computes the partial derivative of the linearization wrt the contact time (given by ee and idx).
         * Computes this partial only for the linearization at the given time and about the given input and state
         * @param dA
         * @param dB
         * @param dC
         * @param state
         * @param inputs
         * @param time
         * @param ee
         * @param idx
         */
//        void ComputeLinearizationPartialWrtContactTimes(matrix_t& dA,
//                                                        matrix_t& dB,
//                                                        vector_t& dC,
//                                                        const vector_t& state,
//                                                        const Inputs& inputs,
//                                                        double time,
//                                                        int ee, int idx);

    protected:
    private:

        vector_t ComputePinocchioVelocities(const vector_t& state, const vector_t& joint_vels, const matrix_t& CMMbinv,
                                            const matrix_t& CMMj) const;

        // TODO: Should potentially templatize
        void ConvertMPCStateToPinocchioState(const vector_t& state, Eigen::Ref<vector_t> q_pin) const;
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

//        vector_t ref_state_;

        std::map<std::string, int> frame_map_;
        std::vector<std::string> frames_;

        // integrator to be used
        std::unique_ptr<RKIntegrator> integrator_;

        const Eigen::Vector3d GRAVITY;

        // Pre-dynamically allocated values
        matrix_t Ac_, Bc_;
        vector_t Cc_, Cc2_;
        vector_t xdot_;
        vector_t q_pin_;
    };
}


#endif //BILEVEL_GAIT_GEN_CENTROIDAL_MODEL_H
