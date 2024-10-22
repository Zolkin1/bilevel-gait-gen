//
// Created by zolkin on 1/18/24.
//

#ifndef BILEVEL_GAIT_GEN_MODEL_H
#define BILEVEL_GAIT_GEN_MODEL_H

#include <string>
#include <vector>
#include <Eigen/Core>
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"

//#include "rk_integrator.h"
#include "qp/qp_data.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;
    using vector_3t = Eigen::Vector3d;
    using vector6_t = Eigen::Vector<double, 6>;
    using matrix_3t = Eigen::Matrix3Xd;
    using matrix_33t = Eigen::Matrix3d;

    // TODO: Have the model hold which constraints it works with
    class Model {
    public:
        static int constexpr FLOATING_BASE_OFFSET = 7;      // 7 for quaternion, 6 for euler
        static int constexpr FLOATING_VEL_OFFSET = 6;
        static int constexpr MOMENTUM_OFFSET = 6;
        static int constexpr POS_VARS = 3;

        Model(const std::string& robot_urdf, const std::vector<std::string>& frames,
              int discretization_steps, double dt, bool uses_joints, const std::vector<Constraints>& constraints);

        Model(const Model& other);

        Model& operator=(const Model& model);

        // TODO: Make it so that this function takes in a trajectory instead of an input
//        virtual void GetLinearDiscreteDynamics(const vector_t& state, const vector_t& ref_state, const Inputs& input,
//                                       double time, matrix_t& A, matrix_t& B, matrix_t& C) = 0;

        int GetNumEndEffectors() const;

        std::vector<int> GetContactFrames() const;

        double GetMass() const;

        bool UsesJoints() const;

        virtual int GetNumJoints() const;

        const std::vector<Constraints>& GetApplicableConstraints() const;

        virtual Eigen::Vector3d GetCOMPosition(const vector_t& state) const = 0;

        virtual int GetNumTangentStates() const = 0;
        virtual int GetNumManifoldStates() const = 0;

        virtual vector_t ConvertManifoldStateToTangentState(const vector_t& state, const vector_t& ref_state) const = 0;
        virtual vector_t ConvertTangentStateToManifoldState(const vector_t& state, const vector_t& ref_state) const = 0;

        int GetFullModelConfigSpace() const;

        virtual std::vector<Eigen::Vector3d> GetEndEffectorLocations(const vector_t& q) = 0;

        matrix_33t GetEEJacobian(int end_effector, const vector_t& q);
        matrix_33t GetEEJacobianDeriv(int end_effector, const vector_t& q, const vector_t& v);

        matrix_33t GetBaseRotationMatrix(const vector_t& q);

        matrix_33t GetOperationalSpaceInertia(int end_effector, const vector_t& q);
//        matrix_t GetCoriolisMat(const vector_t& q, const vector_t& v);
//        vector_t GetGravityVec(const vector_t& q);
        vector_t GetNonlinearEffects(const vector_t& q, const vector_t& v);

    protected:
        virtual void ConvertMPCStateToPinocchioState(const vector_t& state, Eigen::Ref<vector_t> q_pin) const = 0;

        bool uses_joints_; // was const


        // pinocchio model
        pinocchio::Model pin_model_;

        // pinocchio data
        std::unique_ptr<pinocchio::Data> pin_data_;

        double robot_mass_;

        int num_ee_; // was const

        int discretization_steps_;

        std::map<std::string, int> frame_map_;
        std::vector<std::string> frames_;

        // integrator to be used
//        std::unique_ptr<RKIntegrator> integrator_;

        const Eigen::Vector3d GRAVITY;

        std::vector<Constraints> constraints_; // was const
    private:
        void CreateFrameMap(const std::vector<std::string>& frames);

    };

} // mpc


#endif //BILEVEL_GAIT_GEN_MODEL_H
