//
// Created by zolkin on 1/18/24.
//

#ifndef BILEVEL_GAIT_GEN_SINGLE_RIGID_BODY_MODEL_H
#define BILEVEL_GAIT_GEN_SINGLE_RIGID_BODY_MODEL_H

#include <string>
#include <vector>

#include "model.h"
#include "trajectory.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using vector_3t = Eigen::Vector3d;
    using matrix_t =  Eigen::MatrixXd;
    using vector6_t = Eigen::Vector<double, 6>;
    using matrix_3t = Eigen::Matrix3Xd;
    using matrix_33t = Eigen::Matrix3d;

    // TODO: How to calculate the rotational inertia
    class SingleRigidBodyModel : public Model {
    public:
        SingleRigidBodyModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                             int discretization_steps, double dt);

        void GetLinearDiscreteDynamics(const vector_t& state, const vector_t& ref_state, const Trajectory& traj,
                                       double time, matrix_t& A, matrix_t& B, vector_t& C);

//        vector_t GetDiscreteDynamics(const vector_t& state, const Trajectory& traj, double time,
//                                     const vector_t& ref_state);

        void ComputeLinearizationPartialWrtContactTimes(matrix_t& dA,
                                                        matrix_t& dB,
                                                        vector_t& dC,
                                                        const vector_t& state,
                                                        const Trajectory& traj,
                                                        double time,
                                                        int ee, int idx);

        vector_t CalcDynamics(const vector_t& state, const Trajectory& traj, double time,
                              const vector_t& ref_state);

        vector_3t GetCOMPosition(const vector_t& state) const override;

        static vector_3t ConvertManifoldToTangentQuat(const Eigen::Vector4d& state,
                                               const Eigen::Vector4d& ref_state) ;

        int GetNumTangentStates() const override;
        int GetNumManifoldStates() const override;
        vector_t ConvertManifoldStateToTangentState(const vector_t& state, const vector_t& ref_state) const override;
        vector_t ConvertTangentStateToManifoldState(const vector_t& state, const vector_t& ref_state) const override;

    protected:
        void ConvertMPCStateToPinocchioState(const vector_t& state, Eigen::Ref<vector_t> q_pin) const override;

        const int num_tangent_states_;
        const int num_manifold_states_;

        matrix_33t Ir_;   // Rotational inertia
        matrix_33t Ir_inv_;

       static constexpr int QUAT_SIZE = 4;
       static constexpr int QUAT_START = 6;
    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_SINGLE_RIGID_BODY_MODEL_H
