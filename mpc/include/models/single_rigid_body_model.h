//
// Created by zolkin on 1/18/24.
//

#ifndef BILEVEL_GAIT_GEN_SINGLE_RIGID_BODY_MODEL_H
#define BILEVEL_GAIT_GEN_SINGLE_RIGID_BODY_MODEL_H

#include <string>
#include <vector>

#include "models/model.h"
#include "trajectory.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using vector_3t = Eigen::Vector3d;

    // TODO: Make in the Model namespace
    using tan_state_t = Eigen::Vector<double, 12>;
    using man_state_t = Eigen::Vector<double, 13>;

    using matrix_t =  Eigen::MatrixXd;
    using vector6_t = Eigen::Vector<double, 6>;
    using matrix_3t = Eigen::Matrix3Xd;
    using matrix_33t = Eigen::Matrix3d;

    // TODO: How to calculate the rotational inertia
    class SingleRigidBodyModel : public Model {
    public:
        SingleRigidBodyModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                             int discretization_steps, double dt, const vector_t& nom_state);

        void GetLinearDynamics(const vector_t& state,
                               const vector_t& ref_state,
                               const Trajectory& traj,
                               double dt,
                               double time,
                               matrix_t& A, matrix_t& B, vector_t& C, vector_t& C2);

//        vector_t GetDiscreteDynamics(const vector_t& state, const Trajectory& traj, double time,
//                                     const vector_t& ref_state);

        tan_state_t CalcDynamics(const vector_t& state, const Trajectory& traj, double time,
                              const vector_t& ref_state);

        vector_3t GetCOMPosition(const vector_t& state) const override;

        static vector_3t ConvertManifoldToTangentQuat(const Eigen::Vector4d& state,
                                               const Eigen::Vector4d& ref_state) ;

        int GetNumTangentStates() const override;
        int GetNumManifoldStates() const override;
        vector_t ConvertManifoldStateToTangentState(const vector_t& state, const vector_t& ref_state) const override;
        vector_t ConvertTangentStateToManifoldState(const vector_t& state, const vector_t& ref_state) const override;

        vector_3t GetCOMToHip(int end_effector) const;

        vector_t InverseKinematics(const man_state_t& state, const std::vector<vector_3t>& end_effector_location,
                                   const vector_t& state_guess, const vector_t& joint_limits_ub,
                                   const vector_t& joint_limits_lb);

//        vector_t VelocityInverseKinematics(const vector_t& config,
//                                           const std::vector<vector_3t>& end_effector_velocity,
//                                           const vector_t& vel_guess);

        std::vector<vector_3t> GetEndEffectorLocations(const vector_t& q) override;

        void ComputeLinearizationPartialWrtContactTimes(matrix_t& dA,
                                                        matrix_t& dB,
                                                        vector_t& dC,
                                                        const vector_t& state,
                                                        const Trajectory& traj,
                                                        double time,
                                                        int end_effector,
                                                        int contact_time_idx);

        static constexpr int QUAT_SIZE = 4;
        static constexpr int QUAT_START = 6;
        static constexpr int ORIENTATION_START = 6;
        static constexpr int ANG_VEL_START = 9;
        static constexpr int LIN_MOM_START = 3;
        static constexpr int POS_START = 0;
    protected:
        void ConvertMPCStateToPinocchioState(const vector_t& state, Eigen::Ref<vector_t> q_pin) const override;

        pinocchio::SE3 GetSE3Error(int joint_id, const pinocchio::SE3& des_state);

        void ComputeJacobianForIK(const vector_t& q, const pinocchio::SE3& error,
                                  int id, Eigen::Matrix<double, 6, Eigen::Dynamic>& J,
                                  bool is_frame);

        const int num_tangent_states_;
        const int num_manifold_states_;

        matrix_33t Ir_;   // Rotational inertia
        matrix_33t Ir_inv_;
    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_SINGLE_RIGID_BODY_MODEL_H
