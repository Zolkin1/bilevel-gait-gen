//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_MPC_H
#define BILEVEL_GAIT_GEN_MPC_H

#include <memory>

#include <Eigen/Core>

#include "centroidal_model.h"
#include "qp_interface.h"
#include "trajectory.h"
#include "qp_data.h"

namespace mpc {
    /**
     * This is the MPC class.
     * TODO: Eventually make this multi threaded and check for cache friendliness
     */

    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;

    struct MPCInfo {
        int num_nodes;
        double time_horizon;
        int num_qp_iterations;
        int num_contacts;
        double friction_coef;
        vector_t vel_bounds;
        vector_t joint_bounds;
        std::vector<std::string> ee_frames;
        int discretization_steps;
        int num_switches;
        double integrator_dt;

        MPCInfo();
        MPCInfo(const MPCInfo& info);
    };

    enum Gaits {
        Trot = 0,
        Amble = 1,
        Static_Walk = 2
    };

    class MPC final {
    public:
        MPC(const MPCInfo& info, const std::string& robot_urdf);

        Trajectory Solve(const vector_t& centroidal_state, double init_time);

        void SetWarmStartTrajectory(const Trajectory& trajectory);

        // TODO: Support gauss newton on the cost. For now just accept quadratic cost
        void SetQuadraticCostTerm(const matrix_t& Q);

        void SetLinearCostTerm(const vector_t& w);

        void SetQuadraticFinalCost(const matrix_t& Phi);

        void SetLinearFinalCost(const vector_t& w);

        void AddQuadraticTrackingCost(const vector_t& state_des, const matrix_t& Q);

        // TODO: Support gauss newton on the cost. For now just accept quadratic cost
        void AddHessianApproxCost(const vector_t& state, double time, int node);

        void AddGradientCost(const vector_t& state, double time, int node);

        void AddFinalCost();

        /**
         * Creates a default switching time vector for use in initialization
         */
        static std::vector<std::vector<double>> CreateDefaultSwitchingTimes(int num_switches, int num_ee, double horizon);

        void SetDefaultGaitTrajectory(Gaits gait, int num_polys, const std::array<std::array<double, 3>, 4>& ee_pos);

        void SetStateTrajectoryWarmStart(const std::vector<vector_t>& states);

        // Gets the first target config. Assumes solve has already been called.
        // Currently, does not return any momentum info
        vector_t GetTargetConfig() const;

        // Gets first target velocity. Gets the target body velocities and joints
        vector_t GetTargetVelocity() const;

        const CentroidalModel& GetModel() const;

    protected:
    private:
        // ---------------- Private Member Functions ---------------- //
        // Assumes flat ground and constant coef of friction
        // Since we have flat ground and a constant coefficient of friction, we can just make one pyramid
        void SetFrictionPyramid();

        void ResetQPMats();

        void AddDynamicsConstraints(const vector_t& state, double time, int node);

        void AddFKConstraints(const vector_t& state, double time, int node);

        void AddForceConstraints();

        void AddPositivityConstraints(double time, int node);

        void AddBoxConstraints(const vector_t& state, double time, int node);

        void AddGroundIntersectConstraints();

        void AddFrictionConeConstraints(const vector_t& state, double time, int node);

        void AddInequalityConstraints(const vector_t& state, double time, int node);

        int GetForceSplineStartIdx() const;

        int GetVelocityIndex(int node) const;

        int GetJointIndex(int node) const;

        int GetPosSplineStartIdx() const;

        Trajectory ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const;

        void UpdateQPSizes();

        // Temp functions
        void PrintDynamicsConstraints() const;
        void PrintEqualityConstraints() const;
        void PrintInequalityConstraints() const;

        // ---------------- Member Variables ---------------- //
        // Centroidal model
        CentroidalModel model_;

        QPData data_;

        // MPC info
        const MPCInfo info_;

        // TODO: make const
        int num_states_;    // number of states in the MPC model, not in the underlying pinocchio model
        int num_joints_;
        int num_ee_;

        // force parameterization

        // position parameterization

        // friction pyramid
        Eigen::Matrix<double, 4, 3> friction_pyramid_;

        // QP Interface
        std::unique_ptr<QPInterface> qp_solver;

        // previous trajectory
        Trajectory prev_traj_;

        // Cost terms
        matrix_t Q_;
        vector_t w_;
        matrix_t Phi_;
        vector_t Phi_w_;

        // index helpers
        int cone_constraint_start_;
        int box_constraint_start_;
        int positive_force_start_;
        int ground_positive_start_;
        int foot_on_ground_start_;


        // constants
        static int constexpr POS_VARS = 3;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_MPC_H
