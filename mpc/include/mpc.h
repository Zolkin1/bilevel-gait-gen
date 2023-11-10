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

        MPCInfo();
        MPCInfo(const MPCInfo& info);
    };

    // mpc data, stored as a struct to be cache friendly
    // - gradients and hessians
    // - constant terms
    // - warm start
    // TODO: is this really the struct I want?
    struct MPCData {
        matrix_t dynamics_constraints;
        vector_t dynamics_constants;

        matrix_t equality_constraints;
        vector_t equality_constants;

        matrix_t inequality_constraints;
        vector_t inequality_constants_ub;
        vector_t inequality_constants_lb;

        matrix_t cost_quadratic;
        vector_t cost_linear;
    };

    class MPC final {
    public:
        MPC(const MPCInfo& info, const std::string& robot_urdf);

        Trajectory Solve(const vector_t& centroidal_state);

        void SetWarmStartTrajectory(const Trajectory& trajectory);

        matrix_t GetCostHessianApprox(const vector_t& state, const vector_t& input);

        vector_t GetCostGradient(const vector_t& state, const vector_t& input);
    protected:
    private:
        // ---------------- Private Member Functions ---------------- //
        // Assumes flat ground and constant coef of friction
        // Since we have flat ground and a constant coefficient of friction, we can just make one pyramid
        void SetFrictionPyramid();

        /**
         * Creates a default switching time vector for use in initialization
         */
        static std::vector<std::vector<double>> CreateDefaultSwitchingTimes(int num_switches, int num_ee, double horizon);

        // ---------------- Member Variables ---------------- //
        // Centroidal model
        CentroidalModel model_;

        MPCData data_;

        // MPC info
        const MPCInfo info_;

        int dynamics_constraints_;
        int equality_constraints_;
        int inequality_constraints_;
        int decision_vars_;

        int num_states_;    // number of states in the MPC model, not in the underlying pinocchio model
        int num_inputs_;    // number of inputs in the MPC model, now in the underlying pinocchio model
        int num_joints_;
        int num_ee_;

        int force_start_idx_;
        int pos_start_idx_;
        int vel_start_idx;
        int com_start_idx_;
        int config_start_idx_;

        // force parameterization

        // position parameterization

        // friction pyramid
        Eigen::Matrix<double, 4, 3> friction_pyramid_;

        // QP Interface
        //QPInterface qp_solver;

        // previous trajectory
        Trajectory prev_traj_;

        // constants
        static int constexpr POS_VARS = 3;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_MPC_H
