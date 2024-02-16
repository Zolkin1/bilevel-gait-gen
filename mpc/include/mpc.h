//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_MPC_H
#define BILEVEL_GAIT_GEN_MPC_H

#include <memory>

#include <Eigen/Core>

#include "models/centroidal_model.h"
#include "qp/osqp_interface.h"
#include "trajectory.h"
#include "qp/qp_data.h"
#include "controller.h"
#include "timer.h"
#include "gait_optimizer.h"
#include "models/single_rigid_body_model.h"
#include "rk_integrator.h"
#include "qp/clarabel_interface.h"

namespace mpc {
    /**
     * This is the MPC class.
     * TODO: Eventually make this multi threaded and check for cache friendliness
     */

    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;

    enum MPCVerbosityLevel {
        Nothing = 0,
        Timing = 1,
        Optimization = 2,
        All = 3
    };

    struct MPCInfo {
        int num_nodes;
        int num_qp_iterations;
        int num_contacts;
        double friction_coef;
        vector_t vel_bounds;
        vector_t joint_bounds_lb;
        vector_t joint_bounds_ub;
        std::vector<std::string> ee_frames;
        int discretization_steps;
        int num_switches;
        double integrator_dt;
        double force_bound;
        double swing_height;
        double foot_offset;
        vector_t nom_state;
        Eigen::Vector2d ee_box_size;
        int real_time_iters;
        MPCVerbosityLevel verbose;

        MPCInfo();
        MPCInfo(const MPCInfo& info);
    };

    enum Gaits {
        Trot = 0,
        Amble = 1,
        Static_Walk = 2
    };

    class MPC {
    public:
        MPC(const MPCInfo& info, const std::string& robot_urdf);

        // TODO: Figure out how to not explicitly pass the ee location
        Trajectory CreateInitialRun(const vector_t& state, const std::vector<vector_3t>& ee_start_locations);

        Trajectory GetRealTimeUpdate(const vector_t& state, double init_time,
                                     const std::vector<vector_3t>& ee_start_locations,
                                     bool high_quality);

        virtual Trajectory Solve(const vector_t& state, double init_time,
                                 const std::vector<vector_3t>& ee_start_locations) = 0;

        void SetWarmStartTrajectory(const Trajectory& trajectory);

        // TODO: Support gauss newton on the cost. For now just accept quadratic cost
        void SetQuadraticCostTerm(const matrix_t& Q);

        void SetLinearCostTerm(const vector_t& w);

        void SetQuadraticFinalCost(const matrix_t& Phi);

        void SetLinearFinalCost(const vector_t& w);

        void AddQuadraticTrackingCost(const vector_t& state_des, const matrix_t& Q);

        /**
         * Creates a default switching time vector for use in initialization
         */
        static std::vector<std::vector<double>> CreateDefaultSwitchingTimes(int num_switches, int num_ee, double horizon);

        void SetDefaultGaitTrajectory(Gaits gait, int num_polys, const std::vector<vector_3t>& ee_pos);

        void SetStateTrajectoryWarmStart(const std::vector<vector_t>& states);

        // Gets the first target config. Assumes solve has already been called.
        // Currently, does not return any momentum info
        vector_t GetTargetConfig(double time) const;

        const Model* GetModel() const;

        void AddForceCost(double weight);

        void PrintStats();

        // TODO: Do I really need this?
        controller::Contact GetDesiredContacts(double time) const;

//        virtual vector_t GetFullTargetState(double time) = 0;

        Trajectory GetTrajectory() const;

        vector_t GetForceTarget(double time) const;

        int GetNode(double time) const;

        int GetNumDecisionVars() const;

        int GetNumConstraints() const;

        vector_t Getdx();

        bool ComputeDerivativeTerms();

        bool GetQPPartials(QPPartials& partials) const;

        /**
         * Computes how each term in the QP changes wrt the impact time given by ee and idx
         * @param partials
         * @param ee
         * @param idx
         */
//        virtual bool ComputeParamPartials(const Trajectory& traj, QPPartials& partials, int ee, int idx) = 0;

        vector_t GetQPSolution() const;

        virtual std::vector<std::vector<Eigen::Vector3d>> CreateVizData() = 0;

        void UpdateContactTimes(std::vector<time_v>& contact_times);

        const QPData& GetQPData() const;

    protected:
        // ---------------- Protected Member Functions ---------------- //
        // Assumes flat ground and constant coef of friction
        // Since we have flat ground and a constant coefficient of friction, we can just make one pyramid
        void SetFrictionPyramid();

        virtual void AddDynamicsConstraints(const vector_t& state) = 0;

        void AddFKConstraints(const vector_t& state);

        void AddFKConstraintsProjection(const vector_t& state);

        void AddBoxConstraints();

        void AddForceBoxConstraints();

        void AddForceBoxConstraintPartials(utils::SparseMatrixBuilder& builder, int contact_idx,
                                           int start_idx, int ee);

        void AddFrictionConeConstraints();

        void AddFrictionConeConstraintPartials(utils::SparseMatrixBuilder& builder, int contact_idx,
                                               int start_idx, int ee);

        virtual int GetForceSplineStartIdx() const = 0;

        virtual int GetPosSplineStartIdx() const = 0;

        Trajectory ConvertQPSolToTrajectory(const vector_t& qp_sol, const vector_t& init_state) const;

        virtual vector_t ConvertTrajToQPVec(const Trajectory& traj) const = 0;

        double LineSearch(const vector_t& direction, const vector_t& init_state);

        double GetMeritValue(const vector_t& x, double mu, const vector_t& init_state);

        double GetMeritValue(const Trajectory& traj, double mu, const vector_t& init_state);

        double GetCostValue(const vector_t& x) const;

        vector_t GetEqualityConstraintValues(const Trajectory& traj, const vector_t& init_state);

        // Gets the time at a node
        double GetTime(int node) const;

        double GetMeritGradient(const vector_t& x, const vector_t& p, double mu, const vector_t& init_state);

        void RecordStats(double alpha, const vector_t& direction, const SolveQuality& solve_type,
                         const vector_t& ref_state, double solve_time, double cost);

        int GetNodeIntersectMutableForces() const;

        virtual void SetInitQPSizes() = 0;

        void UpdateQPSizes();

        // TODO: Support gauss newton on the cost. For now just accept quadratic cost
        void AddHessianApproxCost();

        void AddGradientCost();

        void AddFinalCost();

        int UpdateNumInputs();

        void AddDiagonalCost();

        // ---------------- Member Variables ---------------- //
        // Centroidal model
        SingleRigidBodyModel model_;

        QPData data_;

        // MPC info
        const MPCInfo info_;

        const int num_states_;    // number of states in the MPC model, not in the underlying pinocchio model
        const int num_ee_;
        int num_inputs_;

        // friction pyramid
        Eigen::Matrix<double, 4, 3> friction_pyramid_;

        // QP Interface
        std::unique_ptr<ClarabelInterface> qp_solver;

        // previous trajectory
        Trajectory prev_traj_;

        // Cost terms
        matrix_t Q_;
        vector_t w_;
        matrix_t Phi_;
        vector_t Phi_w_;

        matrix_t Q_forces_;

        vector_t prev_qp_sol;

        double init_time_;

        double force_weight_;

        // constants
        static int constexpr POS_VARS = 3;


        int run_num_;

        vector_t line_search_res_;

        vector_t prev_dual_sol_;

        double mu_;

        int num_ineq_fk_;

        std::vector<double> equality_constraint_violations_;
        std::vector<double> step_norm_;
        std::vector<double> alpha_;
        std::vector<double> cost_result_;
        std::vector<double> merit_result_;
        std::vector<double> merit_directional_deriv_;
        std::vector<SolveQuality> solve_type_;
        std::vector<vector_t> ref_state_;
        std::vector<double> solve_time_;
        std::vector<double> cost_;
        std::vector<int> contact_sched_change_;

        bool in_real_time_;

        int constraint_idx_;

        matrix_t A_, B_;
        vector_t C_, C2_;

        const bool constraint_projection_;

        const bool using_clarabel_;

        RKIntegrator integrator_;
    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_MPC_H
