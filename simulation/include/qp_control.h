//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_QP_CONTROL_H
#define BILEVEL_GAIT_GEN_QP_CONTROL_H

#include <Eigen/Core>

#include "OsqpEigen/OsqpEigen.h"

#include "controller.h"

/*
 * This controller should work as follows:
 *  --- Before running ---
 *   - Add all the constraints
 *  --- During run ---
 *   - Evaluate QP
 *   - Update constraints (only the ones that change over time)
 *
 *   The exact structure of this may depends on the QP solver.
 *   I will start by using OSQP. May change later.
 *
 *   I will start by just solving one QP with only a few constraints:
 *   - Dynamics
 *   - Holonomic (?)
 *   - Friction cone
 *   - Stationary contacts
 *   - Torque limits
 *
 *  Eventually I can add the kinematic CBF constraint
 *
 *   The costs will be:
 *   - Leg tracking
 *   - Torso acceleration
 *   - Force tracking
 *
 *   decision variables: [qddot, lambda] \in R^{nv + 3*ncontacts}
 *   constraints:
 *   - dynamics \in R^{FLOATING_VEL}
 *   - friction cone (pyramid)  \in R^{4 * ncontacts}
 *   - stationary contacts \in R^{3 * ncontacs}
 *   - torque limits \in R^{njoints}
 *
 *   total constraints: FLOATING_VEL + 4*ncontacts + 3*ncontacts + njoints
 */

namespace simulator {
    class QPControl : public Controller {
    public:
        QPControl(double control_rate, std::string robot_urdf, const std::string& foot_type, int nv,
                  const Eigen::VectorXd& torque_bounds, double friction_coef,
                  std::vector<double> base_pos_gains,
                  std::vector<double> base_ang_gains,
                  std::vector<double> joint_gains,
                  double leg_weight,
                  double torso_weight,
                  double force_weight);

        void InitSolver(const mjModel* model, const mjData* data) override;

        std::vector<mjtNum> ComputeControlAction(const mjModel* model, const mjData* data) override;

        void SetBasePosGains(double kv, double kp);

        void SetBaseAngleGains(double kv, double kp);

        /**
         * Note: assumes the same kp and kv for each joint.
         * @param kv
         * @param kp
         */
        void SetJointGains(double kv, double kp);

    protected:
    private:
        /**
         * Adds all the constraints to the QP
         */
        void AddConstraintsToQP();

        // Functions to add in all the necessary constraints
        /**
         * Computes all the dynamics used by the constraints and cost.
         * All terms are then stored in the pinocchio data struct if possible.
         * Constraint jacobian and its time derivative are stored in the associated local variable.
         * @param q generalized position vector
         * @param v generalized velocity vector
         */
        void ComputeDynamicsTerms(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Eigen::VectorXd& a);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         * When using generalized accelerations and contact forces, I only need to enforce the dynamics in the
         * dimensions where there is no torque input. When every joint is actuated then this just becomes the
         * floating base coordinates.
         * @param v generalized velocity
         */
        void AddDynamicsConstraints(const Eigen::VectorXd& v);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         * @param v generalized velocity
         */
        void AddContactMotionConstraints(const Eigen::VectorXd& v);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         */
        void AddTorqueConstraints(const Eigen::VectorXd& v);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         * Note: Assumes standing on flat ground
         */
        void AddFrictionConeConstraints();

        // Functions to add the necessary costs
        /**
         * Tracks a desired joint acceleration, velocity, and angle. Uses a PD controller.
         * @param q
         * @param v
         */
        void AddLegTrackingCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

        /**
         * Adds the cost for tracking the desired base (torso) position and orientation
         * @param q generalized configuration vector
         * @param v generalized velocity vector
         */
        void AddTorsoCost(const Eigen::VectorXd& q, const Eigen::VectorXd& v);


        /**
         * Tries to track a pre-defined force. Currently there are no PD gains.
         */
        void AddForceTrackingCost();

        void UpdateConstraintsAndCost(const mjData* data);
        /**
         * Recovers the control inputs from the qp_solution. Involves a inverse dynamics calculation
         * @param qp_sol solution to the QP
         * @return the control inputs
         */
        void RecoverControlInputs(const Eigen::VectorXd& qp_sol, const Eigen::VectorXd& v, std::vector<mjtNum>& control);

        Eigen::MatrixXd GetConstraintJacobian(const Eigen::VectorXd& q);
        Eigen::MatrixXd GetConstraintJacobianDerivative(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

        // QP Solver
        OsqpEigen::Solver qp_solver_;

        Eigen::MatrixXd P_;      // quadratic term
        Eigen::VectorXd w_;      // linear term

        Eigen::MatrixXd A_;      // constraint matrix
        Eigen::VectorXd lb_;     // lower bound
        Eigen::VectorXd ub_;     // upper bound

        double torso_tracking_weight_;
        double leg_tracking_weight_;
        double force_tracking_weight_;

        int num_vel_;
        int num_actuators_;
        int num_constraints_;
        int num_decision_vars_;

        int prev_num_contacts_;

        Eigen::VectorXd torque_bounds_;

        Eigen::MatrixXd Js_;
        Eigen::MatrixXd Jsdot_;

        // Floating base position gains
        double kv_pos_;
        double kp_pos_;

        // Floating base angle gains
        double kv_ang_;
        double kp_ang_;

        // Joint gains
        double kv_joint_;
        double kp_joint_;

        double friction_coef_;

        Eigen::VectorXd force_target_;
    };
}   // simulator


#endif //BILEVEL_GAIT_GEN_QP_CONTROL_H
