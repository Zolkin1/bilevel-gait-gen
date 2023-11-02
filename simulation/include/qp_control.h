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
 *   - dynamics \in R^{nv}
 *   - friction cone (pyramid)  \in R^{4 * ncontacts}
 *   - stationary contacts \in R^{3 * ncontacs}
 *   - torque limits \in R^{njoints}
 *
 *   total constraints: nv + 4*ncontacts + 3*ncontacts + njoints
 */

// TODO: Is this a dense problem?

namespace simulator {
    class QPControl : public Controller {
    public:
        QPControl(double control_rate, std::string robot_urdf, const std::string& foot_type, int nv,
                  const Eigen::VectorXd& torque_bounds);

        void InitSolver(const mjData* data) override;

        std::vector<mjtNum> ComputeControlAction(const mjModel* model, const mjData* data) override;

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
        void ComputeDynamicsTerms(const Eigen::VectorXd& q, const Eigen::VectorXd& v);
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
         */
        void AddFrictionConeConstraints();

        // Functions to add the necessary costs
        void AddLegTrackingCost();
        void AddTorsoCost();
        void AddForceTrackingCost();

        void UpdateConstraintsAndCost(const mjData* data);
        /**
         * Recovers the control inputs from the qp_solution. Involves a inverse dynamics calculation
         * @param qp_sol solution to the QP
         * @return the control inputs
         */
        std::vector<double> RecoverControlInputs(const Eigen::VectorXd& qp_sol);

        Eigen::MatrixXd GetConstraintJacobian(const Eigen::VectorXd& q);
        Eigen::MatrixXd GetConstraintJacobianDerivative(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

        // Hold the underlying QP problem and any other needed data structures
        OsqpEigen::Solver qp_solver_;

        Eigen::MatrixXd P_;      // quadratic term
        Eigen::VectorXd w_;      // linear term

        Eigen::MatrixXd A_;      // constraint matrix
        Eigen::VectorXd lb_;     // lower bound
        Eigen::VectorXd ub_;     // upper bound

        // Other member variables
        int num_vel_;
        int num_actuators_;
        int num_constraints_;
        int num_decision_vars_;

        Eigen::VectorXd torque_bounds_;

        Eigen::MatrixXd Js_;
        Eigen::MatrixXd Jsdot_;
    };
}   // simulator


#endif //BILEVEL_GAIT_GEN_QP_CONTROL_H
