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

namespace controller {
    class QPControl : public Controller {
        using vector_t = Eigen::VectorXd;
    public:
        QPControl(double control_rate, std::string robot_urdf, const std::string& foot_type, int nv,
                  const Eigen::VectorXd& torque_bounds, double friction_coef,
                  const std::vector<double>& base_pos_gains,
                  const std::vector<double>& base_ang_gains,
                  const vector_t& kp_joint_gains,
                  const vector_t& kd_joint_gains,
                  double leg_weight,
                  double torso_weight,
                  double force_weight,
                  int num_contacts);

        Eigen::VectorXd ComputeControlAction(const Eigen::VectorXd& q,
                                             const Eigen::VectorXd& v,
                                             const Eigen::VectorXd& a,
                                             const Contact& contact,
                                             double time) override;

        void SetBasePosGains(double kv, double kp);

        void SetBaseAngleGains(double kv, double kp);

        /**
         * Note: assumes the same kp and kv for each joint.
         * @param kv
         * @param kp
         */
        void SetJointGains(const vector_t& kv, const vector_t& kp);

        void UpdateDesiredContacts(const Contact& contact) override;

        void UpdateForceTargets(const Eigen::VectorXd& force);

    protected:
    private:
        /**
         * Computes all the dynamics used by the constraints and cost.
         * All terms are then stored in the pinocchio data struct if possible.
         * Constraint jacobian and its time derivative are stored in the associated local variable.
         * @param q generalized position vector
         * @param v generalized velocity vector
         */
        void ComputeDynamicsTerms(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Contact& contact);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         * When using generalized accelerations and contact forces, I only need to enforce the dynamics in the
         * dimensions where there is no torque input. When every joint is actuated then this just becomes the
         * floating base coordinates.
         * @param v generalized velocity
         */
        void AddDynamicsConstraints(const Eigen::VectorXd& v, const Contact& contact);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         * @param v generalized velocity
         */
        void AddContactMotionConstraints(const Eigen::VectorXd& v, int num_contacts);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         */
        void AddTorqueConstraints(const Eigen::VectorXd& v, const Contact& contact);
        /**
         * Note: Assumes ComputeDynamicsTerms has already been called.
         * Note: Assumes standing on flat ground
         */
        void AddFrictionConeConstraints(const Contact& contact);

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
        void AddForceTrackingCost(const Contact& contact);

        void UpdateConstraintsAndCost(const Eigen::VectorXd& q,
                                      const Eigen::VectorXd& v,
                                      const Eigen::VectorXd& a,
                                      const Contact& contact);
        /**
         * Recovers the control inputs from the qp_solution. Involves a inverse dynamics calculation
         * @param qp_sol solution to the QP
         * @return the control inputs
         */
        void RecoverControlInputs(const Eigen::VectorXd& qp_sol, const Eigen::VectorXd& v,
                                  Eigen::VectorXd& control, const Contact& contact);

        Eigen::MatrixXd GetConstraintJacobian(const Eigen::VectorXd& q, const Contact& contact);
        Eigen::MatrixXd GetConstraintJacobianDerivative(const Eigen::VectorXd& q, const Eigen::VectorXd& v,
                                                             const Contact& contact);

        int GetNumBothContacts(const Contact& contact1, const Contact& contact2) const;

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
        int num_constraints_{};
        int num_decision_vars_{};

        int prev_num_contacts_;

        Eigen::VectorXd torque_bounds_;

        Eigen::MatrixXd Js_;
        Eigen::MatrixXd Jsdot_;

        // Floating base position gains
        double kv_pos_{};
        double kp_pos_{};

        // Floating base angle gains
        double kv_ang_{};
        double kp_ang_{};

        // Joint gains
        Eigen::VectorXd kv_joint_;
        Eigen::VectorXd kp_joint_;

        double friction_coef_;

        Eigen::VectorXd force_target_;

        Contact des_contact_;
    };
}   // controller


#endif //BILEVEL_GAIT_GEN_QP_CONTROL_H
