//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_INPUTS_H
#define BILEVEL_GAIT_GEN_INPUTS_H

#include <Eigen/Core>

#include "spline.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;

    class Inputs {
    public:
        Inputs(const std::vector<std::vector<double>>& switching_times, int num_joints, int num_nodes, double node_dt);

        Inputs(const Inputs& input);

        Inputs operator=(const Inputs& input);

        /**
         * Gets the force given by the spline.
         * @param end_effector end effector to use
         * @param time to query the spline at
         * @return Force on that end effector at that time.
         */
        Eigen::Vector3d GetForce(int end_effector, double time) const;

        Eigen::Vector3d GetPosition(int end_effector, double time) const;

        /**
         * Gets the joint vels at that time
         * @param time to query the input at
         * @return vector of joint velocities at that time.
         */
        vector_t GetVels(double time) const;

        /**
         *
         * @return const ref to the force splines for all the end effectors
         */
        const std::vector<std::array<Spline, 3>>& GetForces() const;

        /**
         *
         * @return const ref to the position splines for all the end effectors
         */
        const std::vector<std::array<Spline, 3>>& GetPositions() const;

        /**
         *
         * @return Node time spacing
         */
        double GetNodeDt() const;

        /**
         * Sets the force splines for a given end effector over the whole time period.
         * @param end_effector the end effector to change
         * @param force the force splines for each dimension (x,y,z)
         */
        void SetEndEffectorForce(int end_effector, const std::array<Spline, 3>& force);

        void SetEndEffectorPosition(int end_effector, const std::array<Spline, 3>& position);

        /**
         * Sets the joint velocity input at a given time. Note that a ZOH is used and thus this sets the vels for the
         * while period.
         * @param vels velocities for each of the joints
         * @param time time at which to set the joint velocities
         */
        void SetJointVels(const vector_t& vels, double time);

        // TODO: do better
        void SetJointVelsNoTime(const vector_t& vels, int index);

        /**
         *
         * @return the total number of scalar parameters used to describe the input
         */
        int GetNumInputs() const;

        /**
         *
         * @return the number of spline used to describe the contact forces acting on the robot.
         */
        int GetNumForces() const;

        int GetNumPositions() const;

        /**
         * Converts the inputs into a vector to be used with matrix multiplication.
         * Has the following order:
         * [f_{c1,p1}^{x}(1), ... , f_{c1, p1}^{x}(4), f_{c1,p2}^{x}(1), ... , f_{c1, p2}^{x}(4), ... ]
         * = [f_{c1,p1}^{x}, f_{c1,p2}^{x}, f_{c1,p3}^{x}, ... ]
         * = [f_{c1}^{x}, f_{c1}^{y}, f_{c1}^{z}, f_{c2}^{x}, f_{c2}^{y}, f_{c2}^{z}, ... ]
         * i.e. all the forces for ee1, (x,y,z) then for ee2 (x,y,z), etc...
         * Note: returns polynomial variables for the spline, NOT the value of the spline.
         * @param time
         * @return vectorized input.
         */
        vector_t GetInputVector(double time) const;

        /**
         * Gets all the joint velocities. i.e. joint velocities at all times.
         * @return vector of joint velocities.
         */
        std::vector<vector_t> GetAllVels() const;

        /**
         * Gets the index of the polynomial in the force vector for a given time.
         * @param time to query the spline at
         * @return Index of the polynomial at the given time.
         */
        int GetForcePolyIdx(double time) const;

        void UpdateForce(int end_effector, int coord, const std::vector<std::array<double, Spline::POLY_ORDER>>& vars);

        /**
         * Gets the index of the end point of the variables that define the polynomial at that point in time
         * Also returns the number of variables that define the polynomial.
         * Does NOT return the index into the decision variable vector.
         * @param input
         * @param end_effector
         * @param time
         * @param coord
         * @return
         */
        std::pair<int, int> GetForceSplineIndex(int end_effector, double time, int coord) const;

        int GetTotalForceSplineVars() const;

        void UpdateForcePoly(int end_effector, int coord, int idx, const vector_t& vars);

        int GetTotalForceConstants() const;

        const vector_t& GetVel(int idx) const;

        /**
         * Gets the number of parameters describing the non-constant force values
         * (i.e. not the derivative values)
         * @return
         */
        int GetNumForceValsZ() const;

        // Returns a vector form to be multiplied by the B matrix
        vector_t AsQPVector(double time) const;

        void AddPolys(double time);

        void RemoveUnusedPolys(double time);

        void SetInitTime(double time);

        double GetInitTime() const;

        void SetForceSpline(int ee, int coord, const Spline& spline);

    protected:
    private:
        void UpdateForceSplineVarsCount();

        std::vector<std::array<Spline, 3>> forces_;    // We need a spline for each coordinate of each end effector
        // TODO: remove positions
        std::vector<std::array<Spline, 3>> positions_;
        std::vector<vector_t> joint_vels_;      // These are discrete, ZOH velocities.
        double node_dt_;

//        static int constexpr NUM_POLY = 2;
        static int constexpr POS_VARS = 3;

        int force_spline_vars_;

        double init_time_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_INPUTS_H
