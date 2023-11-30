//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_TRAJECTORY_H
#define BILEVEL_GAIT_GEN_TRAJECTORY_H

#include <Eigen/Core>

#include "inputs.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;

    class Trajectory {
        static int constexpr POS_VARS = 3;

    public:
        Trajectory(int len, int state_size, int num_joints,
                   const std::vector<std::vector<double>>& switching_times, double node_dt);

//        Trajectory(const Trajectory& traj);

        std::vector<vector_t> GetStates() const;

        const Inputs& GetInputs() const;

        const std::vector<std::array<Spline, 3>>& GetPositions() const;

        void UpdatePosition(int end_effector, int coord, const std::vector<std::array<double, Spline::POLY_ORDER>>& vars);

        void UpdateForce(int end_effector, int coord, const std::vector<std::array<double, Spline::POLY_ORDER>>& vars);

        /**
         * Resets all states_ and inputs to 0
         */
        void Reset();

        void SetState(int idx, const vector_t& state);
        void SetInput(const Inputs& input);
        void SetInputVels(int idx, const vector_t& joint_vels);

        int GetTotalPosSplineVars() const;

        void UpdateForceSpline(int end_effector, int coord, const vector_t& vars);

        void UpdatePositionSpline(int end_effector, int coord, const vector_t& vars);

        /**
         * Gets the index into the position splines. i.e. if you stack all the spline variables for the position
         * vars into a vector then query it.
         * @param end_effector
         * @param time
         * @param coord
         * @return
         */
        std::pair<int, int> GetPositionSplineIndex(int end_effector, double time, int coord) const;

        void SetPositionsForAllTime(int ee, const std::array<double, POS_VARS>& ee_pos);

        void PrintTrajectoryToFile(const std::string& file_name) const;

        int GetTotalPosConstants() const;

        double GetTotalTime() const;

        // Add an additional spline knot point
        void AddPolys(double final_time);

        void RemoveUnusedPolys(double init_time);

        void SetInitTime(double time);

        void SetEndEffectorSplines(int ee, const Spline& spline);

    protected:
    private:
        void UpdatePosSplineVarsCount();

        std::vector<vector_t> states_;
        Inputs inputs_;
        std::vector<std::array<Spline, 3>> end_effector_pos_;

        int pos_spline_vars_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_TRAJECTORY_H
