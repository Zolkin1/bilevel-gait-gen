//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_TRAJECTORY_H
#define BILEVEL_GAIT_GEN_TRAJECTORY_H

#include <Eigen/Core>

#include "inputs.h"
#include "centroidal_model.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;

    class Trajectory {
        static int constexpr POS_VARS = 3;

    public:
        Trajectory(int len, int state_size, int num_joints,
                   const std::vector<std::vector<double>>& switching_times, double node_dt,
                   double swing_height, double foot_offset);

        Trajectory& operator=(const Trajectory& traj);

//        Trajectory(const Trajectory& traj);

        std::vector<vector_t> GetStates() const;

        const Inputs& GetInputs() const;

        const std::vector<std::array<Spline, 3>>& GetPositions() const;

//        void UpdatePosition(int end_effector, int coord, const std::vector<std::array<double, Spline::POLY_ORDER>>& vars);

//        void UpdateForce(int end_effector, int coord, const std::vector<std::array<double, Spline::POLY_ORDER>>& vars);

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

        int GetTotalPosConstantsZ() const;

        double GetTotalTime() const;

        // Add an additional spline knot point
        void AddPolys(double final_time);

        void RemoveUnusedPolys(double init_time);

        void SetInitTime(double time);

        void SetEndEffectorSplines(int ee, const Spline& force_spline,
                                   const Spline& pos_spline);

        vector_t ConvertToQPVector() const;

        vector_t PositionAsQPVector() const;

        vector_t GetState(int node) const;

        Eigen::Vector3d GetPosition(int ee, double time) const;

        std::vector<bool> GetContacts(double time) const;

        bool PosIsConstant(int ee, int coord, double time) const;

        int GetTotalPosNonConstantZ() const;

        bool IsSplineMutable(int ee, int coord) const;

        void UpdateFullVelocity(int node, const vector_t& vel);

        vector_t GetFullVelocity(int node);

        vector_t GetAcc(int node, double dt);

        std::vector<std::vector<Eigen::Vector3d>> CreateVizData(const CentroidalModel& model);

        int GetNumContactNodes(int ee) const;

        std::vector<std::vector<double>> GetContactTimes() const;

    protected:
    private:
        void UpdatePosSplineVarsCount();

        void SetSwingPosZ();

        void UpdateContactTimes();


        std::vector<vector_t> states_;
        Inputs inputs_;
        std::vector<std::array<Spline, 3>> end_effector_pos_;
        std::vector<std::array<bool, 3>> mut_flags_;

        int pos_spline_vars_;

        double swing_height_;
        double foot_offset_;
        std::vector<vector_t> full_velocities_; // TODO: Do this cleaner

        std::vector<std::vector<Eigen::Vector3d>> fk_traj_;

        std::vector<std::vector<double>> contact_times_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_TRAJECTORY_H
