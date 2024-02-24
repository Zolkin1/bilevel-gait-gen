//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_TRAJECTORY_H
#define BILEVEL_GAIT_GEN_TRAJECTORY_H

#include <Eigen/Core>

#include "spline/spline.h"
#include "controller.h"
#include "spline/end_effector_splines.h"
//#include "model.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;
    using vector_3t = Eigen::Vector3d;

    class Trajectory {
        static int constexpr POS_VARS = 3;

    public:
        enum SplineTypes {
            Position,
            Force
        };

        Trajectory(int len, int state_size, bool using_joints,
                   const std::vector<std::vector<double>>& switching_times, double node_dt,
                   double swing_height, double foot_offset);

//        Trajectory& operator=(const Trajectory& traj);

        std::vector<vector_t> GetStates() const;

        /**
         * Resets all states_ and inputs to 0
         */
        void Reset();

        void SetState(int idx, const vector_t& state);
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

        std::pair<int, int> GetForceSplineIndex(int ee, double time,  int coord) const;

        void SetPositionsForAllTime(int ee, const vector_3t& ee_pos);

        void PrintTrajectoryToFile(const std::string& file_name) const;

        // Add an additional spline knot point
        void AddPolys(double final_time);

        void RemoveUnusedPolys(double init_time);

        void SetInitTime(double time);

        vector_t GetState(int node) const;

        std::vector<bool> GetContacts(double time) const;

        void UpdateFullVelocity(int node, const vector_t& vel);
        void UpdateFullConfig(int node, const vector_t& q);

        vector_t GetFullVelocity(int node) const;
        vector_t GetFullConfig(int node) const;

        int GetNumContactNodes(int ee) const;

        std::vector<time_v> GetContactTimes() const;

        int GetTotalForceSplineVars() const; // TODO: impelement

        bool IsForceMutable(int ee, double time) const;
//        bool IsForceMutable(int ee, int coord, int idx) const;

        vector_t GetSplineLin(const SplineTypes& spline_type, int end_effector, int coord, double time) const;

        int GetTotalPolyVars(const SplineTypes& spline_type, int end_effector, int coord);

        Eigen::Vector3d GetForce(int end_effector, double time) const;

        Eigen::Vector3d GetEndEffectorLocation(int end_effector, double time) const;

        double GetTime(int node) const;

        // Should give pos_spline_vars + force_spline_vars + all vells + all states
        int GetTotalVariables() const;    // TODO: Consider if we have joint information

        vector_t SplinesAsVec() const;

        vector_t SplinePartialsAsVec(int contact_idx) const;

        int GetNode(double time) const;

        controller::Contact GetDesiredContacts(double time) const;

        vector_3t GetForcePartialWrtContactTime(int end_effector, double time, int contact_idx) const;

        vector_3t GetPositionPartialWrtContactTime(int end_effector, double time, int contact_idx) const;

        vector_t GetForceCoefPartialsWrtContactTime(int end_effector, int coord, double time, int contact_idx) const;

        vector_t GetForceCoefPartialsWrtContactTime(int end_effector, int coord, double time, int contact_idx, double dtwdth) const;

        vector_t GetPositionCoefPartialsWrtContactTime(int end_effector, int coord, double time, int contact_idx) const;

        void UpdateContactTimes(std::vector<time_v>& contact_times);

        double GetNextContactTime(int ee, double time) const;

        void SetEEInContact(int ee, double time);

        double GetCurrentSwingTime(int ee) const;

    protected:
    private:
        void UpdateSplineVarsCount();

        void SetSwingPosZ();

        void UpdateContactTimes();

        void UpdateForceSplineVarsCount();

        bool using_joints_;

        std::vector<vector_t> states_;
//        Inputs inputs_;
//        std::vector<std::array<Spline, 3>> end_effector_pos_;
//        std::vector<std::array<Spline, 3>> forces_;
//        std::vector<std::array<bool, 3>> mut_flags_;

        std::vector<EndEffectorSplines> ee_splines_;

        int pos_spline_vars_;
        int force_spline_vars_;

        double swing_height_;
        double foot_offset_;
//        std::vector<vector_t> full_velocities_; // TODO: Do this cleaner

        std::vector<std::vector<Eigen::Vector3d>> fk_traj_;

//        std::vector<std::vector<double>> contact_times_;

        std::array<vector_t, 51> full_config_;
        std::array<vector_t, 51> full_velocity_;

        double init_time_;

        double node_dt_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_TRAJECTORY_H
