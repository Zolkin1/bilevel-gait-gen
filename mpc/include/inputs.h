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

        Eigen::Vector3d GetForce(int end_effector, double time) const;

        Eigen::Vector3d GetPosition(int end_effector, double time) const;

        vector_t GetVels(double time) const;

        const std::vector<std::array<Spline, 3>>& GetForces() const;

        const std::vector<std::array<Spline, 3>>& GetPositions() const;

        double GetNodeDt() const;

        void SetEndEffectorForce(int end_effector, const std::array<Spline, 3>& force);

        void SetEndEffectorPosition(int end_effector, const std::array<Spline, 3>& position);

        void SetJointVels(const vector_t& vels, double time);

        int GetNumInputs() const;
        int GetNumForces() const;
        int GetNumPositions() const;

        vector_t GetInputVector(double time) const;

        std::vector<vector_t> GetAllVels() const;

    protected:
    private:
        std::vector<std::array<Spline, 3>> forces_;    // We need a spline for each coordinate of each end effector
        std::vector<std::array<Spline, 3>> positions_;
        std::vector<vector_t> joint_vels_;      // These are discrete, ZOH velocities.
        double node_dt_;

        static int constexpr NUM_POLY = 2;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_INPUTS_H
