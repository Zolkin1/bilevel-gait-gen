//
// Created by zolkin on 1/19/24.
//

#ifndef BILEVEL_GAIT_GEN_MPC_SINGLE_RIGID_BODY_H
#define BILEVEL_GAIT_GEN_MPC_SINGLE_RIGID_BODY_H

#include "mpc.h"

namespace mpc {
    class MPCSingleRigidBody : public MPC {
    public:
        MPCSingleRigidBody(const MPCInfo& info, const std::string& robot_urdf);

        Trajectory Solve(const vector_t& state, double init_time) override;

        vector_t ConvertTrajToQPVec(const Trajectory& traj) const override;

        std::vector<std::vector<Eigen::Vector3d>> CreateVizData() override;

        vector_t GetFullTargetState(double time, const vector_t& state_guess);

    protected:
        void AddDynamicsConstraints(const vector_t& state) override;

        void AddEELocationConstraints();

        int GetForceSplineStartIdx() const override;
        int GetPosSplineStartIdx() const override;

        void SetInitQPSizes() override;

        void InitalizeQPData();

        int num_run_;

    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_MPC_SINGLE_RIGID_BODY_H