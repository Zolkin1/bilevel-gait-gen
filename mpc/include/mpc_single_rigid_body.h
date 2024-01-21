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

        virtual vector_t ConvertTrajToQPVec(const Trajectory& traj) const override;

        virtual std::vector<std::vector<Eigen::Vector3d>> CreateVizData();
    protected:
        void AddDynamicsConstraints(const vector_t& state) override;

        int GetForceSplineStartIdx() const override;
        int GetPosSplineStartIdx() const override;

        void SetInitQPSizes() override;

        void InitalizeQPData();

    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_MPC_SINGLE_RIGID_BODY_H
