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

        MPCSingleRigidBody(const MPCSingleRigidBody& other);

        Trajectory Solve(const vector_t& state, double init_time,
                         const std::vector<vector_3t>& ee_start_locations) override;

        vector_t ConvertTrajToQPVec(const Trajectory& traj) const override;

        std::vector<std::vector<Eigen::Vector3d>> CreateVizData() override;

        vector_t GetFullTargetState(double time, const vector_t& state_guess);

        std::vector<Eigen::Vector2d> GetEEBoxCenter();

        void UpdateWholeBodyTrajectory(Trajectory& traj);

        // TODO: Merge and/or clean up these functions
        bool ComputeParamPartialsOSQP(const Trajectory& traj, QPPartials& partials, int ee, int idx);

        bool ComputeParamPartialsClarabel(const Trajectory& traj, QPPartials& partials, int ee, int idx);

        double GetCost() const;

        double GetModifiedCost(int num_nodes) const;

        SingleRigidBodyModel GetModelCopy() const;

    protected:
        void AddDynamicsConstraints(const vector_t& state) override;

        void AddEELocationConstraints();

        void AddTDPositionConstraints();

        void AddTDPositionConstraintPartial(utils::SparseMatrixBuilder& builder, vector_t& b, int contact_idx, int eq_idx, int ee);

        void AddRaibertHeuristic();

        void AddRaibertPartials(utils::SparseMatrixBuilder& builder, int contact_idx, int eq_idx, int ee);

        void AddEEStartConstraints(const std::vector<vector_3t>& ee_start_locations);

        int GetForceSplineStartIdx() const override;
        int GetPosSplineStartIdx() const override;

        void SetInitQPSizes() override;

        void InitalizeQPData();

        void AddEEPosCost();

        void IncreaseEEBox();

        void DecreaseEEBox();

        int num_run_;

        static constexpr int EE_NODE_START = 4; // TODO: Tune

        vector_2t ee_bounds_;

    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_MPC_SINGLE_RIGID_BODY_H
