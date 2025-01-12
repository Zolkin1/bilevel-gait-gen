//
// Created by zolkin on 11/29/23.
//

#ifndef BILEVEL_GAIT_GEN_MPC_CONTROLLER_H
#define BILEVEL_GAIT_GEN_MPC_CONTROLLER_H

#include <thread>
#include <queue>
#include <mutex>

#include "mpc_single_rigid_body.h"
#include "qp_control.h"
#include "gait_optimizer.h"

namespace controller {
    /**
     * Class that uses MPC as a feedback controller
     */
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;

    class MPCController : public Controller {
    public:
        MPCController(double control_rate, std::string robot_urdf, const std::string& foot_type, int nv,
                  const Eigen::VectorXd& torque_bounds, double friction_coef,
                  const std::vector<double>& base_pos_gains,
                  const std::vector<double>& base_ang_gains,
                  const vector_t& kp_joint_gains,
                  const vector_t& kd_joint_gains,
                  double leg_weight,
                  double torso_weight,
                  double force_weight,
                  mpc::MPCInfo info,
                  const std::vector<vector_t>& warm_start_states,
                  const vector_t& state_des,
                  int num_polys,
                  const matrix_t& Q,
                  int gait_opt_freq,
                  const std::string& log_file);

        vector_t ComputeControlAction(const vector_t& q,
                                      const vector_t& v,
                                      const vector_t& a,
                                      const Contact& contact,
                                      double time) override;
                                      //double time) override;

        void InitSolver(const vector_t& full_body_state, const vector_t& mpc_state) override;

        std::vector<std::vector<Eigen::Vector3d>> GetTrajViz() override;

        void UpdateTrajViz();

        std::vector<Eigen::Vector2d> GetEEBoxCenter() override;

        const mpc::MPCSingleRigidBody& GetMPC() const;

    protected:
    private:
        vector_t ReconstructState(const vector_t& q, const vector_t& v, const vector_t& a) const;

        void MPCUpdate();

        void FullBodyTrajUpdate(mpc::Trajectory& traj);

        void GetTargetsFromTraj(const mpc::Trajectory& traj, double time);

        bool GaitOpt(double cost_red, double time, const std::vector<Eigen::Vector3d>& ee_locations);

        vector_t ComputeGroundForceTorques(int end_effector);

        vector_t ComputeSwingTorques(int end_effector);

        vector_t NormalizeQuat(const vector_t& state);

        void PrintContactTimes() const;

        double prev_time_;

        QPControl qp_controller_;
        mpc::MPCSingleRigidBody mpc_;
        mpc::GaitOptimizer gait_opt_;

        int gait_opt_freq_;

        bool computed_;

        std::thread mpc_computations_;

        vector_t state_;
        double time_;
        vector_t q_des_;
        vector_t v_des_;
        vector_t a_des_;
        std::vector<mpc::vector_3t> force_des_;
        mpc::Trajectory traj_;

        std::mutex state_time_mut_;
        std::mutex mpc_res_mut_;
        std::mutex traj_viz_mut_;
        std::mutex sync_mut_;
        std::mutex model_mut_;

        std::vector<std::vector<Eigen::Vector3d>> fk_traj_;

        std::vector<mpc::vector_3t> ee_locations_;

        mpc::GaitOptimizer gait_optimizer_;

        mpc::SingleRigidBodyModel model_;

        mpc::MPCInfo info_;

        int num_polys_;

        vector_t full_body_state_;

        vector_t kp_joints_;
        vector_t kv_joints_;

        int run_num;

        Contact contact_;

        std::ofstream log_file_;
    };
} // controller


#endif //BILEVEL_GAIT_GEN_MPC_CONTROLLER_H
