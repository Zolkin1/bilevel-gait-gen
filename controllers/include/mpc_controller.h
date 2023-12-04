//
// Created by zolkin on 11/29/23.
//

#ifndef BILEVEL_GAIT_GEN_MPC_CONTROLLER_H
#define BILEVEL_GAIT_GEN_MPC_CONTROLLER_H

#include "mpc.h"
#include "qp_control.h"

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
                  const std::vector<double>& joint_gains,
                  double leg_weight,
                  double torso_weight,
                  double force_weight,
                  mpc::MPCInfo info,
                  const std::vector<vector_t>& warm_start_states,
                  const vector_t& state_des);

        vector_t ComputeControlAction(const vector_t& q,
                                      const vector_t& v,
                                      const vector_t& a,
                                      const Contact& contact,
                                      double time) override;
                                      //double time) override;

    protected:
    private:
        vector_t ReconstructState(const vector_t& q, const vector_t& v, const vector_t& a) const;

        double prev_time_;

        QPControl qp_controller_;
        mpc::MPC mpc_;
    };
} // controller


#endif //BILEVEL_GAIT_GEN_MPC_CONTROLLER_H
