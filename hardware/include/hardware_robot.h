//
// Created by zolkin on 2/25/24.
//

#ifndef BILEVEL_GAIT_GEN_HARDWARE_ROBOT_H
#define BILEVEL_GAIT_GEN_HARDWARE_ROBOT_H

#include <Eigen/Core>

#include "mpc_controller.h"

#include "../unitree_lib/loop.h"
#include "../unitree_lib/udp.h"
#include "../unitree_lib/comm.h"

namespace hardware {
    using vector_t = Eigen::VectorXd;

    class HardwareRobot {
    public:
        HardwareRobot(const vector_t& init_config, const vector_t& init_vel,
                      std::unique_ptr<controller::MPCController>& controller,
                      int robot_id);

        void StartControlLoop();
    protected:
    private:
        void ControlCallback();

        vector_t ComputeControlAction(double time, const vector_t& q, const vector_t& v);

        bool VerifyControlAction(const vector_t& control);

        void SendUDP();
        void RecieveUDP();

        vector_t ConvertHardwareJointsToPinocchio(const vector_t& q);
        vector_t ConvertHardwareConfigToPinocchio(const vector_t& q);
        vector_t ConvertHardwareVelToPinocchio(const vector_t& v);

        vector_t ConvertPinocchioConfigToHardware(const vector_t& q);
        vector_t ConvertPinocchioJointsToHardware(const vector_t& q);
        vector_t ConvertPinocchioVelToHardware(const vector_t& v);

        std::unique_ptr<controller::MPCController> controller_;

        ofstream log_file_;

        vector_t init_config_;
        vector_t init_vel_;

        static constexpr int POS_VARS = 3;
        static constexpr int FLOATING_BASE_OFFSET = 7;
        static constexpr int FLOATING_VEL_OFFSET = 6;

        static constexpr int NUM_EE = 4;

        static constexpr int NUM_CONFIG = 19;
        static constexpr int NUM_VEL_LIKE = 18;
        static constexpr int NUM_INPUTS = 12;

        int robot_id_;

        UNITREE_LEGGED_SDK::UDP udp;
        UNITREE_LEGGED_SDK::LowCmd cmd = {0};
        UNITREE_LEGGED_SDK::LowState state = {0};
        UNITREE_LEGGED_SDK::LoopFunc control_loop;
        const double dt = 0.002;

        const double motor_kp = 20;
        const double motor_kv = 7;

        int packet_recieved_;
    };
} // hardware


#endif //BILEVEL_GAIT_GEN_HARDWARE_ROBOT_H
