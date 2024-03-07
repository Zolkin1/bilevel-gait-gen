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
#include "../unitree_lib/safety.h"

#include "optitrack-stream-client/client_interface.h"

namespace hardware {
    using vector_t = Eigen::VectorXd;

    struct LPFData {
        vector_t current_sample;
        vector_t prev_sample;
        vector_t prev_output;

        void SetSize(int size);
        void SetZero();
    };

    struct JointGains {
        double hip_kp;
        double hip_kv;

        double thigh_kp;
        double thigh_kv;

        double calf_kp;
        double calf_kv;
    };

    class HardwareRobot {
    public:
        enum State {
            Hold,
            Stand,
            MPC,
            Testing
        };

        HardwareRobot(const vector_t& init_config, const vector_t& init_vel,
                      const vector_t& init_mpc_state,
                      std::unique_ptr<controller::MPCController>& controller,
                      int robot_id, const JointGains& gains, double optitrack_rate);

        void ControlCallback();

        void SendUDP();
        void RecieveUDP();

        void ChangeState(State new_state);
        void SetNominalStandingState(const vector_t& q);

    protected:
    private:

        void OptiTrackMonitor();

        void ComputeCOMStateEstimate(vector_t& q, vector_t& v, vector_t& a, const UNITREE_LEGGED_SDK::LowState& state);

        static std::string StateToString(const State& robot_state);

        vector_t ComputeControlAction(double time, const vector_t& q, const vector_t& v);

        bool VerifyControlAction(const vector_t& control);

        vector_t LPF(const LPFData& data, double fpass, double fsample);

        static void AssignConfigToMotors(const vector_t& q, UNITREE_LEGGED_SDK::LowCmd& cmd);
        static void AssignVelToMotors(const vector_t& v, UNITREE_LEGGED_SDK::LowCmd& cmd);
        static void AssignTorqueToMotors(const vector_t& tau, UNITREE_LEGGED_SDK::LowCmd& cmd);
        static void RecoverStateFromMotors(vector_t& q, vector_t& v, vector_t& a,
                                           const UNITREE_LEGGED_SDK::LowState& state);

        void AssignMPCGains(const controller::Contact& contact);

        vector_t ConvertHardwareJointsToPinocchio(const vector_t& q);
        vector_t ConvertHardwareConfigToPinocchio(const vector_t& q);
        vector_t ConvertHardwareVelToPinocchio(const vector_t& v);

        vector_t ConvertPinocchioConfigToHardware(const vector_t& q);
        vector_t ConvertPinocchioJointsToHardware(const vector_t& q);
        vector_t ConvertPinocchioVelToHardware(const vector_t& v);

        std::unique_ptr<controller::MPCController> controller_;

        ofstream log_file_;
        ofstream state_log_file_;
        ofstream optitrack_log_file_;

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
        UNITREE_LEGGED_SDK::Safety safe;

        std::thread optitrack_client_;
        std::mutex optitrack_mut_;
        Eigen::Vector<double, 7> opti_data_, prev_opti_data_;
        stream_client::ClientInterface client_interface_;

        State robot_state_;

        const double dt = 0.002;

        double motor_kp_; //20;
        double motor_kv_; //7;

        const int state_record_pattern = 1;

        double standing_time_;
        double standing_start_;
        vector_t standing_start_config_;
        vector_t hold_state_;

        int packet_recieved_;
        vector_t standing_nominal_;
        double prev_time_;
        double optitrack_rate_;

        double mpc_time_offset_;
        Eigen::Vector2d mpc_offsets_;
        bool in_mpc_;

        LPFData v_com;
        LPFData a_com;
        LPFData v_joints;
        LPFData grf_lpf;

        JointGains gains_;

        std::chrono::high_resolution_clock::time_point init_time_;

        static constexpr int ZEROING_SAMPLES = 20;

        const double gravity_offset_;
    };
} // hardware


#endif //BILEVEL_GAIT_GEN_HARDWARE_ROBOT_H
