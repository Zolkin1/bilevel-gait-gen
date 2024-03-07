//
// Created by zolkin on 2/25/24.
//

#include "include/hardware_robot.h"
#include "unitree_lib/quadruped.h"

namespace hardware {
    void LPFData::SetZero() {
        current_sample.setZero();
        prev_sample.setZero();
        prev_output.setZero();
    }

    void LPFData::SetSize(int size) {
        current_sample.resize(size);
        prev_sample.resize(size);
        prev_output.resize(size);
    }

    HardwareRobot::HardwareRobot(const hardware::vector_t& init_config, const hardware::vector_t& init_vel,
                                 const hardware::vector_t& init_mpc_state,
                                 std::unique_ptr<controller::MPCController>& controller, int robot_id,
                                 const JointGains& gains, double optitrack_rate) :
                                 udp(UNITREE_LEGGED_SDK::LOWLEVEL),
                                 safe(UNITREE_LEGGED_SDK::LeggedType::A1),
//                                 8080, UNITREE_LEGGED_SDK::UDP_SERVER_IP_BASIC,
//                                     8007, sizeof(UNITREE_LEGGED_SDK::LowCmd),
//                                     sizeof(UNITREE_LEGGED_SDK::LowState)),
                                 controller_(std::move(controller)),
                                 gravity_offset_(9.8536) { // TODO: make this offset not hard coded
        init_time_ = std::chrono::high_resolution_clock::now();

        log_file_.open("hardware_log.txt");
        auto time_now = std::chrono::system_clock::now();
        std::time_t time1_now = std::chrono::system_clock::to_time_t(time_now);

        log_file_ << "Initializing hardware... " << std::ctime(&time1_now) << std::endl;

        state_log_file_.open("hardware_state_log.txt");
        optitrack_log_file_.open("optitrack_client_data.txt");

        init_config_ = init_config;
        init_vel_ = init_vel;

        robot_id_ = robot_id;

        packet_recieved_ = 0;

//        udp.SwitchLevel(UNITREE_LEGGED_SDK::LOWLEVEL);
        udp.InitCmdData(cmd);

        robot_state_ = Hold;
        standing_time_ = 0.5; // 2 seconds to get to the standing position
        prev_time_ = 0;

        hold_state_.resize(NUM_INPUTS + FLOATING_BASE_OFFSET);

        gains_ = gains;

        // These are used when NOT using the MPC
        motor_kp_ = 30;
        motor_kv_ = 7;

        mpc_time_offset_ = 0;
        mpc_offsets_.setZero();
        in_mpc_ = false;

        // Low pass filter data init
        v_com.SetSize(FLOATING_VEL_OFFSET);
        v_com.SetZero();

        v_joints.SetSize(NUM_INPUTS);
        v_joints.SetZero();

        a_com.SetSize(FLOATING_VEL_OFFSET);
        a_com.SetZero();

        grf_lpf.SetSize(NUM_EE);
        grf_lpf.SetZero();


        // Initialize the controller and MPC
        controller_->InitSolver(init_config, init_mpc_state);

        optitrack_rate_ = optitrack_rate;

        // Start the optitrack monitor
        optitrack_client_ = std::thread(&HardwareRobot::OptiTrackMonitor, this);
    }

    void HardwareRobot::SendUDP() {
        udp.Send();
    }

    void HardwareRobot::RecieveUDP() {
        udp.Recv();
    }

    void HardwareRobot::ControlCallback() {
        udp.GetRecv(state);     // Get the state

        using namespace UNITREE_LEGGED_SDK;

        const double time_s = static_cast<double>(state.tick) / 1e3;

        if (time_s >= prev_time_) {
            auto time_now = std::chrono::system_clock::now();
            std::time_t time1_now = std::chrono::system_clock::to_time_t(time_now);

            log_file_ << "[Communication] udp packet received at quad time " << time_s << " and global time: " <<
                               std::ctime(&time1_now) << std::endl;
            log_file_ << "[Communication] " << packet_recieved_ << " packet received." << std::endl;

            log_file_ << "[State] robot state is: " << StateToString(robot_state_) << std::endl;

            vector_t q, v, a;
            q.resize(NUM_INPUTS + FLOATING_BASE_OFFSET);
            v.resize(NUM_INPUTS + FLOATING_VEL_OFFSET);
            a.resize(NUM_INPUTS + FLOATING_VEL_OFFSET);
            RecoverStateFromMotors(q, v, a, state);

            // TODO: Joint accelerations are always 0

            ComputeCOMStateEstimate(q, v, a, state);

            // Checking for un-updated optitrack packets
            if (v.head<FLOATING_VEL_OFFSET>() == vector_t::Zero(FLOATING_VEL_OFFSET)) {
                v.head<FLOATING_VEL_OFFSET>() = v_com.current_sample;
            }

            Eigen::Vector4d grf;
            grf(0) = state.footForce[FL_];
            grf(1) = state.footForce[FR_];
            grf(2) = state.footForce[RL_];
            grf(3) = state.footForce[RR_];

            // ----- Low Pass Filter ----- //
            a_com.prev_sample = a_com.current_sample;
            a_com.current_sample = a.head<FLOATING_VEL_OFFSET>();

            v_com.prev_sample = v_com.current_sample;
            v_com.current_sample = v.head<FLOATING_VEL_OFFSET>(); //v.head<FLOATING_VEL_OFFSET>();

            v_joints.prev_sample = v_joints.current_sample;
            v_joints.current_sample = v.tail<NUM_INPUTS>();

            grf_lpf.prev_sample = grf_lpf.current_sample;
            grf_lpf.current_sample = grf;

            a.head<FLOATING_VEL_OFFSET>() = LPF(a_com, 15, 1000);
            a_com.prev_output = a.head<FLOATING_VEL_OFFSET>();

            // Note: increased to 100 from 50!
            v.head<FLOATING_VEL_OFFSET>() = LPF(v_com, 10, 240);
            v_com.prev_output = v.head<FLOATING_VEL_OFFSET>(); //v.head<FLOATING_VEL_OFFSET>();

            // Note: increased from 2 to 10!
            v.tail<NUM_INPUTS>() = LPF(v_joints, 200, 1000);
            v_joints.prev_output = v.tail<NUM_INPUTS>();

            grf = LPF(grf_lpf, 50, 1000);
            grf_lpf.prev_output = grf;
            // ---------------------------- //

            a(2) = a(2) + gravity_offset_;


            if (!(packet_recieved_ % state_record_pattern)) {
                state_log_file_ << std::setprecision(9) << time_s << " " << q.transpose() << " " << v.transpose() << " " << a.transpose()
                                << " " << grf.transpose() << std::endl;
            }


//            if (state.robotID == robot_id_) {
            if (robot_state_ == Stand) {
                if (standing_start_ == 0) {
                    standing_start_config_ = q.tail<NUM_INPUTS>();
                    standing_start_ = time_s;
                }

                // Get to the mpc standing state -- won't move the feet. Might need to figure that one out.
                double ratio = std::min(1.0, ((time_s - standing_start_) / standing_time_));
                vector_t q_des =
                        ratio * (init_config_.tail<NUM_INPUTS>() - standing_start_config_) + standing_start_config_;

                AssignConfigToMotors(q_des, cmd);
                AssignVelToMotors(init_vel_.tail<NUM_INPUTS>(), cmd);
                AssignTorqueToMotors(vector_t::Constant(NUM_INPUTS, 0), cmd);

                for (int i = 0; i < NUM_INPUTS; i++) {
                    cmd.motorCmd[i].Kp = motor_kp_;
                    cmd.motorCmd[i].Kd = motor_kv_;
                }

            } else if (robot_state_ == Hold) {

                AssignConfigToMotors(q.tail<NUM_INPUTS>(), cmd);
                AssignVelToMotors(vector_t::Constant(NUM_INPUTS, 0), cmd);
                AssignTorqueToMotors(vector_t::Constant(NUM_INPUTS, 0), cmd);

                for (int i = 0; i < NUM_INPUTS; i++) {
                    cmd.motorCmd[i].Kp = motor_kp_;
                    cmd.motorCmd[i].Kd = motor_kv_;
                }

            } else if (robot_state_ == MPC) {

                if (!in_mpc_) {
                    // Gather initial offsets
                    mpc_offsets_ = q.head<2>(); // Making x and y at zero, but keeping whatever z offset we have
//                    mpc_offsets_(2) -= 0.3; // TODO: Check this
                    mpc_time_offset_ = time_s;
                    in_mpc_ = true;
                    log_file_ << "[Robot] mpc offsets: " << mpc_offsets_.transpose() << std::endl;
                    log_file_ << "[Robot] time offset: " << time_s << std::endl;

                    // TODO: Will want to add Some form of MPC reset to handle the time if I ever want to enter
                    //  MPC twice in the same run. But this only makes sense if I can:
                    //   - Change the initial condition on MPC
                    //   - Handle the new time
                    //   - Adjust the target pose
                }

                // Assuming we are close the "init mpc state" value in the constructor

                controller::Contact contact;
                contact.in_contact_.resize(NUM_EE);
                contact.contact_frames_.resize(NUM_EE);

                // Get contact status -  hardcoded for a1 and pinocchio urdf
                contact.in_contact_.at(0) = grf(0) > 0;
                contact.contact_frames_.at(0) = 14;

                contact.in_contact_.at(1) = grf(1) > 0;
                contact.contact_frames_.at(1) = 24;

                contact.in_contact_.at(2) = grf(2) > 0;
                contact.contact_frames_.at(2) = 34;

                contact.in_contact_.at(3) = grf(3) > 0;
                contact.contact_frames_.at(3) = 44;

                log_file_ << "[Control] q: " << q.transpose() << ", v: " << v.transpose() << ", a: "
                          << a.transpose()
                          << std::endl;

                // Apply the offsets
                q.head<2>() -= mpc_offsets_;

//                v.setZero(); // TODO: Remove
//                a.setZero(); // TODO: Remove

                vector_t control_action = controller_->ComputeControlAction(q, v, a, contact,
                                                                            time_s - mpc_time_offset_);

                if (VerifyControlAction(control_action)) {
                    AssignConfigToMotors(control_action.head<NUM_INPUTS>(), cmd);
                    AssignVelToMotors(control_action.segment<NUM_INPUTS>(NUM_INPUTS), cmd);
//                    AssignTorqueToMotors(control_action.tail<NUM_INPUTS>(), cmd);
                    AssignTorqueToMotors(vector_t::Zero(NUM_INPUTS), cmd);

                    // Assign joint kp and kv
                    AssignMPCGains(contact);
                } else {
                    log_file_ << "[Control] Reverting to a default state." << std::endl;
                    std::cerr << "[Control] [Robot] MPC Control action rejected. Returning to standing state." << std::endl;
                    AssignConfigToMotors(init_config_.tail<NUM_INPUTS>(), cmd);
                    AssignVelToMotors(init_vel_.tail<NUM_INPUTS>(), cmd);
                    AssignTorqueToMotors(vector_t::Constant(NUM_INPUTS, 0), cmd);

                    for (int i = 0; i < NUM_INPUTS; i++) {
                        cmd.motorCmd[i].Kp = motor_kp_;
                        cmd.motorCmd[i].Kd = motor_kv_;
                    }

                    // To be safe, move back to standing
                    ChangeState(Stand);
                }
            }

            if (robot_state_ != Stand) {
                standing_start_ = 0;
            }

            if (robot_state_ != Hold) {
                hold_state_(2) = -1;
            }

            if (robot_state_ != MPC) {
                in_mpc_ = false;
            }

            safe.PowerProtect(cmd, state, 7);

            udp.SetSend(cmd);       // Send the command
            log_file_ << "[Communication] command sent." << std::endl;
//            }

            packet_recieved_++;

            prev_time_ = time_s;
        }
    }

    bool HardwareRobot::VerifyControlAction(const vector_t& control) {
        for (int i = 0; i < control.size(); i++) {
            if (std::abs(control(i)) > 200) {
                log_file_ << "[Control] Control action at index " << i << " is too large." << std::endl;
                return false;
            }
        }

        return true;
    }

    std::string HardwareRobot::StateToString(const hardware::HardwareRobot::State& robot_state) {
        switch (robot_state) {
            case Stand:
                return "Stand";
            case Hold:
                return "Hold";
            case MPC:
                return "MPC";
            case Testing:
                return "Testing";
            default:
                return "Invalid State";
        }
    }

    void HardwareRobot::ChangeState(hardware::HardwareRobot::State new_state) {
        robot_state_ = new_state;
        log_file_ << "[State] Robot changed state to " << StateToString(new_state) << std::endl;
    }

    void HardwareRobot::OptiTrackMonitor() {
        while (true) {
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            optitrack_mut_.lock();
            prev_opti_data_ = opti_data_;
            client_interface_.ReadOptiTrackData(opti_data_);
            optitrack_mut_.unlock();

            optitrack_log_file_ << std::chrono::duration_cast<std::chrono::milliseconds>(
                    t1 - init_time_).count() << " " << opti_data_.transpose() << std::endl;

            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t2 - t1);

            const long time_left = 2000000 - std::ceil(time_span.count() * 1e9);

            timespec *remaining;
            const timespec request = {0, time_left};
            nanosleep(&request, remaining);
        }
    }

    void HardwareRobot::ComputeCOMStateEstimate(hardware::vector_t& q, hardware::vector_t& v, hardware::vector_t& a,
                                                const UNITREE_LEGGED_SDK::LowState& state) {
        optitrack_mut_.lock();
        q.head<7>() = opti_data_;
        v.head<3>() = (opti_data_.head<3>() - prev_opti_data_.head<3>())/optitrack_rate_;
        optitrack_mut_.unlock();

//        q.segment<4>(3) = -q.segment<4>(3);

        // Get the quaternion
//            q(POS_VARS) = state.imu.quaternion[1];
//            q(POS_VARS + 1) = state.imu.quaternion[2];
//            q(POS_VARS + 2) = state.imu.quaternion[3];
//            q(POS_VARS + 3) = state.imu.quaternion[4];

        // TODO: Somehow verify
        v(POS_VARS) = state.imu.gyroscope[0];
        v(POS_VARS + 1) = state.imu.gyroscope[1];
        v(POS_VARS + 2) = state.imu.gyroscope[2];

        a(0) = state.imu.accelerometer[0];
        a(1) = state.imu.accelerometer[1];
        a(2) = -state.imu.accelerometer[2]; // Nominally gives +9.81

        // Note: assuming angular acceleration is very small, so keeping it at 0
    }

    // Hard coded for the urdf to the a1
    void HardwareRobot::AssignConfigToMotors(const hardware::vector_t& q, UNITREE_LEGGED_SDK::LowCmd& cmd) {
        using namespace UNITREE_LEGGED_SDK;

        if (q.size() != NUM_INPUTS) {
            throw std::runtime_error("Invalid config size.");
        }

        // FL Leg
        cmd.motorCmd[FL_0].q = q(0);
        cmd.motorCmd[FL_1].q = q(1);
        cmd.motorCmd[FL_2].q = q(2);

        // FR Leg
        cmd.motorCmd[FR_0].q = q(3);
        cmd.motorCmd[FR_1].q = q(4);
        cmd.motorCmd[FR_2].q = q(5);

        // RL Leg
        cmd.motorCmd[RL_0].q = q(6);
        cmd.motorCmd[RL_1].q = q(7);
        cmd.motorCmd[RL_2].q = q(8);

        // RR Leg
        cmd.motorCmd[RR_0].q = q(9);
        cmd.motorCmd[RR_1].q = q(10);
        cmd.motorCmd[RR_2].q = q(11);
    }

    // Hard coded for the urdf to the a1
    void HardwareRobot::AssignVelToMotors(const hardware::vector_t& v, UNITREE_LEGGED_SDK::LowCmd& cmd) {
        using namespace UNITREE_LEGGED_SDK;

        if (v.size() != NUM_INPUTS) {
            throw std::runtime_error("Invalid velocity size.");
        }

        // FL Leg
        cmd.motorCmd[FL_0].dq = v(0);
        cmd.motorCmd[FL_1].dq = v(1);
        cmd.motorCmd[FL_2].dq = v(2);

        // FR Leg
        cmd.motorCmd[FR_0].dq = v(3);
        cmd.motorCmd[FR_1].dq = v(4);
        cmd.motorCmd[FR_2].dq = v(5);

        // RL Leg
        cmd.motorCmd[RL_0].dq = v(6);
        cmd.motorCmd[RL_1].dq = v(7);
        cmd.motorCmd[RL_2].dq = v(8);

        // RR Leg
        cmd.motorCmd[RR_0].dq = v(9);
        cmd.motorCmd[RR_1].dq = v(10);
        cmd.motorCmd[RR_2].dq = v(11);
    }

    void HardwareRobot::AssignTorqueToMotors(const vector_t& tau, UNITREE_LEGGED_SDK::LowCmd& cmd) {
        using namespace UNITREE_LEGGED_SDK;

        if (tau.size() != NUM_INPUTS) {
            throw std::runtime_error("Invalid velocity size.");
        }

        // FL Leg
        cmd.motorCmd[FL_0].tau = tau(0);
        cmd.motorCmd[FL_1].tau = tau(1);
        cmd.motorCmd[FL_2].tau = tau(2);

        // FR Leg
        cmd.motorCmd[FR_0].tau = tau(3);
        cmd.motorCmd[FR_1].tau = tau(4);
        cmd.motorCmd[FR_2].tau = tau(5);

        // RL Leg
        cmd.motorCmd[RL_0].tau = tau(6);
        cmd.motorCmd[RL_1].tau = tau(7);
        cmd.motorCmd[RL_2].tau = tau(8);

        // RR Leg
        cmd.motorCmd[RR_0].tau = tau(9);
        cmd.motorCmd[RR_1].tau = tau(10);
        cmd.motorCmd[RR_2].tau = tau(11);
    }


    void HardwareRobot::RecoverStateFromMotors(vector_t& q, vector_t& v, vector_t& a, const UNITREE_LEGGED_SDK::LowState& state) {
        using namespace UNITREE_LEGGED_SDK;

        // FL Leg
        q(0 + FLOATING_BASE_OFFSET) = state.motorState[FL_0].q;
        q(1 + FLOATING_BASE_OFFSET) = state.motorState[FL_1].q;
        q(2 + FLOATING_BASE_OFFSET) = state.motorState[FL_2].q;

        v(0 + FLOATING_VEL_OFFSET) = state.motorState[FL_0].dq;
        v(1 + FLOATING_VEL_OFFSET) = state.motorState[FL_1].dq;
        v(2 + FLOATING_VEL_OFFSET) = state.motorState[FL_2].dq;

        a(0 + FLOATING_VEL_OFFSET) = state.motorState[FL_0].ddq;
        a(1 + FLOATING_VEL_OFFSET) = state.motorState[FL_1].ddq;
        a(2 + FLOATING_VEL_OFFSET) = state.motorState[FL_2].ddq;


        // FR Leg
        q(3 + FLOATING_BASE_OFFSET) = state.motorState[FR_0].q;
        q(4 + FLOATING_BASE_OFFSET) = state.motorState[FR_1].q;
        q(5 + FLOATING_BASE_OFFSET) = state.motorState[FR_2].q;

        v(3 + FLOATING_VEL_OFFSET) = state.motorState[FR_0].dq;
        v(4 + FLOATING_VEL_OFFSET) = state.motorState[FR_1].dq;
        v(5 + FLOATING_VEL_OFFSET) = state.motorState[FR_2].dq;

        a(3 + FLOATING_VEL_OFFSET) = state.motorState[FR_0].ddq;
        a(4 + FLOATING_VEL_OFFSET) = state.motorState[FR_1].ddq;
        a(5 + FLOATING_VEL_OFFSET) = state.motorState[FR_2].ddq;


        // RL Leg
        q(6 + FLOATING_BASE_OFFSET) = state.motorState[RL_0].q;
        q(7 + FLOATING_BASE_OFFSET) = state.motorState[RL_1].q;
        q(8 + FLOATING_BASE_OFFSET) = state.motorState[RL_2].q;

        v(6 + FLOATING_VEL_OFFSET) = state.motorState[RL_0].dq;
        v(7 + FLOATING_VEL_OFFSET) = state.motorState[RL_1].dq;
        v(8 + FLOATING_VEL_OFFSET) = state.motorState[RL_2].dq;

        a(6 + FLOATING_VEL_OFFSET) = state.motorState[RL_0].ddq;
        a(7 + FLOATING_VEL_OFFSET) = state.motorState[RL_1].ddq;
        a(8 + FLOATING_VEL_OFFSET) = state.motorState[RL_2].ddq;


        // RR Leg
        q(9 + FLOATING_BASE_OFFSET) = state.motorState[RR_0].q;
        q(10 + FLOATING_BASE_OFFSET) = state.motorState[RR_1].q;
        q(11 + FLOATING_BASE_OFFSET) = state.motorState[RR_2].q;

        v(9 + FLOATING_VEL_OFFSET) = state.motorState[RR_0].dq;
        v(10 + FLOATING_VEL_OFFSET) = state.motorState[RR_1].dq;
        v(11 + FLOATING_VEL_OFFSET) = state.motorState[RR_2].dq;

        a(9 + FLOATING_VEL_OFFSET) = state.motorState[RR_0].ddq;
        a(10 + FLOATING_VEL_OFFSET) = state.motorState[RR_1].ddq;
        a(11 + FLOATING_VEL_OFFSET) = state.motorState[RR_2].ddq;
    }

    vector_t HardwareRobot::LPF(const LPFData& data, double fpass, double fsample) {
        const double x = std::tan(3.14159*fpass/fsample);
        vector_t output = (x*(data.current_sample + data.prev_sample) + (1-x)*data.prev_output)/(x+1);

        return output;
    }

    void HardwareRobot::AssignMPCGains(const controller::Contact& contact) {
        int thigh_idx = 1;
        int calf_idx = 2;

        int ee = 0;
        for (int hip_idx = 0; hip_idx < NUM_INPUTS; hip_idx += 3) {
//            if (!contact.in_contact_.at(ee)) {

                // TODO: Tune
                cmd.motorCmd[hip_idx].Kp = gains_.hip_kp;
                cmd.motorCmd[hip_idx].Kd = gains_.hip_kv;

//            std::cout << calf_idx << ", ";

                cmd.motorCmd[thigh_idx].Kp = gains_.thigh_kp;
                cmd.motorCmd[thigh_idx].Kd = gains_.thigh_kv;

                cmd.motorCmd[calf_idx].Kp = gains_.calf_kp;
                cmd.motorCmd[calf_idx].Kd = gains_.calf_kv;

                cmd.motorCmd[hip_idx].mode = 0x0A;
                cmd.motorCmd[thigh_idx].mode = 0x0A;
                cmd.motorCmd[calf_idx].mode = 0x0A;
//            } else {
//                cmd.motorCmd[hip_idx].Kp = 0;
//                cmd.motorCmd[hip_idx].Kd = 0;
//
//                cmd.motorCmd[thigh_idx].Kp = 0;
//                cmd.motorCmd[thigh_idx].Kd = 0;
//
//                cmd.motorCmd[calf_idx].Kp = 0;
//                cmd.motorCmd[calf_idx].Kd = 0;
//
//                cmd.motorCmd[hip_idx].mode = 0x0A;
//                cmd.motorCmd[thigh_idx].mode = 0x0A;
//                cmd.motorCmd[calf_idx].mode = 0x0A;
//            }

            thigh_idx += 3;
            calf_idx += 3;
        }

//        std::cout << "calf gains: " << gains_.thigh_kp << ", " << gains_.thigh_kv << std::endl;

//        std::cout << std::endl;

//        std::cout << cmd.motorCmd[UNI].mode << std::endl;
//        std::cout << cmd.motorCmd[5].mode << std::endl;

//        cmd.motorCmd[UNITREE_LEGGED_SDK::FL_2].Kp = 30;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::FL_2].Kd = 20;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::FR_2].Kp = 30; //gains_.calf_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::FR_2].Kd = 20; //gains_.calf_kv;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_2].Kp = 30;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_2].Kd = 20;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_2].Kp = 15;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_2].Kd = 5;


//        cmd.motorCmd[7].Kp = 10;
//        cmd.motorCmd[7].Kd = 5;

//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_0].Kp = gains_.hip_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_0].Kd = gains_.hip_kv;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_1].Kp = gains_.thigh_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_1].Kd = gains_.thigh_kv;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_2].Kp = gains_.calf_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RL_2].Kd = gains_.calf_kv;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_0].Kp = gains_.hip_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_0].Kd = gains_.hip_kv;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_1].Kp = gains_.thigh_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_1].Kd = gains_.thigh_kv;
//
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_2].Kp = gains_.calf_kp;
//        cmd.motorCmd[UNITREE_LEGGED_SDK::RR_2].Kd = gains_.calf_kv;


    }

} // hardware