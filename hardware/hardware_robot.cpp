//
// Created by zolkin on 2/25/24.
//

#include "include/hardware_robot.h"

namespace hardware {
    HardwareRobot::HardwareRobot(const hardware::vector_t& init_config, const hardware::vector_t& init_vel,
                                 std::unique_ptr<controller::MPCController>& controller, int robot_id) :
                                 udp(UNITREE_LEGGED_SDK::LOWLEVEL),
                                 controller_(std::move(controller)),
                                 control_loop("control_loop", dt, boost::bind(&HardwareRobot::ControlCallback, this)) {
        log_file_.open("hardware_log.txt");
        log_file_ << "Initializing hardware." << std::endl;

        init_config_ = init_config;
        init_vel_ = init_vel;

        robot_id_ = robot_id;

        packet_recieved_ = 0;


    }

    void HardwareRobot::StartControlLoop() {
        log_file_ << "Creating call back loop." << std::endl;
        control_loop.start();
    }

    void HardwareRobot::SendUDP() {
        udp.Send();
    }

    void HardwareRobot::RecieveUDP() {
        udp.Recv();
    }

    void HardwareRobot::ControlCallback() {
        udp.GetRecv(state);     // Get the state

        if (state.robotID == robot_id_) {

            const double time_s = static_cast<double>(state.tick) / 1e6;

            log_file_ << "[Communication] udp packet received at time " << time_s << std::endl;
            log_file_ << "[Communication] " << packet_recieved_ << " packet recieved." << std::endl;
            packet_recieved_++;

            vector_t q(NUM_CONFIG);
            vector_t v(NUM_VEL_LIKE);
            vector_t a(NUM_VEL_LIKE);
            for (int i = 0; i < NUM_INPUTS; i++) {
                // TODO: Figure out which motors to use/the mapping
                q(FLOATING_BASE_OFFSET + i) = state.motorState[i].q;
                v(FLOATING_VEL_OFFSET + i) = state.motorState[i].dq;
                a(FLOATING_VEL_OFFSET + i) = state.motorState[i].ddq;
            }

            // Get the quaternion
            q(POS_VARS) = state.imu.quaternion[1];
            q(POS_VARS + 1) = state.imu.quaternion[2];
            q(POS_VARS + 2) = state.imu.quaternion[3];
            q(POS_VARS + 3) = state.imu.quaternion[4];

            // TODO: Filter
            v(POS_VARS) = state.imu.gyroscope[0];
            v(POS_VARS + 1) = state.imu.gyroscope[1];
            v(POS_VARS + 2) = state.imu.gyroscope[2];

            // TODO: Filter
            a(POS_VARS) = state.imu.accelerometer[0];
            a(POS_VARS + 1) = state.imu.accelerometer[1];
            a(POS_VARS + 2) = state.imu.accelerometer[2];

            // TODO: Check/set the motor mode

            // TODO: Check which sensor measurement to use
            controller::Contact contact;
            contact.in_contact_.resize(NUM_EE);
            contact.contact_frames_.resize(NUM_EE);
            for (int i = 0; i < NUM_EE; i++) {
                contact.contact_frames_.at(i) = (10 * i) + 4; // TODO: Get correct foot to frame mapping
                if (state.footForce[i] > 0) { // TODO: Tune this value
                    contact.in_contact_.at(i) = true;
                } else {
                    contact.in_contact_.at(i) = false;
                }
            }

            log_file_ << "[Control] q: " << q.transpose() << ", v: " << v.transpose() << ", a: " << a.transpose()
                      << std::endl;

            vector_t control_action = controller_->ComputeControlAction(q, v, a, contact, time_s);

            if (VerifyControlAction(control_action)) {
                for (int i = 0; i < NUM_INPUTS; i++) {
                    // TODO: Figure out proper mapping and figure out mode
                    cmd.motorCmd[i].q = control_action(i);
                    cmd.motorCmd[i].dq = control_action(i + NUM_INPUTS);
                    cmd.motorCmd[i].tau = control_action(i + 2 * NUM_INPUTS);
                    cmd.motorCmd[i].Kp = motor_kp; // TODO: tune
                    cmd.motorCmd[i].Kd = motor_kv; // TODO: tune
                }
            } else {
                log_file_ << "[Control] Reverting to a default state." << std::endl;
                for (int i = 0; i < NUM_INPUTS; i++) {
                    // TODO: Figure out proper mapping and figure out mode
                    cmd.motorCmd[i].q = init_config_(i);
                    cmd.motorCmd[i].dq = init_vel_(i);
                    cmd.motorCmd[i].tau = 0;
                    cmd.motorCmd[i].Kp = motor_kp; // TODO: tune
                    cmd.motorCmd[i].Kd = motor_kv; // TODO: tune
                }
            }


            udp.SetSend(cmd);       // Send the command
            log_file_ << "[Communication] command sent." << std::endl;
        }
    }

    bool HardwareRobot::VerifyControlAction(const vector_t& control) {
        for (int i = 0; i < control.size(); i++) {
            if (std::abs(control(i)) > 1000) {
                log_file_ << "[Control] Control action at index " << i << " is too large." << std::endl;
                return false;
            }
        }

        return true;
    }

} // hardware