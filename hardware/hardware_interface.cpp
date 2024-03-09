//
// Created by zolkin on 2/25/24.
//

#include "hardware_robot.h"
#include "config_parser.h"

using namespace hardware;

int main() {
    std::string config_file("/home/zolkin/AmberLab/bilevel-gait-gen/hardware/hardware_a1_config.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    const double viz_rate = config.ParseNumber<double>("viz_rate");

    mpc::MPCInfo info;
    info.discretization_steps = config.ParseNumber<double>("discretization_steps");
    info.num_nodes = config.ParseNumber<int>("num_nodes");
    info.num_qp_iterations = config.ParseNumber<int>("num_qp");
    info.friction_coef = config.ParseNumber<double>("friction_coef");
    info.vel_bounds = config.ParseEigenVector("vel_bounds");
    info.joint_bounds_lb = config.ParseEigenVector("joint_bounds_lb");
    info.joint_bounds_ub = config.ParseEigenVector("joint_bounds_ub");
    info.ee_frames = config.ParseStdVector<std::string>("collision_frames");
    info.num_switches = config.ParseNumber<int>("num_switches");
    info.integrator_dt = config.ParseNumber<double>("integrator_dt");
    info.num_contacts = info.ee_frames.size();
    info.force_bound = config.ParseNumber<double>("force_bound");
    info.swing_height = config.ParseNumber<double>("swing_height");
    info.foot_offset = config.ParseNumber<double>("foot_offset");
    info.nom_state = config.ParseEigenVector("init_config");
    info.ee_box_size = config.ParseEigenVector("ee_box_size");
    info.real_time_iters = config.ParseNumber<int>("run_time_iterations");
    switch (config.ParseNumber<int>("mpc_verbosity")) {
        case 0:
            info.verbose = mpc::Nothing;
            break;
        case 1:
            info.verbose = mpc::Timing;
            break;
        case 2:
            info.verbose = mpc::Optimization;
            break;
        case 3:
            info.verbose = mpc::All;
            break;
        default:
            throw std::runtime_error("Not a valid verbosity level for MPC.");
    }


    // Read in the inital config and parse it for MPC.
    const vector_t standing = config.ParseEigenVector("init_config");
    const vector_t init_state = config.ParseEigenVector("srb_init");

    // Create the warm start
    const std::vector<vector_t> warm_start(info.num_nodes+1, init_state);

    // Create the goal state
    vector_t mpc_des_state = config.ParseEigenVector("srb_target");
//    mpc_des_state.head<2>() << config.ParseNumber<double>("x_des"), config.ParseNumber<double>("y_des");
//    mpc_des_state.segment<2>(3) << config.ParseNumber<double>("xdot_des"), config.ParseNumber<double>("ydot_des");

    // TODO: Get this from hardware!
    // Inital guess end effector positions
    std::array<std::array<double, 3>, 4> ee_pos{};
    ee_pos.at(0) = {0.1526, 0.12523, 0.011089};
    ee_pos.at(1) = {0.1526, -0.12523, 0.011089};
    ee_pos.at(2) = {-0.208321844, 0.1363286, 0.01444};
    ee_pos.at(3) = {-0.208321844, -0.1363286, 0.01444};

    std::vector<Eigen::Vector3d> ee_locations(4);
    for (int i = 0; i < ee_locations.size(); i++) {
        for (int j = 0; j < 3; j++) {
            ee_locations.at(i)(j) = ee_pos.at(i).at(j);
        }
    }

//    mpc::MPCSingleRigidBody mpc = CreateMPC(info, config, warm_start, mpc_des_state, init_state, ee_locations);

    // Create the MPC controller (only used here for the visualizer)
    std::unique_ptr<controller::MPCController> mpc_controller;
    mpc_controller = std::make_unique<controller::MPCController>(config.ParseNumber<double>("control_rate"),
                                                                 config.ParseString("robot_urdf"),
                                                                 config.ParseString("foot_type"),
                                                                 config.ParseEigenVector("init_vel").size(),
                                                                 config.ParseEigenVector("torque_bounds"),
                                                                 config.ParseNumber<double>("friction_coef"),
                                                                 config.ParseStdVector<double>("base_pos_gains"),
                                                                 config.ParseStdVector<double>("base_ang_gains"),
                                                                 config.ParseEigenVector("kp_joint_gains"),
                                                                 config.ParseEigenVector("kd_joint_gains"),
                                                                 config.ParseNumber<double>("leg_tracking_weight"),
                                                                 config.ParseNumber<double>("torso_tracking_weight"),
                                                                 config.ParseNumber<double>("force_tracking_weight"),
                                                                 info,
                                                                 warm_start,
                                                                 mpc_des_state,
                                                                 config.ParseNumber<int>("num_polys"),
                                                                 config.ParseEigenVector("Q_srbd_diag").asDiagonal(),
                                                                 config.ParseNumber<int>("gait_opt_freq"),
                                                                 config.ParseString("log_file"));

    const vector_t init_vel = config.ParseEigenVector("init_vel");

    JointGains gains{};
    gains.hip_kp = config.ParseNumber<double>("hip_joint_kp");
    gains.hip_kv = config.ParseNumber<double>("hip_joint_kv");

    gains.thigh_kp = config.ParseNumber<double>("thigh_joint_kp");
    gains.thigh_kv = config.ParseNumber<double>("thigh_joint_kv");

    gains.calf_kp = config.ParseNumber<double>("calf_joint_kp");
    gains.calf_kv = config.ParseNumber<double>("calf_joint_kv");

    HardwareRobot robot(standing, init_vel,
                        config.ParseEigenVector("srb_init"),
                        mpc_controller, 2,
                        gains,1.0/240.0);

    robot.ChangeState(hardware::HardwareRobot::Hold);

//    robot.StartControlLoop();

    const double dt = 0.0005;
    UNITREE_LEGGED_SDK::LoopFunc control_loop("control_loop", dt, boost::bind(&HardwareRobot::ControlCallback, &robot));
    UNITREE_LEGGED_SDK::LoopFunc udp_send("udp_send", dt, 3, boost::bind(&HardwareRobot::SendUDP, &robot));
    UNITREE_LEGGED_SDK::LoopFunc udp_recv("udp_recv", dt, 3, boost::bind(&HardwareRobot::RecieveUDP, &robot));

    udp_send.start();
    udp_recv.start();
    control_loop.start();


    while (true) {
        std::cout << "Enter the next robot state: [H] Hold, [S] Stand, [M] MPC, [T] Testing." << std::endl;
        string next_state;
        std::cin >> next_state;

        if (next_state == "H") {
            robot.ChangeState(hardware::HardwareRobot::Hold);
            std::cout << "Robot entering Hold mode." << std::endl;
        } else if (next_state == "S") {
            robot.ChangeState(hardware::HardwareRobot::Stand);
            std::cout << "Robot entering Stand mode." << std::endl;
        } else if (next_state == "M") {
            robot.ChangeState(hardware::HardwareRobot::MPC);
            std::cout << "Robot entering MPC mode." << std::endl;
        } else if (next_state == "T") {
            robot.ChangeState(hardware::HardwareRobot::Testing);
            std::cout << "Robot entering Testing mode." << std::endl;
        } else {
            std::cout << "Invalid robot state." << std::endl;
        }

        sleep(0.01);
//        robot.ControlCallback();
    }
}
