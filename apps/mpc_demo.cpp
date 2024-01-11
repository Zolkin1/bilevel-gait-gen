//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"
#include "mpc.h"
#include "inputs.h"
#include "spline.h"
#include "mpc_controller.h"
#include "visualization.h"
#include "simulation_robot.h"

using vector_t = Eigen::VectorXd;
using matrix_t = Eigen::MatrixXd;

Eigen::Vector4d ConvertMujocoQuatToPinocchioQuat(const Eigen::Vector4d& quat) {
    Eigen::Vector4d pin_quat;
    pin_quat(3) = quat(0);
    pin_quat(0) = quat(1);
    pin_quat(1) = quat(2);
    pin_quat(2) = quat(3);

    return pin_quat;
}

int main() {
    // TODO: when I make the file local via the cmake file then I only get it copied over after a compilation
    std::string config_file("/home/zolkin/AmberLab/bilevel-gait-gen/apps/a1_configuration.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

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

    mpc::MPC mpc(info, config.ParseString("robot_urdf"));


    // Read in the inital config and parse it for MPC.
    vector_t standing = config.ParseEigenVector("init_config");
    vector_t init_state = vector_t::Zero(6 + standing.size());
    init_state.tail(standing.size()) = standing;

    // Create the warm start
    std::vector<vector_t> warm_start(info.num_nodes+1, init_state);

    // Create the goal state
    vector_t mpc_des_state = init_state;
    mpc_des_state.segment<2>(6) << config.ParseNumber<double>("x_des"), config.ParseNumber<double>("y_des");

    // Inital guess end effector positions
    std::array<std::array<double, 3>, 4> ee_pos{};
    ee_pos.at(0) = {0.2, 0.2, 0};
    ee_pos.at(1) = {0.2, -0.2, 0};
    ee_pos.at(2) = {-0.2, 0.2, 0};
    ee_pos.at(3) = {-0.2, -0.2, 0};

    // Set warm starts and defaults
    mpc.SetDefaultGaitTrajectory(mpc::Gaits::Trot, config.ParseNumber<int>("num_polys"), ee_pos);
    mpc.SetStateTrajectoryWarmStart(warm_start);

    // Create weights
    matrix_t Q = 2*matrix_t::Identity(24, 24);
    Q.topLeftCorner<6,6>() = matrix_t::Zero(6, 6);
    Q(0,0) = 0.0;
    Q(1,1) = 0.0;
    Q(2,2) = 10;
    Q(6,6) = 300; //300;
    Q(7,7) = 300; //300;
    Q(8,8) = 750; //450;
    Q(9,9) = 200;
    Q(10,10) = 200;
    Q(11,11) = 200;

    Q(12,12) = 600;
    Q(15,15) = 600;
    Q(18,18) = 600;
    Q(21,21) = 600;

    Q(13,13) = 50;
    Q(16,16) = 50;
    Q(19,19) = 50;
    Q(22,22) = 50;

    Q(14,14) = 10;
    Q(17,17) = 10;
    Q(20,20) = 10;
    Q(23,23) = 10;

    // Desried state in the lie algebra
    const vector_t des_alg = mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(mpc_des_state, init_state);
    std::cout << des_alg << std::endl;

    // Add in costs
    mpc.AddQuadraticTrackingCost(des_alg, Q);
    mpc.AddForceCost(0.001);  // Note: NEED to adjust this based on the number of nodes otherwise it is out-weighed
    mpc.SetQuadraticFinalCost(5000*Q);
    mpc.SetLinearFinalCost(-5000*Q*des_alg);

    // Create the MPC controller (only used here for the visualizer)
    std::unique_ptr<controller::Controller> mpc_controller;
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
                                                                 config.ParseNumber<int>("num_polys"));

    // Make the robot for visualization
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);

    mpc.CreateInitialRun(init_state);
    mpc.PrintStats();

    // Visualize results
    simulation::Visualizer viz(config.ParseString("robot_xml"));
    robot->SetSimModel(viz.GetModel());
    for (int i = 0; i < info.num_nodes + 200; i++) {
        vector_t temp_state = mpc.GetFullTargetState(i*info.integrator_dt);
//        temp_state(14) = temp_state(14) + 0.01;
//        temp_state(17) = temp_state(17) - 0.01;
//        temp_state(21) = temp_state(21) - 0.04;
//        temp_state(22) = temp_state(22) + 0.01;
        viz.UpdateState(robot->ConvertPinocchioConfigToMujoco(mpc.GetTargetConfig(i*info.integrator_dt)));
        viz.GetTrajViz(mpc.GetTrajectory().CreateVizData(mpc.GetModel()));
        viz.UpdateViz(config.ParseNumber<double>("viz_rate"));
        mpc.GetRealTimeUpdate(config.ParseNumber<int>("run_time_iterations"), temp_state, i*info.integrator_dt);
//        mpc.GetRealTimeUpdate(config.ParseNumber<int>("run_time_iterations"), temp_state, i*info.integrator_dt);
    }

    // Print the final trajectory to a file for viewing
    mpc::Trajectory traj = mpc.GetTrajectory();
    traj.PrintTrajectoryToFile("demo_final_traj.txt");

    mpc.PrintStats();

}