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

    mpc::MPC mpc(info, config.ParseString("robot_urdf"));


    // Read in the inital config and parse it for MPC.
    vector_t standing = config.ParseEigenVector("init_config");
    standing.segment<4>(3) = ConvertMujocoQuatToPinocchioQuat(standing.segment<4>(3));
    vector_t init_state = vector_t::Zero(6 + standing.size());
    init_state.tail(standing.size()) = standing;

    // Create the warm start
    std::vector<vector_t> warm_start(info.num_nodes+1, init_state);

    // Create the goal state
    vector_t mpc_des_state = init_state;
    mpc_des_state.segment<2>(6) << 1, 0;

    // Inital guess end effector positions
    std::array<std::array<double, 3>, 4> ee_pos{};
    ee_pos.at(0) = {0.2, 0.2, 0};
    ee_pos.at(1) = {0.2, -0.2, 0};
    ee_pos.at(2) = {-0.2, 0.2, 0};
    ee_pos.at(3) = {-0.2, -0.2, 0};

    // Set warm starts and defaults
    mpc.SetDefaultGaitTrajectory(mpc::Gaits::Trot, 2, ee_pos);
    mpc.SetStateTrajectoryWarmStart(warm_start);

    // Create weights
    // TODO: Can probably tune this to get better performance
    matrix_t Q = 20*matrix_t::Identity(24, 24);
    Q.topLeftCorner<6,6>() = matrix_t::Zero(6,6);
    Q(6,6) = 300; //30;
    Q(7,7) = 300; //30;
    Q(8,8) = 450; //10;
    Q(9,9) = 50;
    Q(10,10) = 50;
    Q(11,11) = 50;
    Q(12,12) = 50;

    // TODO: Weight shoulder joints to keep the legs in line

    // Desried state in the lie algebra
    const vector_t des_alg = mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(mpc_des_state, init_state);
    std::cout << des_alg << std::endl;

    // Add in costs
    mpc.AddQuadraticTrackingCost(des_alg, Q);
    mpc.AddForceCost(0.01);  // Note: NEED to adjust this based on the number of nodes otherwise it is out-weighed
    mpc.SetQuadraticFinalCost(50*Q);
    mpc.SetLinearFinalCost(-50*Q*des_alg);

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
                                                                 config.ParseStdVector<double>("joint_gains"),
                                                                 config.ParseNumber<double>("leg_tracking_weight"),
                                                                 config.ParseNumber<double>("torso_tracking_weight"),
                                                                 config.ParseNumber<double>("force_tracking_weight"),
                                                                 info,
                                                                 warm_start,
                                                                 mpc_des_state);

    // Make the robot for visualization
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);

    for (int i = 0; i < 7; i++) {
        mpc.Solve(init_state, 0);
    }
    mpc.PrintStats();

    // Visualize results
    simulation::Visualizer viz(config.ParseString("robot_xml"));
    robot->SetSimModel(viz.GetModel());
    for (int i = 0; i < info.num_nodes+1; i++) {
        vector_t temp_state = mpc.GetFullTargetState(i*info.integrator_dt);
        viz.UpdateState(robot->ConvertPinocchioConfigToMujoco(mpc.GetTargetConfig(i*info.integrator_dt)));
        viz.UpdateViz(config.ParseNumber<double>("viz_rate"));
//        mpc.Solve(temp_state, i*info.integrator_dt);
    }

    // Print the final trajectory to a file for viewing
    mpc::Trajectory traj = mpc.GetTrajectory();
    traj.PrintTrajectoryToFile("demo_final_traj.txt");

}