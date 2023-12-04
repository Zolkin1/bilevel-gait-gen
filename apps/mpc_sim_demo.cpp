//
// Created by zolkin on 11/30/23.
//

#include "config_parser.h"
#include "mpc.h"
#include "inputs.h"
#include "spline.h"
#include "mpc_controller.h"
#include "simulator.h"
#include "simulation_robot.h"


using vector_t = Eigen::VectorXd;

Eigen::Vector4d ConvertMujocoQuatToPinocchioQuat(const Eigen::Vector4d& quat) {
    Eigen::Vector4d pin_quat;
    pin_quat(3) = quat(0);
    pin_quat(0) = quat(1);
    pin_quat(1) = quat(2);
    pin_quat(2) = quat(3);

    return pin_quat;
}

int main() {
    std::string config_file("/home/zolkin/AmberLab/bilevel-gait-gen/apps/a1_configuration.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    mpc::MPCInfo info;
    info.discretization_steps = config.ParseNumber<double>("discretization_steps");
    info.num_nodes = config.ParseNumber<int>("num_nodes");
    info.time_horizon = config.ParseNumber<double>("time_horizon");
    info.num_qp_iterations = config.ParseNumber<int>("num_qp");
    info.friction_coef = config.ParseNumber<double>("friction_coef");
    info.vel_bounds = config.ParseEigenVector("vel_bounds");
    info.joint_bounds = config.ParseEigenVector("joint_bounds");
    info.ee_frames = config.ParseStdVector<std::string>("collision_frames");
    info.num_switches = config.ParseNumber<int>("num_switches");
    info.integrator_dt = config.ParseNumber<double>("integrator_dt");
    info.num_contacts = info.ee_frames.size();


    vector_t standing = config.ParseEigenVector("standing_config");
    standing.segment<4>(3) = ConvertMujocoQuatToPinocchioQuat(standing.segment<4>(3));
    vector_t state = vector_t::Zero(6 + standing.size());
    state.tail(standing.size()) = standing;

    std::vector<vector_t> warm_start_states;
    for (int i = 0; i < info.num_nodes+1; i++) {
        warm_start_states.push_back(state);
    }

    // TODO: Won't work until the MPC is more robust (numerically)
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
                                              warm_start_states,
                                              state);

    // Make the robot for simulation
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);

    // Create a simulation object and populate the robot with mujoco data
    simulator::Simulator sim(robot);

    // Set up controller solver
    robot->InitController();
    robot->DefineContacts(config.ParseStringVector("collision_frames"),
                          config.ParseStdVector<int>("collision_bodies"));

    // Set the robot's initial condition
    Eigen::VectorXd init_config = config.ParseEigenVector("init_config");
    Eigen::VectorXd init_vel = config.ParseEigenVector("init_vel");
    robot->SetInitialCondition(init_config, init_vel);
//    robot->UpdateTargetConfig(config.ParseEigenVector("standing_config"));
//    robot->UpdateTargetVel(config.ParseEigenVector("standing_vel"));

    // Setup sim
    sim.SetupSimulator(robot);

    // Run sim
    sim.RunSimulator(robot);


}