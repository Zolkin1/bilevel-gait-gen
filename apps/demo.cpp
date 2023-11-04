//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

/*
 * Demo file. This should show how to setup the controller and connect it to a simulation.
 * I should be able to take in a urdf file (or maybe a mujoco file) and create the model, create the controller,
 * create the simulation and run it all.
 */

#include <iostream>
#include <string>
#include <Eigen/Core>

#include "config_parser.h"
#include "simulator.h"
#include "simulation_robot.h"
#include "pd_grav_comp.h"
#include "qp_control.h"


int main(int argc, char* argv[]) {
    // Load in the config
    if (argc != 2) {
        std::cerr << "Wrong number of input arguments. Expecting one argument for the name of the yaml file." << std::endl;
        return 1;
    }

    // Needs the absolute path for now
    std::string config_file(argv[1]);
    mpc::utils::ConfigParser config = mpc::utils::ConfigParser(config_file);

    // Make the low level controller
    std::unique_ptr<simulator::Controller> controller;
    if(config.ParseString("controller_type") == "PD_GRAV_COMP") {
        controller = std::make_unique<simulator::PDGravComp>(config.ParseNumber("control_rate"),
                                                            config.ParseString("robot_urdf"),
                                                            config.ParseString("foot_type"),
                                                            config.ParseEigenVector("standing_config"),
                                                            config.ParseEigenVector("standing_vel"));
    } else if (config.ParseString("controller_type") == "QP_CONTROL") {
        controller = std::make_unique<simulator::QPControl>(config.ParseNumber("control_rate"),
                                                            config.ParseString("robot_urdf"),
                                                            config.ParseString("foot_type"),
                                                            config.ParseEigenVector("init_vel").size(),
                                                            config.ParseEigenVector("torque_bounds"),
                                                            config.ParseNumber("friction_coef"),
                                                            config.ParseStdVector<double>("base_pos_gains"),
                                                            config.ParseStdVector<double>("base_ang_gains"),
                                                            config.ParseStdVector<double>("joint_gains"),
                                                            config.ParseNumber("leg_tracking_weight"),
                                                            config.ParseNumber("torso_tracking_weight"),
                                                            config.ParseNumber("force_tracking_weight"));

    } else {
        throw std::runtime_error("Invalid controller type in the yaml file.");
    }


    controller->DefineContacts(config.ParseStringVector("collision_frames"),
                                    config.ParseStdVector<int>("collision_bodies"));

    controller->PrintConfigNames();

    // TODO: Restructure mujoco and controller interface

    // Make the robot for simulation
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, controller);

    robot->UpdateTargetConfig(config.ParseEigenVector("standing_config"));
    robot->UpdateTargetVel(config.ParseEigenVector("standing_vel"));

    // Set the robot's initial condition
    // TODO: check the quaternion is normalized
    Eigen::VectorXd init_config = config.ParseEigenVector("init_config");
    Eigen::VectorXd init_vel = config.ParseEigenVector("init_vel");
    robot->SetInitialCondition(init_config, init_vel);

    // Get the standing configuration
    Eigen::VectorXd standing_config = config.ParseEigenVector("standing_config");

    // Everything in the robot needs to be setup by here because it is consumed by the simulator
    simulator::Simulator sim(robot);

    sim.SetupSimulator();

    sim.RunSimulator();


    // Create a model to be used by MPC (this will let the user specify anything not in the file)
    // this class will wrap pinocchio and provide everything the MPC needs

    // Create the controller
    // - Define cost functions and constraints
    // - Specify how the constraints will be treated

    // The simulation should take in a SimulationRobot class object. This object
    // should have the xml file for Mujoco, the MPC controller, and the low level controller.
    // We can probably assume that the PD control is dealt with in the Mujoco XML. This object also
    // deals with all the internal connections between the MPC and robot etc...

    // Initialize the simulation

    // Run the controller and retrieve diagnostics

    return 0;
}
