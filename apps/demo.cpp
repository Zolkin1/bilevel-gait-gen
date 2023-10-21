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

#include "yaml-cpp/yaml.h"

#include "simulator.h"
#include "simulation_robot.h"

int main(int argc, char* argv[]) {
    // Load in the config
    if (argc != 2) {
        std::cerr << "Wrong number of input arguments. Expecting one argument for the name of the yaml file." << std::endl;
        return 1;
    }

    // Needs the absolute path for now
    std::string config_file(argv[1]);
    YAML::Node config = YAML::LoadFile(config_file);

    // Make the low level controller
    std::unique_ptr<simulator::Controller> controller = std::make_unique<simulator::Controller>();

    // Make the robot for simulation
    std::string robot_file = config["robot_xml"].as<std::string>();
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, controller);

    // Set the robot's initial condition
    // TODO: parse better and check the quaternion is normalized
    std::vector<double> temp_config = config["init_config"].as<std::vector<double>>();
    std::vector<double> temp_vel = config["init_vel"].as<std::vector<double>>();
    Eigen::VectorXd init_config(temp_config.size());
    Eigen::VectorXd init_vel(temp_vel.size());
    for (int i = 0; i < temp_config.size(); i++) {
        init_config(i) = temp_config.at(i);
    }

    for (int i = 0; i < temp_vel.size(); i++) {
        init_vel(i) = temp_vel.at(i);
    }


    robot->SetInitialCondition(init_config, init_vel);

    // Everything in the robot needs to be setup by here because it is consumed by the simulator
    simulator::Simulator sim(robot, 1.0/500.0);

    sim.SetupSimulator();

    sim.RunSimulator();

    // Read in robot file (command line argument)

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
