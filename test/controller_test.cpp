//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include "config_parser.h"
#include "simulator.h"
#include "simulation_robot.h"
#include "pd_grav_comp.h"

TEST_CASE("configuration parser", "[utils]") {
    std::string config_file = "/home/zach/AmberLab/bilevel-gait-generation/bilevel-gait-gen/test/test_config.yaml";
    mpc::utils::ConfigParser config = mpc::utils::ConfigParser(config_file);

    REQUIRE(config.ParseString("test_string") == "hello");
    REQUIRE(config.ParseNumber("test_num") == 4);
    REQUIRE(config.ParseStringVector("test_string_vec").at(0) == "hello");
    REQUIRE(config.ParseStringVector("test_string_vec").at(1) == "world");

    REQUIRE(config.ParseEigenVector("test_number_vec")(0) == 1);
    REQUIRE(config.ParseEigenVector("test_number_vec")(1) == 2);
    REQUIRE(config.ParseEigenVector("test_number_vec")(2) == 3);
}

TEST_CASE("Basic simulator functions", "[simulator]") {
    // Hard code the yaml
    std::string config_file("/home/zach/AmberLab/bilevel-gait-generation/bilevel-gait-gen/test/a1_configuration.yaml");
    mpc::utils::ConfigParser config = mpc::utils::ConfigParser(config_file);

    // Make the low level controller
    std::unique_ptr<simulator::Controller> controller =
            std::make_unique<simulator::PDGravComp>(config.ParseNumber("control_rate"),
                                                    config.ParseString("robot_urdf"),
                                                    config.ParseString("foot_type"),
                                                    config.ParseEigenVector("standing_config"),
                                                    config.ParseEigenVector("standing_vel"));   // number of inputs

    controller->DefineContacts(config.ParseStringVector("collision_frames"),
                               config.ParseIntVector("collision_bodies"));

    // Check the controller was create correctly
    REQUIRE(controller->GetNumInputs() == 12);
    REQUIRE(controller->GetRate() == 0.0025);

    // Make the robot for simulation
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, controller);

    // Set the robot's initial condition
    Eigen::VectorXd init_config = config.ParseEigenVector("init_config");
    Eigen::VectorXd init_vel = config.ParseEigenVector("init_vel");
    robot->SetInitialCondition(init_config, init_vel);

    // Get the standing configuration
    Eigen::VectorXd standing_config = config.ParseEigenVector("standing_config");

    // Everything in the robot needs to be setup by here because it is consumed by the simulator
    simulator::Simulator sim(robot);

    sim.SetupSimulator();

    // Check the mujoco interface worked properly
    REQUIRE(simulator::model_ != nullptr);
    REQUIRE(simulator::data_ != nullptr);

    //sim.RunSimulator();


}