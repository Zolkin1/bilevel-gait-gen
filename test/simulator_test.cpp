//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

// TODO: write a controller unit test

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include "config_parser.h"
#include "simulator.h"
#include "simulation_robot.h"
#include "pd_grav_comp.h"

TEST_CASE("configuration parser", "[utils]") {
    std::string config_file = "test_config.yaml";
    utils::ConfigParser config = utils::ConfigParser(config_file);

    REQUIRE(config.ParseString("test_string") == "hello");
    REQUIRE(config.ParseNumber<int>("test_num") == 4);
    REQUIRE(config.ParseStringVector("test_string_vec").at(0) == "hello");
    REQUIRE(config.ParseStringVector("test_string_vec").at(1) == "world");

    REQUIRE(config.ParseEigenVector("test_number_vec")(0) == 1);
    REQUIRE(config.ParseEigenVector("test_number_vec")(1) == 2);
    REQUIRE(config.ParseEigenVector("test_number_vec")(2) == 3);
}

TEST_CASE("Basic simulator/interface functions", "[simulator]") {
    // Hard code the yaml
    std::string config_file("a1_configuration_test.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    // Make the low level controller
    std::unique_ptr<controller::Controller> controller =
            std::make_unique<controller::PDGravComp>(config.ParseNumber<double>("control_rate"),
                                                    config.ParseString("robot_urdf"),
                                                    config.ParseString("foot_type"),
                                                    config.ParseEigenVector("standing_config"),
                                                    config.ParseEigenVector("standing_vel"));   // number of inputs


    // Check the controller was create correctly
    int num_inputs = controller->GetNumInputs();
    REQUIRE(num_inputs == 12);
    REQUIRE(controller->GetRate() == 0.0025);

    // Make the robot for simulation
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, controller);


    // Everything in the robot needs to be setup by here because it is consumed by the simulator
    simulator::Simulator sim(robot);
    REQUIRE(simulator::model_ != nullptr);      // Check the mujoco model was made properly

    robot->InitController();
    robot->DefineContacts(config.ParseStringVector("collision_frames"),
                          config.ParseStdVector<int>("collision_bodies"));

    // ----------- check initial conditions ----------- //
    Eigen::VectorXd init_config = config.ParseEigenVector("init_config");
    Eigen::VectorXd init_vel = config.ParseEigenVector("init_vel");
    robot->SetInitialCondition(init_config, init_vel);
    REQUIRE(robot->GetInitConfig() == init_config);
    REQUIRE(robot->GetInitVelocities() == init_vel);


    // Set targets
    robot->UpdateTargetConfig(config.ParseEigenVector("standing_config"));
    robot->UpdateTargetVel(config.ParseEigenVector("standing_vel"));

    // ----------- check joint map ----------- //
    Eigen::VectorXd pin_con = robot->ConvertMujocoVecConfigToPinocchio(init_config);
    Eigen::VectorXd pin_vel = robot->ConvertMujocoVecVelLikeToPinocchio(init_vel);

    REQUIRE(robot->ConvertPinocchioJointToMujoco(pin_con.tail(num_inputs)) == init_config.tail(num_inputs));
    REQUIRE(robot->ConvertPinocchioVelToMujoco(pin_vel) == init_vel);

    // Check the data is created
    sim.SetupSimulator(robot);
    REQUIRE(simulator::data_ != nullptr);

    // ----------- check joint map through mujoco structures ----------- //
    for (int i = 0; i < init_config.size(); i++) {
        simulator::data_->qpos[i] = init_config(i);
    }
    Eigen::VectorXd pin_con2 = robot->ConvertMujocoConfigToPinocchio(simulator::data_).tail(num_inputs);
    REQUIRE(robot->ConvertPinocchioJointToMujoco(pin_con2) == init_config.tail(num_inputs));

    for (int i = 0; i < init_vel.size(); i++) {
        simulator::data_->qvel[i] = init_vel(i);
    }
    Eigen::VectorXd pin_vel2 = robot->ConvertMujocoVelToPinocchio(simulator::data_).tail(num_inputs);
    REQUIRE(robot->ConvertPinocchioJointToMujoco(pin_vel2) == init_vel.tail(num_inputs));
}