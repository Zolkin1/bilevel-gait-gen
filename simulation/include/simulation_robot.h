//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H
#define BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H

#include <string>
#include <Eigen/Core>
#include <memory>

#include "mujoco.h"

#include "controller.h"

/**
 * SimulationRobot class.
 * Defines a robot to be used by the simulator. Holds all the information needed by the simulator.
 */
namespace simulator {

    class SimulationRobot {
    public:
        SimulationRobot(const std::string& robot_xml_path, std::unique_ptr<Controller>& controller);

        void SetInitialCondition(const Eigen::VectorXd& initial_condition, const Eigen::VectorXd& initial_vel);

        std::string GetRobotXMLFile() const;

        Eigen::VectorXd GetConfig() const;

        Eigen::VectorXd GetVelocities() const;

        /**
         * Interface with Mujoco to provide the current control action.
         */
        void GetControlAction(const mjModel* model, mjData* data);
    private:
        std::string robot_xml_path_;
        std::unique_ptr<Controller> low_level_controller_;

        Eigen::VectorXd initial_config_;
        Eigen::VectorXd initial_vel_;
    };
} // simulator


#endif //BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H